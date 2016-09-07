/*
 * METACORTEX
 * Copyright 2011-2013 Richard Leggett
 *
 * Based on code from CORTEX
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 *
 * CORTEX project contacts:
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * CORTEX Development team:
 *       R. Ramirez-Gonzalez (Ricardo.Ramirez-Gonzalez@bbsrc.ac.uk)
 *       R. Leggett (richard@leggettnet.org.uk)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "global.h"
#include "binary_kmer.h"
#include "flags.h"
#include "element.h"
#include "seq.h"
#include "open_hash/hash_table.h"
#include "file_reader.h"
#include "dB_graph.h"
#include "logger.h"
#include "graph_tools.h"
#include "graph_formats.h"
#include "node_queue.h"
#include "coverage_walk.h"
#include "perfect_path.h"
#include "metagraphs.h"
#include "cleaning.h"

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
int grow_graph_from_node(dBNode* start_node, dBNode** best_node, dBGraph* graph, Queue* graph_queue)
{
    Queue* nodes_to_walk;
    dBNode* node;
    int orientation;
    int depth;
    int current_graph_size = 0;
    int best_coverage = 0;
    int best_edges = 0;

    *best_node = 0;

    // Nucleotide iterator, used to walk all possible paths from a node
    void walk_if_exists(Nucleotide n) {
        //if (debug) printf("Trying nucleotide %i\n", n);

        // If there is an edge in any colour for this nucleotide...
        if (db_node_edge_exist_any_colour(node, n, orientation)) {

            //if (debug) printf("  Edge exists\n");

            // Get first node along this edge and check we've not already visited it...
            Orientation next_orientation;
            Nucleotide reverse_nucleotide;
            dBNode * next_node;
            next_node = db_graph_get_next_node(node, orientation, &next_orientation, n, &reverse_nucleotide, graph);
            if (!next_node) {
                log_and_screen_printf("Error: Something went wrong with db_graph_get_next_node\n");
                exit(-1);
            }

            // If not already visited the first node, walk it...
            if (!db_node_check_flag_visited(next_node)) {
                pathStep first_step;
                Path * new_path;
                dBNode* end_node;
                int i = 0;

                // Get path
                first_step.node = node;
                first_step.orientation = orientation;
                first_step.label = n;
                new_path = path_new(MAX_EXPLORE_NODES, graph->kmer_size);
                if (!new_path) {
                    log_and_screen_printf("ERROR: Not enough memory to allocate new path.\n");
                    exit(-1);
                }

                db_graph_get_perfect_path_with_first_edge_all_colours(&first_step, &db_node_action_do_nothing, new_path, graph);

                // Add end node to list of nodes to visit
                end_node = new_path->nodes[new_path->length-1];
                if (!db_node_check_flag_visited(end_node)) {
                    if (!db_node_is_blunt_end_all_colours(end_node, new_path->orientations[new_path->length-1])) {
                        if (queue_push_node(nodes_to_walk, end_node, depth+1) == NULL) {
                            log_and_screen_printf("Queue too large. Ending.\n");
                            exit(1);
                        }
                    }
                }

                // Now go through all nodes, look for best and mark all as visited
                for (i=0; i<new_path->length; i++) {
                  //  MARTIN : new_path established where?
                    if (!db_node_check_flag_visited(new_path->nodes[i])) {
                        int this_coverage = element_get_coverage_all_colours(new_path->nodes[i]);
                        int this_edges = db_node_edges_count_all_colours(new_path->nodes[i], forward) + db_node_edges_count_all_colours(new_path->nodes[i], reverse);

                        if ((best_node == 0) ||
                            (this_coverage > best_coverage) ||
                            ((this_coverage == best_coverage) && (this_edges < best_edges)))
                        {
                            best_coverage = this_coverage;
                            best_edges = this_edges;
                            *best_node = new_path->nodes[i];
                        }

                        db_node_action_set_flag_visited(new_path->nodes[i]);
                        queue_push(graph_queue, new_path->nodes[i]);
                        current_graph_size++;
                    }
                }

                // Clean up
                path_destroy(new_path);
            }
        }
    }

    // Start a queue of nodes to walk
    //log_and_screen_printf("Allocating %d Mb to store queue information (max %d nodes, when full each node could be %d)...\n", ((METACORTEX_QUEUE_SIZE * sizeof(QueueItem*)) / 1024) / 1024, METACORTEX_QUEUE_SIZE, sizeof(QueueItem));
    nodes_to_walk = queue_new(METACORTEX_QUEUE_SIZE);
    if (!nodes_to_walk) {
        log_and_screen_printf("Couldn't get memory for node queue.\n");
        exit(-1);
    }

    // Add start node to list of nodes to visit
    if (queue_push_node(nodes_to_walk, start_node, 0) == NULL) {
        log_and_screen_printf("Queue too large. Ending.\n");
        exit(-1);
    }

    if (!db_node_check_flag_visited(start_node)) {
        db_node_action_set_flag_visited(start_node);
        current_graph_size++;
    }

    // Now keep visiting nodes and walking paths
    while (nodes_to_walk->number_of_items > 0) {
        // Take top node from list
        node = queue_pop_node(nodes_to_walk, &depth);

        // Look at all paths out from here
        orientation = forward;
        nucleotide_iterator(&walk_if_exists);
        orientation = reverse;
        nucleotide_iterator(&walk_if_exists);
    }

    queue_free(nodes_to_walk);

    // If we didn't find a start node, presumably this is a singleton?
    if (*best_node == 0) {
        //log_printf("Note: didn't find a best node, setting to start node\n");
        *best_node = start_node;
    }

    return current_graph_size;
}

void metacortex_find_subgraphs(dBGraph* graph, char* consensus_contigs_filename, int min_subgraph_kmers, int min_contig_length, boolean multiple_subgraph_contigs)
{
    FILE* fp_analysis;
    FILE* fp_contigs;
    Queue* graph_queue;
    Path *path_fwd = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);
    Path *path_rev = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);
    Path *final_path = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);
    char seq[256];
    char analysis_filename[MAX_EXPLORE_PATH_LENGTH];
    long int total_nodes = 0;
    int n_seeds = 0;
    int counter= 0;

    graph_queue = queue_new(METACORTEX_QUEUE_SIZE);
    if (!graph_queue) {
        log_and_screen_printf("Couldn't get memory for graph queue.\n");
        exit(-1);
    }

    sprintf(analysis_filename, "%s.analysis", consensus_contigs_filename);
    log_and_screen_printf("Running metacortex subgraph analysis...\n");
    log_and_screen_printf("          Contig file: %s\n", consensus_contigs_filename);
    log_and_screen_printf("        Analysis file: %s\n", analysis_filename);
    log_and_screen_printf("Minimum subgraph size: %i\n", min_subgraph_kmers);
    log_and_screen_printf("Minimum contig length: %i\n", min_contig_length);

    /* Initialise temporaray path array buffers */
    path_array_initialise_buffers(graph->kmer_size);

    /* Open the analysis file */
    fp_analysis = fopen(analysis_filename, "w");
    if (!fp_analysis) {
        log_and_screen_printf("ERROR: Can't open analysis file.\n");
        exit(-1);
    }

    /* Open consensus contigs file */
    fp_contigs = fopen(consensus_contigs_filename, "w");
    if (!fp_contigs) {
        log_and_screen_printf("ERROR: Can't open contig file.\n");
        exit(-1);
    }

    /* For each node, if it's not pruned or visited, try and grow a graph */
    void explore_node(dBNode * node) {
        if (node == NULL) {
            log_and_screen_printf("Error: NULL node passed to explore_node.\n");
            exit(-1);
        }

        if (db_node_check_for_any_flag(node, PRUNED | VISITED) == false) {
            dBNode* seed_node;
            int nodes_in_graph;

            /* Grow graph from this node, returning the 'best' (highest coverage) node to store as seed point */
            log_printf("Growing graph from node\n");
            graph_queue->number_of_items = 0;
            nodes_in_graph = grow_graph_from_node(node, &seed_node, graph, graph_queue);
            total_nodes += nodes_in_graph;

            if (seed_node == NULL) {
                printf("ERROR: Seed node is NULL, nodes in graph is %d\n", nodes_in_graph);
            } else {
                /* Write data to analysis file */
                binary_kmer_to_seq(&(node->kmer), graph->kmer_size, seq);
                fprintf(fp_analysis, "%i\t%i\t%ld\t%s\t", n_seeds, nodes_in_graph, total_nodes, seq);
                binary_kmer_to_seq(&(seed_node->kmer), graph->kmer_size, seq);
                fprintf(fp_analysis, "%s\n", seq);

                /* Enough nodes to bother with? If so, get consensus contig */
                if (nodes_in_graph >= min_subgraph_kmers) {
                    dBNode* queue_node;
                    int pi;

                    binary_kmer_to_seq(&(seed_node->kmer), graph->kmer_size, seq);
                    coverage_walk_get_path(seed_node, forward, NULL, graph, path_fwd);
                    coverage_walk_get_path(seed_node, reverse, NULL, graph, path_rev);
                    path_reverse(path_fwd, final_path);
                    path_append(final_path, path_rev);
                    final_path->id = counter;
                    if (final_path->length >= (min_contig_length - graph->kmer_size)) {
                        log_printf("Write path of size %d\n", final_path->length);
                        log_printf("graph size\t%i\n",nodes_in_graph);
                        path_to_fasta(final_path, fp_contigs);
                    } else {
                        log_printf("Didn't write path of size %d\n", final_path->length);
                    }

                    if (multiple_subgraph_contigs) {
                        /* Now clear visited flags for subgraph */
                        while (graph_queue->number_of_items > 0) {
                            queue_node = (dBNode*)queue_pop(graph_queue);
                            db_node_action_unset_flag(queue_node, VISITED);
                        }

                        /* Now disconnect path from other nodes and mark path as visited, so it's not visited again */
                        for (pi=0; pi<final_path->length; pi++) {
                            cleaning_prune_db_node(final_path->nodes[pi], graph);
                            db_node_action_set_flag(final_path->nodes[pi], VISITED);
                        }
                    }

                    /* Reset paths */
                    path_reset(path_fwd);
                    //perfect_path_get_path(seed_node, forward, &db_node_action_do_nothing, graph, path_fwd);
                    path_reset(path_rev);
                    path_reset(final_path);
                } else {
                    log_printf("  Number of nodes (%i) too small. Not outputting contig.\n", nodes_in_graph);
                }

                counter++;
            }
        }
    }

    /* Traverse each node... */
    db_graph_reset_flags(graph);
    log_and_screen_printf("Finding subgraphs the new way!...\n");
    hash_table_traverse(&explore_node, graph);
    log_and_screen_printf("Finished. Total nodes: %ld\n", total_nodes);

    /* Close files */
    fclose(fp_contigs);
    fclose(fp_analysis);
}
