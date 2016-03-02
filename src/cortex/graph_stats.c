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
#include "graph_stats.h"
#include "cleaning.h"
#include "metagraphs.h"

#define COVERAGE_BINS 150
#define COVERAGE_BIN_SIZE 1
#define MAX_BRANCHES 5

typedef struct {
    int total_size;
    int branch_nodes;
} GraphInfo;

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/

 int grow_graph_from_node_stats(dBNode* start_node, dBNode** best_node, dBGraph* graph, Queue* graph_queue, GraphInfo* nodes_in_graph)
 {
     Queue* nodes_to_walk;
     dBNode* node;
     int orientation;
     int depth;
     int best_coverage = 0;
     int best_edges = 0;

     *best_node = 0;

     // Nucleotide iterator, used to walk all possible paths from a node
     void walk_if_exists(Nucleotide n) {
         //if (debug) printf("Trying nucleotide %i\n", n);
   			int end_orientation;

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
           			 end_orientation = new_path->orientations[new_path->length - 1];
                 if (!db_node_check_flag_visited(end_node)) {
                     if (!db_node_is_blunt_end_all_colours(end_node, new_path->orientations[new_path->length-1])) {
                         if (queue_push_node(nodes_to_walk, end_node, depth+1) == NULL) {
                             log_and_screen_printf("Queue too large. Ending.\n");
                             exit(1);
                         }
                     }
                 }


                // check nodes in path now
                // only really need to check final node as it's a perfect path
                // is it blunt? has it been seen before?
                // things I@m dropping for now - loop detection, branch + loop, catching too large a bubble
                if (db_node_is_blunt_end_all_colours(end_node, end_orientation)) {
                // DO NOTHING WITH THIS
          				//db_graph_check_and_add_path(merged_path, patharray);
               }
               if (db_node_check_flag_visited(end_node)) {
                 // need to count back from here to original branching point?
               }

                 // Now go through all nodes, look for best and mark all as visited
                 for (i=0; i<new_path->length; i++) {
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

                        if (db_node_check_for_any_flag(new_path->nodes[i], BRANCH_NODE_FORWARD | BRANCH_NODE_REVERSE | X_NODE)){
                          nodes_in_graph->branch_nodes++;
                        }

                         db_node_action_set_flag_visited(new_path->nodes[i]);
                         queue_push(graph_queue, new_path->nodes[i]);
                         nodes_in_graph->total_size++;
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

     if (db_node_check_flag_visited(start_node)) {
         db_node_action_set_flag_visited(start_node);
         nodes_in_graph->total_size++;
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

     return 0;
 }


// ----------------------------------------------------------------------
// Work through graph, count cov, X, Y nodes
// ----------------------------------------------------------------------
void find_subgraph_stats(dBGraph * graph, char* consensus_contigs_filename)
{
  FILE* fp_analysis;
	//int branchNodes = 0;
	int X_Nodes_tot[8];
	int Y_Nodes_rev[4];
	int Y_Nodes_for[4];
  long int Contig_Branches[MAX_BRANCHES];
	int X_Nodes = 0;
	int Y_Nodes = 0;
  char* seq = calloc(256, 1);
  long int total_nodes = 0;
  GraphInfo* nodes_in_graph;
  // array to bin coverage 0-5, 5-10, 10-15..95-100
  long int Coverage_Dist[COVERAGE_BINS]; // will this work?
  //int COVERAGE_CEILING = (COVERAGE_BINS-1) * COVERAGE_BIN_SIZE;
  char analysis_filename[strlen(consensus_contigs_filename) + 10];

  Queue* graph_queue;
  int i;
  // Initialise Coverage_Dist  int i;
  for(i=0;i<MAX_BRANCHES;i++){
    Contig_Branches[i]=0;
  }
  for(i=0;i<COVERAGE_BINS;i++){
    Coverage_Dist[i]=0;
  }
  for(i=0;i<4;i++){
    X_Nodes_tot[i*2]=0;
    X_Nodes_tot[(i*2)+1]=0;
    Y_Nodes_for[i]=0;
    Y_Nodes_rev[i]=0;
  } // should I be using pointers and malloc/calloc instead?

  /* Open the analysis file */
  sprintf(analysis_filename, "%s.analysis", consensus_contigs_filename);
  fp_analysis = fopen(analysis_filename, "w");
  if (!fp_analysis) {
      log_and_screen_printf("ERROR: Can't open analysis file.\n");
      exit(-1);
  }

	db_graph_reset_flags(graph);

	// Hash table iterator to label nodes
	void identify_branch_nodes(dBNode * node) {
		//if (!db_node_check_flag_visited(node)) {
    int this_coverage = element_get_coverage_all_colours(node);
    int edges_forward= db_node_edges_count_all_colours(node, forward);
    int edges_reverse = db_node_edges_count_all_colours(node, reverse);
    if (this_coverage<=0) {
        log_and_screen_printf("Error: Coverage is <1 in the graph?\n");
        exit(-1);
    }
    this_coverage = ((this_coverage-1) / COVERAGE_BIN_SIZE);
    if(this_coverage>COVERAGE_BINS){
      this_coverage = COVERAGE_BINS-1;
    }

    Coverage_Dist[this_coverage]++;
    total_nodes++;

		// Look for Y shape branch forward orientation
		// The nodes at the top of the Y should contain different colours
		if (edges_forward > 1
		    && edges_reverse == 1) {
			db_node_action_set_flag(node, BRANCH_NODE_FORWARD);
      Y_Nodes_for[edges_forward - 1]++;
      Y_Nodes++;
		}
		// Look for Y shape branch reverse orientation
		if (edges_reverse > 1
		    && edges_forward == 1) {
			db_node_action_set_flag(node, BRANCH_NODE_REVERSE);
      Y_Nodes_rev[edges_reverse - 1]++;
      Y_Nodes++;
		}
		// Look for X-shaped branch
		if (edges_reverse > 1
		    && edges_forward > 1) {
			db_node_action_set_flag(node, X_NODE);
      X_Nodes_tot[edges_reverse + edges_forward + 1]++;
      X_Nodes++;
		}
  } // identify_branch_nodes()


  graph_queue = queue_new(METACORTEX_QUEUE_SIZE);
  if (!graph_queue) {
      log_and_screen_printf("Couldn't get memory for graph queue.\n");
      exit(-1);
  }
    /* Initialise temporaray path array buffers */
  path_array_initialise_buffers(graph->kmer_size);

	// Hash table iterator to walk nodes, looking for branches
  void explore_node(dBNode * node) {
    if(db_node_check_for_any_flag(node, PRUNED | VISITED) == false){
      dBNode* seed_node;
      nodes_in_graph->total_size = 0;
      nodes_in_graph->branch_nodes = 0;
      // Grow graph from this node, returning the 'best' (highest coverage) node to store as seed point
      log_printf("Growing graph from node\n");
      graph_queue->number_of_items = 0;

      // now with a subgraph, walk the graph looking for bubbles. check bubblefind.c for ideas.
      grow_graph_from_node_stats(node, &seed_node, graph, graph_queue, nodes_in_graph);
      if (seed_node == NULL) {
        printf("ERROR: Seed node is NULL, nodes in graph is %d\n", nodes_in_graph->total_size);
      } else if (nodes_in_graph->total_size) {
        // print out the size of the current subgraph
        log_printf("graph size\t%i\n",nodes_in_graph->total_size);
        fprintf(fp_analysis, "%i\t%i",nodes_in_graph->branch_nodes,nodes_in_graph->total_size);
        if (nodes_in_graph->branch_nodes){
          binary_kmer_to_seq(&(seed_node->kmer), graph->kmer_size, seq);
          fprintf(fp_analysis, "\t%s\n", seq);
        }
        else{
          fprintf(fp_analysis, "\n");
        }
      } else {
        // catch graph size of zero? Not sure why this happens - grow-graph must be failing
        log_printf("graph size of zero?\n");
      }
      if (nodes_in_graph->branch_nodes>(MAX_BRANCHES-1)){
        nodes_in_graph->branch_nodes=MAX_BRANCHES-1;
      }
      Contig_Branches[nodes_in_graph->branch_nodes]++;
    }
  } // explore_node


  // check each node in the graph, FLAG X&Y nodes (mark all nodes as visited)
	hash_table_traverse(&identify_branch_nodes, graph);

  // Output graph wide stats (X/Y node numbers)
  fprintf(fp_analysis, "#total\t%li\n#X-nodes\t%i\n#Y-nodes\t%i\n\n",total_nodes, X_Nodes, Y_Nodes);
  fprintf(fp_analysis, "#total\t%li\n\t\n\t#X\t#Y-FOR\t#Y-REV\n",total_nodes);
  for(i=0;i<4;i++){
    fprintf(fp_analysis, "%i\t%i\t%i\t%i\n",i, X_Nodes_tot[i], Y_Nodes_for[i], Y_Nodes_rev[i]);
  }
  for(i=4;i<8;i++){
    fprintf(fp_analysis, "%i\t%i\n",i, X_Nodes_tot[i]);
  }
  // first line for stats output file
  fprintf(fp_analysis, "\n#Subgraph sizes\n");

  // second travesal - build subgraphs out.
	hash_table_traverse(&explore_node, graph);

  // Output graph wide stats (coverage)
  fprintf(fp_analysis, "\n#Complexity_dist (# X/Y nodes)\t---\n");
  for(i=0;i<MAX_BRANCHES;i++){
    fprintf(fp_analysis, "%i\t%li\n",i, Contig_Branches[i]);
  }

  fprintf(fp_analysis, "\n#Coverage_dist\t---\n");
  for(i=0;i<(COVERAGE_BINS-1);i++){
    fprintf(fp_analysis, "#>%i<=%i\t%li\n",i*COVERAGE_BIN_SIZE, (i + 1)*COVERAGE_BIN_SIZE, Coverage_Dist[i]);
  }
  fprintf(fp_analysis, "#>=%i   \t%li\n",(COVERAGE_BINS-1)*COVERAGE_BIN_SIZE, Coverage_Dist[i]);
  fclose(fp_analysis);


  	db_graph_reset_flags(graph);
}
