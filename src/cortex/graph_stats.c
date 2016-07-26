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
#include <math.h>
#include <unistd.h>
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

#define COVERAGE_BINS 10
#define COVERAGE_BIN_SIZE 1
#define MAX_BRANCHES 5
#define GRAPH_LOG10_LIMIT 10 // little hacky to do this here, because it needs to match size of subgraph_dist in graph_stats.h


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
     int best_edges[5];
     int i;
     *best_node = 0;
     for(i=0; i<5; i++){
       best_edges[i]=0;
     }

     // Nucleotide iterator, used to walk all possible paths from a node
     void walk_if_exists(Nucleotide n) {
         //if (debug) printf("Trying nucleotide %i\n", n);
   			int end_orientation;

         // If there is an edge in any colour for this nucleotide...
         if (db_node_edge_exist_any_colour(node, n, orientation)) {


             log_and_screen_printf("\nNEW NODE\n");  // DEBUG BUBBLE BUG

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
                  log_and_screen_printf("\n\tUNVISITED, WALKING\n");  // DEBUG BUBBLE BUG
                 pathStep first_step;
                 Path * new_path;
                 dBNode* end_node;
                 i = 0;

                 // Get path
                 first_step.node = node;
                 first_step.orientation = orientation;
                 first_step.label = n;
                 new_path = path_new(MAX_EXPLORE_NODES, graph->kmer_size);
                 if (!new_path) {
                     log_and_screen_printf("ERROR: Not enough memory to allocate new path.\n");
                     exit(-1);
                 }


                log_and_screen_printf("\n\tGETTING PERFECT PATH...\n");  // DEBUG BUBBLE BUG
                 db_graph_get_perfect_path_with_first_edge_all_colours(&first_step, &db_node_action_do_nothing, new_path, graph);
                log_and_screen_printf("\n\t\t...GOT PERFECT PATH\n");  // DEBUG BUBBLE BUG

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
                // things I'm dropping for now - loop detection, branch + loop, catching too large a bubble
                if (db_node_is_blunt_end_all_colours(end_node, end_orientation)) {
                  // DO NOTHING WITH THIS
                  //db_graph_check_and_add_path(merged_path, patharray);
                }
                if (db_node_check_flag_visited(end_node)) {
                   // need to count back from here to original branching point?
                   log_and_screen_printf("\nBUBBLE FOUND, path length\t%i\n", new_path->length);
                   // length of path here? not perfect - if bubble structure is complex, it will only report on the most immediate perfect path size.
                }

               // Now go through all nodes, look for best and mark all as visited

              log_and_screen_printf("\n\tCHECKING PERFECT PATH length %d\n", new_path->length);  // DEBUG BUBBLE BUG
               for (i=0; i<new_path->length; i++) {
                   log_and_screen_printf("\t\tnode %d\n", i);  // DEBUG BUBBLE BUG
                   if (!db_node_check_flag_visited(new_path->nodes[i])) {
                       int this_coverage = element_get_coverage_all_colours(new_path->nodes[i]);
                       int this_FOR_edges = db_node_edges_count_all_colours(new_path->nodes[i], forward);
                       int this_REV_edges = db_node_edges_count_all_colours(new_path->nodes[i], reverse);

                       // add node degrees to 2D array of all degrees in subgraph
                       nodes_in_graph->node_degree[this_FOR_edges][this_REV_edges]++;

                       // if this is the new best node update the other bests
                       if ((best_node[0] == 0) ||
                           (this_coverage > nodes_in_graph->best_coverage[0]) ||
                           ((this_coverage == nodes_in_graph->best_coverage[0]) && ((this_FOR_edges + this_REV_edges) < best_edges[0])))
                       {
                           best_edges[0] = (this_FOR_edges + this_REV_edges);
                           *best_node = new_path->nodes[i];
                       }

                       if (this_coverage>nodes_in_graph->highest_cov){
                          nodes_in_graph->highest_cov=this_coverage;
                          binary_kmer_assignment_operator(nodes_in_graph->current_kmer,new_path->nodes[i]->kmer);
                          binary_kmer_assignment_operator(nodes_in_graph->highest_cov_in_subgraph,nodes_in_graph->current_kmer);
                       }

                       // if this is better than the lowest 'good' node (top five coverage)
                       // NOTE: should top five be an option?
                       if ((best_node == 0) ||
                           (this_coverage > nodes_in_graph->best_coverage[4]) ||
                           ((this_coverage == nodes_in_graph->best_coverage[4]) && ((this_FOR_edges + this_REV_edges) > best_edges[4])))
                       {
                          // sort algorithm - because I sort as I build array, no need to make more than one pass
                          int temp_cov=0;
                          binary_kmer_initialise_to_zero(&(nodes_in_graph->temp_kmer));
                          // yes, this is the same as above.
                          binary_kmer_assignment_operator(nodes_in_graph->current_kmer,new_path->nodes[i]->kmer);
                          // seed_node->kmer

                          int j=0;
                          while(this_coverage){
                            log_and_screen_printf("\n\t\t\tj %d\t%d\t%d\n", j, this_coverage, nodes_in_graph->best_coverage[j]);  // DEBUG BUBBLE BUG
                            if(j>4){
                              this_coverage=0;
                            }
                            else if (this_coverage>nodes_in_graph->best_coverage[j]||
                            ((this_coverage == nodes_in_graph->best_coverage[j]) && ((this_FOR_edges + this_REV_edges) > best_edges[j]))){
                              temp_cov = nodes_in_graph->best_coverage[j];
                              nodes_in_graph->best_coverage[j] = this_coverage;
                              this_coverage=temp_cov;

                              // recycle temp_cov for one line
                              temp_cov=best_edges[j];
                              best_edges[j] = (this_FOR_edges + this_REV_edges);
                              // set the current edge count for the rest of the loop
                              this_FOR_edges=temp_cov;
                              this_REV_edges=0;
                              temp_cov=0;

                              binary_kmer_assignment_operator(nodes_in_graph->temp_kmer,nodes_in_graph->kmer[j]);
                              binary_kmer_assignment_operator(nodes_in_graph->kmer[j],nodes_in_graph->current_kmer);
                              binary_kmer_assignment_operator(nodes_in_graph->current_kmer,nodes_in_graph->temp_kmer);
                              j++;
                            }
                            else{
                              j++;
                            }
                          }
                      }

                      if (db_node_check_for_any_flag(new_path->nodes[i], BRANCH_NODE_FORWARD)){
                        nodes_in_graph->branch_nodes++;
                      }
                      else if (db_node_check_for_any_flag(new_path->nodes[i], BRANCH_NODE_REVERSE)){
                        nodes_in_graph->branch_nodes++;
                      }
                      else if (db_node_check_for_any_flag(new_path->nodes[i], X_NODE)){
                        nodes_in_graph->branch_nodes++;
                      }

                       db_node_action_set_flag_visited(new_path->nodes[i]);
                       queue_push(graph_queue, new_path->nodes[i]);
                       nodes_in_graph->total_size++;
                   }
               } // path->length loop
               // Clean up
               path_destroy(new_path);
             }
         }
     } // walk_if_exists

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
void find_subgraph_stats(dBGraph * graph, char* consensus_contigs_filename, int min_subgraph_kmers)
{
  FILE* fp_analysis;
  FILE* fp_analysis_DIGEST;
  FILE* fp_degrees;
  long int Contig_Branches[MAX_BRANCHES];
  char* seq = calloc(256, 1);
  long int total_nodes = 0;
  int i;  int j; float percentage;


  GraphInfo* nodes_in_graph = calloc(1,sizeof(GraphInfo));
  // need a small function for initialising this?
  nodes_in_graph->largest_subgraph = 0;
  nodes_in_graph->num_subgraphs = 0;
  nodes_in_graph->num_subgraphs_2k = 0;
  for(i=0;i<GRAPH_LOG10_LIMIT;i++){
    nodes_in_graph->subgraph_dist[i]=0;
  }
  for (i=0; i<5; i++) {
    nodes_in_graph->best_coverage[i]=0;
  }
  nodes_in_graph->branch_nodes_total=0;

  // array to bin coverage 0-5, 5-10, 10-15..95-100
  long int Coverage_Dist[COVERAGE_BINS*COVERAGE_BIN_SIZE]; // will this work?
  char analysis_filename[strlen(consensus_contigs_filename) + 10];
  char degrees_filename[strlen(consensus_contigs_filename) + 10];

  Queue* graph_queue;
  for(i=0;i<MAX_BRANCHES;i++){
    Contig_Branches[i]=0;
  }
  // Initialise Coverage_Dist  int i;
  for(i=0;i<(COVERAGE_BINS*COVERAGE_BIN_SIZE);i++){
    Coverage_Dist[i]=0;
  }

  /* Open the analysis file */
  sprintf(analysis_filename, "%s.analysis", consensus_contigs_filename);
  fp_analysis = fopen(analysis_filename, "w");
  if (!fp_analysis) {
      log_and_screen_printf("ERROR: Can't open analysis file.\n");
      exit(-1);
  }

  /* Open the analysis file */
  sprintf(analysis_filename, "%s.digest", consensus_contigs_filename);
  fp_analysis_DIGEST = fopen(analysis_filename, "w");
  if (!fp_analysis_DIGEST) {
      log_and_screen_printf("ERROR: Can't open analysis (DIGEST) file.\n");
      exit(-1);
  }


  /* Open the sugraph degree file */
  sprintf(degrees_filename, "%s.degrees", consensus_contigs_filename);
  fp_degrees = fopen(degrees_filename, "w");
  if (!fp_degrees) {
      log_and_screen_printf("ERROR: Can't open degrees file.\n");
      exit(-1);
  }

  // header line for degrees file
  for(i=0;i<5;i++){
    for(j=0;j<5;j++){
      fprintf(fp_degrees,"for[%d]rev[%d]\t", i, j);
    }
  }
  fprintf(fp_degrees,"total\n");

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
    this_coverage = (this_coverage-1);

    if(this_coverage>COVERAGE_BINS*COVERAGE_BIN_SIZE-1){
      this_coverage = COVERAGE_BINS*COVERAGE_BIN_SIZE-1;
    }
      //log_and_screen_printf("WTFCOVDIST\t%i\n", this_coverage);

    Coverage_Dist[this_coverage]++;
    total_nodes++;

		// Look for Y shape branch forward orientation
		// The nodes at the top of the Y should contain different colours
		if (edges_forward > 1
		    && edges_reverse == 1) {
			db_node_action_set_flag(node, BRANCH_NODE_FORWARD);
		}
		// Look for Y shape branch reverse orientation
		if (edges_reverse > 1
		    && edges_forward == 1) {
			db_node_action_set_flag(node, BRANCH_NODE_REVERSE);
		}
		// Look for X-shaped branch
		if (edges_reverse > 1
		    && edges_forward > 1) {
			db_node_action_set_flag(node, X_NODE);
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
      nodes_in_graph->end_nodes = 0;
      nodes_in_graph->highest_cov = 0;
      for(i=0;i<5;i++){
        for(j=0;j<5;j++){
          nodes_in_graph->node_degree[i][j]=0;
        }
      }

      // Grow graph from this node, returning the 'best' (highest coverage) node to store as seed point
      log_printf("Growing graph from node\n");
      graph_queue->number_of_items = 0;

      // now with a subgraph, walk the graph looking counting degrees by graph and overal
      grow_graph_from_node_stats(node, &seed_node, graph, graph_queue, nodes_in_graph);
      if (seed_node == NULL) {
        printf("ERROR: Seed node is NULL, nodes in graph is %d\n", nodes_in_graph->total_size);
      } else if (nodes_in_graph->total_size) {
        // print out the size of the current subgraph
        log_printf("graph size\t%i\n",nodes_in_graph->total_size);
        fprintf(fp_analysis, "%i\t%i\t",nodes_in_graph->branch_nodes,nodes_in_graph->total_size);
        binary_kmer_to_seq(&nodes_in_graph->highest_cov_in_subgraph, graph->kmer_size, seq);
        fprintf(fp_analysis, "%s\n", seq);

        // update graph wide stats
        print_degree_stats(nodes_in_graph, fp_degrees);
        if(nodes_in_graph->total_size>nodes_in_graph->largest_subgraph){
          nodes_in_graph->largest_subgraph=nodes_in_graph->total_size;
        }
        nodes_in_graph->branch_nodes_total=nodes_in_graph->branch_nodes_total+nodes_in_graph->branch_nodes;
        nodes_in_graph->num_subgraphs++;
        i=log10(nodes_in_graph->total_size);
        if(i>=GRAPH_LOG10_LIMIT){
          i=GRAPH_LOG10_LIMIT-1;
        }
        nodes_in_graph->subgraph_dist[i]++;
        if(nodes_in_graph->total_size>2000){
          nodes_in_graph->num_subgraphs_2k++;
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

  // first line for stats output file
  fprintf(fp_analysis, "\n#Subgraph sizes\n");

  // second travesal - build subgraphs out.
	hash_table_traverse(&explore_node, graph);
  fclose(fp_analysis);

  // Output graph wide stats (coverage)
  fprintf(fp_analysis_DIGEST, "\n#Complexity_dist of total graph (# X/Y nodes)\t---\n");
  for(i=0;i<MAX_BRANCHES;i++){
    fprintf(fp_analysis_DIGEST, "%i\t%li\n",i, Contig_Branches[i]);
  }

  fprintf(fp_analysis_DIGEST, "\n#Coverage_dist\t---\n");

  // first two lines are for 1, 2-4 cov. after that stick revert to cov bin size
  fprintf(fp_analysis_DIGEST, "#1\t%li\n", Coverage_Dist[0]);
  fprintf(fp_analysis_DIGEST, "#2-4\t%i\n", sum_array(Coverage_Dist,1,3));
  if(COVERAGE_BIN_SIZE>4){
    fprintf(fp_analysis_DIGEST, "#>4<=%i\t%i\n", COVERAGE_BIN_SIZE, sum_array(Coverage_Dist,4,COVERAGE_BIN_SIZE-1));
  }
  fprintf(fp_analysis_DIGEST, "\n");

 // now repeat the coverage output, but for every bin
  for(i=0;i<(COVERAGE_BINS*COVERAGE_BIN_SIZE-1);i+=COVERAGE_BIN_SIZE){
    if(COVERAGE_BIN_SIZE>1){
      fprintf(fp_analysis_DIGEST, "#>%i<=%i\t%i\n",(i)*COVERAGE_BIN_SIZE, (i + 1)*COVERAGE_BIN_SIZE, sum_array(Coverage_Dist, i*COVERAGE_BIN_SIZE, (i + 1)*COVERAGE_BIN_SIZE));
    }
    else{
      fprintf(fp_analysis_DIGEST, "#%i\t%li\n",i+1, Coverage_Dist[i]);
    }
  }
  fprintf(fp_analysis_DIGEST, "#>=%i   \t%li\n",(COVERAGE_BINS)*COVERAGE_BIN_SIZE, Coverage_Dist[(COVERAGE_BINS*COVERAGE_BIN_SIZE)-1]);

  // kmer figures
  fprintf(fp_analysis_DIGEST, "\n#kmers\nunique\t%lli\ttotal\t%lli\t",graph->unique_kmers,graph->loaded_kmers);
  percentage=(float)(graph->unique_kmers)/(float)(graph->loaded_kmers);
  fprintf(fp_analysis_DIGEST, "%%_of_total\t%.2f\n", percentage*100);
  percentage=(float)(nodes_in_graph->largest_subgraph)/(float)(graph->unique_kmers);
  fprintf(fp_analysis_DIGEST, "\n#subgraphs\nlargest_subgraph\t%i\t%%_of_total\t%.2f\n",nodes_in_graph->largest_subgraph, percentage*100);
  fprintf(fp_analysis_DIGEST, "num_subgraphs\t%i\tnum_subgraphs>2k\t%i\n", nodes_in_graph->num_subgraphs, nodes_in_graph->num_subgraphs_2k);
  fprintf(fp_analysis_DIGEST, "num_subgraphs_per_E^6kmers\t%f\n", (float)(nodes_in_graph->num_subgraphs)/(float)(graph->unique_kmers));
  fprintf(fp_analysis_DIGEST, "branches\t%i\tper_1000_nodes\t%.2f\n\n#graph size dist:\n", nodes_in_graph->branch_nodes_total, (float)(nodes_in_graph->branch_nodes_total)/1000.0);
  for(i=0;i<GRAPH_LOG10_LIMIT;i++){
    fprintf(fp_analysis_DIGEST, "<=E^%i\t%i\n",i,nodes_in_graph->subgraph_dist[i]);
  }

  fprintf(fp_analysis_DIGEST, "\n#highest coverage kmers\n");
  for(i=0;i<5;i++){
    // seq never re-initialised
    binary_kmer_to_seq(&nodes_in_graph->kmer[i], graph->kmer_size, seq);
    fprintf(fp_analysis_DIGEST, "%d\t%s\n", nodes_in_graph->best_coverage[i], seq);
  }
  fprintf(fp_analysis_DIGEST, "\n");

  fclose(fp_analysis_DIGEST);


  // run R script to produce figures for report
  // will this work? initialising 'cmd' like this?
  printf("\nPATH : %s\n", getenv("PATH"));
    //char cmd = printf("Rscript %s %s", <path_to_src>/degree_plots.R, degrees_filename);
    //system(cmd);  // potential problems with this apparently? is permissions are an initialiseAlignmentSummaryFile

  char command[1024];
  char cwd[1024];
  char r_script_path='/home/aylingm/grimoire/metacortex/';
   if (getcwd(cwd, sizeof(cwd)) != NULL){
     sprintf(command, "Rscript %sdegree_plots.R %s/%s", r_script_path, cwd, degrees_filename);
     //sprintf(command, "Rscript /home/aylingm/testgenome/degree_plots.R /home/aylingm/testgenome/patronol/contigs_46_100k.1_clean.fa.degrees");
     //sprintf(command, "Rscript /home/aylingm/grimoire/metacortex/src/degree_plots.R %s/%s", cwd, degrees_filename);
     log_and_screen_printf("\n%s\n", command);
     system(command);
     log_and_screen_printf("\n");
   }
   else{
     log_and_screen_printf("CWD command reutrned NULL\n");
   }

  // exec("Rscript <path_to_src>/degree_plots.R degrees_filename")

  db_graph_reset_flags(graph);
}

void print_degree_stats(GraphInfo * nodes_in_graph, FILE* fp_degrees){
  int i;  int j;
  int total_nodes=nodes_in_graph->total_size;

  for(i=0;i<5;i++){
    for(j=0;j<5;j++){
      fprintf(fp_degrees,"%f\t",  ((float) nodes_in_graph->node_degree[i][j]) / (float) total_nodes);
    }
  }
  fprintf(fp_degrees,"%d\n", total_nodes);
}

int sum_array(long int * array, int first, int last){
  int sum = 0;
  int i;

  for(i=first; i<=last; i++){
    sum += array[i];
  }
  return sum;
}
