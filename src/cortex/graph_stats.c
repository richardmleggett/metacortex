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

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/

// ----------------------------------------------------------------------
// Work through graph, count cov, X, Y nodes
// ----------------------------------------------------------------------
void metacortex_get_stats(dBGraph * graph, char* consensus_contigs_filename)
{
  FILE* fp_analysis;
	//int branchNodes = 0;
	int X_Nodes_tot[8];
	int Y_Nodes_rev[4];
	int Y_Nodes_for[4];
	int X_Nodes = 0;
	int Y_Nodes = 0;
  long int total_nodes = 0;
  // array to bin coverage 0-5, 5-10, 10-15..95-100
  long int Coverage_Dist[COVERAGE_BINS]; // will this work?
  //int COVERAGE_CEILING = (COVERAGE_BINS-1) * COVERAGE_BIN_SIZE;
  char analysis_filename[strlen(consensus_contigs_filename) + 10];

  Queue* graph_queue;
  int i;
  // Initialise Coverage_Dist  int i;
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

  graph_queue = queue_new(METACORTEX_QUEUE_SIZE);
  if (!graph_queue) {
      log_and_screen_printf("Couldn't get memory for graph queue.\n");
      exit(-1);
  }
  /* Initialise temporaray path array buffers */
  path_array_initialise_buffers(graph->kmer_size);

	// Hash table iterator to label nodes
	void get_node_stats(dBNode * node) {
		if (db_node_check_flag_not_pruned(node)) {
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
		}

    if (db_node_check_for_any_flag(node, PRUNED | VISITED) == false) {
      dBNode* seed_node;
      int nodes_in_graph;
      /* Grow graph from this node, returning the 'best' (highest coverage) node to store as seed point */
      log_printf("Growing graph from node\n");
      graph_queue->number_of_items = 0;
      nodes_in_graph = grow_graph_from_node(node, &seed_node, graph, graph_queue);
      if (seed_node == NULL) {
          printf("ERROR: Seed node is NULL, nodes in graph is %d\n", nodes_in_graph);
      } else {
        // something here holding nodes in graph - array seems too much? could be unwieldy
        // inflexible too. Can't pop onto array like in perl.
        log_printf("graph size\t%i\n",nodes_in_graph);
      }
    }
	}

  // check each node in the graph
	hash_table_traverse(&get_node_stats, graph);


  fprintf(fp_analysis, "#total\t%li\n#X-nodes\t%i\n#Y-nodes\t%i\n\n",total_nodes, X_Nodes, Y_Nodes);
  fprintf(fp_analysis, "#total\t%li\n\t\n\t#X\t#Y-FOR\t#Y-REV\n",total_nodes);
  for(i=0;i<4;i++){
    fprintf(fp_analysis, "%i\t%i\t%i\t%i\n",i, X_Nodes_tot[i], Y_Nodes_for[i], Y_Nodes_rev[i]);
  }
  for(i=4;i<8;i++){
    fprintf(fp_analysis, "%i\t%i\n",i, X_Nodes_tot[i]);
  }
  fprintf(fp_analysis, "#Coverage_dist\t---\n");
  for(i=0;i<(COVERAGE_BINS-1);i++){
    fprintf(fp_analysis, "#>%i<=%i\t%li\n",i*COVERAGE_BIN_SIZE, (i + 1)*COVERAGE_BIN_SIZE, Coverage_Dist[i]);
  }
  fprintf(fp_analysis, "#>=%i   \t%li\n",i*(COVERAGE_BIN_SIZE-1), Coverage_Dist[i]);
  fclose(fp_analysis);
}
