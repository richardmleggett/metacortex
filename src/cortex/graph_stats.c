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

#define COVERAGE_BINS 20
#define COVERAGE_BIN_SIZE 5

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/

// ----------------------------------------------------------------------
// Work through graph, count cov, X, Y nodes
// ----------------------------------------------------------------------
void metacortex_get_stats(dBGraph * db_graph, char* consensus_contigs_filename)
{
  FILE* fp_analysis;
	//int branchNodes = 0;
	int X_Nodes = 0;
	int Y_Nodes = 0;
  long int total_nodes = 0;
  // array to bin coverage 0-5, 5-10, 10-15..95-100
  long int Coverage_Dist[COVERAGE_BINS]; // will this work?
  int COVERAGE_CEILING = (COVERAGE_BINS-1) * COVERAGE_BIN_SIZE;
  char analysis_filename[strlen(consensus_contigs_filename) + 10];

  // Initialise Coverage_Dist
  int i;
  for(i=0;i<COVERAGE_BINS;i++){
    Coverage_Dist[i]=0;
  }

  /* Open the analysis file */
  sprintf(analysis_filename, "%s.analysis", consensus_contigs_filename);
  fp_analysis = fopen(analysis_filename, "w");
  if (!fp_analysis) {
      log_and_screen_printf("ERROR: Can't open analysis file.\n");
      exit(-1);
  }

	db_graph_reset_flags(db_graph);

	// Hash table iterator to label nodes
	void find_end_nodes(dBNode * node) {
		if (db_node_check_flag_not_pruned(node)) {
      int this_coverage = element_get_coverage_all_colours(node);
      // assert (this_coverage>0) ?
      if(this_coverage>COVERAGE_CEILING){
        this_coverage = COVERAGE_BINS-1;
      }
      else if(this_coverage>0){
        this_coverage = ((this_coverage-1) / COVERAGE_BIN_SIZE)+1;
      }
      else {
        // Coverage is zero?
      }
      Coverage_Dist[this_coverage]++;
      total_nodes++;

			// Look for Y shape branch forward orientation
			// The nodes at the top of the Y should contain different colours
			if (db_node_edges_count_all_colours(node, forward) > 1
			    && db_node_edges_count_all_colours(node, reverse) == 1) {
				db_node_action_set_flag(node, BRANCH_NODE_FORWARD);
        Y_Nodes++;
			}
			// Look for Y shape branch reverse orientation
			if (db_node_edges_count_all_colours(node, reverse) > 1
			    && db_node_edges_count_all_colours(node, forward) == 1) {
				db_node_action_set_flag(node, BRANCH_NODE_REVERSE);
        Y_Nodes++;
			}
			// Look for X-shaped branch
			if (db_node_edges_count_all_colours(node, reverse) > 1
			    && db_node_edges_count_all_colours(node, forward) > 1) {
				db_node_action_set_flag(node, X_NODE);
        X_Nodes++;
			}
		}
	}

  // check each node in the graph
	hash_table_traverse(&find_end_nodes, db_graph);

  fprintf(fp_analysis, "#total\t%li\n#X-nodes\t%i\n#Y-nodes\t%i\n#Coverage_dist\t---\n",total_nodes, X_Nodes, Y_Nodes);
  for(i=0;i<COVERAGE_BINS;i++){
    fprintf(fp_analysis, "#%i-%i\t%li\n",i*COVERAGE_BIN_SIZE, (i+1)*COVERAGE_BIN_SIZE, Coverage_Dist[i]);
  }
  fclose(fp_analysis);
}
