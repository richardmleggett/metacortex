/************************************************************************
 *
 * This file is part of MetaCortex
 *
 * Authors:
 *     Richard M. Leggett (richard.leggett@earlham.ac.uk) and
 *     Martin Ayling (martin.ayling@earlham.ac.uk)
 *
 * MetaCortex is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MetaCortex is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MetaCortex.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************/

/************************************************************************
 * report_output.h
 ************************************************************************/

#include "dB_graph.h"
#include "graph_stats.h"

void writeLaTeXHeader(FILE* fp_latex, char* consensus_contigs_filename);

/*void writeLaTeXreport(FILE* fp_latex, int MAX_BRANCHES, int COVERAGE_BIN_SIZE,
 int GRAPH_LOG10_LIMIT, int NUM_BEST_NODES, long int Contig_Branches,
 long int Coverage_Dist, dBGraph * graph, GraphInfo * nodes_in_graph); */

void writeLaTeXreport(FILE* fp_latex, int MAX_BRANCHES, int COVERAGE_BIN_SIZE,
                      int COVERAGE_BINS, int GRAPH_LOG10_LIMIT, int NUM_BEST_NODES, long int * Contig_Branches,
                      long int * Coverage_Dist, dBGraph * graph, GraphInfo * nodes_in_graph);

void writeLaTeXreport_to_log_and_screen(int MAX_BRANCHES, int COVERAGE_BIN_SIZE,
                                        int COVERAGE_BINS, int GRAPH_LOG10_LIMIT, int NUM_BEST_NODES, long int * Contig_Branches,
                                        long int * Coverage_Dist, dBGraph * graph, GraphInfo * nodes_in_graph);

int sum_array(long int * array, int first, int last);
