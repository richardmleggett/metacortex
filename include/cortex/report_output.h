/*
 *
 * CORTEX project contacts:
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * Development team:
 *       R. Ramirez-Gonzalez (Ricardo.Ramirez-Gonzalez@bbsrc.ac.uk)
 *       R. Leggett (richard@leggettnet.org.uk)
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 * This file is to declare the analysis functions. Stuff that is not necesarily part of the main assembly, but that can be useful to analyze the graph.
 */


 #include "dB_graph.h"
 #include "graph_stats.h"

void writeLaTeXHeader(FILE* fp_latex, char* consensus_contigs_filename);

/*void writeLaTeXreport(FILE* fp_latex, int MAX_BRANCHES, int COVERAGE_BIN_SIZE,
 int GRAPH_LOG10_LIMIT, int NUM_BEST_NODES, long int Contig_Branches,
 long int Coverage_Dist, dBGraph * graph, GraphInfo * nodes_in_graph); */

void writeLaTeXreport(FILE* fp_latex, int MAX_BRANCHES, int COVERAGE_BIN_SIZE, \
 int COVERAGE_BINS, int GRAPH_LOG10_LIMIT, int NUM_BEST_NODES, long int * Contig_Branches, \
 long int * Coverage_Dist, dBGraph * graph, GraphInfo * nodes_in_graph);

int sum_array(long int * array, int first, int last);
