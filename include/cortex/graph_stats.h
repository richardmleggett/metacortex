/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
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
 */

#include <binary_kmer.h>

#ifndef GRAPHINFO_H
#define GRAPHINFO_H
typedef struct {
    int total_size;
    int largest_subgraph;
    int subgraph_dist[10];
    int num_subgraphs;
    int num_subgraphs_2k;
    int branch_nodes;
    int branch_nodes_total;
    int end_nodes;
    int node_degree[5][5];
    int best_coverage[5];
    int simple_bubbles;
    BinaryKmer kmer[5];
    BinaryKmer temp_kmer;
    BinaryKmer current_kmer;
    BinaryKmer highest_cov_in_subgraph;
    int highest_cov;
} GraphInfo;
#endif /* GRAPHINFO_H */


void find_subgraph_stats(dBGraph* graph, char* consensus_contigs_filename, int min_subgraph_kmers, int max_node_edges, float delta_coverage);

void print_degree_stats(GraphInfo * nodes_in_graph, FILE* fp_degrees);
