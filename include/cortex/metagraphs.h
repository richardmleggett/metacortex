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

#define MAX_SEEDS 200000000
#define METACORTEX_QUEUE_SIZE 10000000
#define MAX_EXPLORE_PATH_LENGTH 200000
#define MAX_EXPLORE_NODES 200

typedef struct {
     dBNode* seed_node;
     int graph_size;
} SubGraphInfo;

/*typedef struct {
	int max_size;
	int number_of_items;
	int head;
	int tail;
	void** items;
} Queue;
*/
//int grow_graph_from_node(dBNode* start_node, dBNode** best_node, dBGraph* graph, Queue* graph_queue);

void metacortex_find_subgraphs(dBGraph* graph, char* consensus_contigs_filename, int min_subgraph_kmers, int min_contig_length);