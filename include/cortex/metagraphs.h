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
 * metagraphs.h
 ************************************************************************/

#define MAX_SEEDS 200000000
#define METACORTEX_QUEUE_SIZE 20000000 // 10000000
#define MAX_EXPLORE_PATH_LENGTH 2000000
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

void metacortex_find_subgraphs(dBGraph* graph, char* consensus_contigs_filename, int min_subgraph_kmers, int min_contig_length, boolean multiple_subgraph_contigs);
