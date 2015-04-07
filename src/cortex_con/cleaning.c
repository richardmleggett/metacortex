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
 
#include <stdio.h>
#include <stdint.h>
#ifdef THREADS
#include <pthread.h>
#endif
#include <binary_kmer.h>
#include <flags.h>
#include <element.h>
#include <dB_graph.h>
#include <path.h>
#include <assert.h>
#include <perfect_path.h>
#include <cleaning.h>
#include <logger.h>

void cleaning_prune_db_node(dBNode * node, dBGraph * db_graph)
{

	pathStep ps, next_ps, rev_ps;
	ps.node = node;
    if(db_node_check_flag_not_pruned(node)){
        db_graph->pruned_kmers++;
        db_node_action_set_flag_pruned(node);
    }
    Nucleotide n;
    int i;
    
    for(i=0; i <2; i++){
        ps.orientation = i == 0 ? forward:reverse;
        for(n=0; n < 4; n++){
       		ps.label = n;
       		if (db_node_edge_exist_any_colour(ps.node, ps.label, ps.orientation)) {
       			db_graph_get_next_step(&ps, &next_ps, &rev_ps,
       					       db_graph);
       			if (rev_ps.node != NULL) {
       				db_node_reset_edge_all_colours(rev_ps.node,
       							       rev_ps.orientation,
       							       rev_ps.label);
       			}
       		}
       	}
    }
    
	db_node_reset_edges_all_colours(node);
	db_graph->calculated = false;
}


int cleaning_prune_low_coverage_path(int min_cov, int limit, dBGraph * db_graph)
{
    log_and_screen_printf("Removing paths with coverage %i or lower...\n", min_cov);
    
	db_graph_reset_flags(db_graph);
	
	int removed_nodes = 0;
	int count_paths = 0;
	int lone_reads = 0;
	
	PathArray *pa = path_array_get_from_buffer_with_size(4);
	boolean prune = true;
	
	void prune_step_node(pathStep * step) {
		removed_nodes++;
		cleaning_prune_db_node(step->node, db_graph);
	}
	void has_low_coverage_node(pathStep * ps) {
		
		if (element_get_coverage_all_colours(ps->node) > min_cov) {
			prune &= false;

		}
	}
		
	void check_and_remove_low_cov_path(Path * p){
		
		    prune = true;
			if (p->length > 0) {
				
				path_inner_iterator(&has_low_coverage_node,
						    p);
				if (prune) {
					pathStep first, last;
					path_get_last_step(&last, p);
					path_get_step_at_index(0, &first, p);
					if(element_get_coverage_all_colours(last.node) > min_cov || 
					   element_get_coverage_all_colours(first.node) > min_cov){
						
						path_inner_iterator(&prune_step_node, p);
						count_paths++;	
						
					}else{
						path_iterator(&prune_step_node,
									  p);
						lone_reads++;
						count_paths++;			
						
					}
					
				}
			}
		
	}
	
	void remove_low_cov_path(dBNode * node) {
		
		if(!db_node_check_for_any_flag(node, PRUNED  | VISITED)){
		
			
			 if(element_get_coverage_all_colours(node) <= min_cov){
				perfect_path_get_path_with_callback(node, undefined, &db_node_action_set_flag_visited, &check_and_remove_low_cov_path, db_graph);
			}
		}


	}

	hash_table_traverse(&remove_low_cov_path, db_graph);

	path_array_free_from_buffer(pa);

	return removed_nodes;
}

int cleaning_remove_bubbles(int limit, dBGraph * db_graph)
{

	db_graph_reset_flags(db_graph);

	int removed_nodes = 0;
	int count_bubbles = 0;

	
	

	void prune_step_node(pathStep * step) {
		removed_nodes++;
		cleaning_prune_db_node(step->node, db_graph);
	}
	
	void visit_step_node(pathStep * step) {
		db_node_action_set_flag_visited(step->node);
	}

	void remove_bubble_with_orientation(dBNode * node,
					    Orientation orientation) {
		int i, j;
		pathStep psi, psj, psi1, psj1;
		PathArray *pa = path_array_get_from_buffer_with_size(4);
		perfect_path_get_all_paths_from(node, orientation, pa, limit, db_graph);
		
		for (i = 0; i < 4; i++) {
			path_iterator(&visit_step_node, pa->paths[i]);
		}
		for (i = 0; i < 3; i++) {

			if (pa->paths[i]->length > 2) {
							
				path_get_last_step(&psi, pa->paths[i]);
				for (j = i + 1; j < 4; j++) {
                    
                    
                    
                    //TODO: change to ">2"
					if (pa->paths[j]->length > 2/* &&
					     (pa->paths[j]->length == pa->paths[i]->length)*/){
                    
                            
                        path_get_last_step(&psj, pa->paths[j]);

						if (path_step_equals_without_label(&psi, &psj)) {
							path_get_step_at_index(1, &psi1, pa->paths[i]);
							path_get_step_at_index(1, &psj1, pa->paths[j]);
							count_bubbles++;
							if (element_get_coverage_all_colours(psj1.node) < element_get_coverage_all_colours(psi1.node)) {
								path_inner_iterator
								    (&prune_step_node,
								     pa->
								     paths[j]);
							} else {
								path_inner_iterator
								    (&prune_step_node,
								     pa->
								     paths[i]);
							}

						}
					}

				}
			}

		}
		path_array_free_from_buffer(pa);
	}

	void remove_bubble(dBNode * node) {
		if(!db_node_check_for_any_flag(node, PRUNED | VISITED)){
			Orientation o = forward;
			if (db_node_edges_count(node, o) > 1) {
				remove_bubble_with_orientation(node, o);
			}
			o = reverse;
			if (db_node_edges_count(node, o) > 1) {
				remove_bubble_with_orientation(node, o);
			}
		}
	}

	hash_table_traverse(&remove_bubble, db_graph);
	log_and_screen_printf("%'i bubbles removed\n", count_bubbles);

	return removed_nodes;
}

int cleaning_remove_low_coverage(int coverage, dBGraph * db_graph)
{
	int pruned = 0;
	void remove_low_cov(dBNode * node) {
		if (element_get_coverage_all_colours(node) <= coverage) {
			cleaning_prune_db_node(node, db_graph);
			pruned++;
		}
	}
	hash_table_traverse(&remove_low_cov, db_graph);
	log_and_screen_printf("Removed %'i low coverage nodes\n", pruned);
	return pruned;
}

// remove any node under x coverage (and the arrows)
// Node coverage is viewed as sum of ALL colours
// Signature kept to be backward compatible
boolean db_graph_db_node_prune_low_coverage(dBNode * node, int coverage,
		                        		    void (*node_action) (dBNode * node),
				                            dBGraph * db_graph)
{
	return cleaning_remove_low_coverage(coverage, db_graph) > 0;
}

// Remove spurious edges, eg.
// 
//   ----> C ----> B ---->
//                 ^
//          _______/ 
//         /
// E ----> A ----> D ---->
//
// Edge between A & B can be removed if:
//     coverage(B) = 1+coverage(C)
// and coverage(A) = 1+coverage(D)
// and connection exists E->A
int cleaning_remove_spurious_links(int max_difference, int min_coverage, dBGraph * db_graph)
{    
    int links_removed = 0;
    
    void check_and_remove_edge(dBNode * node_A) {
        int o;
        
        // Check we have minimum coverage at node A
        if (element_get_coverage_all_colours(node_A) < min_coverage) {
            return;
        }
                
        // Check both orientations
        for (o=0; o<2; o++) {
            int orientation = o == 0 ? forward:reverse;
            int n_edges_fwd = db_node_edges_count_all_colours(node_A, orientation);
            int n_edges_rev = db_node_edges_count_all_colours(node_A, opposite_orientation(orientation));
            
            // Must be more than 1 edge forward and at least one edge reverse
            if ((n_edges_fwd > 1) && (n_edges_rev > 0)) {            
                // We need to do a pairwise comparison of edges, so we have n1 and n2 to represent
                // pairs of labels. Each pair ends up evaluated twice - once assuming n1 leads to
                // D, the other assuming n1 leads to B.
                // RHRG: Changed the code to use the nucleotide iterator, to make t his compatible with the SOLiD code
                
                // Nucleotide n1, n2;
                //for (n1=Adenine; n1<=Thymine; n1++) {
                void outer_loop(Nucleotide n1){
                    if (db_node_edge_exist_any_colour(node_A, n1, orientation)) {
                        //for (n2=Adenine; n2<=Thymine; n2++) {
                        void inner_loop(Nucleotide n2){
                            if ((n2 != n1) && (db_node_edge_exist_any_colour(node_A, n2, orientation))) {
                                pathStep current_step, next_step, reverse_step, reverse_step_B;
                                dBNode *node_B, *node_C, *node_D;
                                int coverage_difference;
                                
                                // We get passed node A, so find node D            
                                current_step.node = node_A;
                                current_step.label = n1;
                                current_step.orientation = orientation;
                                db_graph_get_next_step(&current_step, &next_step, &reverse_step, db_graph);
                                node_D = next_step.node;
                                if (!node_D) {
                                    fprintf(stderr, "[cleaning_remove_spurious_links] Something went wrong, couldn't find node D\n");
                                    exit(1);
                                }
                                
                                // If this is what we want it to be, coverage(A) - coverage(D) should be > 0 and <= max_difference
                                coverage_difference = element_get_coverage_all_colours(node_A) - element_get_coverage_all_colours(node_D);
                                if ((coverage_difference > 0) && (coverage_difference <= max_difference)) {                                    
                                    // Now find node B
                                    current_step.node = node_A;
                                    current_step.orientation = orientation;
                                    current_step.label = n2;
                                    db_graph_get_next_step(&current_step, &next_step, &reverse_step_B, db_graph);
                                    node_B = next_step.node;
                                    if (!node_B) {
                                        fprintf(stderr, "[cleaning_remove_spurious_links] Something went wrong, couldn't find node B\n");
                                        exit(1);
                                    }
                                    
                                    // Node B coverage should be above the minimum too
                                    if (element_get_coverage_all_colours(node_B) >= min_coverage) {
                                        dBNode* neighbours[4];
                                        Nucleotide labels[4];
                                        int num_edges;
                                        int i;
                                        
                                        for (i=0; i<4; i++) {
                                            neighbours[i] = 0;
                                            labels[i] = Undefined;
                                        }                                        
                                        
                                        // From B, find node C. For now, we'll only take the simple case where there
                                        // is only 1 reverse edge from B (2 including A). TODO: Sum coverage of reverse edges???
                                        num_edges = db_graph_get_neighbouring_nodes_all_colours(next_step.node, opposite_orientation(next_step.orientation), neighbours, labels, db_graph);
                                        if (num_edges == 2) {
                                            // Find the edge that doesn't lead back to A...
                                            if ((neighbours[0] == node_A) && (neighbours[1] != node_A)) {
                                                node_C = neighbours[1];
                                            } else if ((neighbours[0] != node_A) && (neighbours[1] == node_A)) {
                                                node_C = neighbours[0];
                                            } else {
                                                char label[db_graph->kmer_size + 1];
                                                fprintf(stderr, "[cleaning_remove_spurious_links] Something went wrong, couldn't find node C\n");
                                                binary_kmer_to_seq(&node_A->kmer, db_graph->kmer_size, label);
                                                printf("  Node A was %s at coverage %d\n", label, element_get_coverage_all_colours(node_A));
                                                binary_kmer_to_seq(&node_B->kmer, db_graph->kmer_size, label);
                                                printf("  Node B was %s at coverage %d\n", label, element_get_coverage_all_colours(node_B));
                                                binary_kmer_to_seq(&(neighbours[0]->kmer), db_graph->kmer_size, label);
                                                printf("  Neighbour 0 was %s at coverage %d\n", label, element_get_coverage_all_colours(neighbours[0]));
                                                binary_kmer_to_seq(&(neighbours[1]->kmer), db_graph->kmer_size, label);
                                                printf("  Neighbour 1 was %s at coverage %d\n", label, element_get_coverage_all_colours(neighbours[1]));                                                
                                                exit(1);
                                            }
                                            
                                            // If this is what we're looking for, coverage(B) - coverage(C) should be > 0 and <= max_difference
                                            coverage_difference = element_get_coverage_all_colours(node_B) - element_get_coverage_all_colours(node_C);
                                            if ((coverage_difference > 0) && (coverage_difference <= max_difference)) {
                                                // That's it, we've found a link we can remove
                                                char label[db_graph->kmer_size + 1];
                                                printf("Removed link...\n");
                                                binary_kmer_to_seq(&node_A->kmer, db_graph->kmer_size, label);
                                                printf("  Node A was %s at coverage %d\n", label, element_get_coverage_all_colours(node_A));
                                                binary_kmer_to_seq(&node_B->kmer, db_graph->kmer_size, label);
                                                printf("  orientation was %s\n", orientation == forward ? "Fwd":"Rev");
                                                printf("  n1 (to D) was %c and n2 (to B) was %c\n", binary_nucleotide_to_char(n1), binary_nucleotide_to_char(n2));
                                                printf("  Node B was %s at coverage %d\n", label, element_get_coverage_all_colours(node_B));
                                                binary_kmer_to_seq(&node_C->kmer, db_graph->kmer_size, label);
                                                printf("  Node C was %s at coverage %d\n", label, element_get_coverage_all_colours(node_C));
                                                binary_kmer_to_seq(&node_D->kmer, db_graph->kmer_size, label);
                                                printf("  Node D was %s at coverage %d\n", label, element_get_coverage_all_colours(node_D));

                                                db_node_reset_edge_all_colours(node_A, orientation, n2);
                                                db_node_reset_edge_all_colours(node_B, reverse_step_B.orientation, reverse_step_B.label);
                                                links_removed++;
                                            }                                        
                                        }
                                    }
                                }                                
                            }
                        }
                        nucleotide_iterator(&inner_loop);
                    }
                }
                nucleotide_iterator(&outer_loop);
            }
        }
    }
    
    hash_table_traverse(&check_and_remove_edge, db_graph);
    
    return links_removed;
}

static int tip_clip_length;
static dBGraph * ct_db_graph;

static void ct_prune_step_node(pathStep * step) {
    cleaning_prune_db_node(step->node, ct_db_graph);
    
}


static PathEnd is_tip(int length, Path * p){
    pathStep first ;
    pathStep last;
    path_get_step_at_index(0, &first, p);
    path_get_last_step(&last, p);
    
    
    boolean first_blunt = db_node_edges_count_all_colours(first.node, opposite_orientation(first.orientation)) == 0;
    boolean first_blunt_rev = db_node_edges_count_all_colours(first.node, first.orientation) == 0;
    
    boolean first_branch = db_node_edges_count_all_colours(first.node, opposite_orientation(first.orientation)) > 1;
    boolean first_join = db_node_edges_count_all_colours(first.node, first.orientation) > 1;
    
    
    boolean last_blunt = db_node_edges_count_all_colours(last.node, last.orientation) == 0;
    boolean last_blunt_rev = db_node_edges_count_all_colours(last.node, opposite_orientation( last.orientation)) == 0;
    boolean last_branch = db_node_edges_count_all_colours(last.node, last.orientation)> 1;
    boolean last_join = db_node_edges_count_all_colours(last.node, opposite_orientation(last.orientation)) > 1;
    
 
    
        
    
    PathEnd pe = NONE;
    
    if(p->length < length){
        
        if(first_blunt && !first_branch && !first_join){
            if((last_join || last_branch) && !last_blunt && !last_blunt_rev){
                pe = FIRST;
            }
        }
        if(last_blunt && !last_branch && !last_join){
            if((first_join || first_branch) && !first_blunt  && !first_blunt_rev ){
                pe = LAST;
            }
        }
        
        if (first_blunt && last_blunt) {
            if(first_join || first_branch){
                pe = LAST;
            }else if(last_join || last_branch){
                pe = FIRST;
            }
            
        }

        if(first.node == last.node){
            pe = NONE;
        }
    }
    
   /* if(pe == NONE){
        path_to_fasta(p, stdout);
        printf("NO TIP");
    }*/
    return pe;
}


static void check_tip(Path * p){
    pathStep first ;
    pathStep last;
    path_get_step_at_index(0, &first, p);
    path_get_last_step(&last, p);
    
    PathEnd pe = is_tip(tip_clip_length, p);
    
    if (pe != NONE) {
        if(pe == FIRST){
            db_node_action_set_flag(last.node, TIP_START);
        }else if(pe == LAST){
            db_node_action_set_flag(first.node, TIP_START);
        }
    }else{
        double avg; 
        int min, max;
        path_get_statistics(&avg, &min, &max, p);
       /* if(avg > 1 && path_get_length(p) < tip_clip_length){
            path_to_fasta(p, stdout);
            printf("NOT TIP\n");
        }*/
    }
}



static void flag_tips(dBNode * node ){ 
   if (!db_node_check_for_any_flag(node, PRUNED | TIP_START )) {
       Orientation o = undefined;
       boolean go = false;
       Nucleotide n = Undefined;
       
       if(db_node_is_blunt_end(node, forward) && db_node_has_precisely_one_edge_all_colours(node, reverse, &n)){
           o = reverse;
           go = true;
       }else if(db_node_is_blunt_end(node, reverse) && db_node_has_precisely_one_edge_all_colours(node, forward, &n)){
           o = forward;
           go = true;
       }
       if(go){
           
           perfect_path_get_path_with_callback(node, o, &db_node_action_do_nothing, &check_tip, ct_db_graph);
       }
   }
}

static int prune_tip(Path * p){
    pathStep first ;
    pathStep last;
    path_get_step_at_index(0, &first, p);
    path_get_last_step(&last, p);
    
    PathEnd pe = is_tip(tip_clip_length, p);
    
    if (pe == NONE) {
        printf("NOT TIP!\n");
        path_to_fasta(p, stdout);
    }
    assert(pe != NONE);
    
    
    
    if (pe == FIRST) {
        cleaning_prune_db_node(first.node, ct_db_graph);
    }
    if (pe == LAST) {
        cleaning_prune_db_node(last.node, ct_db_graph);
    }
    
    path_inner_iterator(&ct_prune_step_node, p);
    
    
    return path_get_length(p) -1;

}

static void remove_flaged_tips(dBNode * node, void * args ){
    if (db_node_check_for_any_flag(node, TIP_START ) && !db_node_check_for_any_flag(node, PRUNED  )) {
        
        clip_tip_vars * ctv = (clip_tip_vars * ) args;
        ctv->marked_tips++;
        PathArray * paf = path_array_get_from_buffer_with_size(4);
        PathArray * par = path_array_get_from_buffer_with_size(4);
        
        //We get all the pats from the tip. 
        int fwd_found = perfect_path_get_all_paths_from(node, forward, paf, MAX_PATH_LENGTH, ctv->db_graph);
        int rev_found = perfect_path_get_all_paths_from(node, reverse, par, MAX_PATH_LENGTH, ctv->db_graph);
        
        Path * shortest = NULL, *tmp = NULL;
        //Path *sh1 = NULL, *sh2 = NULL;
        int i, shortest_l = ctv->tip_length +1, temp_l;
        
        if (fwd_found > 0) {
            for (i = 0; i < 4; i++) {
                tmp = paf->paths[i];
                if(tmp != NULL && path_get_length(tmp)>0){
                                
                    if(is_tip(tip_clip_length,tmp) != NONE){
                        temp_l =  path_get_length(tmp);
                        if(temp_l < shortest_l  ){
                            shortest = tmp;
                            shortest_l = temp_l;
                            //sh1 = tmp;
                        }
                    }
                }
            }
        }
        
        if(shortest!=NULL){
            int pruned =  prune_tip(shortest);
            if(pruned > 0){
                ctv->tips_nodes_removed += pruned;
                ctv->tips_removed++;
            }
        }
        shortest = NULL;
        if(rev_found > 0){
            for (i = 0; i < 4; i++) {
                tmp = par->paths[i];
                if (tmp != NULL && path_get_length(tmp)>0) {
                    
                    
                    if(is_tip(tip_clip_length,tmp) != NONE){
                        temp_l = path_get_length(tmp);
                        if(temp_l < shortest_l  ){
                            shortest = tmp;
                            shortest_l = temp_l;
                            //sh2 = tmp;
                        }
                    }
                }
               
            }
        }
        
        

        
        if(shortest!=NULL){
            int pruned =  prune_tip(shortest);
            if(pruned > 0){
                ctv->tips_nodes_removed += pruned;
                ctv->tips_removed++;
            }
        }
        
        path_array_free_from_buffer(par);
        path_array_free_from_buffer(paf);
    }
}

int cleaning_remove_tips(int max_length, int max_it, dBGraph * db_graph){
    
    int threads =  1;
#ifdef THREADS
    threads = db_graph->number_of_threads;
#endif
    clip_tip_vars ** ctv = calloc(threads, sizeof(clip_tip_vars *));
    
    int i;
    for(i = 0; i < threads; i++){
        ctv[i] = malloc(sizeof(clip_tip_vars));
        ctv[i]->tips_removed = 0;
        ctv[i]->tips_nodes_removed = 0;
        ctv[i]->db_graph = db_graph;
        ctv[i]->tip_length = max_length;
        ctv[i]->marked_tips = 0;
    }
   
    tip_clip_length = max_length;
    ct_db_graph = db_graph;

    //int tip_count = 0;
    int tips_removed = 0;
    int nodes_removed = 0;
    int marked = 0;
    int tc_it = 1;
    long long prunned_before = 0;
    tips_removed = 0;
    //nodes_removed = 0;
    int current_tips = 0, current_nodes = 0, current_marked = 0;
    
    do {
        prunned_before = db_graph->pruned_kmers;
        //tip_count = tips_removed;
        log_and_screen_printf("clip_tip iteration: %'d\n",tc_it++);
#ifdef THREADS
        hash_table_threaded_traverse(&db_node_action_clear_flags, ct_db_graph);
        hash_table_threaded_traverse(&flag_tips, ct_db_graph);
        hash_table_threaded_traverse_with_args(&removed_flaged_tips, (void **)ctv , db_graph);
#else
        hash_table_traverse(&db_node_action_clear_flags,ct_db_graph);
        hash_table_traverse(&flag_tips,ct_db_graph);
        hash_table_traverse_with_args(&remove_flaged_tips, (void **)ctv , db_graph);
#endif
        current_tips = 0; current_nodes = 0;    current_marked = 0;
        for(i = 0; i < threads; i++){  
            tips_removed += ctv[i]->tips_removed;
            nodes_removed += ctv[i]->tips_nodes_removed ;
            marked += ctv[i]->marked_tips;
            current_tips += ctv[i]->tips_removed;
            current_nodes += ctv[i]->tips_nodes_removed ;
            current_marked += ctv[i]->marked_tips;
            
        }
        log_and_screen_printf("clip_tip removed nodes %'lld, (%'lld tips, %'lld marked)\n", current_nodes, current_tips, current_marked);
        for(i = 0; i < threads; i++){  
            
            ctv[i]->tips_removed = 0;
            ctv[i]->tips_nodes_removed = 0 ;

            ctv[i]->marked_tips = 0;
            //free(ctv[i]);
        }
        hash_table_print_stats(db_graph);
    }while(--max_it > 0 &&  current_tips != 0 && prunned_before != db_graph->pruned_kmers);
    ct_db_graph = NULL;
    for(i = 0; i < threads; i++){  
        free(ctv[i]);
    }
    free(ctv);
    
    log_and_screen_printf("Total tips removed %d nodes from %d tips in %d iterations ", nodes_removed, tips_removed, tc_it);
    return  tips_removed;
}

