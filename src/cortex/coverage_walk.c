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

#include <string.h>
#include <stdint.h>
#ifdef THREADS
#include <pthread.h>
#endif
#include <assert.h>

#include <global.h>
#include <flags.h>
#include <binary_kmer.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <path.h>
#include <perfect_path.h>
#include <logger.h>
#include <metacortex.h>
#include "coverage_walk.h"

//int debugme = 0;

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
static void coverage_walk_pre_step_action(pathStep * ps)
{
    if (ps->orientation == forward) {
        db_node_action_set_flag(ps->node, VISITED_FORWARD);
    } else {
        db_node_action_set_flag(ps->node, VISITED_REVERSE);
    }

}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
static void coverage_walk_post_step_action(pathStep * ps)
{
    if (ps->orientation == forward) {
        db_node_action_unset_flag(ps->node, VISITED_FORWARD);
    } else {
        db_node_action_unset_flag(ps->node, VISITED_REVERSE);
    }
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
Nucleotide coverage_walk_get_best_label(dBNode* node, Orientation orientation, dBGraph* db_graph)
{
    Nucleotide label = Undefined;
    int highest_coverage = 0;

    void check_edge(Nucleotide nucleotide) {
        //if (debugme) log_printf("  Trying nucleotide %c\n", binary_nucleotide_to_char(nucleotide));
        if (db_node_edge_exist_any_colour(node, nucleotide, orientation)) {
            pathStep step, reverse_step, next_step;
            int coverage;

            step.node = node;
            step.label = nucleotide;
            step.orientation = orientation;
            db_graph_get_next_step(&step, &next_step, &reverse_step, db_graph);
            coverage = element_get_coverage_all_colours(next_step.node);

            //if (debugme) log_printf("  Coverage %i\n", coverage);
            if (coverage > highest_coverage) {
                label = nucleotide;
                highest_coverage = coverage;
            }
        }
    }

    nucleotide_iterator(&check_edge);

    return label;
}




/*----------------------------------------------------------------------*
 * Function:  coverage_walk_get_best_label_bubble                       *
 * Purpose:   find first_step->label and first_step->alt_label					*
 *         		walk bubble if encountered. Will return a single					*
 *         		label if one clear alternative, rather than      					*
 *         		multiple possible options     					                  *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
Nucleotide coverage_walk_get_best_label_bubble(pathStep * step, dBNode* node, Orientation orientation, dBGraph* db_graph)
{
    int highest_coverage = 0;
    int highest_coverage_length = 0;
    int all_coverages[4] = {0, 0, 0, 0};	// initialises all elements to zero
    int all_lengths[4] = {0, 0, 0, 0};	// initialises all elements to zero
    int bubble_edge = -1;
    dBNode * nodes[4];	// legal? pointing to other nodes that already exist
    Path* paths[4];
    int i;

    // Clear path array
    for (i=0; i<4; i++) {
        paths[i] = 0;
    }

    step->label=Undefined;

    // check a node for edges
    // if more than one, walk each path
    // are final (branching) nodes the same? then this is a bubble.
    // is sum(coverage) for both paths equal to or highest coverage?
    // keep highest coverage as label, and next highest as alt_label
    // NOTE: only two paths allowed here. Extend later

    // check to see if a simple bubble occurs from this branch point - one that
    //   does rejoin, and has no further branching in between
    void bubble_check()
    {
        int i, j;

        for (i=0; i<4; i++){
            for (j=i+1; j<4; j++){
                if ((all_coverages[i] > 0) &&
                    (all_coverages[j] > 0) &&
                    (nodes[i] == nodes[j]))
                {
                    log_printf("BUBBLE FOUND IN COVERAGE WALK\n");
                    // Only take the bubble route if sum of paths is better than alternative path
                    if ((all_coverages[i] + all_coverages[j]) > highest_coverage) {
                        if (all_coverages[i] > all_coverages[j]) {
                            bubble_edge = i;
                        } else {
                            bubble_edge = j;
                        }
                        highest_coverage = all_coverages[i]+all_coverages[j];
                    } else {
                        // leave highest coverage as it is, ignore bubble
                    }
                } else {
                    // do nothing
                } // i&j exists
            } // j
        } // i
        // DOES THIS WORK WITH BUBBLES THAT DON'T HAVE HIGHEST COVERAGE?
        if (bubble_edge>-1) {
            char seq[1024];

            step->label = bubble_edge;
            db_node_action_set_flag(step->node, POLYMORPHISM);
            db_node_action_set_flag(paths[bubble_edge]->nodes[paths[bubble_edge]->length - 1], POLYMORPHISM);

            binary_kmer_to_seq(&(step->node->kmer), db_graph->kmer_size, seq);
            log_printf("Polymorphism start node %s\n", seq);

            binary_kmer_to_seq(&(paths[bubble_edge]->nodes[paths[bubble_edge]->length - 1]->kmer), db_graph->kmer_size, seq);
            log_printf("Polymorphism end node %s\n", seq);
        }
    }


    void check_edge(Nucleotide nucleotide) {
        if (db_node_edge_exist_any_colour(node, nucleotide, orientation)) {
            pathStep current_step, reverse_step, next_step;
            double avg_coverage;
            int min_coverage;
            int max_coverage;
            int MAX_BRANCH_LENGTH=(db_graph->kmer_size)*2;

            current_step.node = node;
            current_step.label = nucleotide;
            current_step.orientation = orientation;
            db_graph_get_next_step(&current_step, &next_step, &reverse_step, db_graph);

            paths[nucleotide] = path_new(MAX_BRANCH_LENGTH, db_graph->kmer_size);
            db_graph_get_perfect_path_with_first_edge_all_colours(&current_step, &db_node_action_do_nothing, paths[nucleotide], db_graph);
            path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, paths[nucleotide]);
            if(min_coverage<coverage_thresh){
              // destroy path
            }
            all_coverages[nucleotide] = avg_coverage;
            all_lengths[nucleotide] = paths[nucleotide]->length;

            nodes[nucleotide] = paths[nucleotide]->nodes[paths[nucleotide]->length-1];
            // Add end node to list of nodes to visit
        }
    }

    // check for single best edge first on coverage, with length of path breaking ties
    void check_coverages(Nucleotide nucleotide) {
      if (all_coverages[nucleotide] > 1){
        if ((all_coverages[nucleotide] > highest_coverage) ||
            ((all_coverages[nucleotide] == highest_coverage) &&
                (all_lengths[nucleotide] > highest_coverage_length))){
            highest_coverage=all_coverages[nucleotide];
            highest_coverage_length=all_lengths[nucleotide];
            step->label=nucleotide;
        }
      }
    }

    nucleotide_iterator(&check_edge);
    nucleotide_iterator(&check_coverages); // check coverages as individual edges first
    bubble_check();

    for (i=0; i<4; i++) {
        if (paths[i] != 0) path_destroy(paths[i]);
    }

    return 0;
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
pathStep* coverage_walk_get_first_label(pathStep * first_step, dBGraph * db_graph)
{
    //char seq[1024];
    //debugme = 1;
    first_step->label = Undefined;

    if (db_node_edges_count_all_colours(first_step->node, first_step->orientation)==1){ // for simple, non-branching nodes
        first_step->label = coverage_walk_get_best_label(first_step->node, first_step->orientation, db_graph);
    } else {
        coverage_walk_get_best_label_bubble(first_step, first_step->node, first_step->orientation, db_graph);
    }

    //debugme = 0;
    //binary_kmer_to_seq(&(first_step->node->kmer), db_graph->kmer_size, seq);
    //log_printf("  First kmer %s first label %c\n", seq, binary_nucleotide_to_char(first_step->label));
    return first_step;
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
static pathStep *coverage_walk_get_next_step(pathStep * current_step, pathStep * next_step, pathStep * reverse_step, dBGraph * db_graph)
{
    db_graph_get_next_step(current_step, next_step, reverse_step, db_graph);
    assert(next_step != NULL);
    next_step->label = Undefined;

    if (next_step->node != NULL) {
        if (db_node_edges_count_all_colours(next_step->node, next_step->orientation) == 1) {
            next_step->label = coverage_walk_get_best_label(next_step->node, next_step->orientation, db_graph);
        } else if (db_node_edges_count_all_colours(next_step->node, next_step->orientation) > 1) {
            coverage_walk_get_best_label_bubble(next_step, next_step->node, next_step->orientation, db_graph);
        }
    }

    return next_step;
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
static boolean coverage_walk_continue_traversing(pathStep * current_step,
                                                 pathStep * next_step,
                                                 pathStep * reverse_step,
                                                 Path * temp_path,
                                                 dBGraph * db_graph)
{
    pathStep first;

    boolean cont;
    cont = current_step->label != Undefined;

    /* We don't do these checks for the first node - in case it's a Y node */
    if (temp_path->length > 1) {
        /* Check for a cycle - as this is a perfect path, we only need to check the first node. If we come
         back in at one of the other nodes, then it will result in two edges in one orientation */
        path_get_step_at_index(0, &first, temp_path);
        if (path_step_equals_without_label(&first, current_step)) {
            //char seq[1024];
            //binary_kmer_to_seq(&(current_step->node->kmer), db_graph->kmer_size, seq);
            //log_printf("  Stopped for cycle at %s\n", seq);
            path_add_stop_reason(LAST, PATH_FLAG_IS_CYCLE, temp_path);
            cont = false;
        }

        /* Check for visited flag */
        if (db_node_check_for_any_flag(next_step->node, next_step->orientation == forward? VISITED_FORWARD:VISITED_REVERSE)) {
            cont = false;
        }

        /* Now check for one or more edges moving forward */
        if (db_node_edges_count_all_colours(current_step->node, current_step->orientation) == 0) {
            //char seq[1024];
            //binary_kmer_to_seq(&(current_step->node->kmer), db_graph->kmer_size, seq);
            //printf("  Stopped for blunt end at %s\n", seq);
            path_add_stop_reason(LAST, PATH_FLAG_STOP_BLUNT_END, temp_path);
            cont = false;
        }

        /* Check path has space */
        if (!path_has_space(temp_path)) {
            //char seq[1024];
            //binary_kmer_to_seq(&(current_step->node->kmer), db_graph->kmer_size, seq);
            //log_printf("  Stopped for longer than buffer at %s\n", seq);
            path_add_stop_reason(LAST, PATH_FLAG_LONGER_THAN_BUFFER, temp_path);
            cont = false;
        }
    }

    return cont;
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
WalkingFunctions * coverage_walk_get_funtions(WalkingFunctions *walking_functions)
{
    perfect_path_get_funtions(walking_functions);

    // Which to over-rule?
    walking_functions->continue_traversing = &coverage_walk_continue_traversing;
    walking_functions->get_next_step = &coverage_walk_get_next_step;
    walking_functions->pre_step_action = &coverage_walk_pre_step_action;
    walking_functions->post_step_action =&coverage_walk_post_step_action;

    return walking_functions;
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
int coverage_walk_get_path_with_callback(dBNode * node, Orientation orientation,
                                         void (*node_action) (dBNode * node),
                                         void (*path_action) (Path * path),
                                         dBGraph * db_graph)
{
    // Get walking functions
    WalkingFunctions wf;
    coverage_walk_get_funtions(&wf);

    // Setup first step
    pathStep first;
    first.node = node;
    first.orientation = orientation;
    first.label = Undefined;
    wf.get_starting_step = &coverage_walk_get_first_label;

    // Setup step action to include passed in node action
    //void (*action) (pathStep * step);
    //action = wf.step_action;

    void local_step_action(pathStep * ps) {
        //action(ps);
        node_action(ps->node);
        return;
    }
    if (node_action != NULL) {
        wf.step_action = &local_step_action;
    }

    // Setup path action to include passed in path action
    void (*action_path) (Path * p);

    action_path = wf.output_callback;

    void local_path_action(Path * p) {
        action_path(p);
        path_action(p);
        return;
    }
    wf.output_callback = local_path_action;

    // Get a buffer for this path
    Path *path = path_get_buffer_path();

    // Do the walk
    int ret = db_graph_generic_walk(&first, path, &wf, db_graph);

    // Free buffer
    path_free_buffer_path(path);

    return ret;
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
int coverage_walk_get_path(dBNode * node, Orientation orientation, void (*node_action) (dBNode * node), dBGraph * db_graph, Path * path)
{
    void copy_path(Path * p) {
        path_copy(path, p);
    }

    log_printf("\t[NEW-NODE in coverage walk]\n");
    coverage_walk_get_path_with_callback(node, orientation,	node_action, &copy_path, db_graph);

    return path_get_edges_count(path);
}




///////////////////////////////////////////////////////////////////////////////////////

void execute_path_step_callbacks_MinCov(pathStep * p, PathStepActionCallbackArray * callbacks){
    int i;
    for (i = 0; i < callbacks->used; i++) {
        void  (* f)() = callbacks->callback[i];
        if (callbacks->args[i] == NULL) {
            f(p);
        }else{
            f(p, callbacks->args[i]);
        }
    }
}


void execute_path_callbacks_MinCov(Path * p, PathCallbackArray * callbacks){
    int i;
    for (i = 0; i < callbacks->used; i++) {
        void  (* f)() = callbacks->callback[i];
        if (callbacks->args[i] == NULL) {
            f(p);
        }else{
            f(p, callbacks->args[i]);
        }
    }
}

void execute_node_callbacks_MinCov(dBNode * n, NodeActionCallbackArray * callbacks){
    int i;
    for (i = 0; i < callbacks->used; i++) {
        void  (* f)() = callbacks->callback[i];
        if (callbacks->args[i] == NULL) {
            f(n);
        }else{
            f(n, callbacks->args[i]);
        }
    }
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
pathStep* min_coverage_walk_get_first_label(pathStep * first_step, dBGraph * db_graph, int coverage_thresh)
{
    //char seq[1024];
    //debugme = 1;
    first_step->label = Undefined;

    if (db_node_edges_count_all_colours(first_step->node, first_step->orientation)==1){ // for simple, non-branching nodes
        first_step->label = coverage_walk_get_best_label(first_step->node, first_step->orientation, db_graph);
        if(coverage_thresh>){
          first_step->label = Undefined;
        }
    } else {
        min_coverage_walk_get_best_label_bubble(first_step, first_step->node, first_step->orientation, db_graph, coverage_thresh);
    }

    //debugme = 0;
    //binary_kmer_to_seq(&(first_step->node->kmer), db_graph->kmer_size, seq);
    //log_printf("  First kmer %s first label %c\n", seq, binary_nucleotide_to_char(first_step->label));
    return first_step;
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
boolean min_coverage_walk_continue_traversing(pathStep * current_step,
                                                 pathStep * next_step,
                                                 pathStep * reverse_step,
                                                 Path * temp_path,
                                                 dBGraph * db_graph,
                                                 int coverage_thresh)
{
    pathStep first;

    boolean cont;
    cont = current_step->label != Undefined;

    /* We don't do these checks for the first node - in case it's a Y node */
    if (temp_path->length > 1) {
        /* Check for a cycle - as this is a perfect path, we only need to check the first node. If we come
         back in at one of the other nodes, then it will result in two edges in one orientation */
        path_get_step_at_index(0, &first, temp_path);
        if (path_step_equals_without_label(&first, current_step)) {
            path_add_stop_reason(LAST, PATH_FLAG_IS_CYCLE, temp_path);
            cont = false;
        }

        /* Check for visited flag */
        if (db_node_check_for_any_flag(next_step->node, next_step->orientation == forward? VISITED_FORWARD:VISITED_REVERSE)) {
            cont = false;
        }

        /* Now check for one or more edges moving forward */
        if (db_node_edges_count_all_colours(current_step->node, current_step->orientation) == 0) {
            path_add_stop_reason(LAST, PATH_FLAG_STOP_BLUNT_END, temp_path);
            cont = false;
        }

        /* Check path has space */
        if (!path_has_space(temp_path)) {
            path_add_stop_reason(LAST, PATH_FLAG_LONGER_THAN_BUFFER, temp_path);
            cont = false;
        }

        /* Check next step has high enough coverage */
        if (coverage_thresh > element_get_coverage_all_colours(next_step->node)) {
            path_add_stop_reason(LAST, PATH_FLAG_STOP_BLUNT_END, temp_path); // not ideal reason
            cont = false;
        }
    }

    return cont;
}



int db_graph_MinCov_walk(pathStep * first_step, Path * path, WalkingFunctions * functions, dBGraph * db_graph, int coverage_threshold)
{
	//dBNode * node = first_step->node;
	//sanity check
	if (first_step == NULL) {
		printf("[db_graph_generic_walk] can't pass a null node\n");
        assert(0);
		exit(1);
	}

	if (path == NULL) {
		printf("[db_graph_generic_walk] can't pass a null path\n");
        assert(0);
		exit(1);
	}

	if (functions == NULL) {
		printf("[db_graph_generic_walk] can't pass  null functions\n");
        assert(0);
		exit(1);
	}

	if (db_graph == NULL) {
		printf("[db_graph_generic_walk] can't pass  null db_graph\n");
        assert(0);
		exit(1);
	}

	pathStep current_step, next_step, rev_step;

	path_step_assign(&next_step, first_step);
    next_step.path = path; //This way, on all the assignments we keep the pointer to the path that was sent to the function originally.
	functions->get_starting_step(&next_step, db_graph);

	//boolean  try = true;
	int count = 0;
    boolean walked;
    boolean added;

	do {
		walked = false;
        do {
            added = path_add_node(&next_step, path);
            path_step_assign(&current_step, &next_step);
            //try = false;
            functions->pre_step_action(&current_step);

            added = false;
            if (current_step.label != Undefined) {
                functions->get_next_step(&current_step, &next_step, &rev_step, db_graph);
            }

            if (added == false) {
                // Do something?
            }

            functions->step_action(&current_step);
            execute_path_step_callbacks_MinCov(&current_step,&functions->step_actions);
            execute_node_callbacks_MinCov(current_step.node, &functions->node_callbacks);
        }
        while (min_coverage_walk_continue_traversing(&current_step, &next_step, &rev_step, path,db_graph, coverage_threshold));

		if (path->length > 0) {
			walked = true;
		}

		if (walked) {
			count++;
			functions->output_callback(path);
            execute_path_callbacks_MinCov(path, &functions->path_callbacks);

		}
		do {
			pathStep ps;

			path_get_last_step(&ps, path);
			if(ps.node != NULL){
                functions->post_step_action(&ps);
                path_remove_last(path);
            }
		}while (functions->continue_backwards(path, db_graph));
        if(path_get_length(path) > 0){//This is to enable the backtracking.
            path_get_last_step(&current_step, path);//The continue_backwards should fix the last step
            functions->get_next_step(&current_step,&next_step, &rev_step,db_graph);

        }
	}while (path_get_length(path) > 0);
	return count;
}



/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
int min_coverage_walk_get_path_with_callback(dBNode * node, Orientation orientation,
                                         void (*node_action) (dBNode * node),
                                         void (*path_action) (Path * path),
                                         dBGraph * db_graph, int coverage_thresh)
{
    // Get walking functions
    WalkingFunctions wf;
    coverage_walk_get_funtions(&wf);

    // Setup first step
    pathStep first;
    first.node = node;
    first.orientation = orientation;
    first.label = Undefined;
    wf.get_starting_step = &min_coverage_walk_get_first_label;

    // Setup step action to include passed in node action
    //void (*action) (pathStep * step);
    //action = wf.step_action;

    void local_step_action(pathStep * ps) {
        //action(ps);
        node_action(ps->node);
        return;
    }
    if (node_action != NULL) {
        wf.step_action = &local_step_action;
    }

    // Setup path action to include passed in path action
    void (*action_path) (Path * p);

    action_path = wf.output_callback;

    void local_path_action(Path * p) {
        action_path(p);
        path_action(p);
        return;
    }
    wf.output_callback = local_path_action;

    // Get a buffer for this path
    Path *path = path_get_buffer_path();

    // Do the walk
    int ret = db_graph_MinCov_walk(&first, path, &wf, db_graph, coverage_thresh);

    // Free buffer
    path_free_buffer_path(path);

    return ret;
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
int min_coverage_walk_get_path(dBNode * node, Orientation orientation, void (*node_action) (dBNode * node), dBGraph * db_graph, Path * path, int coverage_thresh)
{
    void copy_path(Path * p) {
        path_copy(path, p);
    }

    log_printf("\t[NEW-NODE in coverage walk]\n");
    min_coverage_walk_get_path_with_callback(node, orientation,	node_action, &copy_path, db_graph, coverage_thresh);

    return path_get_edges_count(path);
}
