/*
This file is a part of ``RRT(*)'', an incremental 
sampling-based optimal motion planning library.
Copyright (c) 2010 Sertac Karaman <sertac@mit.edu>, Emilio Frazzoli <frazzoli@mit.edu>

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
*/

/*
RRT* algorithm is explained in the following paper:
http://arxiv.org/abs/1005.0416
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include <glib.h>

#include "opttree.h"


// Allocate memory for a new node
node_t* opttree_new_node (opttree_t *self) {
    node_t *node_new = (node_t *) malloc (sizeof (node_t));
    node_new->parent = NULL;
    node_new->children = NULL;
    node_new->state = optsystem_new_state (self->optsys); 
    node_new->reaches_target = FALSE;
    node_new->distance_from_root = -1.0;
    node_new->distance_from_parent = -1.0;
    node_new->traj_from_parent = NULL;
    node_new->inputs_from_parent = NULL;
    return node_new;
}


// Allocate memory for a new node, but does not allocate memory for the state
node_t* opttree_new_node_no_state () {
    node_t *node_new = (node_t *) malloc (sizeof (node_t));
    node_new->parent = NULL;
    node_new->children = NULL;
    node_new->state = NULL;
    node_new->reaches_target = FALSE;
    node_new->distance_from_root = -1.0;
    node_new->distance_from_parent = -1.0;
    node_new->traj_from_parent = NULL;
    node_new->inputs_from_parent = NULL;
    return node_new;
}


// Free node 
int opttree_free_node (opttree_t *self, node_t *node) {
    
    optsystem_free_state (self->optsys, node->state);
    g_slist_free (node->children);
    
    GSList *traj_ptr = node->traj_from_parent;
    while (traj_ptr) {
        optsystem_free_state (self->optsys, (state_t *)(traj_ptr->data));
        traj_ptr = g_slist_next (traj_ptr);
    }
    g_slist_free (node->traj_from_parent);

    GSList *inputs_ptr = node->inputs_from_parent;
    while (inputs_ptr) {
        optsystem_free_input (self->optsys, (input_t*)(inputs_ptr->data));
        inputs_ptr = g_slist_next (inputs_ptr);
    }
    g_slist_free (node->inputs_from_parent);

    free (node);
    
    return 1;
}


// Find the nearest neighbor in the tree
node_t* opttree_find_nearest_neighbor (opttree_t *self, state_t *state_from) {

    node_t *min_node = NULL;

    kdres_t *kdres = kd_nearest (self->kdtree, optsystem_get_state_key (self->optsys, state_from));
    if (kd_res_end (kdres))  {
        printf ("ERROR: No nearest neighbors\n");
        exit(1);
    }
    min_node = kd_res_item_data (kdres);
    kd_res_free (kdres);
    
    return min_node;
}


// Recursively update the distances
void opttree_update_distance_from_root (opttree_t *self, node_t *node_parent, node_t *node_curr) {

    node_curr->distance_from_root = node_parent->distance_from_root + node_curr->distance_from_parent;

    // Check for reachability of the target, if so update the lower_bound
    if (optsystem_is_reaching_target (self->optsys, node_curr->state)) {
        if (node_curr->distance_from_root < self->lower_bound) {
            self->lower_bound = node_curr->distance_from_root;
            self->lower_bound_node = node_curr;
        }
        node_curr->reaches_target = 1;
    }
    else 
        node_curr->reaches_target = 0;
    
    GSList *node_child_list = node_curr->children;
    while (node_child_list){
        opttree_update_distance_from_root (self, node_curr, (node_t *)(node_child_list->data));
        node_child_list = g_slist_next (node_child_list);
    }
}


// Adds a given trajectory to the graph and returns a pointer to the the last node
node_t *opttree_add_traj_to_graph (opttree_t *self, node_t *node_start, node_t *node_end, 
                                   GSList *trajectory, int num_node_states, int *node_states, GSList *inputs, double *radius) {

    node_t *node_prev = node_start;   // This variable maintains the node last added to the graph to update parents

    int trajectory_count = 0;
    int node_states_count = 0;

    GSList *subtraj_curr = NULL;
    GSList *inputs_curr = NULL;

    GSList *inputs_ptr = inputs; 
    GSList *trajectory_ptr = trajectory;
    while (trajectory_ptr) {

        state_t *state_curr = (state_t *) (trajectory_ptr->data);

        // Check whether this the end of the trajectory
        if (!g_slist_next (trajectory_ptr)) {
            // This must be the last node in the traj.

            node_t *node_curr;
            if (node_end) {    // If node_end is given, then modify node end accordingly
                
                node_curr = node_end;
                node_t *node_parent = node_curr->parent;
                node_parent->children = g_slist_remove (node_parent->children, (gpointer) node_curr);
                // Free traj_from_parent
                GSList *traj_from_parent_ptr = node_curr->traj_from_parent;
                while (traj_from_parent_ptr) {
                    optsystem_free_state (self->optsys, (state_t *) (traj_from_parent_ptr->data));
                    traj_from_parent_ptr = g_slist_next (traj_from_parent_ptr);
                }
                g_slist_free (node_curr->traj_from_parent);
                node_curr->traj_from_parent = NULL;
                // Free input_from_parent
                GSList *inputs_from_parent_ptr = node_curr->inputs_from_parent;
                while (inputs_from_parent_ptr) {
                    optsystem_free_input (self->optsys, (input_t *) (inputs_from_parent_ptr->data));
                    inputs_from_parent_ptr = g_slist_next (inputs_from_parent_ptr);
                }
                g_slist_free (node_curr->inputs_from_parent);
                node_curr->inputs_from_parent = NULL;
                
                // Free this state
                optsystem_free_state (self->optsys, state_curr);
                
            }
            else {   // If node_end is not given, then insert this state into the graph 
                node_curr = opttree_new_node_no_state ();
                node_curr->children = NULL;
                node_curr->state = state_curr;
                node_curr->reaches_target = optsystem_is_reaching_target (self->optsys, node_curr->state);
                node_curr->bowl_radius = *radius;

                kd_insert (self->kdtree, optsystem_get_state_key (self->optsys, node_curr->state), node_curr);

                // Add node to the graph
                self->list_nodes = g_slist_prepend (self->list_nodes, (gpointer) (node_curr));
                self->num_nodes++;
                
            }
            node_curr->parent = node_prev;
            if (inputs) {
                input_t *input_this = (input_t *)(inputs_ptr->data);
                inputs_curr = g_slist_prepend (inputs_curr, input_this);
            }
            node_curr->inputs_from_parent = g_slist_reverse (inputs_curr);
            node_curr->traj_from_parent = g_slist_reverse (subtraj_curr); 
            node_curr->distance_from_parent = optsystem_evaluate_distance_for_cost (self->optsys,
                                                                                    node_curr->inputs_from_parent);
            opttree_update_distance_from_root (self, node_prev, node_curr);

            // Add this node to the children of the previous node
            node_prev->children = g_slist_prepend (node_prev->children, node_curr);

            // Reevaluate reaches_target variables
            if (node_curr->reaches_target) {
                if (! (node_prev->reaches_target) ) {
                    node_t *node_temp = node_prev;
                    while (node_temp )  {
                        node_temp->reaches_target = TRUE;
                        node_temp = node_temp->parent;
                    }
                }
                
                if (node_curr->reaches_target) {
                    if (node_curr->distance_from_root < self->lower_bound) {
                        self->lower_bound = node_curr->distance_from_root;
                        self->lower_bound_node = node_curr;
                    }
                }
            }
            
            // Reset the pointers
            subtraj_curr = NULL;
            node_prev = node_curr;
            
            goto end_iteration;
            
        }
        
        if (node_states) {
            if ( trajectory_count == node_states[node_states_count] ) {
                
                // Create the new node
                node_t *node_curr = opttree_new_node_no_state ();
                node_curr->state =state_curr;
                node_curr->parent = node_prev;
                node_curr->children = NULL;
                node_curr->reaches_target = optsystem_is_reaching_target (self->optsys, node_curr->state);
                if (inputs) {
                    inputs_curr = g_slist_prepend (inputs_curr, inputs_ptr->data);
                }
                node_curr->inputs_from_parent = g_slist_reverse (inputs_curr);
                node_curr->traj_from_parent = g_slist_reverse (subtraj_curr);
                node_curr->distance_from_parent = optsystem_evaluate_distance_for_cost (self->optsys,
                                                                                        node_curr->inputs_from_parent);
                node_curr->distance_from_root = node_prev->distance_from_root + node_curr->distance_from_parent;

                kd_insert (self->kdtree, optsystem_get_state_key (self->optsys, node_curr->state), node_curr);

                // Update the parent node children pointer
                node_prev->children = g_slist_prepend (node_prev->children, node_curr);

                // Reevaluate reaches_target variables
                if (node_curr->reaches_target) {
                    if (! (node_prev->reaches_target) ) {
                        node_t *node_temp = node_prev;
                        while (node_temp )  {
                            node_temp->reaches_target = TRUE;
                            node_temp = node_temp->parent;
                        }
                    }

                    if (node_curr->distance_from_root < self->lower_bound) {
                        self->lower_bound = node_curr->distance_from_root;
                        self->lower_bound_node = node_curr;
                    }
                }

                self->list_nodes = g_slist_prepend (self->list_nodes, node_curr);
                self->num_nodes++;
                
                // Reset the pointers
                subtraj_curr = NULL;
                inputs_curr = NULL;

                node_prev = node_curr;

                // Increase the node_states counter
                node_states_count++;

                goto end_iteration;
                
            }
        }

        // Add current state to the subtrajectory
        subtraj_curr = g_slist_prepend (subtraj_curr, state_curr);
        if (inputs)
            inputs_curr = g_slist_prepend (inputs_curr, inputs_ptr->data);
        
    end_iteration:

        trajectory_ptr = g_slist_next (trajectory_ptr);
        if (inputs) {
            inputs_ptr = g_slist_next (inputs_ptr);
        }
        trajectory_count++;
    }
    
    g_slist_free (trajectory);
    g_slist_free (inputs);

    free (node_states);
            
    return node_prev;
}

// if distance == 1000, it means that the particle collides,
// those particles are not counted in the mean
float calci_stdev(float *d, int size)
{
    float mean = 0;
    for(int i=0; i<size; i++)
    {
        if( *(d+i) != 1000)
            mean += *(d+i);
        else
            size --;
    }
    mean /= size;

    float var = 0;
    for(int i=0; i<size; i++)
    {
        if( *(d+i) != 1000)
            var += ((*(d+i)) - mean)*((*(d+i)) - mean);
    }
    return sqrt(var/size);
}


void create_noisy_state(state_t *curr, double radius)
{
    double t = radius;
    for(int i=0; i<NUM_STATES; i++)
    {
        curr->x[i] += (2*t*random()/(RAND_MAX + 1.0) - t);
    }
}

// Extends a given state towards another state
int optsystem_extend_to (opttree_t *tree, optsystem_t *self, state_t *state_from, state_t *state_towards, 
                         int *fully_extends, GSList **trajectory,
                         int *num_node_states, int **nodes_states, GSList **inputs, double *radius) {

    int discretization_num_steps = 10;
    *radius = 1.0;

    GSList *trajectory_curr = NULL; // Start with an empty trajectory
    GSList *inputs_curr = NULL;

    double dist_x = state_towards->x[0] - state_from->x[0];
    double dist_y = state_towards->x[1] - state_from->x[1];

    double dist = sqrt (dist_x * dist_x + dist_y * dist_y);

    //printf("extending from: %3.4f %3.4f\n", state_from->x[0], state_from->x[1]);
    //printf("extending to: %3.4f %3.4f\n", state_towards->x[0], state_towards->x[1]);
    if (dist < 1.0) 
    {
        if (optsystem_segment_on_obstacle (self, state_from, state_towards, discretization_num_steps) ) 
        {
            *fully_extends = 0;
            return 0;
        }

        int count = 0;
        // try 100 particles from same state
        for(int i=0; i<100; i++)
        {
            if(propagate_to_root(tree, state_from))
                count++;
        }
        //printf("count: %d\n", count);
        if(count > 98)      // extend this
        {
            state_t *state_new = optsystem_new_state (self);
            state_new->x[0] = state_towards->x[0] + NOISE;
            state_new->x[1] = state_towards->x[1] + NOISE;

            if (optsystem_segment_on_obstacle (self, state_from, state_new, discretization_num_steps) ) 
            {
                *fully_extends = 0;
                optsystem_free_state(self, state_new);
                return 0;
            }

            trajectory_curr = g_slist_prepend (trajectory_curr, state_new);        
            //printf("extended to: %3.4f %3.4f\n", state_new->x[0], state_new->x[1]);

            input_t *input_new = optsystem_new_input (self);
            input_new->x[0] = dist;
            inputs_curr = g_slist_prepend (inputs_curr, input_new);
            *fully_extends = 1;
        }
        else
        {
            *fully_extends = 0;
            return 0;
        }
    }
    else 
    { 
        *fully_extends = 0;
        int count = 0;
        
        state_t *state_new = optsystem_new_state (self);  
        for(int i=0; i<100; i++)
        {
            state_new->x[0] = (state_towards->x[0] - state_from->x[0])/dist + state_from->x[0];
            state_new->x[1] = (state_towards->x[1] - state_from->x[1])/dist + state_from->x[1];
            
            if(propagate_to_root(tree, state_from))
                count++;
        }
        optsystem_free_state(self, state_new);
        //printf("count: %d\n", count);
        
        if(count > 98)
        {
            state_t *state_new = optsystem_new_state (self);  
            state_new->x[0] = (state_towards->x[0] - state_from->x[0])/dist + state_from->x[0];
            state_new->x[1] = (state_towards->x[1] - state_from->x[1])/dist + state_from->x[1];
            state_new->x[0] += NOISE;
            state_new->x[1] += NOISE;

            if (optsystem_segment_on_obstacle(self, state_from, state_new, discretization_num_steps)) 
            {
                optsystem_free_state(self, state_new);
                return 0;
            }

            //printf("extended to: %3.4f %3.4f\n", state_new->x[0], state_new->x[1]);
            trajectory_curr = g_slist_prepend (trajectory_curr, state_new);        

            input_t *input_new = optsystem_new_input (self);
            input_new->x[0] = 1.0;
            inputs_curr = g_slist_prepend (inputs_curr, input_new);
        }
        else
        {
            return 0;
        }
    }
    
    //printf("radius: %f\n", *radius);
    *trajectory = trajectory_curr;
    *inputs = inputs_curr;
    *num_node_states = 0;
    *nodes_states = NULL;
    
    return 1;
}


// Extend node towards a given sample
// If node can not be extended , e.g., a collision occurs, then returns FALSE
node_t *opttree_extend_towards_sample (opttree_t *self, node_t *node_from, state_t *state_towards) {

    // Extend the node towards the sample
    int fully_extends = 0;
    GSList *trajectory = NULL;
    int num_node_states = 0;
    int  *node_states = NULL;
    GSList *inputs = NULL;
    double *radius;
    radius = calloc(1, sizeof(double));

    if (optsystem_extend_to (self, self->optsys, node_from->state, state_towards, 
                             &fully_extends, &trajectory, &num_node_states, &node_states, &inputs, radius) == 0) {
        return NULL;
    }
    
    // Add all the nodes in the trajectory to the tree
    node_t *node_curr = opttree_add_traj_to_graph (self, node_from, NULL, trajectory, num_node_states, node_states, inputs, radius);
    
    // Check for reachability of the target, if so update the lower_bound
    if (optsystem_is_reaching_target (self->optsys, node_curr->state)) 
    {
        if (node_curr->distance_from_root < self->lower_bound) 
        {
            self->lower_bound = node_curr->distance_from_root;
            self->lower_bound_node = node_curr;
        }
        node_curr->reaches_target = 1;
    }
    else 
        node_curr->reaches_target = 0;
    
    return node_curr; // Return the last node
}


// Extends a given node towards state_towards and returns the resulting state
//    does not generate a new node or populate the tree
state_t *opttree_extend_towards_sample_no_create_node (opttree_t *self, node_t *node_from, state_t *state_towards) {

    state_t *state = NULL;

    int fully_extends = 0;
    GSList *trajectory = NULL;
    int num_node_states = 0;
    int *node_states = NULL;
    GSList *inputs = NULL;
    double *radius;
    radius = calloc(1, sizeof(double));
    
    if (optsystem_extend_to (self, self->optsys, node_from->state, state_towards, 
                             &fully_extends, &trajectory, &num_node_states, &node_states, &inputs, radius) == 0)
        return NULL;

    {                                      
        // Get the last state in the trajectory
        GSList *state_trajectory = trajectory;
        while (state_trajectory) {
            state_t *state_curr = state_trajectory->data;
            if (!g_slist_next (state_trajectory))
                state = optsystem_clone_state (self->optsys, state_curr);
            state_trajectory = g_slist_next (state_trajectory); 
        }
    }
    
    { // Define some local variables
        GSList *state_trajectory = trajectory;
        while (state_trajectory) {
            state_t *state_curr = (state_t *)(state_trajectory->data);
            optsystem_free_state (self->optsys, state_curr);
            state_trajectory = g_slist_next (state_trajectory); 
        }
        g_slist_free (trajectory);
        GSList *inputs_ptr = inputs;
        while (inputs_ptr) {
            input_t *input_curr = (input_t *)(inputs_ptr->data);
            optsystem_free_input (self->optsys, input_curr);
            inputs_ptr = g_slist_next (inputs_ptr);
        }
        g_slist_free (inputs);
        
        // Free node_states
        free (node_states);
    }

    return state;
}


// Goes through all the nodes in node_list and returns a pointer to the node that 
//    gets to state_towards with minimum cost
node_t* opttree_find_min_node_in_set (opttree_t *self, state_t *state_towards, GSList *list_nodes) {

    node_t *min_node = NULL;
    double min_cost = DBL_MAX;

    GSList *node_curr_list = list_nodes; 
    while (node_curr_list ){
        node_t *node_curr = (node_t *) node_curr_list->data;
        
        
        // Extend this node towards state_towards and evaluate the cost incurred
        //    -- if the node does not extend, then move on             
        int fully_extends = 0;
        GSList *trajectory = NULL;
        int num_node_states = 0;
        int *node_states = NULL;
        GSList *inputs = NULL;
        
        double *radius;
        radius = calloc(1, sizeof(double));

        if (optsystem_extend_to (self, self->optsys, node_curr->state, state_towards, 
                                 &fully_extends, &trajectory, &num_node_states, &node_states, &inputs, radius) != 0) {

            if (fully_extends)
            {
                // Evaluate the total cost
                double total_cost = node_curr->distance_from_root;
                
                double incremental_cost = optsystem_evaluate_distance_for_cost (self->optsys, inputs);
                
                total_cost += incremental_cost;
                
                if (total_cost < min_cost) {
                    min_node = node_curr;
                    min_cost = total_cost;
                }
            }

            // Free the states
            {
                GSList* traj_temp = trajectory;
                while (traj_temp) {
                    state_t *state_this = (state_t *)(traj_temp->data);
                    optsystem_free_state (self->optsys, state_this);
                    traj_temp = g_slist_next (traj_temp);
                }
            }
            // Free the trajectory 
            g_slist_free (trajectory); 
            // Free Inputs
            {
                GSList *inputs_temp = inputs; 
                while (inputs_temp) {
                    input_t *input_this = (input_t *)(inputs_temp->data);
                    optsystem_free_input (self->optsys, input_this);
                    inputs_temp = g_slist_next (inputs_temp);
                }
            }
            // Free the inputs list
            g_slist_free (inputs);
            // Free node_states
            free (node_states);
            
        }
        node_curr_list = g_slist_next (node_curr_list);
    }
    

    return min_node;
}


// Takes a kdres set and returns a list of nodes
GSList *opttree_kdtree_to_gslist (state_t * state, kdres_t *kdres) {
    
    GSList *node_list = NULL;

    kd_res_rewind (kdres);
    while (!kd_res_end(kdres)) {
        node_t * node_curr = kd_res_item_data (kdres);
        node_list = g_slist_prepend (node_list, node_curr);
        kd_res_next (kdres);
    }
    
    return node_list;
}


// Finds the set of nodes with max cost
GSList *opttree_find_nodes_in_ball (opttree_t *self, state_t *state, double ball_radius) {

    GSList *nodes_in_ball = NULL; 

    kdres_t *kdres = kd_nearest_range (self->kdtree, optsystem_get_state_key (self->optsys, state), ball_radius);
    nodes_in_ball = opttree_kdtree_to_gslist (state, kdres);
    kd_res_free (kdres);
    
    return nodes_in_ball;
}


// Extends the node back to the tree 
int opttree_extend_back_to_tree (opttree_t *self, node_t *node_from, GSList *node_list) {
    

    state_t *state_from = node_from->state;

    GSList *node_curr_list = node_list;    
    while (node_curr_list) {
        node_t *node_curr = (node_t *) (node_curr_list->data);
        
        if (node_curr == node_from){ // If current node is the same node_from, then continue normal operation
            node_curr_list = g_slist_next (node_curr_list);
            continue;
        }

        if (node_curr == self->root) {
            node_curr_list = g_slist_next (node_curr_list);
            continue;
        }

        state_t *state_towards = node_curr->state;


        gboolean free_states = FALSE;

        // Try to extend the state towards the sample
        int fully_extends = 0; 
        GSList *trajectory = NULL;
        int num_node_states = 0;
        int *node_states = NULL;
        GSList *inputs = NULL;
        double *radius;
        radius = calloc(1, sizeof(double));
        if (optsystem_extend_to (self, self->optsys, state_from, state_towards, 
                                 &fully_extends, &trajectory, &num_node_states, &node_states, &inputs, radius) == 0) {
            fully_extends = 0;
        }
        
        if (fully_extends) {   // Consider rewiring the tree
            
            // Calculate the cummulative distance from the root till the end of the trajectory
            double dist = node_from->distance_from_root;
            dist += optsystem_evaluate_distance_for_cost (self->optsys, inputs);

            // Check whether the new branch is shorter than the existing one 
            if (dist < node_curr->distance_from_root) {

                // Add the trajectory to the tree
                node_t *node_new; 
                node_new = opttree_add_traj_to_graph (self, node_from, node_curr, trajectory, num_node_states, node_states, inputs, radius);
            } 
            else {                   // If the new distance from the root is not less than the current one
                free_states = TRUE;  // remember to free up the memory used by the states in the trajectory
            }
        } 
        else{                      // If it does not fully extend
            free_states = TRUE;    // remember to free up the memory used by the states in the trajectory
        }                          //     OBSTACLES CAUSE NO EXTENSIONS, IN WHICH CASE fully_extend = 0
        
        // Free up the memory used by the trajectory if the trajectory is not registered into the graph
        if (free_states) {
            GSList *trajectory_ptr = trajectory;
            while (trajectory_ptr) {
                optsystem_free_state (self->optsys, (state_t *)(trajectory_ptr->data));
                trajectory_ptr = g_slist_next (trajectory_ptr);
            }
            g_slist_free (trajectory);
            GSList *inputs_ptr = inputs;
            while (inputs_ptr) {
                optsystem_free_input (self->optsys, (input_t *)(inputs_ptr->data));
                inputs_ptr = g_slist_next (inputs_ptr);
            }
            g_slist_free (inputs);
            free (node_states);
        }

        node_curr_list = g_slist_next (node_curr_list);
    } 
    return 1;
}


int opttree_iteration (opttree_t *self) {

    // 1. Sample a state
    state_t state_random;
    if (optsystem_sample_state (self->optsys, &state_random) == 0)
        return 0;
    
    /*
    //    Target sampling rule
    if (self->lower_bound_node) 
    {
        if (rand()/(RAND_MAX + 1.0) > self->target_sample_prob_after_first_solution) {
            if (optsystem_sample_state (self->optsys, &state_random) == 0)
                return 0;
        }
        else 
        {
            if (optsystem_sample_target_state (self->optsys, &state_random) == 0)
                return 0;
        }
    }
    */

    // A) RRT* ALGORITHM
    if (self->run_rrtstar) 
    {   
        // A.1. Calculate the ball radius constant
        self->ball_radius_last = self->ball_radius_constant 
            * (sqrt(log(1+(double)(self->num_nodes))/((double)(self->num_nodes))));
        if (self->ball_radius_last >= self->ball_radius_max)
            self->ball_radius_last = self->ball_radius_max;

        // A.2. Find the nearest node
        node_t *nearest_node = opttree_find_nearest_neighbor (self, &state_random);
        if (!nearest_node)
            return 0;
        
        // A.3. Extend the node towards the sampled state -- just computer the extended state
        state_t *extended_state = opttree_extend_towards_sample_no_create_node (self, nearest_node, &state_random);
        if (!extended_state)
            return 0;
        
        // A.4. Compute the set of all close nodes around extended_state
        GSList *close_nodes = NULL;
        close_nodes = opttree_find_nodes_in_ball (self, extended_state, self->ball_radius_last);        

        // A.5. Pick the node to be extended
        node_t *node_from = nearest_node;  // If close_nodes is empty,     then extend nearest_nodes
        if (close_nodes) {                 // If close_nodes is non-empty, then extend the min_node in close_nodes
            node_from = opttree_find_min_node_in_set (self, extended_state, close_nodes);
            if (!node_from)                //   If no close node can be extended, then fall back to the nearest 
                node_from = nearest_node;
        }
        
        // A.6. Extend the appropriate node
        node_t *extended_node = opttree_extend_towards_sample (self, node_from, extended_state);
        if (!extended_node) {
            g_slist_free (close_nodes);
            return 0;
        }

        // A.7. Rewire the tree if possible
        opttree_extend_back_to_tree (self, extended_node, close_nodes);
        
        optsystem_free_state (self->optsys, extended_state);
        
        g_slist_free (close_nodes);
    }
    
    // B) RRT ALGORITHM
    else {
        
        // B.1. Find the nearest node
        node_t *nearest_node = opttree_find_nearest_neighbor (self, &state_random);
        if (!nearest_node)
            return 0;
        
        // B.2. Extend the nearest towards the sample
        node_t *extended_node = opttree_extend_towards_sample (self, nearest_node, &state_random);
        if (!extended_node)
            return 0;
    }

    return 1;
} 


int opttree_reinitialize (opttree_t *self) {

    // Clear the tree
    GSList *node_curr_list  = self->list_nodes; 
    while (node_curr_list) {
        node_t *node_curr = node_curr_list->data;
        
        opttree_free_node (self, node_curr);

        node_curr_list = g_slist_next (node_curr_list); 
    }

    g_slist_free (self->list_nodes);
    
    self->list_nodes = NULL;

    // Reinitialize the basics
    self->lower_bound = DBL_MAX;
    self->lower_bound_node = NULL;

    // Reinitialize the kdtree
    kd_clear (self->kdtree);

    // Initialize the root node
    self->root = opttree_new_node (self);
    optsystem_get_initial_state (self->optsys, self->root->state);
    self->root->distance_from_root = 0.0;
    self->root->distance_from_parent = 0.0;
    kd_insert (self->kdtree, optsystem_get_state_key (self->optsys, self->root->state), self->root);

    // Initialize the list of all nodes
    self->list_nodes = NULL;
    self->list_nodes = g_slist_prepend (self->list_nodes, (gpointer) (self->root));
    self->num_nodes = 1;

    return 1;
}


opttree_t* opttree_create () { 

    opttree_t *self = (opttree_t *) calloc (sizeof (opttree_t), 1);

    // Set up the dynamical system
    self->optsys = (optsystem_t *) calloc (sizeof (optsystem_t), 1);
    optsystem_new_system (self->optsys);

    // Initialize the kdtree
    self->kdtree = kd_create (optsystem_get_num_states(self->optsys));

    // Set the lower bound to infinity
    self->lower_bound = DBL_MAX;
    self->lower_bound_node = NULL;

    // Initialize the root node
    self->root = opttree_new_node (self);
    optsystem_get_initial_state (self->optsys, self->root->state);
    self->root->distance_from_root = 0.0;
    self->root->distance_from_parent = 0.0;

    kd_insert (self->kdtree, optsystem_get_state_key (self->optsys, self->root->state), self->root);

    // Initialize the list of all nodes
    self->list_nodes = NULL;
    self->list_nodes = g_slist_prepend (self->list_nodes, (gpointer) (self->root));
    self->num_nodes = 1;

    // Initialize parameters to default values
    self->run_rrtstar = 1;
    self->ball_radius_constant = 30.0;
    self->ball_radius_max = 1.0;
    self->target_sample_prob_after_first_solution = 0.0;
    self->target_sample_prob_before_first_solution = 0.0;
    return self;
}


// Frees the memory allocated for the nodes of the tree 
int opttree_free_tree (opttree_t *self) {

    GSList *list_node_curr = self->list_nodes;
    while (list_node_curr) {
        opttree_free_node (self, (node_t *)(list_node_curr->data) );
        list_node_curr = g_slist_next (list_node_curr);
    }
    g_slist_free (self->list_nodes);

    return 1;
}


int opttree_destroy (opttree_t *self) {

    optsystem_free_system (self->optsys);

    opttree_free_tree (self);

    return 1;
}


int opttree_set_root_state (opttree_t *self, state_t *state) 
{

    optsystem_set_initial_state (self->optsys, state); 
    optsystem_get_initial_state (self->optsys, self->root->state);
    self->root->bowl_radius = 1.0;

    return 1;
}

// propagate state till root, return 0 if collides or cannot find a node within
// a particular bowl
int propagate_to_root(opttree_t *self, state_t *state)
{
    state_t *state_copy = optsystem_clone_state(self->optsys, state);
    create_noisy_state(state_copy, NOISE_LIM);
    //printf("start prop - state_copy: %3.5f %3.5f\n", state_copy->x[0], state_copy->x[1]);
    
    int obstacle = 0;
    state_t *state_new = optsystem_new_state(self->optsys);
    while(1)
    {
        node_t *nearest = opttree_find_nearest_neighbor(self, state_copy);
        state_t *nearest_s = nearest->state;
        //printf("nearest: %3.5f %3.5f\n", nearest_s->x[0], nearest_s->x[1]);

        if(nearest->distance_from_root == 0)
        {
            //printf("found root\n");
            break;
        }
        else
        {
            node_t *parent = nearest->parent;
            state_t *nearest_p_s = parent->state;
            for(int i=0; i<NUM_STATES; i++)
            {
                float del = (nearest_p_s->x[i] - nearest_s->x[i]);
                state_new->x[i] = state_copy->x[i] + del + NOISE;
            }
            //printf("in prop - state_new: %3.5f %3.5f\n", state_new->x[0], state_new->x[1]);

            if(optsystem_on_obstacle(self->optsys, state_new))
            {
                //printf("here1\n");
                obstacle = 1;
                break;
            }
            else if(optsystem_segment_on_obstacle(self->optsys, state_copy, state_new, 10))
            {
                obstacle = 1;
                break;
            }

            for(int i=0; i<NUM_STATES; i++)
            {
                state_copy->x[i] = state_new->x[i];
                //printf("state_copy[%d]: %3.2f\n", i, state_copy->x[i]);
            }
        }
    }
    optsystem_free_state(self->optsys, state_new);
    optsystem_free_state(self->optsys, state_copy);
    if(obstacle)
        return 0;
    //printf("returning 1\n");
    return 1;
}

// propagate state till parent, return 0 if collides or cannot find a node within
// a particular bowl
int propagate_to_parent(opttree_t *self, state_t *towards, state_t *from, double *curr_radius)
{
    state_t *state_new = optsystem_clone_state(self->optsys, towards);
    create_noisy_state(state_new, *curr_radius);
    state_t *state_copy = optsystem_clone_state(self->optsys, state_new);
    
    node_t *nearest = opttree_find_nearest_neighbor(self, from);
    double to_reach_radius = nearest->bowl_radius;
    //printf("radius: %f\n", *curr_radius);

    //printf("start prop - state_copy: %3.5f %3.5f\n", state_copy->x[0], state_copy->x[1]);
    gboolean obstacle = 0, outside_bowl = 0;
    //printf("nearest: %3.5f %3.5f\n", nearest_s->x[0], nearest_s->x[1]);

    for(int i=0; i<NUM_STATES; i++)
    {
        float del = (from->x[i] - state_new->x[i]);
        state_new->x[i] += (del + NOISE);
    }
    //printf("in prop - state_new: %3.5f %3.5f\n", state_new->x[0], state_new->x[1]);
    double dist_to_from = optsystem_evaluate_distance(self->optsys, state_new, from);

    if(optsystem_segment_on_obstacle(self->optsys, state_copy, state_new, 10))
        obstacle = 1;
    else if(dist_to_from > to_reach_radius)
        outside_bowl = 1;

    if(obstacle || outside_bowl)
    {
        //printf("state new: %f %f\n", state_new->x[0], state_new->x[1]);
        double temp = MAX(optsystem_evaluate_distance(self->optsys, state_new, towards) - 0.05, 0);
        //printf("temp: %f\n", temp);
        *curr_radius = MIN(*curr_radius, temp);
    }
    optsystem_free_state(self->optsys, state_new);
    optsystem_free_state(self->optsys, state_copy);
    if(obstacle || outside_bowl)
        return 0;
    //printf("returning 1\n");
    return 1;
}

