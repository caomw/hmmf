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

#ifndef __OPTSYSTEM_H_
#define __OPTSYSTEM_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <glib.h>


#define NUM_STATES 2       // Number of states
#define NUM_INPUTS 1       // Number of inputs

#define NOISE   (0.3*(rand()/(RAND_MAX + 1.0)))
//#define NOISE   0

// Type definitions 
typedef struct _region_2d_t region_2d_t;
typedef struct _state_t state_t;
typedef struct _input_t input_t;
typedef struct _optsystem_t optsystem_t;

// Initializes the optsystem_t structure
int
optsystem_new_system (optsystem_t *self);

// Frees the memory allocated for elements of the optsystem_t structure
int
optsystem_free_system (optsystem_t *self);

// Allocates memory for a new state and returns a pointer to the new state
state_t*
optsystem_new_state (optsystem_t *self);

// Frees the memory allocated for a state
int
optsystem_free_state (optsystem_t *self, state_t *state);

// Allocates memory for a new input and returns a pointer to the new input
input_t*
optsystem_new_input (optsystem_t *self);

// Frees the memory allocated for an input
int
optsystem_free_input (optsystem_t *self, input_t *input);

// Allocates memory for a new state. Copies a given state into the new state 
//   and returns a pointer to the new state
state_t*
optsystem_clone_state (optsystem_t *self, state_t *state);

// Sets the initial state to a given state_t structure
int 
optsystem_set_initial_state (optsystem_t *self, state_t *state);

// Copies the initial state to the given state_t structure
int 
optsystem_get_initial_state (optsystem_t *self, state_t *state);

// Returns the number of states
int
optsystem_get_num_states (optsystem_t *self);

// Returns a key to the state for the kd-tree implementation. The return 
//   value must be a double array of size equal to the number of states
double*
optsystem_get_state_key (optsystem_t *self, state_t *state);

// Returns a new sample state from the free-space 
int 
optsystem_sample_state (optsystem_t *self, state_t *random_state);

// Returns a new sample state from the target state
int 
optsystem_sample_target_state (optsystem_t *self, state_t *random_state);

// Evaluates the distance between two given states for nearest neighbors computation
double 
optsystem_evaluate_distance (optsystem_t *self, state_t *state_from, state_t *state_to);

// Evaluates the cost of a given trajectory
double 
optsystem_evaluate_distance_for_cost (optsystem_t *self, GSList *inputs);


// Returns true iff the given state reaches the goal
gboolean 
optsystem_is_reaching_target (optsystem_t *self, state_t *state);

// Returns 1 iff (state) is on an obstacle
gboolean optsystem_on_obstacle (optsystem_t *self, state_t *state);

// Returns 1 iff the line connecting (state_initial) and (state_final) lies on an obstacle
int optsystem_segment_on_obstacle (optsystem_t *self, state_t *state_initial, state_t *state_final, int num_steps);




// ===========   Environment interface    =============================================
// Updates the operating region
gboolean 
optsystem_update_operating_region (optsystem_t *self, region_2d_t *operating_region);

// Updates the goal region
gboolean 
optsystem_update_goal_region (optsystem_t *self, region_2d_t *goal_region);

// Updates the list of obstacles
gboolean 
optsystem_update_obstacles (optsystem_t *self, GSList *obstacle_list);

// A 3D region to describe operating and goal regions as well as obstacles
struct _region_2d_t {
    double center[2];
    double size[2];
};
// ====================================================================================


// State structure
struct _state_t {
    double x[NUM_STATES];
};


// Input structure
struct _input_t {
    double x[NUM_INPUTS];
};


// System structure
struct _optsystem_t {

    state_t *initial_state;     // Initial state
    
    region_2d_t operating_region;  // Operating region

    region_2d_t goal_region;       // Goal region

    GSList *obstacle_list;         // A list of obstacles
};


#endif
