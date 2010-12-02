#ifndef __RRT_H__
#define __RRT_H__

#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <list>
using namespace std;

#define NUM_STATES      (2)

class State
{
    public:
        double x[NUM_STATES];
        State()
        {
            ;
        };
        State(double *s)
        {
            for(int i=0; i<NUM_STATES; i++)
            {
                x[i] = (*s);
                s++;                // inc the pointer
            }
        };

        bool operator==(State s)
        {
            for(int i=0; i<NUM_STATES; i++)
            {
                if(s.x[i] != x[i])
                    return false;
            }
            return true;
        }
        bool operator!=(State s)
        {
            for(int i=0; i<NUM_STATES; i++)
            {
                if(s.x[i] != x[i])
                    return true;
            }
            return false;
        }

        void print()
        {
            for(unsigned int i=0; i<NUM_STATES; i++)
            {
                printf("%f ", x[i]);
            }
            cout<<endl;
        }
};

class Node
{
    public:
        State state;
        Node *parent;
        double csrc;        // cost from src
        double cgoal;       // heuristic
        double cparent;     // cost from parent

        Node(State s, Node *p, double c1=0, double c2=0, double c3=0)
        {
            state = s;
            parent = p;
            csrc = c1;
            cgoal = c2;
            cparent = c3;
        };
        Node()
        {
        };
};

typedef struct Box{
    State center;
    double size[NUM_STATES];
}Box;


typedef struct Goal{
    State state;
    double size;
}Goal;

typedef struct kdtree kdtree;
typedef struct kdres kdres;

// extern variables

extern list<Node> tree;        // stores the tree
extern list<Node> path;        // stores path that reaches goal
extern list<Node> optpath;     // stores softened path

extern kdtree *obstree;
extern Goal goal;
extern State robot;
extern Box box;
extern double robot_radius;
extern double MAX_OBS_SIZE;
extern double NUM_OBSTACLES;
extern double obs_rad[100];

// prototypes

double dist(State s1, State s2); 
bool is_inside_goal(State c);
bool is_obstructed(State s);
State sample_state();
Node* nearest(State s);
Node* nearest_kdtree(kdtree *node_tree, State s);
Node extend_rrt(Node *near, State s);
bool does_line_hit(State s1, State s2, double *pos, double rad);
bool can_join_nodes(Node n1, Node n2);
void process_tree_rrt(Node goal_node);
double rrt_plan(unsigned int num_iter);


int extend_rrtstar(kdtree *node_tree, Node *near, State s, Node &returned_node, double curr_min_cost);
double rrtstar_plan(unsigned int num_iter);
void process_tree_rrtstar(Node *goal_node);

#endif
