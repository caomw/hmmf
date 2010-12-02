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

#endif
