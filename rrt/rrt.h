#ifndef __RRT_H__
#define __RRT_H__

#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <list>
using namespace std;

#define NUM_STATES      (2)

class State
{
    public:
        vector<double> x;
        
        State()
        {
            ;
        };
        State(vector<double> init)
        {
            x = init;
        };
        bool operator==(State s)
        {
            if(s.x ==  x)
                return true;
            else
                return false;
        }
        bool operator!=(State s)
        {
            if(s.x == x)
                return false;
            else
                return true;
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
};

typedef struct Box{
    State center;
    double width;
    double height;
}Box;

#endif
