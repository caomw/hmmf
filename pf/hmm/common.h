#ifndef __common_h__
#define __common_h__

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <list>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <sys/time.h>
using namespace std;

#define XMAX        (5.0)
#define XMIN        (-5.0)
#define NUM_DIM     (1)

#define randdim     (XMIN + rand()/(RAND_MAX + 1.0)*(XMAX - XMIN))
#define randf       (rand()/(RAND_MAX + 1.0))

typedef struct state_t{
    float x[NUM_DIM];
}state;

class edge;
class vertex;
class graph;

class vertex{

    public:
        state s;

        vector<edge *> edgein;
        vector<edge *> edgeout;

        vertex() {};
        vertex(state st){
            for(int i=0; i<NUM_DIM; i++)
                s.x[i] = st.x[i];
        }
};

class edge{

    public:
        vertex *from;
        vertex *to;
        float prob;

        edge(vertex *f, vertex *t, float p){
            this->from = f;
            this->to = t;
            this->prob = p;
        }
        edge reverse(){
            return edge(this->to, this->from, this->prob);
        }
};

class graph{
    public:
        vector<vertex *> vlist;

        graph() {};
        void add_vertex(vertex *v){
            vlist.push_back(v);
        }
        ~graph(){
            for(vector<vertex*>::iterator i = vlist.begin(); i != vlist.end(); i++)
                delete *i;
            vlist.clear();
        };
 };

state sample(){
    state s;
    for(int i=0; i<NUM_DIM; i++)
        s.x[i] = randdim;
    return s;
}

float get_msec(){
    struct timeval start;
    gettimeofday(&start, NULL);
    return start.tv_sec*1000 + start.tv_usec/1000.0;
}

float randn(float mean,float var){
    float x1, x2, w = 1.5, y1;
    while (w >= 1.0){
        x1 = 2.0*randf - 1.0;
        x2 = 2.0*randf - 1.0;
        w = x1 * x1 + x2 * x2;
    }
    w = sqrt( (-2.0*log(w)) / w );
    y1 = x1 * w;
    
    return mean + y1*var;
};
#endif
