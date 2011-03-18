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
#include <time.h>

#include "halton.h"
using namespace std;

#define XMAX        (2.0)
#define XMIN        (-2.0)
#define NUM_DIM     (1)

#define randu       (XMIN + rand()/(RAND_MAX + 1.0)*(XMAX - XMIN))
#define randf       (rand()/(RAND_MAX + 1.0))

class edge;
class vertex;
class graph;

typedef struct state_t{
    double x[NUM_DIM];
}state;
typedef struct mini_sample{
    state s;
    vertex *parent;
}minis;

class vertex{

    public:
        state s;
        double t;
        vertex *prev;
        double ravg;
        int num_child;

        vector<edge *> edgein;
        vector<edge *> edgeout;
        
        vertex(state st, double tt){
            for(int i=0; i<NUM_DIM; i++)
                s.x[i] = st.x[i];
            t = tt;
            ravg = 0;
            num_child = 0;
        }
        ~vertex(){};
};

class edge{

    public:
        vertex *from;
        vertex *to;
        double prob;

        edge(vertex *f, vertex *t, double p){
            this->from = f;
            this->to = t;
            this->prob = p;
        }
        edge reverse(){
            return edge(this->to, this->from, this->prob);
        }
        ~edge(){};
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
            {
                vertex *v = *i;
                for(vector<edge *>::iterator j = v->edgein.begin(); j != v->edgein.end(); j++)
                {
                    edge *e = *j;
                    //cout<<"deleting e "<<e<<endl;
                    delete e;
                }
            }
            for(vector<vertex*>::iterator i = vlist.begin(); i != vlist.end(); i++)
            {
                vertex *v = *i;
                //cout<<"deleting v "<<v<<endl;
                delete v;
            }
        };
 };


state sample(){
    state s;
    for(int i=0; i<NUM_DIM; i++)
        s.x[i] = randu;
    return s;
}

state sample_quasi(){
    state s;
    for(int i=0; i<NUM_DIM; i++)
    {
        double r;
        halton(&r);
        s.x[i] = XMIN + (XMAX-XMIN)*r;
    }
    return s;
}

double get_msec(){
    struct timeval start;
    gettimeofday(&start, NULL);
    return start.tv_sec*1000 + start.tv_usec/1000.0;
}

double randn(float mean,float var){
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
