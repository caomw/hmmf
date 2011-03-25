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
#include <algorithm>
#include <utility>

#include "halton.h"
using namespace std;

#define NUM_DIM     (2)

#define randux       (XMIN + rand()/(RAND_MAX + 1.0)*(XMAX - XMIN))
#define randuy       (YMIN + rand()/(RAND_MAX + 1.0)*(YMAX - YMIN))
#define randf       (rand()/(RAND_MAX + 1.0))

#define SQ(x)       ((x)*(x))

class edge;
class vertex;
class graph;

typedef struct state_t{
    double x[NUM_DIM];
}state;

class minis{
    public:
        state s;
        vertex *parent;
    minis(){
        parent = NULL;
    };
    ~minis(){};
};

class vertex{

    public:
        state s;
        // prob of best path that ends up here incorporating obs
        vector<double> prob;
        // alphas
        vector<double> alpha;
        // time till which obs are incorporated
        vector<double> t;
        // parent of the best path
        vector<vertex *> prev;
        double voronoi_area;
        int num_child;
        
        vector<edge *> edgein;
        vector<edge *> edgeout;
        
        vertex(state st){
            for(int i=0; i<NUM_DIM; i++)
                s.x[i] = st.x[i];
            voronoi_area = 0;
            num_child = 0;
        }
        ~vertex(){};
};

class edge{

    public:
        vertex *from;
        vertex *to;
        double prob;
        double delt;

        edge(vertex *f, vertex *t, double p){
            from = f;
            to = t;
            prob = p;
            delt = dt;
        }
        edge reverse(){
            return edge(this->to, this->from, this->prob);
        }
        ~edge(){};
};

class graph{
    public:
        vector<vertex *> vlist;
        int num_vert;

        graph() {
            num_vert = 0;
        };
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
    s.x[0] = randux;
    s.x[1] = randuy;
    return s;
}

/*
state sample_quasi(){
    state s;
    double r[NUM_DIM];
    halton(r);
    for(int i=0; i<NUM_DIM; i++)
    {
        s.x[i] = XMIN + (XMAX-XMIN)*r[i];
    }
    return s;
}
*/

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

double dist(state s1, state s2)
{
    double t = 0;
    for(int i=0; i<NUM_DIM; i++)
        t += SQ(s1.x[i] - s2.x[i]);

    return sqrt(t);
};


// prototypes

#endif
