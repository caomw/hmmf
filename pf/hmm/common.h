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
        double prob;
        // alphas
        vector<double> alpha;
        
        // parent of the best path
        vertex * prev;
        vertex *next;

        edge *best_in;
        edge *best_out;

        double voronoi_area;
        unsigned int num_child;
        
        vector<edge *> edgein;
        vector<edge *> edgeout;
        
        vertex(state st){
            for(int i=0; i<NUM_DIM; i++)
                s.x[i] = st.x[i];
            voronoi_area = 0;
            num_child = 0;

            prev = NULL;
            next = NULL;
        }
        ~vertex(){};
};

class edge{

    public:
        vertex *from;
        vertex *to;
        double prob;
        double delt;

        edge(vertex *f, vertex *t, double time){
            from = f;
            to = t;
            prob = 0;
            delt = time;
        }
        edge reverse(){
            return edge(this->to, this->from, this->prob);
        }
        ~edge(){};
};

class graph{
    public:
        vector<vertex *> vlist;
        unsigned int num_vert;

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

double dist(state s1, state s2)
{
    double t = 0;
    for(int i=0; i<NUM_DIM; i++)
        t += (s1.x[i] - s2.x[i])*(s1.x[i] - s2.x[i]);

    return sqrt(t);
};

/*
state sample(state around_which, double radius){
    
double xmin_tmp = around_which.x[0] - radius;
double ymin_tmp = around_which.x[1] - radius;
#define randux       (xmin_tmp + rand()/(RAND_MAX + 1.0)*(2*radius))
#define randuy       (ymin_tmp + rand()/(RAND_MAX + 1.0)*(2*radius))

    state s;
    s.x[0] = randux;
    s.x[1] = randuy;
#undef randux
#undef randuy
    return s;
}
*/

state sample()
{
#define randut      (TMIN + rand()/(RAND_MAX + 1.0)*(TMAX-TMIN))
#define randux      (XMIN + rand()/(RAND_MAX + 1.0)*(XMAX-XMIN))

    state s;
    s.x[0] = randut;
    s.x[1] = randux;
#undef randut
#undef randux
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

void randn(float mean,float var, double &y1, double &y2){
    float x1, x2, w = 1.5;
    while (w >= 1.0){
        x1 = 2.0*randf - 1.0;
        x2 = 2.0*randf - 1.0;
        w = x1 * x1 + x2 * x2;
    }
    w = sqrt( (-2.0*log(w)) / w );
    y1 = x1 * w;
    y2 = x2 * w;

    y1 =  mean + y1*sqrt(var);
    y2 =  mean + y2*sqrt(var);
};

inline double sq(double x)
{
    return (x)*(x);
}

double normal_val(double mean, double var, double tocalci)
{
    double temp = exp(-0.5*sq(mean-tocalci)/var);
    return 1/sqrt(2*PI*var)*temp;
}

// from wiki -- rejection algorithm
// samples from f(x) = N(-a, var) + N(a, var)
double rand_two_n(float a, float var)
{
    double M = 1e3;
    double s = sqrt(var);
    bool got_one = 0;
    double ret = 10;
    double uniform = 2*a + 10*s;

    while (got_one == 0)
    {
        // get g(x)
        ret = randf*uniform - uniform/2;
        double gx = ret/uniform;

        double u = randf;

        double fx = normal_val(-a, var, ret) + normal_val(a, var, ret);
        if(u < fx/gx/M)
        {
            got_one = 1;
        }
    }
    return ret;
}

double rand_uniform(float a, float spread)
{
    return randf*(2*spread) -spread + a;
}

#endif
