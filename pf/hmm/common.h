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
#include "Logger.h"

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
        bool is_open;

        // parent of the best path
        vertex *prev;
        vertex *next;

        edge *best_in;
        edge *best_out;

        vector<edge *> edgein;
        vector<edge *> edgeout;
        
        vertex(state st){
            for(int i=0; i<NUM_DIM; i++)
                s.x[i] = st.x[i];

            prev = NULL;
            next = NULL;
            prob = MIN_PROB;
            is_open = 0;
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
            prob = 1;
            delt = time;
        }
        edge reverse(){
            return edge(this->to, this->from, this->prob);
        }
        ~edge()
        {
            vector<edge *>::iterator edge_iter1;
            edge_iter1 = find( from->edgeout.begin(), from->edgeout.end(), this);
            from->edgeout.erase( edge_iter1 );

            vector<edge *>::iterator edge_iter2;
            edge_iter2 = find( to->edgein.begin(), to->edgein.end(), this);
            to->edgein.erase( edge_iter2 );
        };
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
        void remove_vertex(vertex *v)
        {
            // remove from vlist
            vector<vertex *>::iterator iter;
            iter = find(vlist.begin(), vlist.end(), v);
            vlist.erase(iter);
            
            // delete all edges from v and its neighbors
            for(vector<edge *>::iterator i = v->edgein.begin(); i != v->edgein.end(); i++)
            {
                edge *etmp = *i;
                /*                
                vector<edge *>::iterator edge_iter;
                edge_iter = find( etmp->from->edgeout.begin(), etmp->from->edgeout.end(), etmp);
                etmp->from->edgeout.erase( edge_iter );
                */
                delete etmp;
            }
            for(vector<edge *>::iterator i = v->edgeout.begin(); i != v->edgeout.end(); i++)
            {
                edge *etmp = *i;
                /*
                vector<edge *>::iterator edge_iter;
                edge_iter = find( etmp->to->edgein.begin(), etmp->to->edgein.end(), etmp);
                etmp->to->edgein.erase( edge_iter );
                */
                delete etmp;
            }
            
            // finally delete v
            delete v;
        }

        ~graph(){
            
            /*
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
            */
        };
 };

double dist(state s1, state s2)
{
    double t = 0;
    for(int i=0; i<NUM_DIM; i++)
        t += (s1.x[i] - s2.x[i])*(s1.x[i] - s2.x[i]);

    return sqrt(t);
};

double get_msec(){
    struct timeval start;
    gettimeofday(&start, NULL);
    return start.tv_sec*1000 + start.tv_usec/1000.0;
}

// mean = 0, var = 1
double randn()
{
    static float x1, x2, w = 1.5;
    static bool which = 0;

    if(which == 1)
    {
        which = 0;
        return x2*w;
    }

    while (w >= 1.0){
        x1 = 2.0*randf - 1.0;
        x2 = 2.0*randf - 1.0;
        w = x1 * x1 + x2 * x2;
    }
    w = sqrt( (-2.0*log(w)) / w );
    
    which = 1;
    return x1 * w;
};

void multivar_normal(double *mean, double *var, double *ret, int dim)
{
    for(int i=0; i < dim; i++)
        ret[i] = mean[i] + sqrt(var[i])*randn();
}

inline double sq(double x)
{
    return (x)*(x);
}

double normal_val(double *mean, double *var, double *tocalci, int dim)
{
    double top = 0;
    double det = 1;
    for(int i=0; i<dim; i++)
    {
        top += sq(mean[i] - tocalci[i])/2/var[i];

        det = det*var[i];
    }
    top = exp(-0.5*top);
    double bottom = 1/pow(2*PI, dim/2.0)/ sqrt( det );
    
    return bottom*top;
}

/*
// from wiki -- rejection algorithm
// samples from f(x) = N(-a, var) + N(a, var)
double rand_two_n( double a, double var)
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
*/

#endif
