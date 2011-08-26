#ifndef __common_h__
#define __common_h__

#include <iostream>
#include <cstdio>
#include <math.h>
#include <cfloat>
#include <list>
#include <vector>
#include <queue>
#include <fstream>
#include <cstdlib>
#include <sys/time.h>
#include <algorithm>
#include <utility>
#include <cassert>
#include <string>

#include "lp_lib.h"
#include "kdtree.h"
using namespace std;

#define RANDF       (rand()/(RAND_MAX+1.0))


float randn();
void multivar_normal(float *mean, float *var, float *ret, int dim);
float sq(float x);
float normal_val(float *mean, float *var, float *tocalci, int dim);

inline
float sq(float x)
{
    return (x)*(x);
}

inline
float normal_val(float *mean, float *var, float *tocalci, int dim)
{
    float top = 0;
    float det = 1;
    for(int i=0; i<dim; i++)
    {
        top += sq(mean[i] - tocalci[i])/2/var[i];

        det = det*var[i];
    }
    top = exp(-0.5*top);
    float to_ret = 1/pow(2*M_PI, dim/2.0)/ sqrt( det ) * top;
    
    if ( to_ret < 1e-20)
        to_ret = 1e-20;
    
    return to_ret;
}

void tic();
void toc();
double get_msec();

#endif

