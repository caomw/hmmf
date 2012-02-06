#ifndef __common_h__
#define __common_h__

#include <iostream>
#include <cstdio>
#include <cmath>
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
#include <iomanip>

#include <Eigen/Dense>
#include "kdtree.h"

using namespace std;
using namespace Eigen;

#define RANDF       (rand()/(RAND_MAX+1.0))


double randn();
void multivar_normal(double *mean, double *var, double *ret, int dim);
double sq(double x);
double normal_val(double *mean, double *var, double *tocalci, int dim);

// mean = 0, var = 1
inline
double randn()
{
    static double x1, x2, w = 1.5;
    static bool which = 0;

    if(which == 1)
    {
        which = 0;
        return x2*w;
    }

    while (w >= 1.0){
        x1 = 2.0*RANDF - 1.0;
        x2 = 2.0*RANDF - 1.0;
        w = x1 * x1 + x2 * x2;
    }
    w = sqrt( (-2.0*log(w)) / w );
    
    which = 1;
    return x1 * w;
};

inline
void multivar_normal(double *mean, double *var, double *ret, int dim)
{
    for(int i=0; i < dim; i++)
        ret[i] = mean[i] + sqrt(var[i])*randn();
}

inline
double sq(double x)
{
    return (x)*(x);
}

inline
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
    double to_ret = (top/pow(2*M_PI, dim/2.0))/ sqrt( det );
   
    if ( to_ret < 1e-30)
        to_ret = 1e-30;

    return to_ret;
}

void tic();
double toc();
double get_msec();

#endif

