#ifndef __common_h__
#define __common_h__

#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string.h>
#include <fstream>
#include <cstdlib>
#include <sys/time.h>
#include <string>
#include <unordered_map>
#include <set>

//#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
//#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "kdtree.h"

using namespace std;
using namespace Eigen;

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
    float to_ret = (top/pow(2*M_PI, dim/2.0))/ sqrt( det );
   
    if ( to_ret < 1e-10)
        to_ret = 1e-10;

    return to_ret;
}

void tic();
double toc();
double get_msec();

#endif

