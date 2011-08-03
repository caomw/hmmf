#include "common.h"

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
        x1 = 2.0*RANDF - 1.0;
        x2 = 2.0*RANDF - 1.0;
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
