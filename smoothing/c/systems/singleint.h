#ifndef __singleint__
#define __singleint__

#define ndim        (2)
#define ndim_obs    (2)

#include "sysutil.h"

double zmin[ndim] = {-1, -1};
double zmax[ndim] = {1, 1};
double zdiff[ndim] = {2, 2};
double init_var[ndim] = {1e-2, 1e-2};
double init_state[ndim] = {0.8, 0.8};
double init_state_real[ndim] = {0.8, 0.8};
double pvar[ndim] = {1e-2, 1e-2};
double ovar[ndim] = {1e-3, 1e-3};
double zero[ndim] = {0, 0};


int drift(double* s, double *ret, double dt=1.0, bool real=false)
{
    ret[0] = -s[0]*dt;
    ret[1] = -s[1]*dt;
    return 0;
}
int diffusion(double* s, double* ret, double dt=1.0, bool real=false)
{
    double var[ndim] ={0};
    var[0] = pvar[0]*dt;
    var[1] = pvar[1]*dt;
    multivar_normal(zero, var, ret, ndim);
    return 0;
}
int get_obs(double* s, double* obs, bool is_clean=false)
{
    for(int i=0; i< ndim; i++)
        obs[i] = 0;
    if(is_clean)
    {
        for(int i=0;i<ndim_obs; i++)
            obs[i] = s[i];
        return 0;
    }
    double noise[ndim] = {0};
    multivar_normal(zero, ovar, noise, ndim_obs); 
    obs[0] = s[0] + noise[0];
    obs[1] = s[1] + noise[1];
    return 0;
}
double holding_time(double* s, double r)
{
    double h = r*(zmax[0] - zmin[0]);
    double ret[ndim] ={0};
    drift(s, ret, 1.0, true);
    return h*h/(pvar[0] + h*norm(ret));
}

int integrate_system(double* curr_state, double dt, bool is_clean=false)
{
    double integration_delta = min(1e-3, dt/2.0);
    double runner_time = 0;
    while(runner_time < dt)
    {
        double next_state_delta[ndim] ={0};
        drift(curr_state, next_state_delta, integration_delta, true);
        double noise[ndim] = {0};
        if(!is_clean)
        {
            diffusion(curr_state, noise, integration_delta, true);
            add(noise, next_state_delta);
        }
        add(next_state_delta, curr_state);
        runner_time = runner_time + integration_delta;
    }
    return 0;
}

#endif
