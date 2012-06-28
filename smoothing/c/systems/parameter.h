#ifndef __parameter__
#define __parameter__

#define ndim (2)
#define ndim_obs (1)

#include "sysutil.h"

float zmin[ndim] = {-1, -1};
float zmax[ndim] = {1, 1};
float init_var[ndim] = {1e-2, 1e-1};
float init_state[ndim] = {0.0, 0.8};
float init_state_real[ndim] = {0.0, 0.5};
float pvar[ndim] = {1e-2, 1e-2};
float ovar[ndim] = {1e-2};
float zero[ndim] = {0};

int drift(float* s, float *ret, float dt=1.0, bool real=false)
{
    float phi = s[1];
    if(real)
        phi = 0.5;
    ret[0] = s[0]*cos(2*M_PI*phi*s[0])*dt;
    ret[1] = 0*dt;
    return 0;
}
int diffusion(float* s, float* ret, float dt=1.0, bool real=false)
{
    float var[ndim] ={0};
    var[0] = pvar[0]*dt;
    var[1] = pvar[1]*dt;
    multivar_normal(zero, var, ret, ndim);
    //cout<<ret[0] <<" "<< ret[1] <<endl;
    if(real)
        ret[1] = 0;
    return 0;
}
int get_obs(float* s, float* obs, bool is_clean = false)
{
    for(int i=0; i< ndim; i++)
        obs[i] = 0;
    if(is_clean)
    {
        obs[0] = s[0];
        return 0;
    }
    float noise[ndim_obs] = {0};
    multivar_normal(zero, ovar, noise, ndim_obs);
    obs[0] = s[0] + noise[0];
    return 0;
}
float holding_time(float* s, float r)
{
    float h = r*(zmax[1] - zmin[1]);
    float ret[ndim] = {0};
    drift(s, ret, 1.0, true);
    return h*h/(max_norm(pvar) + h*norm(ret));
}

int integrate_system(float* curr_state, float dt, bool is_clean=false)
{
    float integration_delta = min(1e-3, dt/2.0);
    float runner_time = 0;
    while(runner_time < dt)
    {
        float next_state_delta[ndim] ={0};
        drift(curr_state, next_state_delta, integration_delta, true);
        float noise[ndim] = {0};
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

int get_kalman_path(vector< vector<float> >& kfp, vector< vector<float> >& observations, float dt)
{
    kfp.clear();
    kfp = vector< vector<float> >(observations.size(), vector<float>(ndim,0));
    return 0;
}
#endif
