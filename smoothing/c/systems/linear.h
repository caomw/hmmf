#ifndef __uhlmann__
#define __uhlmann__
/*
 * x1d = x2
 * x2d = beta - alpha x2^2
 * alphad = noise
 * y = x1
 */

#define ndim        (3)
#define ndim_obs    (1)

float zmin[ndim] = {-1, -2, 0};
float zmax[ndim] = {1, 2, 1};
float init_var[ndim] = {1e-2, 1e-2, 1e-1};
float init_state[ndim] = {0, 0, 0.8};
float init_state_real[ndim] = {0, 0, 0.5};
float pvar[ndim] = {1e-4, 1e-2, 1e-3};
float ovar[ndim_obs] = {1e-3};
float zero[ndim] = {0, 0, 0};

float norm(float* s)
{
    float sum = 0;
    for(int i=0; i<ndim;i++)
        sum = sum + sq(s[i]);
    return sqrt(sum);
}
int drift(float* s, float *ret, float dt=1.0, bool real=false)
{
    float alpha = s[2];
    if(real)
        alpha = 0.5;
    ret[0] = s[0]*dt;
    ret[1] = (beta - alpha*sq(s[1]))*dt;
    ret[2] = 0;
    return 0;
}
int diffusion(float* s, float* ret, float dt=1.0, bool real=false)
{
    float var[ndim] ={0};
    var[0] = pvar[0]*dt;
    var[1] = pvar[1]*dt;
    var[2] = pvar[2]*dt;
    multivar_normal(zero, var, ret, ndim);
    if(real)
        ret[2] = 0;
    return 0;
}
int get_obs(float* s, float* obs)
{
    for(int i=0; i< ndim_obs; i++)
        obs[i] = 0;
    float noise[ndim_obs] = {0};
    multivar_normal(zero, ovar, noise, ndim_obs); 
    obs[0] = s[0] + noise[0];
    return 0;
}
float holding_time(float* s, float r)
{
    float h = r*(zmax[0] - zmin[0]);
    float ret[ndim] ={0};
    drift(s, ret, 1.0, true);
    return h*h/(pvar[0] + h*norm(ret));
}

#endif
