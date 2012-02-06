#ifndef __singleint__
#define __singleint__

#define ndim        (2)
#define ndim_obs    (2)

float zmin[ndim] = {-1, -1};
float zmax[ndim] = {1, 1};
float init_var[ndim] = {1e-2, 1e-2};
float init_state[ndim] = {0.8, 0.8};
float init_state_real[ndim] = {0.8, 0.8};
float pvar[ndim] = {1e-2, 1e-2};
float ovar[ndim] = {1e-3, 1e-3};
float zero[ndim] = {0, 0};

float norm(float* s)
{
    float sum = 0;
    for(int i=0; i<ndim;i++)
        sum = sum + sq(s[i]);
    return sqrt(sum);
}
int drift(float* s, float *ret, float dt=1.0, bool real=false)
{
    ret[0] = -s[0]*dt;
    ret[1] = -s[1]*dt;
    return 0;
}
int diffusion(float* s, float* ret, float dt=1.0, bool real=false)
{
    float var[ndim] ={0};
    var[0] = pvar[0]*dt;
    var[1] = pvar[1]*dt;
    multivar_normal(zero, var, ret, ndim);
    return 0;
}
int get_obs(float* s, float* obs)
{
    for(int i=0; i< ndim; i++)
        obs[i] = 0;
    float noise[ndim] = {0};
    multivar_normal(zero, ovar, noise, ndim_obs); 
    obs[0] = s[0] + noise[0];
    obs[1] = s[1] + noise[1];
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
