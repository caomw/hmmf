#ifndef __parameter_hard__
#define __parameter_hard__

#define ndim (3)
#define ndim_obs (1)

float zmin[ndim] = {-0.2, 0.0, 0.5};
float zmax[ndim] = {0.8, 1.0, 1.5};
float init_var[ndim] = {1e-2, 1e-1, 1e-1};
float init_state[ndim] = {0.0, 0.9, 0.7};
float init_state_real[ndim] = {0.0, 0.5, 1.0};
float pvar[ndim] = {1e-2, 1e-3, 1e-3};
float ovar[ndim] = {1e-3};
float zero[ndim] = {0};

float norm(float* s)
{
    float sum = 0;
    for(int i=0; i<ndim;i++)
        sum = sum + sq(s[i]);
    return sqrt(sum);
}
int drift(float* s, float *ret, float dt=1.0, bool real=false)
{
    float phi = s[1];
    if(real)
        phi = 0.5;
    ret[0] = cos(2*M_PI*phi*s[0])*dt;
    ret[1] = 0*dt;
    ret[2] = 0*dt;
    return 0;
}
int diffusion(float* s, float* ret, float dt=1.0, bool real=false)
{
    float mult = s[2];
    if(real)
        mult = 1.0;

    float var[ndim] ={0};
    var[0] = sq(mult)*pvar[0]*dt;
    var[1] = pvar[1]*dt;
    var[2] = pvar[2]*dt;
    multivar_normal(zero, var, ret, ndim);
    if(real)
    {
        ret[1] = 0;
        ret[2] = 0;
    }
    return 0;
}
int get_obs(float* s, float* obs)
{
    for(int i=0; i< ndim; i++)
        obs[i] = 0;
    float noise[ndim] = {0};
    multivar_normal(zero, ovar, noise, ndim_obs); 
    obs[0] = s[0] + noise[0];
    obs[1] = 0;
    obs[2] = 0;
    return 0;
}
float holding_time(float* s, float r)
{
    float h = r*(zmax[0] - zmin[0]);
    float ret[ndim] ={0};
    drift(s, ret, 1.0, true);
    return h*h/( sq(s[2])*pvar[0] + h*norm(ret));
}

#endif
