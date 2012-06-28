#ifndef __rotation__
#define __rotation__

// state is quaternion (w,x,y,z)
// observations are (x,y,z)
#define ndim        (4)
#define ndim_obs    (3)

float init_var[ndim] = {1e-2, 1e-2, 1e-2, 1e-2};
float init_state[ndim] = {1/1.414, 0, 0, 1/1.414};             // pointing upwards
float init_state_real[ndim] = {1/1.414, 0, 0, 1/1.414};
float pvar[ndim] = {1e-4, 1e-4, 1e-4, 1e-4};
float ovar[3] = {1e-4, 1e-4, 1e-4};
float zero[ndim] = {0};
float nomw[3] = {0.1,0.1,0.1};

float norm(float* s)
{
    float sum = 0;
    for(int i=0; i<ndim;i++)
        sum = sum + sq(s[i]);
    return sqrt(sum);
}
int drift(float* q, float *ret, float dt=1.0, bool real=false)
{
    ret[0] = -0.5*(nomw[0]*q[1] + nomw[1]*q[2] + nomw[2]*q[3])*dt;
    ret[1] = 0.5*(nomw[0]*q[0] + nomw[2]*q[2] - nomw[1]*q[3])*dt;
    ret[2] = 0.5*(nomw[1]*q[0] - nomw[2]*q[1] + nomw[0]*q[3])*dt;
    ret[3] = 0.5*(nomw[2]*q[0] + nomw[1]*q[1] + nomw[0]*q[2])*dt;

    for(int i=0;i<ndim; i++)
        ret[i] = 0;

    return 0;
}
int diffusion(float* q, float* ret, float dt=1.0, bool real=false)
{
    float var[ndim] ={0};
    for(int i=0; i<ndim; i++)
        var[i] = pvar[i]*dt;
    
    multivar_normal(zero, var, ret, ndim);
    
    return 0;
}
int get_obs(float* s, float* obs, bool is_clean=false)
{
    float sphi = 1; // sqrt(1 - s[0]*s[0]);
    for(int i=0; i< ndim_obs; i++)
        obs[i] = 0;
    if(is_clean)
    {
        obs[0] = s[1]/sphi;
        obs[1] = s[2]/sphi;
        obs[2] = s[3]/sphi;
        return 0;
    }
    float noise[3] = {0};
    multivar_normal(zero, ovar, noise, 3); 
    
    obs[0] = s[1]/sphi + noise[0];
    obs[1] = s[2]/sphi + noise[1];
    obs[2] = s[3]/sphi + noise[2];
    float len = 1; //sqrt(obs[0]*obs[0] +obs[1]*obs[1] + obs[2]*obs[2]);
    obs[0] = obs[0]/len;
    obs[1] = obs[1]/len;
    obs[2] = obs[2]/len;

    return 0;
}
float holding_time(float* s, float r)
{
    float h = r*2;
    float ret[ndim] ={0};
    drift(s, ret, 1.0, true);
    return h*h/(pvar[0] + h*norm(ret));
}

#endif
