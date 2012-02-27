#ifndef __rotation__
#define __rotation__

// state is quaternion (w,x,y,z)
// observations are (x,y,z)
#define ndim        (4)
#define ndim_obs    (3)

double init_var[ndim] = {1e-2, 1e-2, 1e-2, 1e-2};
double init_state[ndim] = {1/1.414, 0, 0, 1/1.414};             // pointing upwards
double init_state_real[ndim] = {1/1.414, 0, 0, 1/1.414};
double pvar[ndim] = {1e-4, 1e-4, 1e-4, 1e-4};
double ovar[3] = {1e-4, 1e-4, 1e-4};
double zero[ndim] = {0};
double nomw[3] = {0.1,0.1,0.1};

double norm(double* s)
{
    double sum = 0;
    for(int i=0; i<ndim;i++)
        sum = sum + sq(s[i]);
    return sqrt(sum);
}
int drift(double* q, double *ret, double dt=1.0, bool real=false)
{
    ret[0] = -0.5*(nomw[0]*q[1] + nomw[1]*q[2] + nomw[2]*q[3])*dt;
    ret[1] = 0.5*(nomw[0]*q[0] + nomw[2]*q[2] - nomw[1]*q[3])*dt;
    ret[2] = 0.5*(nomw[1]*q[0] - nomw[2]*q[1] + nomw[0]*q[3])*dt;
    ret[3] = 0.5*(nomw[2]*q[0] + nomw[1]*q[1] + nomw[0]*q[2])*dt;

    for(int i=0;i<ndim; i++)
        ret[i] = 0;

    return 0;
}
int diffusion(double* q, double* ret, double dt=1.0, bool real=false)
{
    double var[ndim] ={0};
    for(int i=0; i<ndim; i++)
        var[i] = pvar[i]*dt;
    
    multivar_normal(zero, var, ret, ndim);
    
    return 0;
}
int get_obs(double* s, double* obs, bool is_clean=false)
{
    double sphi = 1; // sqrt(1 - s[0]*s[0]);
    for(int i=0; i< ndim_obs; i++)
        obs[i] = 0;
    if(is_clean)
    {
        obs[0] = s[1]/sphi;
        obs[1] = s[2]/sphi;
        obs[2] = s[3]/sphi;
        return 0;
    }
    double noise[3] = {0};
    multivar_normal(zero, ovar, noise, 3); 
    
    obs[0] = s[1]/sphi + noise[0];
    obs[1] = s[2]/sphi + noise[1];
    obs[2] = s[3]/sphi + noise[2];
    double len = 1; //sqrt(obs[0]*obs[0] +obs[1]*obs[1] + obs[2]*obs[2]);
    obs[0] = obs[0]/len;
    obs[1] = obs[1]/len;
    obs[2] = obs[2]/len;

    return 0;
}
double holding_time(double* s, double r)
{
    double h = r*2;
    double ret[ndim] ={0};
    drift(s, ret, 1.0, true);
    return h*h/(pvar[0] + h*norm(ret));
}

#endif
