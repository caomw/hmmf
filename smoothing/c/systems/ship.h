#ifndef __ship__
#define __ship__

#define ndim        (4)
#define ndim_obs    (2)

double zmin[ndim] = {-5, -5, -2, -2};
double zmax[ndim] = {5, 5, 2, 2};
double zdiff[ndim] = {10, 10, 4, 4};
double init_var[ndim] = {1e-2, 1e-2, 1e-2, 1e-2};
double init_state[ndim] = {-2, -2, 1, 1};
double init_state_real[ndim] = {-2, -2, 1, 1};
double pvar[ndim] = {1e-1, 1e-1, 1e-1, 1e-1};
double ovar[ndim_obs] = {1e-1, 1e-1};
double zero[ndim] = {0};

double norm(double* s)
{
    double sum = 0;
    for(int i=0; i<ndim;i++)
        sum = sum + sq(s[i]);
    return sqrt(sum);
}
int drift(double* s, double *ret, double dt=1.0, bool real=false)
{
    double f1dt=0, f2dt = 0;
    double dist = sqrt(s[0]*s[0] + s[1]*s[1]);
    if( ( dist >= 9) && ( (s[0]*s[2] + s[1]*s[3]) >= 0) )
    {
        f1dt = -50*dt*s[0]/dist;
        f2dt = -50*dt*s[1]/dist;
    }
    ret[0] = s[2]*dt; 
    ret[1] = s[3]*dt;
    ret[2] = f1dt;
    ret[3] = f2dt;

    return 0;
}
int diffusion(double* s, double* ret, double dt=1.0, bool real=false)
{
    double var[ndim] ={0};
    for(int i=0; i<ndim; i++)
        var[i] = pvar[i]*dt;
    multivar_normal(zero, var, ret, ndim);
    return 0;
}
int get_obs(double* s, double* obs, bool is_clean=false)
{
    double range = sqrt(s[0]*s[0] + s[1]*s[1]);
    double theta = atan2(s[1], s[0]);
    for(int i=0; i< ndim_obs; i++)
        obs[i] = 0;
    if(is_clean)
    {
        obs[0] = range;
        obs[1] = theta;
        return 0;
    }
    double noise[ndim_obs] = {0};
    multivar_normal(zero, ovar, noise, ndim_obs); 
    obs[0] = range + noise[0];
    obs[1] = theta + noise[1];
    return 0;
}
double holding_time(double* s, double r)
{
    double h = r*(zmax[0] - zmin[0]);
    double ret[ndim] ={0};
    drift(s, ret, 1.0, true);
    return h*h/(pvar[0] + h*norm(ret));
}

#endif
