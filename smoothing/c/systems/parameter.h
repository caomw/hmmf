#ifndef __parameter__
#define __parameter__

#define ndim (2)
#define ndim_obs (1)

double zmin[ndim] = {-0.2, 0};
double zmax[ndim] = {1, 1.5};
double init_var[ndim] = {1e-2, 1e-1};
double init_state[ndim] = {0.0, 0.9};
double init_state_real[ndim] = {0.0, 0.5};
double pvar[ndim] = {1e-3, 1e-4};
double ovar[ndim] = {1e-3};
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
    double phi = s[1];
    if(real)
        phi = 0.5;
    ret[0] = cos(2*M_PI*phi*s[0])*dt;
    ret[1] = 0*dt;
    return 0;
}
int diffusion(double* s, double* ret, double dt=1.0, bool real=false)
{
    double var[ndim] ={0};
    var[0] = pvar[0]*dt;
    var[1] = pvar[1]*dt;
    multivar_normal(zero, var, ret, ndim);
    //cout<<ret[0] <<" "<< ret[1] <<endl;
    if(real)
        ret[1] = 0;
    return 0;
}
int get_obs(double* s, double* obs, bool is_clean = false)
{
    for(int i=0; i< ndim; i++)
        obs[i] = 0;
    if(is_clean)
    {
        obs[0] = s[0];
        return 0;
    }
    double noise[ndim_obs] = {0};
    multivar_normal(zero, ovar, noise, ndim_obs);
    obs[0] = s[0] + noise[0];
    return 0;
}
double holding_time(double* s, double r)
{
    double h = r*(zmax[0] - zmin[0]);
    double ret[ndim] = {0};
    drift(s, ret, 1.0, true);
    return h*h/(pvar[0] + h*norm(ret));
}

#endif
