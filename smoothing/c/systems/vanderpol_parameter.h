#ifndef __vanderpol_parameter__
#define __vanderpol_paramter__

#define ndim (3)

double zmin[ndim] = {-2.5, -4.0, 1.5};
double zmax[ndim] = {2.0, 2.0, 2.5};
double init_var[ndim] = {1e-2, 1e-2, 1e-2};
double init_state[ndim] = {0, 1.0, 2.2};
double pvar[ndim] = {1e-2, 1e-2, 1e-4};
double ovar[ndim] = {1e-1, 1e-1};
double zero[ndim] = {0, 0, 0};

double norm(double* s)
{
    double sum = 0;
    for(int i=0; i<ndim;i++)
        sum = sum + sq(s[i]);
    return sqrt(sum);
}
int drift(double* s, double *ret, double dt=1.0, bool real=false)
{
    double mu = s[2];
    if(real)
        mu = 2.0;
    ret[0] = s[1]*dt;
    ret[1] = (-s[0] + mu*s[1]*(1-sq(s[0])))*dt;
    ret[2] = 0*dt;
    return 0;
}
int diffusion(double* s, double* ret, double dt=1.0, bool real=false)
{
    double var[ndim] ={0};
    var[0] = pvar[0]*dt;
    var[1] = pvar[1]*dt;
    var[2] = pvar[2]*dt;
    multivar_normal(zero, var, ret, ndim);
    if(real)
        ret[2] = 0;
    return 0;
}
int get_obs(double* s, double* obs)
{
    for(int i=0; i< ndim; i++)
        obs[i] = 0;
    double noise[ndim] = {0};
    multivar_normal(zero, ovar, noise, 2); 
    obs[0] = s[0] + noise[0];
    obs[1] = s[1] + noise[1];
    return 0;
}
double holding_time(double* s, double r)
{
    double h = r*(zmax[1] - zmin[0]);
    double ret[ndim] ={0};
    drift(s, ret, 1.0, true);
    return h*h/(pvar[0] + h*norm(ret));
}

#endif
