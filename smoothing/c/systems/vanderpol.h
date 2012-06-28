#ifndef __vanderpol__
#define __vanderpol__

#define ndim        (2)
#define ndim_obs    (2)

#include "sysutil.h"

float zmin[ndim] = {-2.5, -4.0};
float zmax[ndim] = {2.0, 2.5};
float init_var[ndim] = {1e-2, 1e-2};
float init_state[ndim] = {0.0, 1.0};
float init_state_real[ndim] = {0.0, 1.0};
float pvar[ndim] = {5*1e-2, 5*1e-2};
float ovar[ndim_obs] = {5*1e-2, 1e-2};
float zero[ndim] = {0};

int drift(float* s, float *ret, float dt=1.0, bool real=false)
{
    float mu = 2.0;
    ret[0] = s[1]*dt;
    ret[1] = (-s[0] + mu*s[1]*(1-sq(s[0])))*dt;
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
int get_obs(float* s, float* obs, bool is_clean=false)
{
    for(int i=0; i< ndim_obs; i++)
        obs[i] = 0;
    if(is_clean)
    {
        for(int i=0; i< ndim_obs; i++)
            obs[i] = s[i];
        return 0;
    }
    float noise[ndim_obs] = {0};
    multivar_normal(zero, ovar, noise, ndim_obs); 
    for(int i=0; i< ndim_obs; i++)
        obs[i] = s[i] + noise[i];
    return 0;
}
float holding_time(float* s, float r)
{
    float h = r*sqrt(sq(zmax[0] - zmin[0])+ sq(zmax[1] - zmin[1]));
    float ret[ndim] ={0};
    drift(s, ret, 1.0, true);
    return h*h/(pvar[0] + h*norm(ret));
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

#if 0
int get_kalman_path(vector< vector<float> >& kfp, vector< vector<float> >& observations, float dt)
{
    kfp.clear();

    Matrix2d Q;
    Q << init_var[0], 0, 0, init_var[1];
    Vector2d curr_state;
    curr_state(0) = init_state[0];
    curr_state(1) = init_state[1];
    
    MatrixXd Rk(1,1);
    Rk(0,0) = ovar[0];
    MatrixXd Cd(1,2);
    Cd(0,0) = 1; Cd(0,1) = 0;

    Matrix2d Wk;
    Wk(0,0) = pvar[0]*dt; Wk(0,1) = 0;
    Wk(1,1) = pvar[1]*dt; Wk(1,0) = 0;
    
    for(unsigned int i=0; i< observations.size(); i++)
    {
        vector<float> next_obs = observations[i];
        MatrixXd noisy_obs(1,1);
        noisy_obs(0,0) = next_obs[0];

        float stmp1[ndim] ={0}, next_state[ndim] = {0};
        stmp1[0] = curr_state(0);
        stmp1[1] = curr_state(1);
        copy(stmp1, next_state);

        integrate_system(next_state, dt, true);
        float clean_obs[ndim_obs] = {0};
        get_obs(next_state, clean_obs, true);
        
        Matrix2d Ad;
        Ad(0,0) = 0; Ad(0,1) = 1; Ad(1,0) = -1 -4*stmp1[0]*stmp1[1]; Ad(1,1) = 2*(1-stmp1[0]*stmp1[0]);
        
        Matrix2d Q_new = Ad * Q * Ad.transpose() + Wk;
        Vector2d Lk = Q_new*Cd.transpose()*(Cd*Q_new*Cd.transpose() + Rk).inverse();
        
        MatrixXd Sk(1,1);
        curr_state(0) = next_state[0];
        curr_state(1) = next_state[1];
        Sk = noisy_obs - Cd*curr_state;
        
        Vector2d estimate = curr_state + Lk*Sk;
        
        MatrixXd covar = (Matrix2d::Identity() - Lk*Cd)*Q_new;

        Q = covar;
        curr_state = estimate;

        vector<float> stmp2(ndim,0);
        stmp2[0] = curr_state(0);
        stmp2[1] = curr_state(1);
        kfp.push_back(stmp2);

    }
    return 0;
}
#endif

#endif
