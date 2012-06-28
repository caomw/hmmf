#ifndef __ship__
#define __ship__

#define ndim        (4)
#define ndim_obs    (2)

#include "sysutil.h"

float zmin[ndim] = {-10, -10, -2, -2};
float zmax[ndim] = {10, 10, 2, 2};
float zdiff[ndim] = {20, 20, 4, 4};
float init_var[ndim] = {1e-2, 1e-2, 1e-2, 1e-2};
float init_state[ndim] = {-2, -2, 1, 1};
float init_state_real[ndim] = {-2, -2, 1, 1};
float pvar[ndim] = {1e-1, 1e-1, 1e-1, 1e-1};
float ovar[ndim] = {1e-2, 1e-2, 1e-6, 1e-6};
float zero[ndim] = {0};

int drift(float* s, float *ret, float dt=1.0, bool real=false)
{
    float f1dt=0, f2dt = 0;
    float dist = sqrt(s[0]*s[0] + s[1]*s[1]);
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
int diffusion(float* s, float* ret, float dt=1.0, bool real=false)
{
    float var[ndim] ={0};
    for(int i=0; i<ndim; i++)
        var[i] = pvar[i]*dt;
    multivar_normal(zero, var, ret, ndim);
    return 0;
}
int get_obs(float* s, float* obs, bool is_clean=false)
{
    float range = sqrt(s[0]*s[0] + s[1]*s[1]);
    float theta = atan2(s[1], s[0]);
    for(int i=0; i< ndim_obs; i++)
        obs[i] = 0;
    if(is_clean)
    {
        obs[0] = range;
        obs[1] = theta;
        return 0;
    }
    float noise[ndim_obs] = {0};
    multivar_normal(zero, ovar, noise, ndim_obs); 
    obs[0] = range + noise[0];
    obs[1] = theta + noise[1];
    return 0;
}
float holding_time(float* s, float r)
{
    float h = r*(zmax[0] - zmin[0]);
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

int get_kalman_path(vector< vector<float> >& kfp, vector< vector<float> >& observations, float dt)
{
    kfp.clear();
    
    Matrix4d Q = Matrix4d::Zero();
    for(int i=0; i<4; i++)
        Q(i,i) = init_var[i];

    Vector4d curr_state;
    for(int j=0; j<ndim; j++)
        curr_state(j) = init_state[j];
    
    Matrix4d Rk = Matrix4d::Zero();
    for(int i=0; i<4; i++)
        Rk(i,i) = ovar[i];

    for(unsigned int i=0; i< observations.size(); i++)
    {
        vector<float>& next_obs = observations[i];
        Vector4d noisy_obs;
        noisy_obs(0) = next_obs[0];
        noisy_obs(1) = next_obs[1];
        noisy_obs(2) = 0;
        noisy_obs(3) = 0;

        float stmp1[ndim] = {0};
        for(int j=0; j<ndim; j++)
            stmp1[j] = curr_state(j);
        
        float next_state[ndim] = {0}, next_state_obs[ndim] = {0};
        copy(stmp1, next_state);
        integrate_system(next_state, dt, true);
        get_obs(next_state, next_state_obs, true);
        
        Vector4d obs_vector;
        obs_vector(0) = next_state_obs[0];
        obs_vector(1) = next_state_obs[1];
        obs_vector(2) = 0;
        obs_vector(3) = 0;

        // Linearize Ad
        Matrix4d Ad = Matrix4d::Zero();
        Ad(0,2) = 1; Ad(1,3) = 1;
        if( ( sqrt(stmp1[0]*stmp1[0] + stmp1[1]*stmp1[1]) >= 9) && \
                ( (stmp1[0]*stmp1[2] + stmp1[1]*stmp1[3]) >= 0) )
        {
            float den = pow(stmp1[0]*stmp1[0] + stmp1[1]*stmp1[1],1.5);
            Ad(2,1) = -50*stmp1[1]*stmp1[1]/den;
            Ad(2,2) = -50*(-stmp1[0]*stmp1[1])/den;
            Ad(3,1) = -50*(-stmp1[0]*stmp1[1])/den;
            Ad(3,2) = -50*stmp1[0]*stmp1[0]/den;

        }
        Ad = Ad*dt;

        // Linearize Cd
        Matrix4d Cd = Matrix4d::Zero();
        float range = sqrt(next_state[0]*next_state[0] + next_state[1]*next_state[1]);
        Cd(0,0) = next_state[0]/range; Cd(0,1) = next_state[1]/range;
        Cd(1,0) = -next_state[1]/range/range;
        Cd(1,1) = next_state[0]/range/range;
        
        Matrix4d Wk = Matrix4d::Zero();
        for(int j=0; j<4; j++)
            Wk(j,j) = pvar[j]*dt;

        Matrix4d Q1 = Ad * Q * Ad.transpose() + Wk;
        Matrix4d Lk = Q1*Cd.transpose()*(Cd*Q1*Cd.transpose() + Rk).inverse();

        Vector4d Sk;
        Sk = noisy_obs - obs_vector;
        for(int j=0; j<ndim; j++)
            curr_state(j) = next_state[j];
        
        Vector4d estimate = curr_state + Lk*Sk;

        Q = (Matrix4d::Identity() - Lk*Cd)*Q1;
        curr_state = estimate;
        
        //cout<< "Q: "<< Q << endl << endl;

        vector<float> stmp2(ndim,0);
        for(int j=0; j<ndim; j++)
            stmp2[j] = curr_state(j);
        
        kfp.push_back(stmp2);
    }
    return 0;
}
#endif
