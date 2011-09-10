#include "ship.h"

System::System()
{
    min_states = new float[NUM_DIM];
    max_states = new float[NUM_DIM];
    obs_noise = new float[NUM_DIM];
    process_noise = new float[NUM_DIM];
    init_var = new float[NUM_DIM];

    for(int i=0; i< NUM_DIM; i++)
    {
        min_states[i] = -10;
        max_states[i] = 10;
        init_state.x[i] = -7;
    }
    init_state.x[2] = 0;
    init_state.x[3] = 0;

    for(int i=0; i< NUM_DIM; i++)
    {
        process_noise[i] = 0.3*0.3;
        obs_noise[i] = 1e-2;
        init_var[i] = 1e-2;
    }
    
    sim_time_delta = 0.1;
}

System::~System()
{
    delete[] min_states;
    delete[] max_states;
    delete[] obs_noise;
    delete[] process_noise;
    delete[] init_var;
}

State System::sample()
{
    State s;
    while(1)
    {
        for(int i=0; i< NUM_DIM; i++)
        {
            s.x[i] = min_states[i] + RANDF*( max_states[i] - min_states[i]);
        }

        if( is_free(s) )
            break;
    }
    return s;
}

bool System::is_free(State &s)
{
    return 1;
}

int System::get_key(State& s, double *key)
{
    for(int i =0; i < NUM_DIM; i++)
    {
        key[i] = (s.x[i] - min_states[i])/(max_states[i] - min_states[i]);
        //assert(key[i] <= 1.1);
    }
    return 0;
}

State System::integrate(State& s, float duration, bool is_clean)
{
    State t;

    float *var = new float[NUM_DIM];
    float *mean = new float[NUM_DIM];
    float *tmp = new float[NUM_DIM];

    for(int i=0; i<NUM_DIM; i++)
    {
        t.x[i] = s.x[i];
    }

    for(int i=0; i<NUM_DIM; i++)
    {   
        var[i] = process_noise[i]*0.01;
        tmp[i] = 0;
        mean[i] = 0;
    }

    float curr_time = 0;
    while(curr_time < duration)
    {
        if( !is_clean)
            multivar_normal( mean, var, tmp, NUM_DIM);
        
        float old_x3 = t.x[2]; float old_x4 = t.x[3];
        
        float f1dt=0, f2dt = 0;
        if( ( sqrt(t.x[0]*t.x[0] + t.x[1]*t.x[1]) >= 9) && ( (t.x[0]*t.x[2] + t.x[1]*t.x[3]) >= 0) )
        {
            f1dt = -50*sim_time_delta*t.x[0]/sqrt(t.x[0]*t.x[0] + t.x[1]*t.x[1]);
            f2dt = -50*sim_time_delta*t.x[1]/sqrt(t.x[0]*t.x[0] + t.x[1]*t.x[1]);
        }
        t.x[2] = t.x[2] + f1dt + tmp[2];
        t.x[3] = t.x[3] + f2dt + tmp[3];

        t.x[0] = t.x[0] + sim_time_delta*(old_x3) + tmp[0]; 
        t.x[1] = t.x[1] + sim_time_delta*(old_x4) + tmp[1];

        curr_time += 0.01;
    }

    delete[] mean;
    delete[] tmp;
    delete[] var;

    return t;
}

void System::get_variance(State& s, float duration, float* var)
{
    for(int i=0; i<NUM_DIM; i++)
    {   
        var[i] = process_noise[i]*duration;
    } 
}

State System::observation(State& s, bool is_clean)
{
    State t;

    float *tmp = new float[NUM_DIM_OBS];
    float *mean = new float[NUM_DIM_OBS];

    if( !is_clean)  
        multivar_normal( mean, obs_noise, tmp, NUM_DIM_OBS);
    else
    {
        for(int i=0; i<NUM_DIM_OBS; i++)
            tmp[i] = 0;
    }
    

    for(int i=0; i<NUM_DIM_OBS; i++)
        t.x[i] = 0;

    for(int i=0; i<NUM_DIM_OBS; i++)
        t.x[i] = s.x[i] + tmp[i];

    delete[] mean;
    delete[] tmp;

    return t;
}

void System::get_kalman_path( vector<State>& obs, vector<float>& obs_times, list<State>& kalman_path)
{
    return;
}

