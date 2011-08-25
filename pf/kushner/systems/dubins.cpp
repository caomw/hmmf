#include "dubins.h"

System::System()
{
    min_states = new float[NUM_DIM];
    min_states[0] = 0; min_states[1] = 0; min_states[2] = 0;
    
    max_states = new float[NUM_DIM];
    max_states[0] = 1; max_states[1] = 1; max_states[2] = M_PI_2;

    obs_noise = new float[NUM_DIM];
    obs_noise[0] = 1e-2; obs_noise[1] = 1e-2; obs_noise[2] = 1e-2;
    
    process_noise = new float[NUM_DIM];
    process_noise[0] = 1e-2; process_noise[1] = 1e-2; process_noise[2] = 1e-2;
    
    init_var = new float[NUM_DIM];
    init_var[0] = 1e-2; init_var[1] = 1e-2; init_var[2] = 1e-2;
    
    init_state.x[0] = 0;
    init_state.x[1] = 0;
    init_state.x[2] = 0;

    sim_time_delta = 0.01;
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

    bool retflag = 0;

    // obs 1
    if( (s[0] >= 0.127) && (s[0] <= 0.26) )
    {
        if( (s[1] >= 0) && (s[1] <= .217) )
            retflag = 0;
        else
            retflag = 1;
    }
    else
        retflag = 1;

    if (retflag == 0)
        return 0;

    // obs 2
    if( (s[0] >= 0.1) && (s[0] <= 0.2) )
    {
        if( (s[1] >= .32) && (s[1] <= .5) )
            retflag = 0;
        else
            retflag = 1;
    }
    else
        retflag = 1;

    return retflag;
}

int System::get_key(State& s, double *key)
{
    for(int i =0; i < NUM_DIM; i++)
    {
        key[i] = (s.x[i] - min_states[i])/(max_states[i] - min_states[i]);
        assert(key[i] <= 1.1);
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
        var[i] = process_noise[i]*duration;
        tmp[i] = 0;
        mean[i] = 0;
    }
    if( !is_clean)  
        multivar_normal( mean, var, tmp, NUM_DIM);
        
    float curr_time = 0;
    while(curr_time < duration)
    {
        t.x[0] += 1.0*cos(t.x[2])*sim_time_delta;
        t.x[1] += 1.0*sin(t.x[2])*sim_time_delta;
        t.x[2] += 1.0*sim_time_delta;

        curr_time += sim_time_delta;
    }

    for(int i=0; i<NUM_DIM; i++)
        t.x[i] += tmp[i];

    delete[] mean;
    delete[] tmp;
    delete[] var;

    return t;
}


State System::observation(State& s, bool is_clean)
{
    State t;

    float *tmp = new float[NUM_DIM];
    float *mean = new float[NUM_DIM];

    if( !is_clean)  
        multivar_normal( mean, obs_noise, tmp, NUM_DIM);

    if(is_clean)
    {
        for(int i=0; i<NUM_DIM; i++)
            tmp[i] = 0;
    }

    for(int i=0; i<NUM_DIM; i++)
        t.x[i] = s.x[i] + tmp[i];

    delete[] mean;
    delete[] tmp;

    return t;
}

