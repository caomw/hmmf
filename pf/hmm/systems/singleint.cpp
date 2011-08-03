#include "singleint.h"

System::System()
{
    min_states = new double[NUM_DIM];
    min_states[0] = 0; min_states[1] = 0; min_states[2] = 0; min_states[3] = 0;
    
    max_states = new double[NUM_DIM];
    max_states[0] = 0.5; max_states[1] = 1; max_states[2] = 1; max_states[3] = 1;

    obs_noise = new double[NUM_DIM-1];
    obs_noise[0] = 1e-4; obs_noise[1] = 1e-4; obs_noise[2] = 1e-4;
    
    process_noise = new double[NUM_DIM-1];
    process_noise[0] = 1e-4; process_noise[1] = 1e-4; process_noise[2] = 1e-4;
    
    init_var = new double[NUM_DIM-1];
    init_var[0] = 1e-4; init_var[1] = 1e-4; init_var[2] = 1e-4;
    
    init_state.x[0] = 0;
    init_state.x[1] = 1;
    init_state.x[2] = 1;
    init_state.x[3] = 1;

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
    //return 1;

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

State System::integrate(State& s, double duration, bool is_clean)
{
    State t;

    double *var = new double[NUM_DIM-1];
    double *mean = new double[NUM_DIM-1];
    double *tmp = new double[NUM_DIM-1];

    for(int i=0; i<NUM_DIM; i++)
    {
        t.x[i] = s.x[i];
    }

    for(int i=0; i<NUM_DIM-1; i++)
    {
        var[i] = process_noise[i]*duration;
        tmp[i] = 0;
        mean[i] = 0;
    }
    if( !is_clean)  
        multivar_normal( mean, var, tmp, NUM_DIM-1);

    t.x[0] = t.x[0] + duration;
    for(int i=1; i<NUM_DIM; i++)
        t.x[i] = exp(-3*duration)*t.x[i] + tmp[i-1];

    delete[] mean;
    delete[] tmp;
    delete[] var;

    return t;
}


State System::observation(State& s, bool is_clean)
{
    State t;

    double *tmp = new double[NUM_DIM-1];
    double *mean = new double[NUM_DIM-1];

    if( !is_clean)  
        multivar_normal( mean, obs_noise, tmp, NUM_DIM-1);

    if(is_clean)
    {
        for(int i=0; i<NUM_DIM-1; i++)
            tmp[i] = 0;
    }

    t.x[0] = s.x[0];                        // time is same
    for(int i=1; i<NUM_DIM; i++)
        t.x[i] = s.x[i] + tmp[i-1];

    delete[] mean;
    delete[] tmp;

    return t;
}

