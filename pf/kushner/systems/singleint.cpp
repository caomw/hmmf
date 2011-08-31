#include "singleint.h"

System::System()
{
    min_states = new float[NUM_DIM];
    max_states = new float[NUM_DIM];
    obs_noise = new float[NUM_DIM];
    process_noise = new float[NUM_DIM];
    init_var = new float[NUM_DIM];

    for(int i=0; i< NUM_DIM; i++)
    {
        min_states[i] = 0;
        max_states[i] = 1.0;
        init_state.x[i] = 1.0;
    }
    
    for(int i=0; i< NUM_DIM; i++)
    {
        process_noise[i] = 1e-3;
        obs_noise[i] = 1e-3;
        init_var[i] = 1e-3;
    }
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
        var[i] = process_noise[i]/6*( exp(6*duration) -1);
        tmp[i] = 0;
        mean[i] = 0;
    }
    if( !is_clean)  
        multivar_normal( mean, var, tmp, NUM_DIM);

    for(int i=0; i<NUM_DIM; i++)
        t.x[i] = exp(-3*duration)*t.x[i] + tmp[i];

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
    else
    {
        for(int i=0; i<NUM_DIM; i++)
            tmp[i] = 0;
    }

    for(int i=0; i<NUM_DIM; i++)
        t.x[i] = s.x[i] + tmp[i-1];

    delete[] mean;
    delete[] tmp;

    return t;
}

void System::get_kalman_path( vector<State>& obs, list<State>& kalman_path)
{
    /*
    kalman_path.clear();

    kalman_path.push_back(init_state);

    float *Q = new float[NUM_DIM];
    for(int i=0; i < NUM_DIM; i++)
        Q[i] = init_var[i];

    State curr_state = init_state;
    for(list<State>::iterator i = obs.begin(); i!= obs.end(); i++)
    {
        State& next_obs = *i;
        float delta_t = next_obs.x[0] - curr_state.x[0];
        
        //cout<<"delta_t: "<< delta_t << endl;
        State next_state = integrate(curr_state, delta_t, true);
        State clean_obs = observation(next_state, true);

        for(int j=0; j < NUM_DIM-1; j++)
        {
            Q[j] = exp(-6*delta_t)*Q[j] + process_noise[j]/6*( exp(6*delta_t) -1);
            float S = next_obs.x[j+1] - clean_obs.x[j+1];
            float L = Q[j]/(Q[j] + obs_noise[j]);

            next_state.x[j+1] += L*S;

            Q[j] = (1 -L)*Q[j];
        }
        curr_state = next_state;

        kalman_path.push_back(curr_state);
    }
    */

    /*
    for(int i= 0; i< max_states[0]/DT; i++)
    {
            // update xkf
            xkf[i].x[0] = x[i].x[0];
            state curr_obs = observation(xkf[i], 1);
            float S = y[i].x[dim] - curr_obs.x[dim];
            float L = Q/(Q + OBS_VAR[dim]);
            xkf[i].x[dim] += L*S;

            // update covar
            Q = (1 - L)*Q;

            // propagate
            xkf[i+1].x[dim] = exp(-3*DT)*(xkf[i].x[dim]);
            Q = exp(-6*DT)*Q + PRO_VAR[dim]/6*(exp(6*DT) - 1);
        }
    }
    */
}

