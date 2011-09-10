#include "vanderpol.h"

System::System()
{
    min_states = new double[NUM_DIM];
    max_states = new double[NUM_DIM];
    obs_noise = new double[NUM_DIM];
    process_noise = new double[NUM_DIM];
    init_var = new double[NUM_DIM];

    for(int i=0; i< NUM_DIM; i++)
    {
        min_states[i] = -6;
        max_states[i] = 6;
    }
    init_state.x[0] = 0;
    init_state.x[1] = 1.0;

    for(int i=0; i< NUM_DIM; i++)
    {
        process_noise[i] = 1e-4;
        obs_noise[i] = 1e-4;
        init_var[i] = 1e-4;
    }
    
    sim_time_delta = 0.005;
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

State System::integrate(State& s, double duration, bool is_clean)
{
    State t;

    double *var = new double[NUM_DIM];
    double *mean = new double[NUM_DIM];
    double *tmp = new double[NUM_DIM];
    
    double delta_t = min(duration, 0.005);

    for(int i=0; i<NUM_DIM; i++)
    {
        t.x[i] = s.x[i];
    }

    for(int i=0; i<NUM_DIM; i++)
    {   
        var[i] = process_noise[i]*delta_t;
        tmp[i] = 0;
        mean[i] = 0;
    }

    double curr_time = 0;
    while(curr_time < duration)
    {
        if( !is_clean)
            multivar_normal( mean, var, tmp, NUM_DIM);
        
        double f1dt=0, f2dt = 0;
        f1dt = t.x[1]*delta_t;
        f2dt = (-t.x[0] + 2.0*t.x[1]*(1 - t.x[0]*t.x[0]) )*delta_t;
    
        t.x[0] = t.x[0] + f1dt + tmp[0];
        t.x[1] = t.x[1] + f2dt + tmp[1];

        curr_time += min(delta_t, duration - curr_time);
    }

    delete[] mean;
    delete[] tmp;
    delete[] var;

    return t;
}

void System::get_variance(State& s, double duration, double* var)
{
    for(int i=0; i<NUM_DIM; i++)
    {   
        var[i] = process_noise[i]*duration;
    } 
}

State System::observation(State& s, bool is_clean)
{
    State t;

    double *tmp = new double[NUM_DIM_OBS];
    double *mean = new double[NUM_DIM_OBS];

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

void System::get_kalman_path( vector<State>& obs, vector<double>& obs_times, list<State>& kalman_path)
{
#if 0
    kalman_path.clear();

    kalman_path.push_back(init_state);

    double *Q = new double[NUM_DIM];
    for(int i=0; i < NUM_DIM; i++)
        Q[i] = init_var[i];

    State curr_state = init_state;
    double prev_time = 0;
    for(unsigned int i=0; i< obs.size(); i++)
    {
        State& next_obs = obs[i];
        double delta_t = obs_times[i] - prev_time;
        
        //cout<<"delta_t: "<< delta_t << endl;
        State next_state = integrate(curr_state, delta_t, true);
        State clean_obs = observation(next_state, true);

        for(int j=0; j < NUM_DIM; j++)
        {
            Q[j] = exp(-6*delta_t)*Q[j] + process_noise[j]/6*( exp(6*delta_t) -1);
            double S = next_obs.x[j] - clean_obs.x[j];
            double L = Q[j]/(Q[j] + obs_noise[j]);

            next_state.x[j] += L*S;

            Q[j] = (1 -L)*Q[j];
        }
        curr_state = next_state;

        kalman_path.push_back(curr_state);
        prev_time = obs_times[i];
    }
#endif
}

