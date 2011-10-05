#include "singleint.h"

System::System()
{
    min_states = new double[NUM_DIM];
    max_states = new double[NUM_DIM];
    obs_noise = new double[NUM_DIM_OBS-1];
    process_noise = new double[NUM_DIM-1];
    init_var = new double[NUM_DIM-1];

    for(int i=0; i< NUM_DIM; i++)
    {
        min_states[i] = -2;
        max_states[i] = 2;
        init_state.x[i] = 1;
    }
    max_states[0] = 1.0;
    min_states[0] = 0;
    init_state.x[0] = 0;

    for(int i=0; i< NUM_DIM-1; i++)
    {
        process_noise[i] = 1e-2;
        init_var[i] = 1e-2;
    }
    for(int i=0; i< NUM_DIM_OBS-1; i++)
    {
        obs_noise[i] = 1e-2;
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
    double border = RANDF;
    double which_time = 0;
    if(border < 0.05)
    {
        which_time = max_states[0];
    }
    while(1)
    {
        for(int i=0; i< NUM_DIM; i++)
        {
            s.x[i] = min_states[i] + RANDF*( max_states[i] - min_states[i]);
        }
        if(border < 0.05)
            s.x[0] = which_time;

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

    double *var = new double[NUM_DIM-1];
    double *mean = new double[NUM_DIM-1];
    double *tmp = new double[NUM_DIM-1];

    for(int i=0; i< NUM_DIM; i++)
    {
        t.x[i] = s.x[i];
    }

    for(int i=0; i<NUM_DIM-1; i++)
    {   
        var[i] = process_noise[i]/6*( exp(6*duration) -1);
        tmp[i] = 0;
        mean[i] = 0;
    }
    if( !is_clean)  
        multivar_normal( mean, var, tmp, NUM_DIM-1);

    for(int i=1; i<NUM_DIM; i++)
        t.x[i] = exp(-3*duration)*t.x[i] + tmp[i-1];
    
    t.x[0] += duration;

    delete[] mean;
    delete[] tmp;
    delete[] var;

    return t;
}


void System::get_fdt(State& s, double duration, double* next)
{
    for(int i=1; i < NUM_DIM; i++)
    {
        next[i] = -3*s.x[i]*duration;
    }
    next[0] = s.x[0] + duration;
}

void System::get_variance(State& s, double duration, double* var)
{
    for(int i=0; i<NUM_DIM-1; i++)
    {   
        var[i] = process_noise[i]/6*( exp(6*duration) -1);
    } 
}
void System::get_obs_variance(State& s, double* var)
{
    for(int i=0; i<NUM_DIM_OBS-1; i++)
    {   
        var[i] = obs_noise[i];
    } 
}

State System::observation(State& s, bool is_clean)
{
    State t;

    double *tmp = new double[NUM_DIM_OBS-1];
    double *mean = new double[NUM_DIM_OBS-1];

    if( !is_clean)  
        multivar_normal( mean, obs_noise, tmp, NUM_DIM_OBS);
    else
    {
        for(int i=0; i<NUM_DIM_OBS-1; i++)
            tmp[i] = 0;
    }

    for(int i=1; i<NUM_DIM_OBS; i++)
        t.x[i] = s.x[i] + tmp[i-1];

    delete[] mean;
    delete[] tmp;

    return t;
}

void System::get_kalman_path( vector<State>& obs, vector<double>& obs_times, list<State>& kalman_path, list<State>& kalman_covar)
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

void System::get_pf_path( vector<State>& obs, vector<double>& obs_times, list<State>& pf_path)
{

}
