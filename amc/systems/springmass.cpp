#include "springmass.h"

System::System()
{
    min_states = new double[NUM_DIM];
    max_states = new double[NUM_DIM];
    obs_noise = new double[NUM_DIM];
    process_noise = new double[NUM_DIM];
    init_var = new double[NUM_DIM];

    for(int i=0; i< NUM_DIM; i++)
    {
        min_states[i] = -12;
        max_states[i] = 12;
    }
    max_states[0] = 5;
    min_states[0] = -5;
    init_state.x[0] = 1.2;
    init_state.x[1] = 1;

    for(int i=0; i< NUM_DIM; i++)
    {
        process_noise[i] = 1;
        obs_noise[i] = 1;
        init_var[i] = 1;
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

    double delta_t = min(duration, 0.01);
    for(int i=0; i<NUM_DIM; i++)
    {
        t.x[i] = s.x[i];
    }

    for(int i=0; i<NUM_DIM; i++)
    {
        tmp[i] = 0;
        mean[i] = 0;
    }
        
    double curr_time = 0;
    while(curr_time < duration)
    {
        if( !is_clean)
        {
            for(int i=0; i<NUM_DIM; i++)
                tmp[i] = 0;
            
            var[0] = process_noise[0]*delta_t;
            var[1] = pow(exp(-t.x[0]*t.x[0]/100),2)*process_noise[1]*delta_t;
            multivar_normal( mean, var, tmp, NUM_DIM);
        }
        
        t.x[0] += (t.x[1]*delta_t + tmp[0]);
        t.x[1] += ((-t.x[0] - pow(t.x[0],3))*delta_t + tmp[1]);

        curr_time += min(delta_t, duration - curr_time);
    }

    delete[] mean;
    delete[] tmp;
    delete[] var;

    return t;
}

void System::get_variance(State& s, double duration, double* var)
{
    var[0] = process_noise[0]*duration;
    var[1] = pow(exp(-s.x[0]*s.x[0]/100),2)*process_noise[1]*duration;
}


State System::observation(State& s, bool is_clean)
{
    State t;

    double *tmp = new double[NUM_DIM_OBS];
    double *mean = new double[NUM_DIM_OBS];

    if( !is_clean)  
        multivar_normal( mean, obs_noise, tmp, NUM_DIM_OBS);

    if(is_clean)
    {
        for(int i=0; i<NUM_DIM_OBS; i++)
            tmp[i] = 0;
    }

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
    
    Matrix3d Q;
    Q << init_var[0], 0, 0, 0, init_var[1], 0, 0, 0, init_var[2];
    Vector3d curr_state;
    curr_state(0) = init_state.x[0];
    curr_state(1) = init_state.x[1];
    curr_state(2) = init_state.x[2];
    
    MatrixXd Rk(3,3);
    Rk << obs_noise[0], 0, 0, 0, obs_noise[1], 0, 0, 0, obs_noise[2];
    Matrix3d Cd = Matrix3d::Identity();

    double prev_time = 0;
    for(unsigned int i=0; i< obs.size(); i++)
    {
        State& next_obs = obs[i];
        Vector3d noisy_obs;
        noisy_obs(0) = next_obs.x[0];
        noisy_obs(1) = next_obs.x[1];
        noisy_obs(2) = next_obs.x[2];

        double delta_t = obs_times[i] - prev_time;
        
        //cout<<"delta_t: "<< delta_t << endl;
        State stmp1;
        stmp1.x[0] = curr_state(0);
        stmp1.x[1] = curr_state(1);
        stmp1.x[2] = curr_state(2);
        State next_state = integrate(stmp1, delta_t, true);
        State clean_obs = observation(next_state, true);
        
        Matrix3d Ad;
        Ad << 0, 0, -sin(stmp1.x[2]), 0, 0, cos(stmp1.x[2]), 0, 0, 0; 
        
        Matrix3d Wk;
        Wk << process_noise[0]*delta_t, 0, 0, 0, process_noise[1]*delta_t, 0, 0, 0, process_noise[2]*delta_t;
        
        Matrix3d Q_new = Ad * Q * Ad.transpose() + Wk;
        Matrix3d Lk = Q_new*Cd.transpose()*(Cd*Q_new*Cd.transpose() + Rk).inverse();
        
        Vector3d Sk;
        curr_state(0) = next_state.x[0];
        curr_state(1) = next_state.x[1];
        curr_state(2) = next_state.x[2];
        Sk = noisy_obs - Cd*curr_state;
        
        Vector3d estimate = curr_state + Lk*Sk;
        
        Matrix3d covar = (Matrix3d::Identity() - Lk*Cd)*Q_new;

        Q = covar;
        curr_state = estimate;

        State stmp2;
        stmp2.x[0] = curr_state(0);
        stmp2.x[1] = curr_state(1);
        stmp2.x[2] = curr_state(2);
        kalman_path.push_back(stmp2);
        prev_time = obs_times[i];
    }

#endif
}

