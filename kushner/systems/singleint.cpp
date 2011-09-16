#include "singleint.h"

System::System()
{
    min_states = new double[NUM_DIM];
    max_states = new double[NUM_DIM];
    obs_noise = new double[NUM_DIM];
    process_noise = new double[NUM_DIM];
    init_var = new double[NUM_DIM];

    for(int i=0; i< NUM_DIM; i++)
    {
        min_states[i] = 0;
        max_states[i] = 1;
        init_state.x[i] = 0.90;
        beacon_state.x[i] = 1.0;
    }
    
    for(int i=0; i< NUM_DIM; i++)
    {
        process_noise[i] = 1e-2;
        obs_noise[i] = 1e-2;
        init_var[i] = 1e-2;
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
    bool free_state = 1;

    // obs 1
    if( (s[0] <= 0.7) && (s[0] >= 0.56) )
    {
        if( (s[1] <= 0.6) && (s[1] >= 0.4) )
            free_state = 0;
    }
    
    return free_state;
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

void System::get_variance(State& s, double duration, double* var)
{
    for(int i=0; i<NUM_DIM; i++)
    {   
        var[i] = process_noise[i]/6*( exp(6*duration) -1);
    } 
}

void System::get_obs_variance(State& s, double* var)
{
    double range = sqrt((s.x[0]-beacon_state.x[0])*(s.x[0]-beacon_state.x[0]) +\
            (s.x[1]-beacon_state.x[1])*(s.x[1]-beacon_state.x[1]) );
    
    for(int i=0; i<NUM_DIM_OBS; i++)
    {
        var[i] = obs_noise[i]; //*pow((1 + range),2);
    }
}

State System::observation(State& s, bool is_clean)
{
    State t;

    double range = sqrt((s.x[0]-beacon_state.x[0])*(s.x[0]-beacon_state.x[0]) +\
            (s.x[1]-beacon_state.x[1])*(s.x[1]-beacon_state.x[1]) );
    
    double *tmp = new double[NUM_DIM_OBS];
    double *mean = new double[NUM_DIM_OBS];
    double *var = new double[NUM_DIM_OBS];

    get_obs_variance(s, var);
    if( !is_clean)  
        multivar_normal( mean, var, tmp, NUM_DIM_OBS);
    else
    {
        for(int i=0; i<NUM_DIM_OBS; i++)
            tmp[i] = 0;
    }
     
    t.x[0] = range + tmp[0];

    delete[] mean;
    delete[] tmp;
    delete[] var;

    return t;
}

void System::get_kalman_path( vector<State>& obs, vector<double>& obs_times, list<State>& kalman_path, list<State>& kalman_covar)
{
#if 1
    kalman_path.clear();
    kalman_covar.clear();

    kalman_path.push_back(init_state);
    kalman_covar.push_back(init_var);

    Matrix2d Q;
    Q << init_var[0], 0, 0, init_var[1];
    Vector2d curr_state;
    curr_state(0) = init_state.x[0];
    curr_state(1) = init_state.x[1];
    
    double prev_time = 0;
    for(unsigned int i=0; i< obs.size(); i++)
    {
        State& next_obs = obs[i];
        MatrixXd noisy_obs(1,1);
        noisy_obs(0,0) = next_obs.x[0];

        double delta_t = obs_times[i] - prev_time;
        
        //cout<<"delta_t: "<< delta_t << endl;
        State stmp1;
        stmp1.x[0] = curr_state(0);
        stmp1.x[1] = curr_state(1);
        State next_state = integrate(stmp1, delta_t, true);
        State clean_obs = observation(next_state, true);
        MatrixXd obs_vector(1,1);
        obs_vector(0,0) = clean_obs.x[0];

        Matrix2d Ad = Matrix2d::Zero();
        Ad(0,0) = -3;
        Ad(1,1) = -3;
        
        Matrix2d Wk;
        Wk(0,0) = process_noise[0]*delta_t; Wk(0,1) = 0;
        Wk(1,1) = process_noise[1]*delta_t; Wk(1,0) = 0;
        
        MatrixXd Rk(1,1);
        double *obs_var = new double[NUM_DIM_OBS];
        get_obs_variance(next_state, obs_var);
        Rk(0,0) = obs_var[0];
        delete[] obs_var;

        MatrixXd Cd(1,2);
        double range = sqrt((next_state.x[0]-beacon_state.x[0])*(next_state.x[0]-beacon_state.x[0]) +\
            (next_state.x[1]-beacon_state.x[1])*(next_state.x[1]-beacon_state.x[1]) );
        Cd(0,0) = (next_state.x[0] - beacon_state.x[0])/range;
        Cd(0,1) = (next_state.x[1] - beacon_state.x[1])/range;
        cout<<"Cd: " << Cd << endl;
        cout<<"range: "<< range << endl;

        Matrix2d Q_new = Ad * Q * Ad.transpose() + Wk;
        cout<< "Q_new: "<< Q_new << endl << endl;
        Vector2d Lk = Q_new*Cd.transpose()*(Cd*Q_new*Cd.transpose() + Rk).inverse();
        cout<< "Lk: "<< Lk << endl << endl;
        
        MatrixXd Sk(1,1);
        curr_state(0) = next_state.x[0];
        curr_state(1) = next_state.x[1];
        Sk = noisy_obs - obs_vector;
        cout<<"Sk: "<< Sk << endl;

        Vector2d estimate = curr_state + Lk*Sk;
        
        MatrixXd covar = (Matrix2d::Identity() - Lk*Cd)*Q_new;

        Q = covar;
        curr_state = estimate;
        
        cout<< "covar: "<< covar << endl << endl;

        State stmp2;
        stmp2.x[0] = curr_state(0);
        stmp2.x[1] = curr_state(1);
        kalman_path.push_back(stmp2);

        State stmp3;
        stmp3.x[0] = Q(0,0);
        stmp3.x[1] = Q(1,1);
        kalman_covar.push_back(stmp3);

        prev_time = obs_times[i];

        getchar();
    }

#endif
}

