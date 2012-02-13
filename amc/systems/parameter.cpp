#include "parameter.h"

System::System()
{
    min_states = new double[NUM_DIM];
    max_states = new double[NUM_DIM];
    obs_noise = new double[NUM_DIM];
    process_noise = new double[NUM_DIM];
    init_var = new double[NUM_DIM];

    for(int i=0; i< NUM_DIM; i++)
    {
        min_states[i] = -0.5;
        max_states[i] = 1.5;
    }
    init_state.x[0] = 0;
    init_state.x[1] = 0.9;

    for(int i=0; i< NUM_DIM; i++)
    {
        process_noise[i] = 1e-2;
        obs_noise[i] = 1e-2;
        init_var[i] = 1e-2;
    }
    process_noise[1] = 1e-4;
    init_var[1] = 1e-1;
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

State System::get_fdt(State& s, double duration)
{
    State stmp;
    stmp.x[0] = cos(2*M_PI*stmp.x[0]*stmp.x[1])*duration;
    stmp.x[1] = 0*duration;
    return stmp;
}

State System::integrate(State& s, double duration, bool is_clean, bool is_propagate)
{
    State t;

    double *var = new double[NUM_DIM];
    double *mean = new double[NUM_DIM];
    double *tmp = new double[NUM_DIM];
    
    double delta_t = min(duration, 0.001);

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
        
        double phi = t.x[1];
        if(is_propagate)
        {
            phi = 0.5;
            tmp[1] = 0;
        }
        double f1dt=0, f2dt = 0;
        f1dt = (cos(2*M_PI*phi*t.x[0]))*delta_t;
        f2dt = 0*delta_t;
    
        t.x[0] = t.x[0] + f1dt + tmp[0];
        t.x[1] = t.x[1] + f2dt + tmp[1];

        curr_time += min(delta_t, duration - curr_time);
    }

    delete[] mean;
    delete[] tmp;
    delete[] var;

    return t;
}





void System::get_kalman_path( vector<State>& obs, vector<double>& obs_times, list<State>& kalman_path, list<State>& kalman_covar)
{

#if 0
    kalman_path.clear();

    kalman_path.push_back(init_state);
    
    Matrix2d Q;
    Q << init_var[0], 0, 0, init_var[1];
    Vector2d curr_state;
    curr_state(0) = init_state.x[0];
    curr_state(1) = init_state.x[1];
    
    MatrixXd Rk(1,1);
    Rk(0,0) = obs_noise[0];
    MatrixXd Cd(1,2);
    Cd(0,0) = 1; Cd(0,1) = 0;

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
        
        Matrix2d Ad;
        Ad(0,0) = (2*M_PI*stmp1.x[1])*(-sin(2*M_PI*stmp1.x[1]*stmp1.x[0]));
        Ad(0,1) = (2*M_PI*stmp1.x[0])*(-sin(2*M_PI*stmp1.x[1]*stmp1.x[0]));
        Ad(1,0) = 0; Ad(1,1) = 0;
        Matrix2d Wk;
        Wk(0,0) = process_noise[0]*delta_t; Wk(0,1) = 0;
        Wk(1,1) = process_noise[1]*delta_t; Wk(1,0) = 0;
        
        Matrix2d Q_new = Ad * Q * Ad.transpose() + Wk;
        Vector2d Lk = Q_new*Cd.transpose()*(Cd*Q_new*Cd.transpose() + Rk).inverse();
        
        MatrixXd Sk(1,1);
        curr_state(0) = next_state.x[0];
        curr_state(1) = next_state.x[1];
        Sk = noisy_obs - Cd*curr_state;
        
        Vector2d estimate = curr_state + Lk*Sk;
        
        MatrixXd covar = (Matrix2d::Identity() - Lk*Cd)*Q_new;

        Q = covar;
        curr_state = estimate;

        State stmp2;
        stmp2.x[0] = curr_state(0);
        stmp2.x[1] = curr_state(1);
        kalman_path.push_back(stmp2);
        prev_time = obs_times[i];
    }

#endif
}

// importance resampling
int pfilter_resample(State* parts, double *weights, int num_particles)
{
    double totweight = 0;
    for(int i=0; i<num_particles; i++)
        totweight += weights[i];

    double* cum = new double[num_particles];
    double curr_tot = 0;
    for(int i=0; i<num_particles; i++)
    {
        weights[i] = weights[i]/totweight;
        curr_tot += weights[i];
        cum[i] = curr_tot;

        //reset to equal weights
        weights[i] = 1.0/num_particles;
    }

    State* newparts = new State[num_particles];
#if 1
    double u = RANDF/(double)num_particles;
    int i = 0;
    for(int j=0; j<num_particles; j++)
    {
        double tocheck = u + j/(double)num_particles;
        while(tocheck > cum[i])
            i++;

        newparts[j] = parts[i];
    }
    for(int j=0; j<num_particles; j++)
        parts[j] = newparts[j]; 
#endif

    delete[] cum;
    delete[] newparts;
    
    // weights are already set above
    return 0;
}

void System::get_pf_path( vector<State>& obs, vector<double>& obs_times, list<State>& pf_path, int num_particles)
{
    pf_path.clear();

    int obs_size = obs.size();

    State parts[num_particles];
    double weights[num_particles];
    for(int i=0; i< num_particles; i++)
    {
        multivar_normal(init_state.x, init_var, parts[i].x, NUM_DIM);
        weights[i] = 1/(double)num_particles;
    }

    //cout<<"stime: "<< stime<<" transition_time: "<< e->transition_time << " etime: "<< etime << endl;

    double curr_time = 0;
    pf_path.push_back(init_state);
    for(int obs_iter=0; obs_iter<obs_size; obs_iter++)
    {
        double delta_t = obs_times[obs_iter] - curr_time;           // for how much time to propagate
        //cout<<"parts: "<< parts[0].x[0]<<" delta_t: "<< delta_t << " obs: "<< obs_states[obs_iter].x[0]<<endl;

        for(int i=0; i< num_particles; i++)
        {
            parts[i] = integrate(parts[i], delta_t, false);
        }

        for(int i=0; i< num_particles; i++)
        {
            State particle_obs = observation(parts[i], true);
            weights[i] = weights[i] * normal_val( obs[obs_iter].x, obs_noise,\
                    particle_obs.x, NUM_DIM_OBS);
        }

        pfilter_resample(parts, weights, num_particles);

        // get mean
        State mean_state;
        double totw = 0;
        for(int i=0; i<NUM_DIM; i++)
        {
            mean_state.x[i] = 0;
        }
        for(int i=0; i< num_particles; i++)
        {
            //cout<<"weights[i]: "<< weights[i] << endl;
            for(int j=0; j< NUM_DIM; j++)
            {
                mean_state.x[j] = mean_state.x[j] + (parts[i].x[j]) * weights[i];
            }
            totw += weights[i];
        }
        for(int j=1; j< NUM_DIM; j++)
        {
            mean_state.x[j] = mean_state.x[j]/totw;
        }

        pf_path.push_back(mean_state);
        curr_time = obs_times[obs_iter];
        
        //cout<<"press key: "; getchar(); cout<<endl;
    }
}
