#include "ship.h"

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
        init_state.x[i] = -11;
        init_state_real.x[i] = -11;
    }
    init_state.x[2] = 0;
    init_state.x[3] = 0;
    init_state_real.x[2] = 0;
    init_state_real.x[3] = 0;

    for(int i=0; i< NUM_DIM; i++)
    {
        process_noise[i] = 1e-1;
        obs_noise[i] = 1e-1;
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

State System::integrate(State& s, double duration, bool is_clean, bool is_real)
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
        
        double old_x3 = t.x[2]; double old_x4 = t.x[3];
        
        double f1dt=0, f2dt = 0;
        if( ( sqrt(t.x[0]*t.x[0] + t.x[1]*t.x[1]) >= 9) && ( (t.x[0]*t.x[2] + t.x[1]*t.x[3]) >= 0) )
        {
            f1dt = -50*delta_t*t.x[0]/sqrt(t.x[0]*t.x[0] + t.x[1]*t.x[1]);
            f2dt = -50*delta_t*t.x[1]/sqrt(t.x[0]*t.x[0] + t.x[1]*t.x[1]);
        }
        t.x[2] = t.x[2] + f1dt + tmp[2];
        t.x[3] = t.x[3] + f2dt + tmp[3];

        t.x[0] = t.x[0] + delta_t*(old_x3) + tmp[0]; 
        t.x[1] = t.x[1] + delta_t*(old_x4) + tmp[1];

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

void System::get_obs_variance(State& s, double* var)
{
    for(int i=0; i<NUM_DIM_OBS; i++)
    {   
        var[i] = obs_noise[i];
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
    
    double range = sqrt(s.x[0]*s.x[0] + s.x[1]*s.x[1]);
    double theta = atan2(s.x[1], s.x[0]);
    t.x[0] = range + tmp[0];
    t.x[1] = theta + tmp[1];

    delete[] mean;
    delete[] tmp;

    return t;
}

void System::get_kalman_path( vector<State>& obs, vector<double>& obs_times, list<State>& kalman_path, list<State>& kalman_covar)
{
#if 1
    kalman_path.clear();

    kalman_path.push_back(init_state);
    
    Matrix4d Q = Matrix4d::Zero();
    for(int i=0; i<4; i++)
        Q(i,i) = init_var[i];

    Vector4d curr_state;
    curr_state(0) = init_state.x[0];
    curr_state(1) = init_state.x[1];
    curr_state(2) = init_state.x[2];
    curr_state(3) = init_state.x[3];
    
    Matrix4d Rk = Matrix4d::Zero();
    for(int i=0; i<4; i++)
        Rk(i,i) = obs_noise[i];

    double prev_time = 0;
    for(unsigned int i=0; i< obs.size(); i++)
    {
        State& next_obs = obs[i];
        Vector4d noisy_obs;
        noisy_obs(0) = next_obs.x[0];
        noisy_obs(1) = next_obs.x[1];
        noisy_obs(2) = next_obs.x[2];
        noisy_obs(3) = next_obs.x[3];

        double delta_t = obs_times[i] - prev_time;
        
        //cout<<"delta_t: "<< delta_t << endl;
        State stmp1;
        stmp1.x[0] = curr_state(0);
        stmp1.x[1] = curr_state(1);
        stmp1.x[2] = curr_state(2);
        stmp1.x[3] = curr_state(3);
        State next_state = integrate(stmp1, delta_t, true);
        State clean_obs = observation(next_state, true);
        Vector4d obs_vector;
        obs_vector(0) = clean_obs.x[0];
        obs_vector(1) = clean_obs.x[1];
        obs_vector(2) = clean_obs.x[2];
        obs_vector(3) = clean_obs.x[3];

        // Linearize Ad
        Matrix4d Ad = Matrix4d::Zero();
        Ad(0,2) = 1; Ad(1,3) = 1;
        if( ( sqrt(stmp1.x[0]*stmp1.x[0] + stmp1.x[1]*stmp1.x[1]) >= 9) && \
                ( (stmp1.x[0]*stmp1.x[2] + stmp1.x[1]*stmp1.x[3]) >= 0) )
        {
            double den = pow(stmp1.x[0]*stmp1.x[0] + stmp1.x[1]*stmp1.x[1],1.5);
            Ad(2,1) = -50*stmp1.x[1]*stmp1.x[1]/den;
            Ad(2,2) = -50*(-stmp1.x[0]*stmp1.x[1])/den;
            Ad(3,1) = -50*(-stmp1.x[0]*stmp1.x[1])/den;
            Ad(3,2) = -50*stmp1.x[0]*stmp1.x[0]/den;

        }
        // Linearize Cd
        Matrix4d Cd = Matrix4d::Zero();
        double range = sqrt(next_state.x[0]*next_state.x[0] + next_state.x[1]*next_state.x[1]);
        Cd(0,0) = next_state.x[0]/range; Cd(0,1) = next_state.x[1]/range;
        Cd(1,0) = -next_state.x[1]/range/range;
        Cd(1,1) = next_state.x[0]/range/range;
        
        Matrix4d Wk = Matrix4d::Zero();
        for(int j=0; j<4; j++)
            Wk(j,j) = process_noise[j]*delta_t;

        Matrix4d Q_new = Ad * Q * Ad.transpose() + Wk;
        Matrix4d Lk = Q_new*Cd.transpose()*(Cd*Q_new*Cd.transpose() + Rk).inverse();
        
        Vector4d Sk;
        curr_state(0) = next_state.x[0];
        curr_state(1) = next_state.x[1];
        curr_state(2) = next_state.x[2];
        curr_state(3) = next_state.x[3];
        Sk = noisy_obs - obs_vector;
        
        Vector4d estimate = curr_state + Lk*Sk;
        
        Matrix4d covar = (Matrix4d::Identity() - Lk*Cd)*Q_new;

        Q = covar;
        curr_state = estimate;
        
        //cout<< "Q: "<< Q << endl << endl;

        State stmp2;
        stmp2.x[0] = curr_state(0);
        stmp2.x[1] = curr_state(1);
        stmp2.x[2] = curr_state(2);
        stmp2.x[3] = curr_state(3);
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

        /*
           for(int i=0; i< num_particles; i++)
           {
           cout<<"weight: "<< weights[i]<< " state: ";
           for(int j=0; j<NUM_DIM; j++)
           {
           cout<<parts[i].x[j]<<" ";
           }
           cout<<endl;
           }
           cout<<"press key: "; getchar(); cout<<endl;
           */


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

