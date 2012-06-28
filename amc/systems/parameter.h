#ifndef __parameter_h__
#define __parameter_h__

#include "../utils/common.h"
#define NUM_DIM         (2)
#define NUM_DIM_OBS     (1)
// no time in this algorithm

class State
{
    public:
        double x[NUM_DIM];

        State()
        {
            for(int i=0; i<NUM_DIM; i++)
                x[i] = 0;
        }
        State(double *val)
        {
            for(int i=0; i<NUM_DIM; i++)
                x[i] = val[i];
        }
        ~State()
        {
        }

        double operator[](int which_dim)
        {
            assert(which_dim < NUM_DIM);
            return x[which_dim];
        }
        
        double norm2()
        {
            double sum = 0;
            for(int i=0; i< NUM_DIM; i++)
                sum += (x[i]*x[i]);

            return sum;
        }
        double norm()
        {
            double sum = 0;
            for(int i=0; i< NUM_DIM; i++)
                sum += (x[i]*x[i]);

            return sqrt(sum);
        }


        double operator*(const State& s1)
        {
            double ret = 0;
            for(int i=0; i< NUM_DIM; i++)
                ret += (this->x[i]*s1.x[i]);

            return ret;
        }
        State operator+(const State& s1)
        {
            State ret;
            for(int i=0; i< NUM_DIM; i++)
                ret.x[i] = (this->x[i] + s1.x[i]);

            return ret;
        }
        State operator-(const State& s1)
        {
            State ret;
            for(int i=0; i< NUM_DIM; i++)
                ret.x[i] = (this->x[i] - s1.x[i]);

            return ret;
        }
        State& operator=(const State &that)
        {
            if(this != &that)
            {
                for(int i=0; i< NUM_DIM; i++)
                    x[i] = that.x[i];
                
                return *this;
            }
            else
                return *this;
        }
};

class System
{
    public:

        double *obs_noise;
        double *process_noise;
        double *init_var;

        double *min_states;
        double *max_states;
        double sim_time_delta;
        
        State init_state;
        State init_state_real;

        System();
        ~System();

        // functions
        
        double get_holding_time(State& s, double gamma, int num_vert)
        {
            double h = gamma * pow( log(num_vert)/(num_vert), 1.0/(double)NUM_DIM);
            double num = h*h;

            num = num*sq(max_states[0] - min_states[0]);
            
            double sqnum = sqrt(num);
            double den = 0;
            for(int i=0; i< NUM_DIM; i++)
                den += process_noise[i];
           
            State f;
            f.x[0] = cos(2*M_PI*s.x[1]*s.x[0]);
            f.x[1] = 0;
            den += (sqnum*f.norm());
            
            return num/(den);
        }
        
        double get_min_holding_time(double gamma, int num_vert)
        {
            return -1;
        }


        int get_key(State& s, double *key);
        bool is_free(State &s);
        State sample();
        State get_fdt(State& s, double duration);
        State integrate(State& s, double duration, bool is_clean, bool is_propagate=false);
        void get_variance(State& s, double duration, double* var)
        {
            for(int i=0; i<NUM_DIM; i++)
            {   
                var[i] = process_noise[i]*duration;
            } 
        }

        void get_obs_variance(State& s, double* var)
        {
            for(int i=0; i<NUM_DIM_OBS; i++)
            {   
                var[i] = obs_noise[i];
            } 
        }

        State observation(State& s, bool is_clean)
        {
            State t;

            double *tmp = new double[NUM_DIM_OBS];
            double *mean = new double[NUM_DIM_OBS];

            for(int i=0; i< NUM_DIM_OBS; i++)
            {
                mean[i] = 0;
                tmp[i] = 0;
            }

            if( !is_clean)
            {
                multivar_normal( mean, obs_noise, tmp, NUM_DIM_OBS);
            }
            else
            {
                for(int i=0; i<NUM_DIM_OBS; i++)
                {
                    tmp[i] = 0;
                }
            }

            for(int i=0; i<NUM_DIM_OBS; i++)
                t.x[i] = 0;

            for(int i=0; i<NUM_DIM_OBS; i++)
            {
                t.x[i] = s.x[i] + tmp[i];
            }

            delete[] mean;
            delete[] tmp;

            return t;
        }
        void get_kalman_path(vector<State>& obs, vector<double>& obs_times, list<State>& kalman_path, list<State>& kalman_covar);
        void get_pf_path( vector<State>& obs, vector<double>& obs_times, list<State>& pf_path, int nparticles);

        double dist(State s1, State s2)
        {
            double t = 0;
            for(int i=0; i<NUM_DIM; i++)
                t += (s1.x[i] - s2.x[i])*(s1.x[i] - s2.x[i]);

            return sqrt(t);
        };

};


#endif
