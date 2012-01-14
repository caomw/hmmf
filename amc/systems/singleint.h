#ifndef __singleint_h__
#define __singleint_h__

#include "../utils/common.h"
#define NUM_DIM         (2)
#define NUM_DIM_OBS     (2)
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

        System();
        ~System();

        // functions
        
        double get_holding_time(State& s, double gamma, int num_vert)
        {
            double h = gamma * pow( log(num_vert)/(num_vert), 1.0/(double)NUM_DIM);
            h = h*(max_states[0] - min_states[0]);
            
            double num = h*h;

            double sqnum = sqrt(num);
            double den = process_noise[0];
            
            State absf = get_fdt(s, 1.0);
            den += (1*sqnum*absf.norm());
            return num/(den);
            
            /*
            return (sqrt(process_noise[0])*sqrt(4*s.norm() + process_noise[0]) + 
                    process_noise[0] + 2*h*s.norm())/2.0/s.norm2();
            */
        }
        
        double get_min_holding_time(double gamma, int num_vert)
        {
            double h = gamma * pow( log(num_vert)/(num_vert), 1.0/(double)NUM_DIM);
            double num = h*h;
            num = num*sq(max_states[0] - min_states[0]);

            double sqnum = sqrt(num);

            double den = 0;
            for(int i=0; i< NUM_DIM; i++)
                den += process_noise[i];
            
            den += (sqnum*max_states[0]);
            
            return num/(den);
        }


        int get_key(State& s, double *key);
        bool is_free(State &s);
        State sample();
        State get_fdt(State& s, double duration);
        State integrate(State& s, double duration, bool is_clean);
        void get_obs_variance(State& s, double* var);
        void get_variance(State& s, double duration, double* var);
        
        State observation(State& s, bool is_clean);

        void get_kalman_path(vector<State>& obs, vector<double>& obs_times, list<State>& kalman_path, list<State>& kalman_covar);
        void get_pf_path( vector<State>& obs, vector<double>& obs_times, list<State>& pf_path, int nparticles);
};



#endif
