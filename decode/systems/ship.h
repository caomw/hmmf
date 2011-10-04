#ifndef __ship_h__
#define __ship_h__

#include "../utils/common.h"
#define NUM_DIM         (4)
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
        
        int num_particles;

        State init_state;

        System();
        ~System();

        // functions
        
        double get_holding_time(State& s, double gamma, int num_vert)
        {
            double h = gamma * pow( log(num_vert)/(num_vert), 1.0/(double)NUM_DIM);
            double num = h*h;

            //for(int i=0; i< NUM_DIM; i++)
            //    num = num*(max_states[i] - min_states[i]);
            
            double sqnum = sqrt(num);
            double den = 0;
            for(int i=0; i< NUM_DIM; i++)
                den += process_noise[i];
           
            State f;
            f.x[0] = s.x[2]; f.x[1] = s.x[3]; f.x[2] = 0; f.x[3] = 0;
            if( ( sqrt(s.x[0]*s.x[0] + s.x[1]*s.x[1]) >= 9) && ( (s.x[0]*s.x[2] + s.x[1]*s.x[3]) >= 0) )
            {
                f.x[2] = -50*s.x[0]/sqrt(s.x[0]*s.x[0] + s.x[1]*s.x[1]);
                f.x[3] = -50*s.x[1]/sqrt(s.x[0]*s.x[0] + s.x[1]*s.x[1]);
            }
            den += (sqnum*f.norm());
            
            return num/(den);
        }
        
        double get_min_holding_time(double gamma, int num_vert)
        {
            double h = gamma * pow( log(num_vert)/(num_vert), 1.0/(double)NUM_DIM);
            double num = h*h;
            //for(int i=0; i< NUM_DIM; i++)
            //    num = num*(max_states[i] - min_states[i]);
            
            double sqnum = sqrt(num);

            double den = 0;
            for(int i=0; i< NUM_DIM; i++)
            {
                den += process_noise[i];
            }
            double max_abs_f = sqrt(max_states[2]*max_states[2] + max_states[3]*max_states[3] + 2*25);
            den += (sqnum*max_abs_f);
            
            return num/(den);
        }


        int get_key(State& s, double *key);
        bool is_free(State &s);
        State sample();
        State integrate(State& s, double duration, bool is_clean);
        void get_variance(State& s, double duration, double* var);
        void get_obs_variance(State& s, double* var);
        
        State observation(State& s, bool is_clean);

        void get_kalman_path(vector<State>& obs, vector<double>& obs_times, list<State>& kalman_path, list<State>& kalman_covar);
        void get_pf_path( vector<State>& obs, vector<double>& obs_times, list<State>& pf_path);
        
        double dist(State s1, State s2)
        {
            double t = 0;
            for(int i=0; i<NUM_DIM; i++)
                t += (s1.x[i] - s2.x[i])*(s1.x[i] - s2.x[i]);

            return sqrt(t);
        };
};



#endif
