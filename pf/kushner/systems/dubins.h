#ifndef __dubins_h__
#define __dubins_h__

#include "../utils/common.h"

#define NUM_DIM         (3)
#define NUM_DIM_OBS     (3)

class State
{
    public:
        float x[NUM_DIM];

        State()
        {
            for(int i=0; i<NUM_DIM; i++)
                x[i] = 0;
        }
        State(float *val)
        {
            for(int i=0; i<NUM_DIM; i++)
                x[i] = val[i];
        }
        ~State()
        {
        }

        float operator[](int which_dim)
        {
            assert(which_dim < NUM_DIM);
            return x[which_dim];
        }
        
        float norm2()
        {
            float sum = 0;
            for(int i=0; i< NUM_DIM; i++)
                sum += (x[i]*x[i]);

            return sum;
        }
        float norm()
        {
            float sum = 0;
            for(int i=0; i< NUM_DIM; i++)
                sum += (x[i]*x[i]);

            return sqrt(sum);
        }


        float operator*(const State& s1)
        {
            float ret = 0;
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

        float *obs_noise;
        float *process_noise;
        float *init_var;

        float *min_states;
        float *max_states;
        float sim_time_delta;

        State init_state;

        System();
        ~System();

        // functions
 
        float get_holding_time(State& s, float gamma, int num_vert)
        {
            float h = gamma * pow( log(num_vert)/(num_vert), 1.0/(float)NUM_DIM);
            float num = h*h;

            float den = 0;
            for(int i=0; i< NUM_DIM; i++)
                den += process_noise[i];
            
            den += (h*sqrt(5));
            // |f| = sqrt(5) always
            
            return num/(den);
        }
        float get_min_holding_time(float gamma, int num_vert)
        {
            float h = gamma * pow( log(num_vert)/(num_vert), 1.0/(float)NUM_DIM);
            float num = h*h;

            float den = 0;
            for(int i=0; i< NUM_DIM; i++)
                den += process_noise[i];
            
            den += (h*sqrt(5));
            // |f| = sqrt(5) always
            
            return num/(den);
        }
       
        int get_key(State& s, double *key);
        bool is_free(State &s);
        State sample();
        State integrate(State& s, float duration, bool is_clean);
        void get_variance(State& s, float duration, float* var);
        State observation(State& s, bool is_clean);
        
        void get_kalman_path(vector<State>& obs, vector<float>& obs_times, list<State>& kalman_path);
};

inline float dist(State& s1, State& s2)
{
    float t = 0;
    for(int i=0; i<NUM_DIM; i++)
        t += (s1.x[i] - s2.x[i])*(s1.x[i] - s2.x[i]);
    
    //cout<<"dist: "<< t << endl;
    return sqrt(t);
};
#endif
