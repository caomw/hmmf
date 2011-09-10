#ifndef __dubins_h__
#define __dubins_h__

#include "../utils/common.h"

#define NUM_DIM         (3)
#define NUM_DIM_OBS     (3)

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
            double num = h*h;

            double den = 0;
            for(int i=0; i< NUM_DIM; i++)
                den += process_noise[i];
            
            den += (h*sqrt(5));
            // |f| = sqrt(5) always
            
            return num/(den);
        }
        double get_min_holding_time(double gamma, int num_vert)
        {
            double h = gamma * pow( log(num_vert)/(num_vert), 1.0/(double)NUM_DIM);
            double num = h*h;

            double den = 0;
            for(int i=0; i< NUM_DIM; i++)
                den += process_noise[i];
            
            den += (h*sqrt(5));
            // |f| = sqrt(5) always
            
            return num/(den);
        }
       
        int get_key(State& s, double *key);
        bool is_free(State &s);
        State sample();
        State integrate(State& s, double duration, bool is_clean);
        void get_variance(State& s, double duration, double* var);
        State observation(State& s, bool is_clean);
        
        void get_kalman_path(vector<State>& obs, vector<double>& obs_times, list<State>& kalman_path);
};

inline double dist(State& s1, State& s2)
{
    double t = 0;
    for(int i=0; i<NUM_DIM; i++)
        t += (s1.x[i] - s2.x[i])*(s1.x[i] - s2.x[i]);
    
    //cout<<"dist: "<< t << endl;
    return sqrt(t);
};
#endif
