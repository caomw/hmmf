#ifndef __singleint_h__
#define __singleint_h__

#include "common.h"

class State
{
    public:
        double x[NUM_DIM];

        State()
        {
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
        
        int get_key(State& s, double *key);
        bool is_free(State &s);
        State sample();
        State integrate(State& s, double duration, bool is_clean);
        State observation(State& s, bool is_clean);
};

inline double dist(State s1, State s2);




#endif
