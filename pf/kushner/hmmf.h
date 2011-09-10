#ifndef __hmmf_h__
#define __hmmf_h__

#include "utils/common.h"
#include "systems/ship.h"

class Edge;
class Vertex;
class Graph;

class Vertex
{
    public:

        State s;
        
        // prob of best path that ends up here incorporating obs
        float prob_best_path;
        float prob_best_path_buffer;
        float obs_update_time;
        
        float holding_time;
        float holding_time_delta;

        // parent of the best path
        Vertex *prev;

        list<Edge *> edges_in;
        list<Edge *> edges_out;
        
        float self_transition_prob;

        Vertex(State& st);
        ~Vertex(){};
        
        friend class System;
};

class Edge{

    public:
        Vertex *from;
        Vertex *to;
    
        list<Edge*>::iterator from_iter;
        list<Edge*>::iterator to_iter;
        list<Edge*>::iterator elist_iter;

        float transition_prob_delta;
        float transition_prob;
        float transition_time;
        
        Edge(Vertex* f, Vertex* t, float prob, float trans_time);

        Edge reverse(){
            return Edge(this->to, this->from, this->transition_prob, this->transition_time);
        }
        ~Edge()
        {
            //cout<<"called destructor"<<endl;

        };
};

class Graph{

    public:
        
        int obs_interval;
        float max_obs_time;
        float delta;

        float gamma, gamma_t;
        struct kdtree *state_tree;
       
        System* system;

        Graph(System& sys);
        ~Graph();
        
        vector<Vertex *> vlist;
        list<Edge *> elist;
        
        unsigned int num_vert;
        list<State> truth;
        int obs_curr_index;
        vector<float> obs_times;
        vector<State> obs;
        list<State> best_path;
        list<State> kalman_path;
        
        // graph sanity check
        list< list<State> > monte_carlo_trajectories;
        list<float> monte_carlo_probabilities;
        list< list<float> > monte_carlo_times;

        // graph functions
        unsigned int get_num_vert(){return num_vert; };

        void remove_vertex(Vertex* v);
        int vertex_delete_edges(Vertex* v);
        void remove_edge(Edge *e);
        
        int insert_into_kdtree(Vertex *v);
        Vertex* nearest_vertex(State s);
        void normalize_edges(Vertex *from);
        
        void print_rrg();
        void plot_graph();
        void plot_trajectory();
        void plot_monte_carlo_trajectories(); 

        // algorithm functions
        
        void iterate();
        Vertex* add_sample();
        bool is_edge_free( Edge *etmp);
        
        int reconnect_edges_neighbors(Vertex* v);
        int connect_edges(Vertex *v);
        int connect_edges_approx(Vertex *v);

        void propagate_system();
        
        void put_init_samples(int howmany);
       
        void normalize_density();
        void propagate_density(Vertex* v);
        void update_density_explicit(Vertex* v);
        void update_density_implicit(Vertex* v);
       
        int calculate_delta();
        int calculate_probabilities_delta(Vertex* v);
        int calculate_probabilities_delta_all();
        float make_holding_time_constant();

        void propagate_viterbi(Vertex* v);
        void update_viterbi(Vertex* v);
        void update_viterbi_neighbors(Vertex* v);
        
        void update_observation_prob(State& yt);
       

        void get_kalman_path();
        
        bool is_everything_normalized();
        int simulate_trajectory_explicit();
        int simulate_trajectory_implicit();
        
        friend class System;
};


#endif
