#ifndef __filter_h__
#define __filter_h__

#include "utils/common.h"
#include "systems/singleint.h"

class Edge;
class Vertex;
class Graph;

class Vertex
{
    public:

        State s;
        
        // prob of best path that ends up here incorporating obs
        double prob_best_path;
        double prob_best_path_buffer;
        double obs_update_time;
         
        double holding_time;
        double holding_time_delta;

        // parent of the best path
        Vertex *prev;

        list<Edge *> edges_in;
        list<Edge *> edges_out;
        
        double self_transition_prob;

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

        double transition_prob_delta;
        double transition_prob;
        double transition_time;
        
        Edge(Vertex* f, Vertex* t, double prob, double trans_time);

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
        double max_obs_time;
        double delta;
        double min_holding_time;
        bool seeding_finished;

        double gamma, gamma_t;
        struct kdtree *state_tree;
       
        System* system;

        Graph(System& sys);
        ~Graph();
        
        vector<Vertex *> vlist;
        list<Edge *> elist;
        
        unsigned int num_vert;
        list<State> truth;
        int obs_curr_index;
        vector<double> obs_times;
        vector<State> obs;
        list<State> best_path;
        list<State> kalman_path;
        list<State> pf_path;
        list<State> kalman_covar;
        
        // graph sanity check
        list< list<State> > monte_carlo_trajectories;
        list< list<State> > actual_trajectories;
        list<double> monte_carlo_probabilities;
        list< list<double> > monte_carlo_times;

        // graph functions
        unsigned int get_num_vert(){return num_vert; };

        void remove_vertex(Vertex* v);
        int vertex_delete_edges(Vertex* v);
        void remove_edge(Edge *e);
        
        int insert_into_kdtree(Vertex *v);
        Vertex* nearest_vertex(State s);
        void normalize_edges(Vertex *from);
 
        double dist(State s1, State s2)
        {
            double t = 0;
            for(int i=0; i<NUM_DIM; i++)
                t = t + (s1.x[i] - s2.x[i])*(s1.x[i] - s2.x[i]);

            return sqrt(t);
        };       
        
        void print_rrg();
        void plot_graph();
        void plot_trajectory();
        void plot_monte_carlo_trajectories(); 
        void analyse_monte_carlo_trajectories();
        void plot_monte_carlo_density(char* filename);

        // algorithm functions
        
        void iterate();
        Vertex* add_sample(bool is_seed=false);
        bool is_edge_free( Edge *etmp);
        
        int reconnect_edges_neighbors(Vertex* v);
        int connect_edges(Vertex *v);
        int connect_edges_approx(Vertex *v);

        void propagate_system();
        
        void put_init_samples(int howmany);
       
        void average_density(Vertex* v);
        void approximate_density(Vertex* v);
        void normalize_density();
        void propagate_density(Vertex* v);
        
        void update_density_implicit_no_obs(Vertex* v);
        void update_density_implicit_no_obs_all();
        bool update_density_implicit(Vertex* v);
        void update_density_implicit_all();
        void buffer_prob_copy();
       
        int calculate_delta();
        int calculate_probabilities_delta(Vertex* v);
        int calculate_probabilities_delta_all();

        void update_observation_prob(State& yt);
       
        void get_kalman_path();
        void get_pf_path(int nparticles);
        
        bool is_everything_normalized();
        int simulate_trajectory_implicit();
        
        friend class System;
};


#endif
