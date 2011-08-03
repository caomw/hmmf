#ifndef __hmmf_h__
#define __hmmf_h__

#include "common.h"
#include "singleint.h"

class Edge;
class Vertex;
class Graph;

class Vertex
{
    public:

        State s;
        
        // prob of best path that ends up here incorporating obs
        double prob_best_path;
        
        bool is_open_queue;
        bool is_closed_queue;

        // parent of the best path
        Vertex *prev;
        Vertex *next;

        Edge *best_in;
        Edge *best_out;

        list<Edge *> edges_in;
        list<Edge *> edges_out;
        
        Vertex(State& st);
        ~Vertex(){};
        
        friend class System;
};

class Edge{

    public:
        Vertex *from;
        Vertex *to;
        double transition_prob;
        double transition_time;
        
        Edge(Vertex* f, Vertex* t, double pro=1.0);

        Edge reverse(){
            return Edge(this->to, this->from, this->transition_prob);
        }
        ~Edge()
        {
            //cout<<"called destructor"<<endl;

        };
};

class Graph{

    private:
        System* system;
        int num_particles; 

        double gamma;
        struct kdtree *state_tree;
        struct kdtree *time_tree; 
        
    public:

        Graph(System& sys);
        ~Graph();
        
        vector<Vertex *> vlist;
        unsigned int num_vert;
        int samples_per_obs;
        int num_observations;
        list<State> truth;
        list<State> obs;
        list<State> best_path;
        list<State> kalman_path;
        
        // graph functions
        unsigned int get_num_vert(){return num_vert; };

        void add_vertex(Vertex *v){
            vlist.push_back(v);
            num_vert++;
        }
        void remove_vertex(Vertex* v);
        void remove_edge(Edge *e);
        
        int insert_into_kdtree(Vertex *v);
        Vertex* nearest_vertex(State s);
        void normalize_edges(Vertex *from);
        
        void print_rrg();
        void plot_graph();
        void plot_trajectory();
        

        // algorithm functions
        
        void iterate();
        void add_sample();
        bool is_edge_free( Edge *etmp);
        void connect_edges(Vertex *v);
        void propagate_system();
        
        void put_init_samples();
        int write_transition_prob(Edge *e);
        int write_observation_prob(Edge *e, State& obs);
        void do_viterbi( Vertex *v );
        void update_observation_prob(State& yt);
        void get_best_path();
         
        friend class System;
};


#endif
