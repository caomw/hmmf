#include "hmmf.h"

State get_mean(Graph& graph)
{
    State mstate;
    float mean[NUM_DIM] = {0};
    float totprob = 0;
    for(unsigned int i=0; i<graph.num_vert; i++)
    {
        for(int j=0; j< NUM_DIM; j++)
            mean[j] += ((graph.vlist[i]->s.x[j]) * (graph.vlist[i]->prob_best_path));

        if(graph.vlist[i]->prob_best_path != graph.vlist[i]->prob_best_path)
        {
            cout<<"found prob_best_path nan" << endl;
        }
        totprob += graph.vlist[i]->prob_best_path;
    }
    cout << "mean: ";
    for(int j=0; j< NUM_DIM; j++)
    {
        mstate.x[j] = mean[j]/totprob;

        if(mstate.x[j] != mstate.x[j])
        {
            cout<<"found a nan: " << mstate.x[j]<<" "<<totprob << endl;
        }

        cout<< mstate.x[j] << " ";
    }
    cout<<endl;

    return mstate;
}

int do_batch()
{
    System sys;
    Graph graph(sys);
    
    double start = get_msec();
   
#if 1
    tic();
    for(int i=0; i < 10000; i++)
    {
        graph.add_sample();
    }
    for(unsigned int i=0; i < graph.num_vert; i++)
    {
        Vertex* v = graph.vlist[i];
        graph.connect_edges_approx(v);
        if(i %1000 == 0)
        {
           cout<<"i: "<< i << endl;
           toc();
        }
    }

    graph.make_holding_time_constant();
    graph.system->sim_time_delta = graph.min_holding_time;

    for(unsigned int i=0; i< graph.num_vert; i++)
    {
        Vertex* v = graph.vlist[i];
        v->prob_best_path = normal_val(graph.system->init_state.x, graph.system->init_var,\
                v->s.x, NUM_DIM);
    }
    graph.normalize_density();
#endif

    graph.propagate_system();
    graph.get_kalman_path();

#if 1
    tic();
    graph.best_path.clear();
    //graph.best_path.push_back(get_mean(graph));
    for(unsigned int i=0; i < graph.obs.size(); i++)
    {
        graph.obs_curr_index = i;
        cout<< "obs_times: " << graph.obs_times[graph.obs_curr_index] << " ";
       
        for(unsigned int j = 0; j< graph.num_vert; j++)
        {
            Vertex* v = graph.vlist[j];
            graph.update_density_explicit(v);
        }
        for(unsigned int j = 0; j< graph.num_vert; j++)
        {
            Vertex* v = graph.vlist[j];
            v->prob_best_path = v->prob_best_path_buffer;
        }
        graph.normalize_density();
        
        cout<<i <<": ";
        graph.best_path.push_back(get_mean(graph));
    }
#endif

#if 0
    cout<<"starting simulation of trajectories" << endl;
    for(int i=0; i< 10; i++)
    {
        graph.simulate_trajectory_explicit();
    }
    graph.plot_monte_carlo_trajectories();
#endif

    graph.plot_trajectory();
    graph.plot_graph();
    cout<<"time: "<< get_msec() - start << " [ms]" << endl;

    return 0;
}

int do_incremental()
{
    System sys;
    Graph graph(sys);
    
    double start = get_msec();

#if 1 
    tic();
    for(int i=0; i < 5000; i++)
    {
        Vertex* v = graph.add_sample();
        graph.connect_edges_approx(v);
        graph.reconnect_edges_neighbors(v);
        if(i %1000 == 0)
        {
           cout<<"i: "<< i << endl;
           toc();
        }
    }
    
    for(unsigned int i=0; i< graph.num_vert; i++)
    {
        Vertex* v = graph.vlist[i];
        v->prob_best_path = normal_val(graph.system->init_state.x, graph.system->init_var,\
                v->s.x, NUM_DIM);
    }
    graph.normalize_density();
#endif

    graph.propagate_system();
    graph.get_kalman_path();

#if 1
    tic();
    graph.seeding_finished = true;
    graph.best_path.clear();
    graph.best_path.push_back(get_mean(graph));
    for(unsigned int i=0; i < 1; i++)
    {
        for(int j=0; j< 10; j++)
        {
            Vertex* v = graph.add_sample();
            graph.connect_edges_approx(v);
            graph.update_density_implicit(v);
            v->prob_best_path = v->prob_best_path_buffer;
            
            graph.reconnect_edges_neighbors(v);
            
        }

        graph.obs_curr_index = i;
        cout<< "obs_times: " << graph.obs_times[graph.obs_curr_index] << " ";
      
        graph.normalize_density();
        for(unsigned int j = 0; j< graph.num_vert; j++)
        {
            Vertex* v = graph.vlist[j];
            graph.update_density_implicit(v);
        }
        for(unsigned int j = 0; j< graph.num_vert; j++)
        {
            Vertex* v = graph.vlist[j];
            v->prob_best_path = v->prob_best_path_buffer;
        }
        graph.normalize_density();
        
        cout<<i <<": ";
        graph.best_path.push_back(get_mean(graph));
    }
#endif

#if 0
    cout<<"starting simulation of trajectories" << endl;
    for(int i=0; i< 100; i++)
    {
        graph.simulate_trajectory_implicit();
    }
    graph.plot_monte_carlo_trajectories();
#endif

    graph.plot_trajectory();
    graph.plot_graph();
    cout<<"time: "<< get_msec() - start << " [ms]" << endl;

    return 0;
}

int main()
{
    srand(0);

    //do_batch();
    do_incremental();

    return 0;
}

