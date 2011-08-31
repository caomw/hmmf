#include "hmmf.h"

State get_mean(Graph& graph)
{
    State mstate;
    float mean[NUM_DIM] = {0};
    float totprob = 0;
    for(int i=0; i<graph.num_vert; i++)
    {
        for(int j=0; j< NUM_DIM; j++)
            mean[j] += ((graph.vlist[i]->s.x[j]) * (graph.vlist[i]->prob_best_path));

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
    for(int i=0; i < graph.num_vert; i++)
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
#endif
#if 0
    tic();
    for(int i=0; i < 10000; i++)
    {
       Vertex* v = graph.add_sample();
       graph.connect_edges_approx(v);
       graph.reconnect_edges_neighbors(v);
       if(i % 1000 == 0)
       {
           cout<<"i: "<< i << endl;
           toc();
       }
    }
#endif

    for(int i=0; i< graph.num_vert; i++)
    {
        Vertex* v = graph.vlist[i];
        v->prob_best_path = normal_val(graph.system->init_state.x, graph.system->init_var,\
                v->s.x, NUM_DIM);
    }
    
    graph.system->sim_time_delta = graph.system->get_min_holding_time(graph.gamma, graph.num_vert);
    graph.propagate_system();
    //graph.get_kalman_path();

#if 1
    tic();
    graph.best_path.clear();
    graph.best_path.push_back(get_mean(graph));
    for(int i=0; i < graph.obs.size(); i++)
    {
        /*
        for(int j=0; j< 200; j++)
        {
            Vertex* v = graph.add_sample();
            graph.connect_edges_approx(v);
            graph.reconnect_edges_neighbors(v);
            //graph.update_viterbi(v);
            graph.update_viterbi_neighbors(v);
        }
        */
        graph.obs_curr_index = i;
        cout<< "obs_times: " << graph.obs_times[graph.obs_curr_index] << " ";
       
        for(int j = 0; j< graph.num_vert; j++)
        {
            Vertex* v = graph.vlist[j];
            graph.update_density(v);
        }
        for(int j = 0; j< graph.num_vert; j++)
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
    for(int i=0; i< 1000; i++)
    {
        graph.simulate_trajectory();
    }
    graph.plot_monte_carlo_trajectories();
#endif

    graph.plot_trajectory();
    graph.plot_graph();
    cout<<"time: "<< get_msec() - start << " [ms]" << endl;
}

int graph_sanity_check()
{
    System sys;
    Graph graph(sys);
    
    double start = get_msec();

    tic();
    
    for(int i=0; i < 1000; i++)
    {
        graph.add_sample();
        if(i % 500 == 0)
        {
            cout<< i << endl;
            toc();
        }
    }
    
    // make the holding time constant
    graph.plot_graph();
    cout<<"added samples: "<< graph.num_vert << endl;

    assert(graph.is_everything_normalized());
    cout<<"graph is correct"<<endl;
    
    for(int i=0; i< 100; i++)
    {
        graph.simulate_trajectory();
        //cout<<i<<endl;
    }
    
    graph.propagate_system();
    graph.get_kalman_path();
    graph.plot_trajectory();
    graph.plot_monte_carlo_trajectories();
    
    cout<<"time: "<< get_msec() - start << " [ms]" << endl;
}

int main()
{
    srand(0);

    //graph_sanity_check();
    do_batch();

    return 0;
}

