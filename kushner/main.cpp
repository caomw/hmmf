#include "hmmf.h"

State get_mean(Graph& graph, bool is_cout=true)
{
    State mstate;
    double mean[NUM_DIM] = {0};
    double totprob = 0;
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
    if(is_cout)
        cout << "mean: ";
    for(int j=0; j< NUM_DIM; j++)
    {
        mstate.x[j] = mean[j]/totprob;

        if(mstate.x[j] != mstate.x[j])
        {
            cout<<"found a nan: " << mstate.x[j]<<" "<<totprob << endl;
        }
        
        if(is_cout)
            cout<< mstate.x[j] << " ";
    }
    if(is_cout)
        cout<<endl;

    return mstate;
}

void get_sq_error(Graph& graph, double& bpe, double& kfe)
{
    bpe = 0;
    kfe = 0;
    list<State>::iterator bpiter = graph.best_path.begin();
    list<State>::iterator kfiter = graph.kalman_path.begin();

    for(list<State>::iterator i = graph.truth.begin(); i!=graph.truth.end(); i++)
    {
        State& s1 = (*i);
        State& s2 = (*bpiter);
        State& s3 = (*kfiter);
        
        bpe += graph.dist(s1,s2)*graph.dist(s1,s2);
        kfe += graph.dist(s1,s3)*graph.dist(s1,s3);

        bpiter++;
        kfiter++;
    }
}

int do_batch(int tot_vert)
{
    System sys;
    Graph graph(sys);
    
    double start = get_msec();
   
#if 1
    tic();
    for(int i=0; i< 10; i++)
        graph.add_sample(true);

    for(int i=0; i < tot_vert; i++)
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
    graph.seeding_finished = true;
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
        cout<< "time: " << graph.obs_times[graph.obs_curr_index] << " ";
       
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

int do_incremental(int tot_vert)
{
    System sys;
    Graph graph(sys);
    
    double start = get_msec();

    graph.propagate_system();
    graph.get_kalman_path();

#if 1 
    tic();
    for(int i=0; i< 100; i++)
    {
        Vertex* v = graph.add_sample(true);
        graph.connect_edges_approx(v);
        graph.reconnect_edges_neighbors(v);
    }

    for(int i=0; i < 1000; i++)
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
    graph.seeding_finished = true;
   
#if 0
    // checking approximation
    cout<<"starting----" << endl;
    get_mean(graph);
    for(int j=0; j< 500; j++)
    {
        Vertex* v = graph.add_sample();
        graph.connect_edges_approx(v);
        
        graph.reconnect_edges_neighbors(v);
        
        graph.approximate_density(v);
        graph.normalize_density();
        
        get_mean(graph);
    }
#endif
#endif

#if 1
    int to_add = tot_vert/graph.obs.size();
    tic();
    graph.best_path.clear();
    graph.best_path.push_back(get_mean(graph, false));
    for(unsigned int i=0; i < graph.obs.size(); i++)
    {
        for(int j=0; j< to_add; j++)
        {
            Vertex* v = graph.add_sample();
            graph.connect_edges_approx(v);
            
            graph.reconnect_edges_neighbors(v);
            
            graph.approximate_density(v);
        }
        
        graph.obs_curr_index = i;
        cout<< "time: " << graph.obs_times[graph.obs_curr_index] << "\t";
        
        for(int j=0; j< (int)(graph.delta/graph.system->sim_time_delta); j++)
        {
            graph.update_density_implicit_no_obs_all();
            graph.normalize_density();
        }
        graph.update_density_implicit_all();
        graph.normalize_density();
        
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

#if 1 
    double bpe, kfe;
    get_sq_error(graph, bpe, kfe);
    cout<<"bpe: "<< bpe <<" kfe: "<< kfe << endl;
#endif
    
    return 0;
}

int do_timing_plot()
{
    int max_runs = 1;

    System sys;
    
    Graph g1(sys);
    g1.propagate_system();
    g1.get_kalman_path();
    double bpe1, kfe1;
    get_sq_error(g1, bpe1, kfe1);
    //cout<<"bpe1: "<<bpe1<<" " << " kfe1: "<< kfe1 << endl;
    
    for(int tot_vert=1000; tot_vert < 10000; tot_vert+= 5000)
    {
        double average_time = 0;
        double average_bpe = 0;
        double average_kfe = 0;
        
        for(int how_many=0; how_many < max_runs; how_many++)
        {
            Graph graph(sys);

            graph.truth = g1.truth;
            graph.obs = g1.obs;
            graph.kalman_path = g1.kalman_path;
            graph.obs_times = g1.obs_times;

            tic();
            for(int i=0; i< 10; i++)
            {
                Vertex* v = graph.add_sample(true);
                graph.connect_edges_approx(v);
                graph.reconnect_edges_neighbors(v);
            }

            for(int i=0; i < 100; i++)
            {
                Vertex* v = graph.add_sample();
                graph.connect_edges_approx(v);
                graph.reconnect_edges_neighbors(v);
            }

            for(unsigned int i=0; i< graph.num_vert; i++)
            {
                Vertex* v = graph.vlist[i];
                v->prob_best_path = normal_val(graph.system->init_state.x, graph.system->init_var,\
                        v->s.x, NUM_DIM);
            }
            graph.normalize_density();
            graph.seeding_finished = true;

#if 1
            int to_add = tot_vert/graph.obs.size();
            graph.best_path.clear();
            graph.best_path.push_back(get_mean(graph, false));
            for(unsigned int i=0; i < graph.obs.size(); i++)
            {
                for(int j=0; j< to_add; j++)
                {
                    Vertex* v = graph.add_sample();
                    graph.connect_edges_approx(v);

                    graph.reconnect_edges_neighbors(v);

                    graph.approximate_density(v);
                }

                graph.obs_curr_index = i;

                graph.update_density_implicit_all();
                graph.normalize_density();

                graph.best_path.push_back(get_mean(graph, false));
            }
#endif

            //graph.plot_trajectory();

            double bpe, kfe;
            get_sq_error(graph, bpe, kfe);

            average_time += toc();
            average_bpe += bpe;
            average_kfe += kfe;
        }

        average_time = average_time/(double)max_runs;
        average_bpe = average_bpe/(double)max_runs;
        average_kfe = average_kfe/(double)max_runs;

        cout<<tot_vert+110<<"\t"<<average_time<<"\t"<< average_bpe <<"\t"<< average_kfe << endl;
    }
    return 0;
}


int main(int argv, char* argc[])
{
    cout.precision(5);
    //do_batch(5000);
    do_incremental(20000);
    
    //do_timing_plot();
    
    return 0;
}

