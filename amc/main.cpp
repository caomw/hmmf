#include "filter.h"

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
        cout << "totprob: "<< totprob<<" mean: ";
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

void get_sq_error(Graph& graph, double& bpe, double& kfe, double& pfe)
{
    bpe = 0;
    kfe = 0;
    pfe = 0;

    list<State>::iterator bpiter = graph.best_path.begin();
    list<State>::iterator kfiter = graph.kalman_path.begin();
    list<State>::iterator pfiter = graph.pf_path.begin();

    for(list<State>::iterator i = graph.truth.begin(); i!=graph.truth.end(); i++)
    {
        State& s1 = (*i);
        State& s2 = (*bpiter);
        State& s3 = (*kfiter);
        State& s4 = (*pfiter);

        bpe = bpe + graph.dist(s1,s2)*graph.dist(s1,s2);
        kfe = kfe + graph.dist(s1,s3)*graph.dist(s1,s3);
        pfe = pfe + graph.dist(s1,s4)*graph.dist(s1,s4);

        bpiter++;
        kfiter++;
        pfiter++;
    }
}

int do_batch(int tot_vert)
{
    System sys;
    Graph graph(sys);

    graph.propagate_system();
    graph.get_kalman_path();
    tic();
    graph.get_pf_path(tot_vert);
    cout<<"pf time: "<< toc() << endl;

#if 1
    tic();
    for(int i=0; i < tot_vert; i++)
    {
        graph.add_sample(false);
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
    for(unsigned int i=0; i< graph.num_vert; i++)
    {
        Vertex* v = graph.vlist[i];
        v->prob_best_path = normal_val(graph.system->init_state.x, graph.system->init_var,\
                v->s.x, NUM_DIM);
    }
    graph.normalize_density();
    graph.seeding_finished = true;

    graph.make_chain_implicit();
     
#endif
    
#if 1
    cout<<"start filter"<<endl;
    tic();
    graph.best_path.clear();
    graph.best_path.push_back(get_mean(graph, false));
    for(unsigned int i=0; i < graph.obs.size(); i++)
    {
        graph.obs_curr_index = i;
        // cout<< "time: " << graph.obs_times[graph.obs_curr_index] << " ";
        graph.update_density_implicit_all(); 
        if(i%10 == 0)
            cout<<i<<endl;
        graph.best_path.push_back(get_mean(graph, false));
    }
    cout<<"filter time: "<< toc() << endl;
#endif

#if 0
    cout<<"starting simulation of trajectories" << endl;
    for(int i=0; i< 100; i++)
    {
        if((i%1000) == 0)
            cout<<i<<endl;
        graph.simulate_trajectory_implicit();
    }
    graph.plot_monte_carlo_trajectories();
#endif

    double bpe, kfe, pfe;
    get_sq_error(graph, bpe, kfe, pfe);
    cout<<"bpe: "<< bpe <<" kfe: "<< kfe << " pfe: "<< pfe << endl;
    
    graph.plot_trajectory();
    //graph.plot_graph();

    return 0;
}

int do_incremental(int tot_vert)
{
    System sys;
    Graph graph(sys);

    double start = get_msec();

    graph.propagate_system();
    graph.get_kalman_path();
    graph.get_pf_path(1000);

#if 1 
    tic();
    for(int i=0; i< 100; i++)
    {
        Vertex* v = graph.add_sample(true);
        graph.connect_edges_approx(v);
        graph.reconnect_edges_neighbors(v);
    }

    for(int i=0; i < tot_vert-100; i++)
    {
        Vertex* v = graph.add_sample();
        graph.connect_edges_approx(v);
        graph.reconnect_edges_neighbors(v);
        if(i %1000 == 0)
        {
            //cout<<"i: "<< i << endl;
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

#if 1
    cout<<"starting filtering"<<endl;
    int to_add = tot_vert/graph.obs.size();
    tic();
    graph.best_path.clear();
    graph.best_path.push_back(get_mean(graph, false));
    for(unsigned int i=0; i < graph.obs.size(); i++)
    {   
        /*
        for(int j=0; j< to_add; j++)
        {
            Vertex* v = graph.add_sample();
            graph.connect_edges_approx(v);

            graph.reconnect_edges_neighbors(v);

            graph.approximate_density(v);
        }
        */
        graph.obs_curr_index = i;
        //cout<< "time: " << graph.obs_times[graph.obs_curr_index] << "\t";

        graph.update_density_implicit_all();
        graph.normalize_density();

        graph.best_path.push_back(get_mean(graph, false));
    }
#endif
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
    double bpe, kfe, pfe;
    get_sq_error(graph, bpe, kfe, pfe);
    cout<<"bpe: "<< bpe <<" kfe: "<< kfe << " pfe: "<< pfe << endl;
#endif

    return 0;
}

int do_error_plot()
{
    ofstream err_out("data/err_out.dat");
    int max_runs = 10;

    System sys;

    for(int tot_vert=10; tot_vert < 200; tot_vert+= 10)
    {
        double average_time_hmm = 0;
        double average_time_pf = 0;
        double average_bpe = 0;
        double average_kfe = 0;
        double average_pfe = 0;

        Graph graph(sys);
        graph.propagate_system();

        for(int how_many=0; how_many < max_runs; how_many++)
        {
            tic();
            for(int i=0; i < tot_vert; i++)
                graph.add_sample();
            for(unsigned int i=0; i< graph.num_vert; i++)
            {
                Vertex* v = graph.vlist[i];
                graph.connect_edges_approx(v);
            }
            //average_time_hmm += toc();

            graph.get_kalman_path();
            tic();
            graph.get_pf_path(tot_vert);
            average_time_pf += toc();
            
            tic();
            for(unsigned int i=0; i< graph.num_vert; i++)
            {
                Vertex* v = graph.vlist[i];
                v->prob_best_path = normal_val(graph.system->init_state.x, graph.system->init_var,\
                        v->s.x, NUM_DIM);
            }
            graph.normalize_density();
            graph.seeding_finished = true;

            graph.best_path.clear();
            graph.best_path.push_back(get_mean(graph, false));
            for(unsigned int i=0; i < graph.obs.size(); i++)
            {
                graph.obs_curr_index = i;
                graph.update_density_implicit_all();
                graph.best_path.push_back(get_mean(graph, false));
            }
            average_time_hmm += toc();

            //graph.plot_trajectory();

            double bpe, kfe, pfe;
            get_sq_error(graph, bpe, kfe, pfe);

            average_bpe += bpe;
            average_kfe += kfe;
            average_pfe += pfe;
        }

        average_time_hmm = average_time_hmm/(double)max_runs;
        average_time_pf = average_time_pf/(double)max_runs;
        average_bpe = average_bpe/(double)max_runs;
        average_kfe = average_kfe/(double)max_runs;
        average_pfe = average_pfe/(double)max_runs;

        cout<<tot_vert<<"\t"<<average_time_hmm<<"\t"<<average_time_pf<<"\t"<< average_bpe <<"\t"<< average_kfe << "\t"<<average_pfe<<endl;
        err_out<<tot_vert<<"\t"<<average_time_hmm<<"\t"<<average_time_pf<<"\t"<< average_bpe <<"\t"<< average_kfe << "\t"<<average_pfe<<endl;
        err_out.flush();
    }
    err_out.close();
    return 0;
}

int do_timing_plot()
{
    ofstream timing_out("data/timing_out.dat");

    System sys;

#if 0
    // batch
    for(int n=20000; n < 75000; n = n+1000)
    {
        Graph graph(sys);
        tic();
        for(int i=0; i < n; i++)
            graph.add_sample();
        
        for(unsigned int i=0; i < graph.num_vert; i++)
        {
            Vertex* v = graph.vlist[i];
            graph.connect_edges_approx(v);
        }
        timing_out<<n<<"\t"<<toc()<<endl;
        timing_out.flush();
        cout<<n<<endl;
    }
#endif

#if 1
    Graph graph(sys);
    tic();
    // incremental construction
    for(int j=0; j< 1000; j++)
    {
        Vertex* v = graph.add_sample();
        graph.connect_edges_approx(v);
        graph.reconnect_edges_neighbors(v);
    }
    graph.seeding_finished = true;
    timing_out<<graph.num_vert<<"\t"<<toc()<<endl;
    for(int j=0; j< 50000; j++)
    {
        Vertex* v = graph.add_sample();
        graph.connect_edges_approx(v);
        graph.reconnect_edges_neighbors(v);

        timing_out<<graph.num_vert<<"\t"<<toc()<<endl;
        timing_out.flush();

        if(graph.num_vert%1000 == 0)
            cout<<graph.num_vert<<endl;
    }
#endif
    timing_out.close();
    return 0;
}

int do_movie(int tot_vert)
{
    System sys;

    Graph graph(sys);
    graph.propagate_system();
    //graph.get_kalman_path();

#if 1
    // incremental construction
    tic();
    for(int i=0; i< 10; i++)
    {
        Vertex* v = graph.add_sample(true);
        graph.connect_edges_approx(v);
        graph.reconnect_edges_neighbors(v);
    }

    for(int i=0; i < 500; i++)
    {
        Vertex* v = graph.add_sample();
        graph.connect_edges_approx(v);
        graph.reconnect_edges_neighbors(v);
    }

    for(int j=0; j< tot_vert-500; j++)
    {
        Vertex* v = graph.add_sample();
        graph.connect_edges_approx(v);

        graph.reconnect_edges_neighbors(v);

        if(graph.num_vert%1000 == 0)
        {
            //cout<<graph.num_vert << endl;
        }
    }
    // normalize density
    for(unsigned int i=0; i< graph.num_vert; i++)
    {
        Vertex* v = graph.vlist[i];
        v->prob_best_path = normal_val(graph.system->init_state.x, graph.system->init_var,\
                v->s.x, NUM_DIM);
    }
    graph.normalize_density();

#if 0
    graph.best_path.clear(); 
    //graph.best_path.push_back(get_mean(graph));
    for(unsigned int i=0; i< graph.obs.size(); i++)
    {
        graph.obs_curr_index = i;
        //cout<< "time: " << graph.obs_times[graph.obs_curr_index] << "\t";

        graph.update_density_implicit_no_obs_all();
        graph.normalize_density();

        graph.best_path.push_back(get_mean(graph, false));
    }
#endif
#if 1
    for(int i=0; i< 10000; i++)
    {
        graph.simulate_trajectory_implicit();
    }
    graph.plot_monte_carlo_density((char*)"data/density.dat");
#endif
#endif

    // output graph
    graph.plot_graph();
    //graph.plot_trajectory();

    return 0;
}

int main(int argc, char* argv[])
{
    //srand(time(NULL));
    cout.precision(5);

    int tot_vert = 1000;
    if (argc > 1)
        tot_vert = atoi(argv[1]);

    do_error_plot();
    // do_batch(tot_vert);
    // do_incremental(tot_vert);

    // do_timing_plot();
    // do_movie(tot_vert);


    return 0;
}

