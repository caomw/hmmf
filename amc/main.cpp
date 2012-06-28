#include "filter.h"

State get_mean(Graph& graph, State& covarstate, bool is_cout=true)
{
    State mstate;
    double mean[NUM_DIM] = {0};
    double covar[NUM_DIM] = {0};
    double totprob = 0;
    for(unsigned int i=0; i<graph.num_vert; i++)
    {
        for(int j=0; j< NUM_DIM; j++)
            mean[j] += ((graph.vlist[i]->s.x[j]) * (graph.vlist[i]->prob_best_path));
        
        for(int j=0; j< NUM_DIM; j++)
            covar[j] += (sq(graph.vlist[i]->s.x[j]) * (graph.vlist[i]->prob_best_path));

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
        covarstate.x[j] = covar[j]/totprob - sq(mstate.x[j]);

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

void get_sq_error(Graph& graph, double& bpe, double& kfe, double& pfe, double max_time=1.0)
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
        
        /*
        bpe = bpe + sq(s1.x[1] - s2.x[1]);
        kfe = kfe + sq(s1.x[1] - s3.x[1]);
        pfe = pfe + sq(s1.x[1] - s4.x[1]);
        */
        bpiter++;
        kfiter++;
        pfiter++;
    }
    bpe = bpe*max_time/(double)graph.truth.size();
    kfe = kfe*max_time/(double)graph.truth.size();
    pfe = pfe*max_time/(double)graph.truth.size();
}

int do_err_convergence()
{
    System sys;
    int inc = 1000;
    for (int n=5000; n< 200000; n=n+inc)
    {
        srand(0);
        Graph graph(sys);
        for(int i=0; i < n; i++)
        {
            graph.add_sample(false);
        }
        for(unsigned int i=0; i < graph.num_vert; i++)
        {
            Vertex* v = graph.vlist[i];
            graph.connect_edges_approx(v);
        }
        graph.calculate_delta();
        graph.system->sim_time_delta = 0.01; // graph.delta;
        graph.calculate_probabilities_delta_all();
        for(int i=0; i< 50000; i++)
            graph.simulate_trajectory_implicit();
        
        double e1=0, e2=0;
        graph.analyse_monte_carlo_trajectories(e1, e2);
        //graph.plot_monte_carlo_trajectories();
    }
    return 0;
}

int do_err_convergence_incremental()
{
    ofstream eout("err.dat");
    System sys;
    Graph graph(sys);
    int inc = 1000;
    for(int n=0; n < 1000; n++)
    {
        Vertex* v = graph.add_sample(false);
        graph.connect_edges_approx(v);
        graph.calculate_probabilities_delta(v);
        graph.reconnect_edges_neighbors(v);
    }
    for (int n=0; n< 1000000; n=n+1)
    {
        Vertex* v = graph.add_sample(false);
        graph.connect_edges_approx(v);
        graph.calculate_probabilities_delta(v);
        graph.reconnect_edges_neighbors(v);
        
        if(n % inc == 0)
        {
            for(int i=0; i< 50000; i++)
                graph.simulate_trajectory_implicit();
            double e1=0, e2=0;
            graph.analyse_monte_carlo_trajectories(e1, e2);
            eout<<graph.num_vert<<" "<<e1<<" "<<e2<<endl;
            eout.flush();
            graph.clear_monte_carlo_trajectories();
        }
    }
    eout.close();
    return 0;
}
int do_batch(int tot_vert)
{
    System sys;
    Graph graph(sys);
    srand(0);
    
#if 1
    tic();
    for(int i=0; i < 100; i++)
    {
        graph.add_sample(true);
    }
    for(int i=0; i < tot_vert-100; i++)
    {
        graph.add_sample(false);
    }
    for(unsigned int i=0; i < graph.num_vert; i++)
    {
        Vertex* v = graph.vlist[i];
        graph.connect_edges_approx(v);
        /*
        if(i %1000 == 0)
        {
            cout<<"i: "<< i << endl;
            toc();
        }
        */
    }
    for(unsigned int i=0; i< graph.num_vert; i++)
    {
        Vertex* v = graph.vlist[i];
        v->prob_best_path = normal_val(graph.system->init_state.x, graph.system->init_var,\
                v->s.x, NUM_DIM);
    }
    graph.normalize_density();
    graph.seeding_finished = true;
    
    graph.calculate_delta();
    cout<<"delta: "<< graph.delta<<endl;
    graph.system->sim_time_delta = graph.delta;
    graph.calculate_probabilities_delta_all();
    
    srand(0);
    graph.propagate_system();
    graph.get_kalman_path();
    tic();
    graph.get_pf_path(1000);
    cout<<"pf time: "<< toc() << endl;
#endif
    
#if 1
    cout<<"start filter"<<endl;
    tic();
    State curr_var;
    State curr_mean;
    graph.best_path.clear();
    for(unsigned int i=0; i < graph.obs.size(); i++)
    {
        curr_mean = get_mean(graph, curr_var, false);
        /*
        for(int j=0; j<NUM_DIM; j++)
            cout<<curr_var.x[j]<<" ";
        cout<<endl;
        */
        graph.best_path.push_back(curr_mean);

        graph.obs_curr_index = i;
        // cout<< "time: " << graph.obs_times[graph.obs_curr_index] << " ";
        graph.update_density_implicit_all(curr_var.norm());
        /*
           if(i%10 == 0)
           cout<<i<<endl;
           */
    }
    curr_mean = get_mean(graph, curr_var, false);
    /*
    for(int j=0; j<NUM_DIM; j++)
        cout<<curr_var.x[j]<<" ";
    cout<<endl;
    */
    graph.best_path.push_back(curr_mean);

    cout<<"filter time: "<< toc() << endl;
#endif

#if 0
    cout<<"starting simulation of trajectories" << endl;
    for(int i=1; i< 1001; i++)
    {
        if((i%1000) == 0)
            cout<<i<<endl;
        graph.simulate_trajectory_implicit();
    }
    graph.analyse_monte_carlo_trajectories();
    //graph.plot_monte_carlo_trajectories();
#endif

    double bpe, kfe, pfe;
    get_sq_error(graph, bpe, kfe, pfe, graph.max_obs_time);
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
    State covar;
    graph.best_path.clear();
    graph.best_path.push_back(get_mean(graph, covar, false));
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

        graph.best_path.push_back(get_mean(graph, covar, false));
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

int do_error_plot(int start_vert, int end_vert=1.0)
{
    ofstream err_out("data/err_out.dat");
    int max_runs = 20;
    
    int svert = start_vert;
    int evert = end_vert;
    if(evert < svert)
        evert = svert + 1;
        
    System sys;

    for(int tot_vert=svert; tot_vert < svert+1; tot_vert+= 100)
    {
        double average_time_hmm = 0;
        double average_time_pf = 0;
        double average_bpe = 0;
        double average_kfe = 0;
        double average_pfe = 0;
       
        srand(0);
        for(int how_many=0; how_many < max_runs; how_many++)
        {
            tic();
            Graph graph(sys);
            for(int i=0; i < 100; i++)
                graph.add_sample(true);
            for(int i=0; i < tot_vert-100; i++)
                graph.add_sample(false);
            for(unsigned int i=0; i< graph.num_vert; i++)
            {
                Vertex* v = graph.vlist[i];
                graph.connect_edges_approx(v);
            }
            graph.calculate_delta();
            graph.system->sim_time_delta = graph.delta;
            graph.calculate_probabilities_delta_all();
            graph.seeding_finished = true;
            average_time_hmm += toc();

            graph.propagate_system();
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

            graph.best_path.clear();
            State covar;
            graph.best_path.push_back(get_mean(graph, covar, false));
            for(unsigned int i=0; i < graph.obs.size(); i++)
            {
                graph.obs_curr_index = i;
                graph.update_density_implicit_all();
                graph.best_path.push_back(get_mean(graph, covar, false));
            }
            average_time_hmm += toc();

            double bpe=0, kfe=0, pfe=0;
            get_sq_error(graph, bpe, kfe, pfe, graph.max_obs_time);
            cout<<how_many<<" "<<bpe<<" "<< kfe<<" "<<pfe<<endl;
            //cout<<how_many<<" "; cout.flush();

            average_bpe += bpe;
            average_kfe += kfe;
            average_pfe += pfe;
        }
        //cout<<endl;
        average_time_hmm = average_time_hmm/(double)max_runs;
        average_time_pf = average_time_pf/(double)max_runs;
        average_bpe = average_bpe/(double)max_runs;
        average_kfe = average_kfe/(double)max_runs;
        average_pfe = average_pfe/(double)max_runs;

        cout<<tot_vert<<" "<<average_time_hmm<<" "<<average_time_pf<<" "<< average_bpe <<" "<< average_kfe << " "<<average_pfe<<endl;
        err_out<<tot_vert<<" "<<average_time_hmm<<" "<<average_time_pf<<" "<< average_bpe <<" "<< average_kfe << " "<<average_pfe<<endl;
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

int main(int argc, char* argv[])
{
    cout.precision(5);
    cout<<setw(10);

    int tot_vert = 1000;
    if (argc > 1)
        tot_vert = atoi(argv[1]);

    // do_err_convergence_incremental();
    // do_err_convergence();
    // do_error_plot(tot_vert);
    do_batch(tot_vert);
    // do_incremental(tot_vert);

    // do_timing_plot();
    // do_movie(tot_vert);


    return 0;
}

