#include "decode.h"

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

        bpe += graph.dist(s1,s2)*graph.dist(s1,s2);
        kfe += graph.dist(s1,s3)*graph.dist(s1,s3);
        pfe += graph.dist(s1,s4)*graph.dist(s1,s4);

        bpiter++;
        kfiter++;
        pfiter++;
    }
}

void get_interpolated_error(Graph& graph, double& bpe)
{
    bpe = 0;
    graph.best_path.reverse();

    list<State>::iterator bp_prev = graph.best_path.begin();
    list<State>::iterator bp_next = graph.best_path.begin();
    list<State>::iterator bp_last_element = graph.best_path.end();
    --bp_last_element;
    bp_next++;
    for(list<State>::iterator i = graph.truth.begin(); i!=graph.truth.end(); i++)
    {
        State& curr = *i;
        State& prev = *bp_prev;
        State& next = *bp_next;
        double inc = 0;
        if( (curr.x[0] >= prev.x[0]) && (curr.x[0] <= next.x[0]))
        {
            double interp = prev.x[1] + (next.x[1] - prev.x[1])/(next.x[0] - prev.x[0])*(curr.x[0] - prev.x[0]);
            inc = sq(interp - curr.x[1]);
            bpe = bpe + inc;
        }
        else if(curr.x[0] > next.x[0])
        {
            bp_prev++;
            bp_next++;
            if(bp_prev == graph.best_path.end())
                bp_prev = bp_last_element;
            if(bp_next == graph.best_path.end())
                bp_next = bp_last_element;
        }
    }
    bpe = bpe/(double)graph.truth.size();
}

int do_batch(int tot_vert)
{
    System sys;
    Graph graph(sys);

    double start = get_msec();
    
    srand(0);
    graph.vlist.reserve(tot_vert);

#if 1 
    tic();
    int begin = 1000;
    int end = 1000;
    // seeding
    for(int i=0; i< begin; i++)
        graph.add_sample(1);
    
    for(unsigned int i=0; i< graph.num_vert; i++)
    {
        Vertex* v = graph.vlist[i];
        v->prob_best_path = normal_val(&(graph.system->init_state.x[1]), graph.system->init_var,\
                &(v->s.x[1]), NUM_DIM-1);
    }
    graph.normalize_density();
    graph.seeding_finished = true;


#if 1
    for(int i=0; i< end; i++)
        graph.add_sample(2);
    
    for(int i=0; i< tot_vert-begin-end; i++)
        graph.add_sample();
    
    for(unsigned int i=0; i< graph.vlist.size(); i++)
    {
        Vertex* v = graph.vlist[i];
        graph.connect_edges_approx(v);
        
        if(i%(graph.vlist.size()/10) == 0)
        {
            cout<<i<<endl;
            toc();
        }
    }
    
    int tries = 1;
    double bpe_avg = 0;
    for(int i=0; i< tries; i++)
    {
        srand(0);
        graph.propagate_system();
        graph.get_best_path();
        double bpe=0;
        get_interpolated_error(graph, bpe);
        //cout<<bpe<<endl;
        bpe_avg += bpe;
    }
    cout<<"bpe: "<< bpe_avg/(double)tries<< endl;
#endif

#endif
     
    graph.plot_trajectory();
    //graph.plot_graph();
    cout<<"time: "<< get_msec() - start << " [ms]" << endl;

    return 0;
}

int do_incremental(int tot_vert)
{
    System sys;
    Graph graph(sys);

    double start = get_msec();
    
    graph.vlist.reserve(tot_vert);
    graph.propagate_system();
    graph.get_kalman_path();
    graph.get_pf_path();

#if 1 
    tic();
    // seeding
    for(int i=0; i< 100; i++)
    {
        graph.add_sample(true);
    }
    for(unsigned int i=0; i< graph.num_vert; i++)
    {
        Vertex* v = graph.vlist[i];
        v->prob_best_path = normal_val(&(graph.system->init_state.x[1]), graph.system->init_var,\
                &(v->s.x[1]), NUM_DIM-1);
    }
    graph.normalize_density();
    graph.seeding_finished = true;


#if 1
    for(int i=0; i< tot_vert-100; i++)
    {
        Vertex* v = graph.add_sample();
        //Vertex* v = graph.vlist[i];
        graph.connect_edges_approx(v);
        graph.reconnect_edges_neighbors(v);
    
        if(i%500 == 0)
        {
            cout<<i<<endl;
            toc();
        }
    }

    graph.get_best_path();

    double kfe, bpe, pfe;
    get_sq_error(graph, bpe, kfe, pfe);
    cout<<"bpe: "<< bpe<<" kfe: " << kfe << endl;
#endif

#endif
     
#if 0
    graph.is_everything_normalized();
    cout<<"starting simulation of trajectories" << endl;
    for(int i=0; i< 10; i++)
    {
        graph.simulate_trajectory();
    }
    graph.plot_monte_carlo_trajectories();
#endif

    graph.plot_trajectory();
    // graph.plot_graph();
    cout<<"time: "<< get_msec() - start << " [ms]" << endl;

#if 0 
    double bpe, kfe, pfe;
    get_sq_error(graph, bpe, kfe, pfe);
    cout<<"bpe: "<< bpe <<" kfe: "<< kfe << " pfe: "<< pfe << endl;
#endif

    return 0;
}

int do_error_plot()
{
    ofstream err_out("data/err.dat");

    int max_runs = 10;

    for(int tot_vert=20000; tot_vert < 50000; tot_vert+= 5000)
    {
        double average_time = 0;
        double average_bpe = 0;
        double stddev_bpe = 0;
        double average_kfe = 0;

        srand(0);
        System sys;
        Graph graph(sys);
        graph.vlist.reserve(tot_vert);

        for(int i=0; i< 1000; i++)
            graph.add_sample(1);

        for(unsigned int i=0; i< graph.num_vert; i++)
        {
            Vertex* v = graph.vlist[i];
            v->prob_best_path = normal_val(&(graph.system->init_state.x[1]), graph.system->init_var,\
                    &(v->s.x[1]), NUM_DIM-1);
        }
        graph.normalize_density();
        graph.seeding_finished = true;

        for(int i=0; i< tot_vert-2000; i++)
            graph.add_sample();
        for(int i=0; i<1000; i++)
            graph.add_sample(2);

        for(unsigned int i=0; i< graph.vlist.size(); i++)
        {
            Vertex* v = graph.vlist[i];
            graph.connect_edges_approx(v);
        }

        srand(0);
        for(int how_many=0; how_many < max_runs; how_many++)
        {
            graph.propagate_system();
            tic();
            graph.get_best_path();
            double kfe=0, bpe=0, pfe=0;
            get_interpolated_error(graph, bpe);
            get_sq_error(graph, pfe, kfe, pfe);
            //cout<<graph.num_vert<<" "<< bpe<<endl;

            average_time += toc();
            average_bpe += bpe;
            stddev_bpe += (bpe*bpe);
            average_kfe += kfe;
        }

        average_time = average_time/(double)max_runs;
        average_bpe = average_bpe/(double)max_runs;
        average_kfe = average_kfe/(double)max_runs;

        if(max_runs > 1)
            stddev_bpe = (stddev_bpe - max_runs*average_bpe*average_bpe)/(max_runs -1.0);

        err_out<<tot_vert<<" "<<average_time<<" "<< average_bpe <<" "<< average_kfe << endl;
        cout<<tot_vert<<"\t"<<average_time<<"\t"<< average_bpe <<"\t"<< average_kfe << endl;
        err_out.flush();
    }
    err_out.close();

    return 0;
}

int main(int argc, char* argv[])
{
    cout.precision(5);

    int tot_vert = 0;
    if (argc > 1)
        tot_vert = atoi(argv[1]);
    else
        tot_vert = 1000;

    do_batch(tot_vert);
    // do_incremental(tot_vert);
    // do_error_plot();

    return 0;
}

