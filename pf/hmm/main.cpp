#include "hmmf.h"



int do_decoding()
{
    System sys;
    Graph graph(sys);
    graph.propagate_system();

    graph.get_kalman_path();
    
    double start = get_msec();

    graph.put_init_samples();

    tic();
    clock_t start_time = clock();
    bool time_finished = false;
    
    //while(time_finished == false)
    //{
    for(int i=0; i < 1000; i++)
    {
        graph.add_sample();
        
        int num_vert = graph.get_num_vert();
        int max_async_updates = 1000*log(num_vert);
        for(int j=0; j< max_async_updates; j++)
        {
            int which = RANDF*num_vert;
            Vertex *v = graph.vlist[which];
            graph.update_viterbi(v);
        }

        if( graph.num_vert % 100 == 0)
        {
            cout<< graph.num_vert <<endl;
            toc();
        }
        
        /*
        clock_t end_time = clock();
        if( (end_time - start_time)/CLOCKS_PER_SEC > 500.0)
            time_finished = true;
        */
    }

    /*
    for(int i=0; i< 1000; i++)
    {
        int max_async_updates = min((double)graph.get_num_vert(), 10*log( graph.get_num_vert()));
        int num_vert = graph.get_num_vert();
        for(int j=0; j< max_async_updates; j++)
        {
            int which = RANDF*num_vert;
            Vertex *v = graph.vlist[which];
            graph.update_viterbi(v);
        }
    }
    */

    graph.get_best_path();
    cout<<"time: "<< get_msec() - start << " [ms]" << endl;

    graph.plot_trajectory();
    graph.plot_graph();

    return 0;
}

int graph_sanity_check()
{
    System sys;
    Graph graph(sys);
    graph.propagate_system();
    graph.get_kalman_path();
    graph.plot_trajectory();
    
    double start = get_msec();

    graph.put_init_samples();
    
    tic();
    for(int i=0; i < 1000; i++)
    {
        graph.add_sample();
        if(i % 100 == 0)
        {
            cout<<i<<endl;
            toc();
        }
    }

    graph.plot_graph();
    cout<<"added samples: "<< graph.num_vert << endl;

    //cout<<"pruned: "<< graph.prune_graph() << endl;
    //graph.print_rrg();
    cout<<"finished putting samples"<<endl;
    assert(graph.is_everything_normalized());
    //cout<<"graph is correct"<<endl;
    
    for(int i=0; i< 100; i++)
    {
        graph.simulate_trajectory();
        //cout<<i<<endl;
    }
    graph.plot_monte_carlo_trajectories();
    
    cout<<"time: "<< get_msec() - start << " [ms]" << endl;
}

int main()
{
    srand(0);

    //graph_sanity_check();
    do_decoding();

    return 0;
}

