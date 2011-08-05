#include "hmmf.h"


double get_msec()
{
    struct timeval start;
    gettimeofday(&start, NULL);
    return start.tv_sec*1000 + start.tv_usec/1000.0;
}

int do_decoding()
{
    System sys;
    Graph graph(sys);
    graph.propagate_system();

    double start = get_msec();

    graph.put_init_samples();

    int count = 0;
    for(list<State>::iterator i= graph.obs.begin(); i != graph.obs.end(); i++)
    {
        State& curr_obs = *i;
        for(int j=0; j< graph.samples_per_obs; j++)
            graph.add_sample();

        graph.update_observation_prob(curr_obs);

        int max_async_updates = min((double)graph.get_num_vert(), 10*log( graph.get_num_vert()));
        int num_vert = graph.get_num_vert();
        for(int j=0; j< max_async_updates; j++)
        {
            int which = RANDF*num_vert;
            Vertex *v = graph.vlist[which];
            graph.update_viterbi(v);
        }

        count++;

        if(count % 10 == 0)
            cout<<count<<endl;
    }
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

    double start = get_msec();

    graph.put_init_samples();
    
    for(int i=0; i < 1000; i++)
    {
        graph.add_sample();
    }
    cout<<"finished putting samples"<<endl;
    assert(graph.is_everything_normalized());
    //cout<<"graph is correct"<<endl;

    for(int i=0; i< 1000; i++)
    {
        graph.simulate_trajectory();
        cout<<i<<endl;
    }

    graph.plot_graph();
    cout<<"time: "<< get_msec() - start << " [ms]" << endl;
}

int main()
{
    srand(0);

    graph_sanity_check();

    return 0;
}

