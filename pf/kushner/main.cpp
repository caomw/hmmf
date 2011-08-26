#include "hmmf.h"

void get_mean(Graph& graph)
{
    float mean[NUM_DIM] = {0};
    for(int i=0; i<graph.num_vert; i++)
    {
        for(int j=0; j< NUM_DIM; j++)
            mean[j] += graph.vlist[i]->s.x[j];
    }
    cout << "mean: ";
    for(int j=0; j< NUM_DIM; j++)
    {
        mean[j] = mean[j]/graph.num_vert;
        cout<< mean[j] << " ";
    }
    cout<<endl;
}

int do_filtering()
{
    System sys;
    Graph graph(sys);
    graph.propagate_system();
    graph.get_kalman_path();
    graph.plot_trajectory();
    
    double start = get_msec();
    
    graph.put_init_samples(100);
    for(int i=0; i< 500; i++)
        graph.add_sample();
    
    cout<<"added samples: "<< graph.num_vert << endl;
    
    int howmany_updates = log(graph.num_vert);

    tic();
    for(int i=0; i < graph.obs.size(); i++)
    {
        graph.obs_curr_index = i;
        cout<< "obs_times: " << graph.obs_times[graph.obs_curr_index] << endl;

        for(int j=0; j < howmany_updates; j++)
        {
            Vertex* v = graph.vlist[RANDF*graph.num_vert];
            graph.propagate_density(v);
        }
        cout<<i <<": ";
        get_mean(graph);
        
    }
    
    graph.plot_graph();
    assert(graph.is_everything_normalized());

    cout<<"time: "<< get_msec() - start << " [ms]" << endl;
    
}

int graph_sanity_check()
{
    System sys;
    Graph graph(sys);
    graph.propagate_system();
    graph.get_kalman_path();
    graph.plot_trajectory();
    
    double start = get_msec();
    

    tic();
    
    graph.put_init_samples(100);
    for(int i=0; i < 1000; i++)
    {
        graph.add_sample();
        if(i % 500 == 0)
        {
            cout<< i << endl;
            toc();
        }
    }

    graph.plot_graph();
    cout<<"added samples: "<< graph.num_vert << endl;

    assert(graph.is_everything_normalized());
    cout<<"graph is correct"<<endl;
    
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
    
    do_filtering();

    return 0;
}

