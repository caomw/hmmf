#include "hmmf.h"

int graph_sanity_check()
{
    System sys;
    Graph graph(sys);
    graph.propagate_system();
    graph.get_kalman_path();
    graph.plot_trajectory();
    
    double start = get_msec();

    tic();
    for(int i=0; i < 1000; i++)
    {
        graph.add_sample();
        if(i % 100 == 0)
        {
            cout<< i << endl;
            toc();
        }
    }

    graph.plot_graph();
    cout<<"added samples: "<< graph.num_vert << endl;

    assert(graph.is_everything_normalized());
    cout<<"graph is correct"<<endl;
    
    for(int i=0; i< 10; i++)
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

    graph_sanity_check();
    
    return 0;
}

