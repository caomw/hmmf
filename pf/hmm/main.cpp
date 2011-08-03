#include "singleint.h"
#include "hmmf.h"


double get_msec()
{
    struct timeval start;
    gettimeofday(&start, NULL);
    return start.tv_sec*1000 + start.tv_usec/1000.0;
}

int main()
{
    srand(0);

    /*
       cout<<"start at: [" << s[0] <<","<< s[1] << "]"<<endl;
       for(int i=0; i< 10; i++)
       {
       State t = sys.integrate(s, 0.1*i, 0);
       State to = sys.observation(t, 0);
       cout<< i<<" at: [" << t[0] <<","<< t[1] << "]"<<endl;
       cout<< i<<" saw: [" << to[0] <<","<< to[1] << "]"<<endl;
       }
       */
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
            graph.do_viterbi(v);
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

