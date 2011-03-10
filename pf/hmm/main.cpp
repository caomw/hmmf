#include "common.h"
#include "gnuplot.cpp"

#define N       1000
#define sobs    0.05
#define spro    0.05

Gnuplot gplt("lines");

state system(state s)
{
    state t;
    t.x[0] = s.x[0]*0.95 + randn(0, spro);
    return t;
}

state obs(state s){
    state t;
    t.x[0] = s.x[0] + randn(0, sobs);
    return t;
}

void gnuplot_init(){
    gplt.reset_all();
    gplt.set_grid();
    gplt<<"set style line 1 lt 1 lw 1 pt 1";
    gplt<<"set term postscript enhanced color";
    gplt<<"set output \"traj.ps\"";
}

void plot(graph g, string options="")
{
    vector<float> vt;
    for(vector<vertex*>::iterator i = g.vlist.begin(); i != g.vlist.end(); i++)
        vt.push_back((*i)->s.x[0]);

    gplt.set_style("lines linestyle 1").plot_x(vt, "graph");
}

int main()
{
    cout.precision(4);
    srand(time(0));
    gnuplot_init();

    graph rrg;
    
    state x0; x0.x[0] = 4.5;
    vector<state> x, y;
    x.push_back(x0);
    y.push_back(obs(x0));
    for(int i=0; i<N; i++){
        x.push_back(system(x.back()));
        y.push_back(obs(x.back()));
    }
     
    vector<float> tmp;
    for(unsigned int i=0; i<y.size(); i++)
    {
        tmp.push_back(y[i].x[0]);
    }
    gplt.set_style("lines linestyle 1").plot_x(tmp, "graph");
   
    for(int i=0; i<2; i++)
    {
        vertex *v = new vertex;
        *v = vertex(sample());
        rrg.add_vertex(v);
    }
    //cout<<"e1: "<<e1.from->s.x[0]<<", "<<e1.to->s.x[0]<<endl;

    return 0;
}
