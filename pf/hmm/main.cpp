#include "common.h"
#include "gnuplot.cpp"
#include "kdtree.h"

#define N       1000
#define sobs    0.05
#define spro    0.05

Gnuplot gplt("lines");
kdtree *state_tree;

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
    //gplt<<"set term postscript eps color";
    //gplt<<"set output \"traj.ps\"";
    gplt<<"set term png";
    gplt<<"set output \"traj.png\"";
}

void plot(graph &g, string options="")
{
    vector<float> vt;
    vector<float> index;
    for(vector<vertex*>::iterator i = g.vlist.begin(); i != g.vlist.end(); i++)
    {
        vertex *tstart = (*i);
        for(vector<edge*>::iterator eo = tstart->edgeout.begin(); eo != tstart->edgeout.end(); eo++)
        {
            vertex *tend = (*eo)->to;
            //draw the edge
            vt.push_back(tstart->s.x[0]);
            index.push_back(tstart->t);
            vt.push_back(tend->s.x[0]);
            index.push_back(tend->t);
        }
    }
    gplt.set_style("lp linestyle 1").plot_xy(vt, index, "graph");
}

int main()
{
    cout.precision(4);
    srand(time(0));
    gnuplot_init();

    vertex v(sample(), 0);

    graph rrg;
    
    state x0; x0.x[0] = 4.5;
    vector<state> x, y;
    x.push_back(x0);
    y.push_back(obs(x0));
    for(int i=0; i<N; i++){
        x.push_back(system(x.back()));
        y.push_back(obs(x.back()));
    }
    
    state_tree = kd_create(NUM_DIM);
    kdres *res;
    for(int i=0; i<20; i++)
    {
        vertex *v = new vertex(sample(), i);
        rrg.add_vertex(v);
        double toput = v->s.x[0];
        kd_insert(state_tree, &toput, v);
        
        double pos;
        res = kd_nearest_range(state_tree, &toput, 3*spro);
        while( !kd_res_end(res))
        {
            vertex *v1 = (vertex *)kd_res_item(res, &pos); 
            //cout<<"found [v, v1] "<<v<<" "<<v1<<endl;
            
            if( v != v1)
            {
                // make edges
                edge *e1 = new edge(v, v1, 0.5);
                edge *e2 = new edge(v1, v, 0.5);

                v->edgeout.push_back(e1);
                v->edgein.push_back(e2);
                v1->edgeout.push_back(e2);
                v1->edgein.push_back(e1);
            }
            kd_res_next(res);
        }
    }
    kd_res_free(res);

    plot(rrg);

    //cout<<"e1: "<<e1.from->s.x[0]<<", "<<e1.to->s.x[0]<<endl;

    kd_free(state_tree);
     
    return 0;
}
