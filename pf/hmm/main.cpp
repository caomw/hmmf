#include "common.h"
#include "gnuplot.cpp"
#include "kdtree.h"

#define M       100
#define N       1000
#define sobs    0.1
#define spro    0.05
#define BETA    (10)
#define GAMMA   (XMAX - XMIN)

// halton
int seed[NUM_DIM];
int base[NUM_DIM];

Gnuplot gplt("lines");
kdtree *state_tree, *mini_tree;
graph rrg;
//vector<minis *> mini_samples;

void halton_init()
{
    seed[0] = 0;
    base[0] = 3;
    halton_dim_num_set (1);
    halton_step_set ( 0);
    halton_seed_set(seed);
    halton_base_set(base);
};

double func_to_int(double x)
{
    return x*x;
    //return 1.0;
}

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
    
    gplt<<"set style line 1 lt 1.2 lc rgb \"red\" pt 6 ps 0.5";// pt 1 pointsize default";
    gplt<<"set style line 2 lt 1.2 lc rgb \"blue\" pt 1 ps 0.2";// pt 2 pointsize default";
    gplt<<"set style line 3 lt 1 lc rgb \"green\"";// pt 4 pointsize default";
    gplt<<"set style line 4 lt 1 lc rgb \"yellow\"";// pt 3 pointsize default";
    gplt<<"set grid";
    gplt<<"set term pdf";
    gplt<<"set output \"traj.pdf\"";
    //gplt<<"set term png";
    //gplt<<"set output \"traj.png\"";
}

void plot_rrg()
{
    vector<float> vt;
    vector<float> index;
    for(vector<vertex*>::iterator i = rrg.vlist.begin(); i != rrg.vlist.end(); i++)
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
    gplt.reset_plot();
    gplt.set_style("linespoints").plot_xy(vt, index, "graph");
}

void plot_traj(vector<state> x, vector<state> y)
{
    vector<float> xf, yf;
    for(int i=0; i< (int)x.size(); i++)
    {
        xf.push_back(x[i].x[0]);
        yf.push_back(y[i].x[0]);
    }
    //gplt<<"set multiplot";
    gplt.set_style("lines ls 1").plot_x(xf, "sys");
    gplt.set_style("points ls 2").plot_x(yf, "obs");
}

double get_edge_prob(vertex *v)
{
    double *tosearch;
    tosearch = v->s.x;
    kdres *res;
    int num_vert = rrg.vlist.size();
    double tmp = 0;

    res = kd_nearest_range(mini_tree, tosearch, GAMMA*log(num_vert)/(float)num_vert);
    while( !kd_res_end(res))
    {
        double pos;
        minis *m = (minis *)kd_res_item(res, &pos);
        if( m->parent == v)
            tmp = tmp+1;

        kd_res_next(res);
    }
    return tmp;
}

void build_graph()
{

    kdres *res;
    for(int i=0; i<M; i++)
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
    for(int j=0; j<M*BETA; j++)
    {
        minis *m;
        m = new (minis);
        m->s = sample_quasi();
        double toput = m->s.x[0];
        kd_insert(mini_tree, &toput, m);
        //mini_samples.push_back(m);

        res = kd_nearest(state_tree, &toput);
        while( !kd_res_end(res))
        {
            double pos = 0;
            vertex *v = (vertex*)kd_res_item(res, &pos);
            m->parent = v;
            
            double tmp = fabs(pos- toput); 
            v->ravg = ((v->ravg)*(v->num_child) + tmp)/(v->num_child + 1.0);
            v->num_child = v->num_child + 1;
            //cout<<tmp<<" "<<v->ravg<<" "<<v->num_child<<endl;

            kd_res_next(res);
        }
    }
    kd_res_free(res);
    
    for(int i=0; i< (int)rrg.vlist.size(); i++)
    {
        vertex *v = rrg.vlist[i];
        for(int j=0; j < (int)v->edgeout.size(); j++)
        {
            edge *e = v->edgeout[j];
            vertex *v1 = e->to;
            e->prob = get_edge_prob(v1);
            //cout<<"e: "<<e<<" "<<e->prob<<endl;
        }
    }
}

void calci_int_samples()
{
    //ofstream vertout("vert.dat");
    kdres *res;
    double intval = 0;
    int num_vert = rrg.vlist.size();

    for(int i=0; i< (int)rrg.vlist.size(); i++)
    {
        double fval = 0;
        vertex *v = rrg.vlist[i];
        res = kd_nearest_range(mini_tree, v->s.x, GAMMA*log(num_vert)/(float)num_vert);
        while( !kd_res_end(res))
        {
            double pos;
            minis *m = (minis *)kd_res_item(res, &pos);
            if( m->parent == v)
            {
                fval += func_to_int(m->s.x[0]);
            }
            kd_res_next(res);
        }
        if(v->num_child != 0)
            fval = fval/((float)v->num_child);
        
        //cout<<v->s.x[0]<<" "<<fval<<" "<<v->ravg<<" "<<v->num_child<<endl;
        //vertout<<v->s.x[0]<<"\t"<<fval<<endl;
        intval += fval;
    }
    cout<<"intval ["<<XMIN<<", "<<XMAX<<"]: "<<intval/M*(XMAX-XMIN)<<endl;
    kd_res_free(res);
    //vertout.close();
}

int main()
{

    cout.precision(4);
    srand(time(0));
    state_tree = kd_create(NUM_DIM);
    mini_tree = kd_create(NUM_DIM);
    halton_init();
    //gnuplot_init();
    
    double ts = get_msec();
    state x0; x0.x[0] = 4.5;
    vector<state> x, y;
    x.push_back(x0);
    y.push_back(obs(x0));
    for(int i=0; i<N; i++){
        x.push_back(system(x.back()));
        y.push_back(obs(x.back()));
    }
    //plot_traj(x, y); 

    build_graph();
    
    calci_int_samples();
    cout<<"dt: "<<get_msec() -ts<<endl;
    
    //plot_rrg();
    
    //mini_samples.clear();
    kd_free(state_tree);
    kd_free(mini_tree); 
    return 0;
}

