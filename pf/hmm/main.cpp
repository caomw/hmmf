#include "common.h"
#include "gnuplot.cpp"
#include "kdtree.h"

#define EXTEND_DIST (0.05)
#define dt      0.1
#define M       100
#define N       50
#define sobs    0.01
#define spro    0.01
#define BETA    (100)
#define GAMMA   (XMAX - XMIN)
#define BOWLR   (GAMMA*sqrt(log(rrg.num_vert)/(float)(rrg.num_vert)))

// halton
int seed[NUM_DIM] = {0, 0};
int base[NUM_DIM] = {2, 3};
int step = 0;

Gnuplot gplt("lines");
kdtree *state_tree, *mini_tree;
graph rrg;

void halton_init()
{
    halton_dim_num_set (NUM_DIM);
    halton_step_set (step);
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
    t.x[0] = s.x[0] + dt*( -0.5*s.x[0] + randn(0, spro));
    t.x[1] = s.x[1] + dt*( -0.9*s.x[1] + randn(0, spro));
    return t;
}

state obs(state s){
    state t;
    t.x[0] = s.x[0] + randn(0, sobs);
    t.x[1] = s.x[1] + randn(0, sobs);
    return t;
}

void gnuplot_init(){
    gplt.reset_all();
    gplt.set_grid();
    
    gplt<<"set style line 1 lt 1.2 lc rgb \"red\" pt 6 ps 0.5";// pt 1 pointsize default";
    gplt<<"set style line 2 lt 1.2 lc rgb \"blue\" pt 1 ps 0.2";// pt 2 pointsize default";
    gplt<<"set style line 3 lt 1 lc rgb \"green\" pt 1 ps 0.2";// pt 4 pointsize default";
    gplt<<"set style line 4 lt 1 lc rgb \"yellow\"";// pt 3 pointsize default";
    gplt<<"set grid";
    gplt<<"set term pdf";
    gplt<<"set output \"traj.pdf\"";
    //gplt<<"set term png";
    //gplt<<"set output \"traj.png\"";
}

void plot_rrg()
{
    vector<float> vt1, vt2;
    for(vector<vertex*>::iterator i = rrg.vlist.begin(); i != rrg.vlist.end(); i++)
    {
        vertex *tstart = (*i);
        for(vector<edge*>::iterator eo = tstart->edgeout.begin(); eo != tstart->edgeout.end(); eo++)
        {
            vertex *tend = (*eo)->to;
            //draw the edge
            vt1.push_back(tstart->s.x[0]);
            vt2.push_back(tstart->s.x[1]);
            vt1.push_back(tend->s.x[0]);
            vt2.push_back(tend->s.x[1]);
        }
    }
    //gplt.reset_plot();
    gplt.set_style("points ls 3").plot_xy(vt1, vt2, "graph");
}

void plot_traj(vector<state> x, vector<state> y)
{
    vector<float> xf1, yf1, xf2, yf2;
    for(int i=0; i< (int)x.size(); i++)
    {
        xf1.push_back(x[i].x[0]);
        yf1.push_back(y[i].x[0]);
        xf2.push_back(x[i].x[1]);
        yf2.push_back(y[i].x[1]);
    }
    //gplt<<"set multiplot";
    gplt.set_style("lines ls 1").plot_xy(xf1, xf2, "sys");
    gplt.set_style("points ls 2").plot_xy(yf1, yf2, "obs");
}

double get_edge_prob(vertex *v)
{
    double *tosearch;
    tosearch = v->s.x;
    kdres *res;
    double hits = 0, tot=0;

    res = kd_nearest_range(mini_tree, tosearch, BOWLR);
    while( !kd_res_end(res))
    {
        double pos;
        minis *m = (minis *)kd_res_item(res, &pos);
        if( m->parent == v)
            hits++;

        tot++;
        kd_res_next(res);
    }
    v->voronoi_area = hits/((float)tot)*BOWLR;
    return hits;
}

void calci_int_samples()
{
    //ofstream vertout("vert.dat");
    kdres *res;
    double intval = 0;

    for(int i=0; i< rrg.num_vert; i++)
    {
        double fval = 0;
        vertex *v = rrg.vlist[i];
        res = kd_nearest_range(mini_tree, v->s.x, BOWLR);
        while( !kd_res_end(res))
        {
            double pos;
            minis *m = (minis *)kd_res_item(res, &pos);
            if( m->parent == v)
                fval += func_to_int(m->s.x[0]);
            
            kd_res_next(res);
        }
        if(v->num_child != 0)
            fval = fval/((float)v->num_child);
        
        //cout<<v->s.x[0]<<" "<<fval<<" "<<v->num_child<<endl;
        //vertout<<v->s.x[0]<<"\t"<<fval<<endl;
        intval += fval;
    }
    cout<<"intval ["<<XMIN<<", "<<XMAX<<"]: "<<intval/M*(XMAX-XMIN)<<endl;
    kd_res_free(res);
    //vertout.close();
}

vertex* nearest_vertex(state s)
{
    kdres *res;
    res = kd_nearest(state_tree, s.x);
    if(kd_res_end(res)){
        cout<<"Error: no nearest"<<endl;
        exit(1);
    }
    vertex *v = (vertex*)kd_res_item_data(res);
    kd_res_free(res);

    return v;
}

void add_mini_samples()
{
    double *pos;
    pos = (double *)malloc(sizeof(state));

    for(int j=0; j<BETA; j++)
    {
        minis *m;
        m = new (minis);
        m->s = sample();
        double *toput = m->s.x;
        kd_insert(mini_tree, toput, m);
        
        vertex *v = nearest_vertex(m->s);
        m->parent = v;
        v->num_child = v->num_child + 1;
    }
    free(pos);
}

void add_major_sample(vertex *v)
{
    double *pos;
    pos = (double *)malloc(sizeof(state));

    if(rrg.num_vert != 0)
    {
        vertex *vnear = nearest_vertex(v->s);
        double len = dist(v->s, vnear->s);
        if( len > EXTEND_DIST)
        {
            // change the state of the vertex to extend state
            for(int i=0; i<NUM_DIM; i++){
                v->s.x[i] = vnear->s.x[i] + EXTEND_DIST*(v->s.x[i] - vnear->s.x[i])/len;
            }
        }
    }
    double *toput = v->s.x;
    kd_insert(state_tree, toput, v);
    rrg.add_vertex(v);
    rrg.num_vert++;

    kdres *res;
    res = kd_nearest_range(state_tree, toput, BOWLR);
    //cout<<"got "<<kd_res_size(res)<<" states"<<endl;
    while( !kd_res_end(res))
    {
        vertex *v1 = (vertex *)kd_res_item(res, pos); 
        
        if( v != v1)
        {
            float t1 = get_edge_prob(v1);
            float t2 = get_edge_prob(v);
            // make edges
            edge *e1 = new edge(v, v1, t1);
            edge *e2 = new edge(v1, v, t2);

            v->edgeout.push_back(e1);
            v->edgein.push_back(e2);
            v1->edgeout.push_back(e2);
            v1->edgein.push_back(e1);
        }
        kd_res_next(res);
    }
    kd_res_free(res);
    free(pos);

    add_mini_samples();
}

int main()
{
    cout.precision(4);
    srand(time(0));
    state_tree = kd_create(NUM_DIM);
    mini_tree = kd_create(NUM_DIM);
    //halton_init();
    gnuplot_init();
    
    double ts = get_msec();
    state x0; x0.x[0] = 1.0, x0.x[1] = -1.0;
    vector<state> x, y;
    x.push_back(x0);
    y.push_back(obs(x0));
    for(int i=0; i<N; i++){
        x.push_back(system(x.back()));
        y.push_back(obs(x.back()));
    }
    plot_traj(x, y);

    // put initial few observations as states
    for(int i=0; i<10; i++)
    {
        vertex *v = new vertex(y[i], 0, 1);
        add_major_sample(v);
    }
    // setup pi_0 here
    int dM = 1;
    for(int i=10; i<N; i++)
    {
        // add the obs as the state
        vertex *v = new vertex(y[i], 0, 1);
        add_major_sample(v);
        
        // add some more states
        for(int j=0; j<dM; j++)
        {
            vertex *v = new vertex(sample(), 0, 1);
            add_major_sample(v);
        }
    }
    plot_rrg();
    
    
    kd_free(state_tree);
    kd_free(mini_tree); 
    return 0;
}

