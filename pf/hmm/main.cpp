#define _GLIBCXX_FULLY_DYNAMIC_STRING   1
#undef _GLIBCXX_DEBUG
#undef _GLIBCXX_DEBUG_PEDANTIC

#define EXTEND_DIST (0.05)
#define dt      0.1
#define M       100
#define N       5
#define sobs    0.01
#define spro    0.01
#define BETA    (100)
#define GAMMA   (XMAX - XMIN)
#define BOWLR   (GAMMA*sqrt(log(rrg.num_vert)/(float)(rrg.num_vert)))

#include "common.h"
#include "kdtree.h"

// halton
int seed[NUM_DIM] = {0, 0};
int base[NUM_DIM] = {2, 3};
int step = 0;

kdtree *state_tree, *mini_tree;
graph rrg;
vector<state> x, y;

void halton_init()
{
    halton_dim_num_set (NUM_DIM);
    halton_step_set (step);
    halton_seed_set(seed);
    halton_base_set(base);
};

state system(state s, int is_clean){
    state t;
    if(!is_clean){
        t.x[0] = s.x[0] + dt*( -0.5*s.x[0] + randn(0, spro));
        t.x[1] = s.x[1] + dt*( -0.9*s.x[1] + randn(0, spro));
    }
    else{
        t.x[0] = s.x[0] + dt*( -0.5*s.x[0]);
        t.x[1] = s.x[1] + dt*( -0.9*s.x[1]);
    }
    return t;
}

state obs(state s, int is_clean){
    state t;
    if(!is_clean){
        t.x[0] = s.x[0] + randn(0, sobs);
        t.x[1] = s.x[1] + randn(0, sobs);
    }
    else{
        t.x[0] = s.x[0];
        t.x[1] = s.x[1];
    }
    return t;
}

/*
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
*/

void plot_rrg()
{
    ofstream rrgout("rrg.dat");

    for(vector<vertex*>::iterator i = rrg.vlist.begin(); i != rrg.vlist.end(); i++)
    {
        vertex *tstart = (*i);
        for(vector<edge*>::iterator eo = tstart->edgeout.begin(); eo != tstart->edgeout.end(); eo++)
        {
            vertex *tend = (*eo)->to;
            
            //draw the edge
            rrgout<<tstart->s.x[0]<<"\t"<<tstart->s.x[1]<<"\t"<<tend->s.x[0]<<"\t"<<tend->s.x[1]<<endl;
        }
    }
    rrgout.close();
}

/*
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
*/

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

double noise_func(state s, state s1, double sigma)
{
    double J = 0;
    for(int i=0; i<NUM_DIM; i++)
        J += -1.0*SQ(s.x[i] - s1.x[i])/2.0/sigma;
    
    double tmp = 1/pow(2*M_PI, NUM_DIM/2.0)/pow(sigma, NUM_DIM/2.0)*exp(J);
    return tmp;
}

inline int edge_index_connected_vertex(vertex *start, vertex *tocheck)
{
    int index = -1;
    for(int i=0; i< (int)start->edgeout.size(); i++)
    {
        if( (start->edgeout[i])->to == tocheck)
            index = i;
    }
    return index;
}

// locate all vertices within bowlr and write prob of that vertex in parent_weights
void update_obs_prob(state yt, vector<vertex*> &nodesinbowl, vector<double> &weights)
{
    nodesinbowl.clear();
    weights.clear();

    double *pos;
    pos = (double *)malloc(sizeof(state));
    double sum = 0;

    kdres *res;
    res = kd_nearest_range(mini_tree, yt.x, BOWLR);
    while( !kd_res_end(res))
    {
        minis *m = (minis *)kd_res_item(res, pos);

        double noise_tmp = noise_func(yt, m->s, sobs);
        vector<vertex *>::iterator iter = find(nodesinbowl.begin(), nodesinbowl.end(), m->parent);
        
        if(nodesinbowl.size() == 0)
        {
            nodesinbowl.push_back(m->parent);
            weights.push_back(noise_tmp);
            sum += noise_tmp;
        }
        else
        {
            if(( iter != nodesinbowl.end()) || (nodesinbowl.back() == (m->parent)))
            {
                weights[iter - nodesinbowl.begin()] = noise_tmp;
                sum += noise_tmp;
            }
            else
            {
                nodesinbowl.push_back(m->parent);
                weights.push_back(noise_tmp);
                sum += noise_tmp;
            }
        }
        kd_res_next(res);
    }
    // normalize
    for(int i=0; i< (int)nodesinbowl.size(); i++)
    {
        weights[i] = weights[i]/sum;
    }

    free(pos);
    kd_res_free(res);
}

double update_edges(vertex *from)
{

    state snew = system(from->s, 1);
    
    double *tosearch;
    tosearch = snew.x;

    kdres *res;
    double hits = 0, tot=0;
    int edge_num = from->edgeout.size();
    double edgeprob[100] = {0};
    int edgehits[100] = {0};

    double *pos;
    pos = (double *)malloc(sizeof(state));
    
    res = kd_nearest_range(mini_tree, tosearch, BOWLR);
    while( !kd_res_end(res))
    {
        minis *m = (minis *)kd_res_item(res, pos);
        
        // rewire mini-samples
        m->parent = nearest_vertex(m->s);
        
        if(m->parent != from){
            int index = edge_index_connected_vertex(from, m->parent);
            if( (index != -1) && (index < edge_num) )
            {
                if(index >= edge_num){
                    cout<<"Aborted: index >= edge_num"<<endl;
                    exit(1);
                }
                edgeprob[index] += noise_func(snew, m->parent->s, spro);
                edgehits[index]++;
            }
        }
        else
            hits++;
        
        kd_res_next(res);
        tot++;
    }
    
    double prob_sum = 0;
    // get average mass for each edge
    for(int i=0; i< edge_num; i++)
    {
        if(edgehits[i] != 0)
            edgeprob[i] = edgeprob[i]/(float)edgehits[i];
        else
            edgeprob[i] = 0;
        prob_sum += edgeprob[i];
    }
    for(int i=0; i< edge_num; i++)
    {
        (from->edgeout[i])->prob = edgeprob[i]/prob_sum;
    }

    kd_res_free(res);
    from->voronoi_area = hits/((float)tot)*BOWLR;       // rough estimate of voronoi region's area
    return hits;
}

void update_viterbi(vertex *v, vector<vertex *> nodesinbowl, vector<double> weights, double obs_time)
{
    double obs_prob = 1.0;

    for(int i=0; i< (int)nodesinbowl.size(); i++)
    {
        if (nodesinbowl[i] == v)
            obs_prob = weights[i];
    }

    // update viterbi
    for(int i=0; i < (int)v->edgein.size(); i++)
    {
        edge *etmp = v->edgein[i];
        vertex *vtmp = etmp->from;
        double prob_tmp = 0;

        prob_tmp = (vtmp->prob) * (etmp->prob) * obs_prob;
        /*
        if( ( (vtmp->t) < obs_time) && (vtmp-> t + etmp->delt >= obs_time))
            prob_tmp = (vtmp->prob) * (etmp->prob) * obs_prob;
        else
            prob_tmp = (vtmp->prob) * (etmp->prob);
        */
        //cout<<"vtmp->prob: "<<vtmp->prob<<" etmp->prob: "<<etmp->prob<<" prob_obs: "<<prob_obs(vtmp, etmp)<<endl;
        //cout<<prob_obs(vtmp, etmp)<<endl;
        if( prob_tmp > (v->prob))
        {
            v->prob = prob_tmp;
            //v->t = (vtmp->t) + (etmp->delt);
            v->prev = vtmp;
        }
    }
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

void add_major_sample(vertex *v, int is_obs)
{
    double *pos;
    pos = (double *)malloc(sizeof(state));

    if(rrg.num_vert != 0)
    {
        if(!is_obs)
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
            // make edges, no probab, update later
            edge *e1 = new edge(v, v1, 0);
            edge *e2 = new edge(v1, v, 0);
            
            // write edges
            v->edgeout.push_back(e1);
            v->edgein.push_back(e2);
            v1->edgeout.push_back(e2);
            v1->edgein.push_back(e1);
        }
        kd_res_next(res);
    }
    kd_res_free(res);
    free(pos);

    update_edges(v);
    
    // other end of this edge is x0, find two edges of x0 (with vertices x1, x2) between which "from" lies
    // update the weights of only those two edges
    // based on? -- re-calci sys_noise_func for edge x0-from wrt x1, x2 & from
    for(int i=0; i< (int)v->edgeout.size(); i++)
    {
        vertex *vtmp = (v->edgeout[i])->to;
        update_edges(vtmp);
    }

    add_mini_samples();
}

int main()
{
    cout.precision(4);
    srand(time(0));
    state_tree = kd_create(NUM_DIM);
    mini_tree = kd_create(NUM_DIM);
    //halton_init();
    
    double ts = get_msec();
    state x0; x0.x[0] = 1.0, x0.x[1] = -1.0;
    x.push_back(x0);
    y.push_back(obs(x0, 0));
    for(int i=0; i<N; i++){
        x.push_back(system(x.back(), 0));
        y.push_back(obs(x.back(), 0));
    }

    vector<vertex *> nodesinbowl;
    vector<double> weights;
    for(int i=0; i<1; i++)
    {
        vertex *v = new vertex(x[i], i, 1.0);
        add_major_sample(v, 1);
        update_obs_prob(y[i], nodesinbowl, weights);
        update_viterbi(v, nodesinbowl, weights, i*dt);
        cout<<"wrote first sample prob: "<<v->prob<<endl;
    }
    
    int dM = 1;
    for(int i=1; i<N; i++)
    {
        vector<vertex*> states_added_in_loop;
        
        // add the obs as the state
        vertex *v = new vertex(y[i], i, 0);
        add_major_sample(v, 1);
        states_added_in_loop.push_back(v);

        // add some more states
        for(int j=0; j<dM; j++)
        {
            vertex *v1 = new vertex(sample(), i+j/(float)dM, 0);
            add_major_sample(v1, 0);
            states_added_in_loop.push_back(v1);
        }

        // update observation prob
        update_obs_prob(y[i], nodesinbowl, weights);
        
        for(int j=0; j< (int)states_added_in_loop.size(); j++)
            update_viterbi(v, nodesinbowl, weights, i*dt);
    }
    plot_rrg(); 
        
    cout<<"dt: "<<get_msec() - ts<<endl; 
    kd_free(state_tree);
    kd_free(mini_tree); 
    return 0;
}

