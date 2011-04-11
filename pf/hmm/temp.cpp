// first dim is time, 1-d system
#define NUM_DIM         (2)

#define TMAX            (1.0)
#define TMIN            (0.0)

#define XMAX            (0.2)
#define XMIN            (-0.2)

#define randf           (rand()/(RAND_MAX + 1.0))

#define EXTEND_DIST     (0.5)
#define dt              (0.125)
#define dM              (100)
#define N               ((int)(TMAX/dt))
#define sobs            (0.01)
#define spro            (0.01)
#define BETA            (100)
#define GAMMA           (1.0)
#define BOWLGAMMA       (GAMMA*pow(log(rrg.num_vert)/(float)(rrg.num_vert), 1/(NUM_DIM)))
#define BOWLR           ( (BOWLGAMMA >= 0.01) ? BOWLGAMMA : 0.01)
#define PI              (3.14156)

#include "common.h"
#include "kdtree.h"

// halton
int seed[NUM_DIM] = {0, 0};
int base[NUM_DIM] = {2, 3};
int step = 0;

kdtree *state_tree, *mini_tree;
graph rrg;
vector<state> x(101);
vector<state> simx;
vector<minis *> mini_samples;

void halton_init()
{
    halton_dim_num_set (NUM_DIM);
    halton_step_set (step);
    halton_seed_set(seed);
    halton_base_set(base);
};

double normal_val(double mean, double var, double tocalci)
{
    double temp = exp(-0.5*sq(mean-tocalci)/var);
    return 1/(2*PI)/sqrt(var)*temp;
}

/*
 * simulate brownian motion
 */
state system(state s, double time){
    
    state t;
    double tmp1=0, tmp2=0;
    
    // cov is time
    randn(0, time, tmp1, tmp2);
    
    t.x[0] = s.x[0] + time;     // add time
    t.x[1] = s.x[1] + tmp1;     // xt+1 = xt + N(0, dt)
    return t;
}

state obs(state s){
    state t;
    double tmp1=0, tmp2=0;
    randn(0, sobs, tmp1, tmp2);
    
    t.x[0] = s.x[0];            // time is same
    t.x[1] = s.x[1] + tmp1;     // add noise to obs
    return t;
}

void plot_rrg()
{
    ofstream rrgout("rrg.dat");
    ofstream rrgpout("rrgp.dat");
    cout<<"write rrg"<<endl;

    for(vector<vertex*>::iterator i = rrg.vlist.begin(); i != rrg.vlist.end(); i++)
    {
        vertex *tstart = (*i);
        /*    
        for(vector<edge*>::iterator eo = tstart->edgeout.begin(); eo != tstart->edgeout.end(); eo++)
        {
            vertex *tend = (*eo)->to;
            edge *etmp = (*eo);
            
            //draw the edge
            rrgout<<tstart->s.x[0]<<"\t"<<tstart->s.x[1]<<"\t"<<tend->s.x[0]<<"\t"<<tend->s.x[1]<<"\t"<<etmp->prob<<"\t"<<etmp->delt<<endl;
        }
        */
        rrgpout<<tstart->s.x[0]<<"\t"<<tstart->s.x[1]<<endl;
    }
    rrgout.close();
    rrgpout.close();
}

void plot_traj()
{
    ofstream traj("traj.dat");
    
    traj<<"system"<<endl;
    cout<<"system"<<endl;
    for(int i=0; i< (int)x.size()-1; i++)
        traj<<x[i].x[0]<<"\t"<<x[i].x[1]<<endl;
    
    traj<<"sim"<<endl;
    cout<<"sim"<<endl;
    for(unsigned int i=0; i< simx.size(); i++)
        traj<<simx[i].x[0]<<"\t"<<simx[i].x[1]<<endl;
    
    traj.close();
}

vertex* nearest_vertex(state s)
{
    kdres *resn;
    resn = kd_nearest(state_tree, s.x);
    if(kd_res_end(resn)){
        cout<<"Error: no nearest"<<endl;
        exit(1);
    }
    vertex *v = (vertex*)kd_res_item_data(resn);
    kd_res_free(resn);

    return v;
}

/*
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

// locate all vertices within bowlr and write prob of that vertex in weights
void update_obs_prob(state yt, vector<vertex*> &nodesinbowl, vector<double> &weights)
{
    nodesinbowl.clear();
    weights.clear();
    vector<int> num_children;

    double *obs_state = yt.x;
    double pos[NUM_DIM] = {0};

    kdres *res;
    res = kd_nearest_range(mini_tree, obs_state, BOWLR);
    while( !kd_res_end(res))
    {
        minis *m = (minis *)kd_res_item(res, pos);

        double noise_tmp = noise_func(yt, m->s, sobs);
        
        if(nodesinbowl.size() == 0)
        {
            nodesinbowl.push_back(m->parent);
            weights.push_back(noise_tmp);
            num_children.push_back(1);
        }
        else
        {
            vector<vertex *>::iterator iter = find(nodesinbowl.begin(), nodesinbowl.end(), m->parent);
            if(( iter != nodesinbowl.end()) || (nodesinbowl.back() == (m->parent)))
            {
                weights[iter - nodesinbowl.begin()] += noise_tmp;
                num_children[iter - nodesinbowl.begin()] += 1;
            }
            else
            {
                num_children.push_back(1);
                nodesinbowl.push_back(m->parent);
                weights.push_back(noise_tmp);
            }
        }
        kd_res_next(res);
    }
    kd_res_free(res);
    
    // normalize
    cout<<"nodesinbowl: "<< nodesinbowl.size() <<endl;
    double weights_sum = 0;
    for(int i=0; i< (int)nodesinbowl.size(); i++)
    {
        vertex *v = nodesinbowl[i];
        weights[i] = (v->voronoi_area)*weights[i]/num_children[i];
        weights_sum += weights[i];
        //cout<<"v: ["<< v->s.x[0] <<", "<< v->s.x[1]<<"]\t"<<weights[i]<<"\tnum_edges: "<< v->edgeout.size() <<"\t"<< noise_func(yt, v->s, sobs)<<endl;        
    }
    for(int i=0; i< (int)nodesinbowl.size(); i++)
    {
        weights[i] = weights[i]/weights_sum;
    }
}
*/

void update_edges(vertex *from)
{
    int edge_num = from->edgeout.size();
    //cout<<"update edges: "<<edge_num<<" "<<from->s.x[0]<<" "<<from->s.x[1]<<endl;
    
    for(int i=0; i< edge_num; i++)
    {
        edge *e = from->edgeout[i];
        double totprob = 0;
        for(int j=0; j< edge_num; j++)
        {
            edge *etmp = from->edgeout[j];
            totprob += normal_val(from->s.x[1], e->delt, etmp->to->s.x[1]);
        }
        e->prob = normal_val(from->s.x[1], e->delt, e->to->s.x[1])/totprob;
    }

    double maxprob = 0;
    edge *bestein = NULL;
    for(unsigned int i=0; i< from->edgein.size(); i++)
    {
        if(from->edgein[i]->prob > maxprob)
        {
            maxprob = from->edgein[i]->prob;
            bestein = from->edgein[i];
        }
    }
    if(bestein != NULL)
        from->best_in = bestein;
    else
        cout<<"best_in is NULL"<<endl;

    maxprob = 0;
    edge *besteout = NULL;
    for(unsigned int i=0; i< from->edgeout.size(); i++)
    {
        if(from->edgeout[i]->prob > maxprob)
        {
            maxprob = from->edgeout[i]->prob;
            besteout = from->edgeout[i];
        }
    }
    
    if(besteout != NULL)
        from->best_out = besteout;
    else
        cout<<"best out is NULL"<<endl;
}

/*
void update_viterbi( const vector<vertex *> &nodesinbowl, const vector<double> &weights, double obs_time)
{
    for(int i=0; i< (int)nodesinbowl.size(); i++)
    {
        vertex *v = nodesinbowl[i];
        double obs_prob = weights[i];
        
        vertex *to_push_prev = NULL;
        double to_push_prob = 0;
        double to_push_alpha = 0;
        double to_push_time = 0;

        int latest_index = -1;
        for(int i=0; i < (int)v->edgein.size(); i++)
        {
            edge *etmp = v->edgein[i];
            vertex *vtmp = etmp->from;
            double prob_tmp = 0;
            
            // find latest time before obs_time in vtmp->prob
            latest_index = -1;
            for(int j = (int)vtmp->t.size() -1; j >= 0; j--)
            {
                //cout<<"obs_time: "<<obs_time<<"\tlooking at: "<< vtmp->s.x[0]<<" "<<vtmp->s.x[1]<<" "<<vtmp->t[j]<<" "<<etmp->prob<<" "<<etmp->delt<<endl;
                if( (vtmp->t[j] + etmp->delt) < obs_time )
                {
                    latest_index = j;
                    //cout<<"found latest index viterbi with t: "<<vtmp->t[j]<<" for obs_time: "<<obs_time<<endl;
                    break;
                }
            }
            
            if( (latest_index != -1) && (vtmp->t.size() != 0))
            {
                prob_tmp = vtmp->prob[latest_index] * (etmp->prob);
                to_push_alpha += prob_tmp;

                if( prob_tmp > (to_push_prob))
                {
                    to_push_prob = prob_tmp;
                    to_push_prev = vtmp;
                    to_push_time = vtmp->t.back() + etmp->delt;
                }
            }
        }
        if( to_push_prob != 0)
        {
            //cout<<"to push prob: "<<to_push_prob<<" "<<obs_prob<<endl;
            v->prob.push_back(to_push_prob * obs_prob);
            v->t.push_back(to_push_time);
            v->prev.push_back(to_push_prev);
            v->alpha.push_back(to_push_alpha * obs_prob);
            //cout<<"pushed "<<v->s.x[0]<<" "<<v->s.x[1]<<endl;
        }
    }
}

void add_mini_samples(state around_which)
{
    for(int j=0; j<BETA; j++)
    {
        minis *m;
        m = new (minis);
        m->s = sample(around_which, BOWLR);
        
        mini_samples.push_back(m);
        kd_insert(mini_tree, m->s.x, m);
        
        vertex *v = nearest_vertex(m->s);
        m->parent = v;
        v->num_child = v->num_child + 1;
    }
}

double time_to_go(vertex *from, vertex* to)
{
    state sf = from->s;
    state st = to->s;
    
    double a = 0, prev_a = 0;
    double del = dt/10;
    double b = del;
    bool min_found = 0;
    while( min_found == 0)
    {
        if( dist(system(sf, b, 1), st) < dist(system(sf, a, 1), st) )
        {
            prev_a = a;
            a = b;
            b += del;
        }
        else
        {
            min_found = 1;
            a = prev_a;
        }
    }
    // min between a and b now, quadratic fit
    // polyfit from mathematica
    double c = (a+b)/2;
    //cout<<"a: "<<a<<" b: "<<b<<" c: "<<c<<endl;
    double t1 = dist( system(sf, a, 1), st);
    double t2 = dist( system(sf, b, 1), st);
    double t3 = dist( system(sf, c, 1), st);
    //cout<<"t1: "<<t1<<" t2: "<<t2<<" t3: "<<t3<<endl;
    double num = t3*(a-b)*(a+b) + t1*(b-c)*(b+c) + t2*(c-a)*(c+a);
    double den = 2*(t3*(a -b) + t1*(b-c) + t2*(-a+c));
    double time = num/den;
    
    if(time > 0)
    {
        //cout<<"min_time: "<<time<<" dist: "<<dist(system(sf, time, 1), st)<<endl;
        return time;
    }
    else
        return -1;
}
*/

void add_major_sample(vertex *v)
{
    double pos[NUM_DIM] = {0};

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
    //cout<<"Adding sample: "<<v->s.x[0]<<" "<<v->s.x[1]<<endl;

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
            double e1t = v1->s.x[1] - v->s.x[1];
            double e2t = v->s.x[1] - v1->s.x[1];

            //cout<<"e1t: "<<e1t<<" e2t: "<<e2t<<endl;
            // make edges, no probab, update later
            // write edges
            if( e1t > 0)
            {
                edge *e1 = new edge(v, v1, e1t);
                v->edgeout.push_back(e1);
                v1->edgein.push_back(e1);
                //cout<<"wrote e: "<<v->s.x[0]<<" "<<v->s.x[1]<<" to "<<v1->s.x[0]<<" "<<v1->s.x[1]<<endl;  
            }
            else if( e2t > 0)
            {
                edge *e2 = new edge(v1, v, e2t);
                v->edgein.push_back(e2);
                v1->edgeout.push_back(e2);
                //cout<<"wrote e: "<<v1->s.x[0]<<" "<<v1->s.x[1]<<" to "<<v->s.x[0]<<" "<<v->s.x[1]<<endl;
            }
        }
        kd_res_next(res);
    }

    // update edges for this vertex
    update_edges(v);

    // update edges for all incoming vertices to this
    for(int k=0; k< (int)v->edgein.size(); k++)
    {
        edge *e = v->edgein[k];
        vertex *vtmp = e->from;
        update_edges(vtmp);
    }
}

int main()
{
    cout.precision(4);
    srand(time(0));
    state_tree = kd_create(NUM_DIM);
    mini_tree = kd_create(NUM_DIM);
    
    state x0; x0.x[0] = 0.0, x0.x[1] = 0.0;
    for(int j=0; j<100; j++)
    {
        x[j].x[1] = 0;
    }
    for(int j=0; j<100; j++)
    {
        vector<state> xt;
        xt.push_back(x0);
        x[0].x[1] = x[0].x[1] + xt.back().x[1];
        x[0].x[0] = xt.back().x[0];

        for(int i=1; i<= 100; i++){
            xt.push_back(system(xt.back(), 0.01));
            x[i].x[1] = x[i].x[1] + xt.back().x[1];
            x[i].x[0] = xt.back().x[0];
        }
    }
    for(int j=0; j < 101; j++)
    {
        x[j].x[1] = x[j].x[1]/100;
    } 

    for(int i=0; i< 100; i++)
    {
        vertex *v = new vertex(sample());
        add_major_sample(v);
        if(i %100 == 0)
            cout<<"i: "<<i<<endl;
    }
    
    // simulate from starting point on the graph
    double time = 0;
    vertex *curr = rrg.vlist[0];
    simx.push_back(curr->s);             // start sim at x0
    
    while(time < 1)
    {
        // move to node with largest prob
        edge *e = curr->best_out;
        curr = e->to;
        time += e->delt;
        simx.push_back(curr->s);
    }
    cout<<"simx size: "<<simx.size()<<endl;

    plot_rrg();
    plot_traj();
    
    kd_free(state_tree);
    kd_free(mini_tree); 
    return 0;
}

