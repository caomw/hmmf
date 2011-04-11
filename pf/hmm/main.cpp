
#define NUM_DIM     (2)

#define XMAX        (1.0)
#define XMIN        (0.6)

#define YMAX        (-0.6)
#define YMIN        (-1.0)

#define randf       (rand()/(RAND_MAX + 1.0))

#define EXTEND_DIST (0.5)
#define dt      (0.10)
#define dM      (20)
#define N       ((int)1.0/dt)
#define sobs    (0.01)
#define spro    (0.01)
#define BETA    (100)
#define GAMMA           (XMAX - XMIN)
#define BOWLGAMMA       (GAMMA*sqrt(log(rrg.num_vert)/(float)(rrg.num_vert)))
#define BOWLR           ( (BOWLGAMMA >= 3*spro) ? BOWLGAMMA : 3*spro)
#define PI              (3.14156)

#include "common.h"
#include "kdtree.h"

// halton
int seed[NUM_DIM] = {0, 0};
int base[NUM_DIM] = {2, 3};
int step = 0;

kdtree *state_tree, *mini_tree;
graph rrg;
vector<state> x, y;
vector<minis *> mini_samples;

void halton_init()
{
    halton_dim_num_set (NUM_DIM);
    halton_step_set (step);
    halton_seed_set(seed);
    halton_base_set(base);
};

// time is assumed to be small, so one step euler integration
state system(state s, double time, int is_clean){
    state t;
    if(!is_clean){
        double tmp1=0, tmp2=0;
        randn(0, spro, tmp1, tmp2);
        t.x[0] = s.x[0] + time*( -0.5*s.x[0] + tmp1);
        t.x[1] = s.x[1] + time*( -0.9*s.x[1] + tmp2);
    }
    else{
        t.x[0] = s.x[0] + time*( -0.5*s.x[0]);
        t.x[1] = s.x[1] + time*( -0.9*s.x[1]);
    }
    return t;
}

state obs(state s, int is_clean){
    state t;
    if(!is_clean){
        double tmp1=0, tmp2=0;
        randn(0, sobs, tmp1, tmp2);
        t.x[0] = s.x[0] + tmp1;
        t.x[1] = s.x[1] + tmp2;
    }
    else{
        t.x[0] = s.x[0];
        t.x[1] = s.x[1];
    }
    return t;
}

void print_rrg()
{
    cout<<endl<<"---printing rrg---"<<endl;
    for(int i=0; i< rrg.num_vert; i++)
    {
        vertex *v = rrg.vlist[i];
        cout<<"\n"<<i<<" vert: "<< v->s.x[0]<<" "<< v->s.x[1]<<" neigh: "<<v->edgeout.size()<<" sizet: "<<v->t.size()<<endl;
        for(int j=0; j< (int)v->t.size(); j++)
        {
            vertex *vprevtmp = v->prev[j];
            if(vprevtmp != NULL)
                cout<<v->t[j]<<"\t"<<v->prob[j]<<"\t"<<"["<< vprevtmp->s.x[0]<<" "<< vprevtmp->s.x[1]<<"]"<<"\t"<<v->alpha[j]<<endl;
        }
        cout<<"edges -- "<<endl;
        if( v->edgeout.size() != 0)
        {
            for(int j=0; j< (int)v->edgeout.size(); j++)
            {
                //cout<<"e: "<<j<<"\t"<<v->edgeout[j]->prob<<"\tto: "<<v->edgeout[j]->to->s.x[0]<<" "<< v->edgeout[j]->to->s.x[1]<<endl;
            }
        }
    }
}

void plot_rrg()
{
    ofstream rrgout("rrg.dat");
    ofstream rrgpout("rrgp.dat");
    ofstream miniout("miniout.dat");

    for(vector<vertex*>::iterator i = rrg.vlist.begin(); i != rrg.vlist.end(); i++)
    {
        vertex *tstart = (*i);
        for(vector<edge*>::iterator eo = tstart->edgeout.begin(); eo != tstart->edgeout.end(); eo++)
        {
            vertex *tend = (*eo)->to;
            
            //draw the edge
            rrgout<<tstart->s.x[0]<<"\t"<<tstart->s.x[1]<<"\t"<<tend->s.x[0]<<"\t"<<tend->s.x[1]<<endl;
        }
        rrgpout<<tstart->s.x[0]<<"\t"<<tstart->s.x[1]<<endl;
    }
    for(int i=0; i< (int)mini_samples.size(); i++)
    {
        miniout<< mini_samples[i]->s.x[0]<<"\t"<< mini_samples[i]->s.x[1]<<endl;
    }
    miniout.close();
    rrgout.close();
    rrgpout.close();
}

void plot_traj()
{
    ofstream traj("traj.dat");

    traj<<"system"<<endl;
    for(int i=0; i< (int)x.size()-1; i++)
        traj<<x[i].x[0]<<"\t"<<x[i].x[1]<<endl;
    
    traj<<"observation"<<endl;
    for(int i=0; i< (int)y.size()-1; i++)
        traj<<y[i].x[0]<<"\t"<<y[i].x[1]<<endl;
    
    /*
    // best traj using alphas
    traj<<"best_path"<<endl;
    for(int time= (N-1); time >= 0; time--)
    {
        double maxpx = 0, maxpy = 0, maxp = 0;
        double avgx=0, avgy =0, totalpha = 0;
        for(int i =0; i<rrg.num_vert; i++)
        {
            vertex *v = rrg.vlist[i];
            if(v->t.size() != 0)
            {
                if( fabs(v->t.back() - time*dt) <= 0.001)
                {
                    // maximum likelihood estimate
                    if( (v->alpha).back() > maxp)
                    {
                        maxp = (v->alpha).back();
                        maxpx = v->s.x[0];
                        maxpy = v->s.x[1];
                        //cout<<"changed max: "<<maxpx<<" "<<maxpy<<endl;    
                    }
                    
                    // weighted mean
                    avgx += (v->s.x[0]) *(v->alpha).back();
                    avgy += (v->s.x[1]) *(v->alpha).back();
                    totalpha += (v->alpha).back();
                    //cout<< v->alpha.back()<<endl;

                    v->t.pop_back();
                    v->prob.pop_back();
                    v->alpha.pop_back();
                }
            }
        }
        if( (totalpha == 0) || (maxp == 0) )
            cout<<"didn't find a single vert with dt: "<<time*dt<<endl;
        
        avgx /= totalpha;
        avgy /= totalpha;
        
        cout<< avgx<<"\t"<< avgy<<endl;
        traj<< avgx<<"\t"<< avgy<<endl;
    }
    */ 
 
    // get best trajectory using deltas

    vertex *best = NULL;
    double max_prob = 0;
    double mptime = 5000;
    int mptime_index = 5000;

    for(int i=0; i < rrg.num_vert; i++)
    {
        vertex *v = rrg.vlist[i];
        if (v->t.size() != 0)
        {
            if(v->t.back() == (N-1)*dt)
            {
                if(v->prob.back() > max_prob)
                {
                    max_prob = v->prob.back();
                    best = v;
                    mptime = v->t.back();
                    mptime_index = v->t.size() - 1;
                }
            }
        }
    }
    if( (mptime_index == 5000) || (best == NULL))
    {
        cout<<"Error no vertex with latest obs"<<endl;
        exit(0);
    }

    traj<<"best_path"<<endl;
    while(mptime != 0)
    {
        //cout<<"mpt: "<< mptime<<" mpti: "<<mptime_index<<endl;
        traj<< best->s.x[0]<<"\t"<< best->s.x[1]<<endl;
        //cout<<"curr_best: "<<best<<"\t"<< best->s.x[0]<<"\t"<< best->s.x[1]<<" neigh: "<< best->edgein.size()<<endl;
        for(int j=0; j< (int)best->t.size(); j++)
        {
            //cout<<best->t[j]<<"\t"<<best->prob[j]<<"\t"<<best->prev[j]<<endl;
        }
        //cout<<"dist between curr and new best: "<< dist(best->s, best->prev[mptime_index]->s)<<endl;
        best = best->prev[mptime_index];
        //cout<<"new_best: "<<best<<" neigh: "<< best->edgein.size()<<endl;
  
        for(int j=0; j< (int)best->t.size(); j++)
        {
            //cout<<best->t[j]<<"\t"<<best->prob[j]<<"\t"<<best->prev[j]<<endl;
        }
        
        bool changed_mptime = 0;
        // get the closest time in prev less than mptime
        for(int i= best->t.size()-1; i >= 0; i--)
        {
            if (best->t[i] < mptime)
            {
                changed_mptime = 1;
                mptime_index = i;
                //cout<<"break with mpti: "<<i<<endl;
                break;
            }
        }
        if(!changed_mptime)
        {
            //cout<<"couldn't change mptime to "<<mptime - dt<<endl;
            return;
        }
        mptime = best->t[mptime_index];
    }
    traj<< x[0].x[0]<<"\t"<< x[0].x[1]<<endl;

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

double noise_func(const state s, const state s1, double sigma)
{
    double J = 0;
    for(int i=0; i<NUM_DIM; i++)
        J += (s.x[i] - s1.x[i])*(s.x[i] - s1.x[i])/2.0/sigma;
    
    double tmp = 1/pow(2*PI, NUM_DIM/2.0)/pow(sigma, NUM_DIM/2.0)*exp(-J);
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

/*
 * update prob of all edges going out from vertex from
 */
double update_edges(vertex *from)
{
    int edge_num = from->edgeout.size();
    double bowlsize = 0;
    for(int i=0; i< edge_num; i++)
    {
        edge *etmp = from->edgeout[i];
        double bowltmp = dist(from->s, etmp->to->s);
        if(bowltmp > bowlsize)
            bowlsize = bowltmp;
    }
    bowlsize *= 2;

    double hits = 0, tot=0;
    vector<double> edgeprob(edge_num, 0.0);
    vector<int> edgehits(edge_num, 0);
    state stmp;
    stmp.x[0] = 0; stmp.x[1] = 0.0;
    vector<state> newstates(edge_num, stmp);
    for(int i=0 ; i< edge_num; i++)
    {
        newstates[i] = system(from->s, from->edgeout[i]->delt, 0);
    }
    
    double pos[NUM_DIM] = {0};
    kdres *res;
    res = kd_nearest_range(mini_tree, from->s.x, bowlsize);
    while( !kd_res_end(res))
    {
        minis *m = (minis *)kd_res_item(res, pos);

        // rewire mini-samples
        m->parent = nearest_vertex(m->s);

        if(m->parent != from)
        {
            int index = edge_index_connected_vertex(from, m->parent);
            if( (index != -1) && (index < edge_num) )
            {
                if(index >= edge_num){
                    cout<<"Aborted: index >= edge_num"<<endl;
                    exit(1);
                }
                edgeprob[index] += noise_func(newstates[index], m->parent->s, spro);
                edgehits[index]++;
            }
        }
        else
        {
            hits++;
        }
        kd_res_next(res);
        tot++;
    }
    kd_res_free(res);
    
    double prob_sum = 0;
    // get average mass for each edge
    for(int j=0; j< edge_num; j++)
    {
        if(edgehits[j] != 0)
        {
            double hitt = edgehits[j];
            edgeprob[j] = ((from->edgeout[j])->to->voronoi_area)* (edgeprob[j])/hitt;
        }
        else
            edgeprob[j] = 0;
        
        prob_sum += edgeprob[j];
    }
    
    // normalize
    for(int j=0; j< edge_num; j++)
    {
        edge *etmp = from->edgeout[j];
        if (prob_sum != 0)
            etmp->prob = edgeprob[j]/prob_sum;
        else
            etmp->prob = 0;
    }
    from->voronoi_area = hits/((float)tot)*BOWLR*BOWLR*PI;       // rough estimate of voronoi region's area
    return hits;
}

/*
 * v is the new vertex just sampled
 * obs_time is the observation's time for which we are updating
 */
void update_viterbi( const vector<vertex *> &nodesinbowl, const vector<double> &weights, double obs_time)
{
    for(int i=0; i< (int)nodesinbowl.size(); i++)
    {
        vertex *v = nodesinbowl[i];
        double obs_prob = weights[i];
        
        vertex *to_push_prev = NULL;
        double to_push_prob = 0;
        double to_push_alpha = 0;

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
                //cout<<"obs_time: "<<obs_time<<"\tlooking at: "<< vtmp->t[j]<<" "<<etmp->prob<<" "<<etmp->delt<<endl;
                if( vtmp->t[j] < obs_time )
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
                }
            }
            /*
            else
            {
                cout<<"latest_index = -1"<<endl;
                exit(0);
            }
            */
        }
        if( to_push_prob != 0)
        {
            //cout<<"to push prob: "<<to_push_prob<<" "<<obs_prob<<endl;
            v->prob.push_back(to_push_prob * obs_prob);
            v->t.push_back(obs_time);
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
    
    double best_delt = 0, best_dist = 100.0;
    int num_iter = 0;
    double deltmax = 1.0;
    double deltmin = dt/10;
    double delt = (deltmax + deltmin)/2.0;
    state stmp = system(sf, delt, 1);
    state stmpmax = system(sf, deltmax, 1);
    state stmpmin = system(sf, deltmin, 1);
    double curr_dist = dist(stmp, st);
    while( (curr_dist > 0.01) && (num_iter <= 20) )
    {
        if( dist(stmp, stmpmax) < dist(stmp, stmpmin) )
        {
            deltmin = delt;
            stmpmin = system(sf, deltmin, 1);
        }
        else
        {
            deltmax = delt;
            stmpmax = system(sf, deltmax, 1);
        }
        delt = (deltmax + deltmin)/2.0;
        num_iter++;

        stmp = system(sf, delt, 1);
        curr_dist = dist(stmp, st);
        if( curr_dist < best_dist)
        {
            best_dist = curr_dist;
            best_delt = delt;
        }
    }
    //cout<<"best_dist: "<<best_dist<<endl;
    return best_delt;
}


void add_major_sample(vertex *v, int is_obs, state around_which)
{
    double pos[NUM_DIM] = {0};
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
            double e1t = time_to_go(v, v1);
            double e2t = time_to_go(v1, v);
            
            //cout<<"e1t: "<<e1t<<" e2t: "<<e2t<<endl;
            // make edges, no probab, update later
            // write edges
            if( e1t > 0)
            {
                edge *e1 = new edge(v, v1, e1t);
                v->edgeout.push_back(e1);
                v1->edgein.push_back(e1);
            }
            if( e2t > 0)
            {
                edge *e2 = new edge(v1, v, e2t);
                v->edgein.push_back(e2);
                v1->edgeout.push_back(e2);
            }
        }
        kd_res_next(res);
    }
    kd_res_free(res);
 
    // add mini samples
    add_mini_samples(around_which);

    update_edges(v);
    
    // other end of this edge is x0, find two edges of x0 (with vertices x1, x2) between which "from" lies
    // update the weights of only those two edges
    // based on? -- re-calci sys_noise_func for edge x0-from wrt x1, x2 & from
    for(int k=0; k< (int)v->edgeout.size(); k++)
    {
        vertex *vtmp = (v->edgeout[k])->to;
        update_edges(vtmp);
    }
}

int main()
{
    cout.precision(4);
    srand(time(0));
    state_tree = kd_create(NUM_DIM);
    mini_tree = kd_create(NUM_DIM);
    
    state x0; x0.x[0] = 1.0, x0.x[1] = -1.0;
    x.push_back(x0);
    y.push_back(obs(x0, 1));
    for(int i=0; i<N; i++){
        x.push_back(system(x.back(), dt, 0));
        y.push_back(obs(x.back(), 0));
    }
    
    vector<vertex *> nodesinbowl;
    vector<double> weights;
    for(int i=0; i<1; i++)
    {
        vertex *v = new vertex(y[i]);
        add_major_sample(v, 1, y[i]);
        add_mini_samples(y[i]);

        // set up initial estimate
        v->prob.push_back(1);
        v->t.push_back(0*dt);
        // note this
        v->prev.push_back(v);
        v->alpha.push_back(1*noise_func(y[i], v->s, sobs));
    }
    //print_rrg();
    //getchar();

    for(int i=1; i<N; i++)
    {
        cout<<"-------------------------------"<<endl;
        double ts = get_msec();
        
        // add some more states
        for(int j=0; j<dM; j++)
        {
            vertex *v1 = new vertex(sample(y[i], BOWLR));
            add_major_sample(v1, 0, y[i]);
            //cout<<"major: "<<rrg.num_vert<<endl;
        }

        //cout<<"updating obs: ["<< y[i].x[0]<<", "<< y[i].x[1] <<"]"<<endl;
        update_obs_prob(y[i], nodesinbowl, weights);
        update_viterbi(nodesinbowl, weights, i*dt);
        ts = get_msec() - ts;
        cout<<"obs: "<<i<<"\t"<<y[i].x[0]<<" "<<y[i].x[1]<<"\tdt: "<<ts<<"\tnum_vert: "<<rrg.num_vert<<"\tbowl: "<<BOWLR<<endl;

        //print_rrg();
        //getchar();
    }
    plot_rrg();
    plot_traj();
    
    for(int i=0; i< (int)mini_samples.size(); i++)
    {
        delete mini_samples[i];
    }
    kd_free(state_tree);
    kd_free(mini_tree); 
    return 0;
}

