// first dim is time, 2-d system
#define NUM_DIM             (3)

#define dM                  (200)
#define GAMMA               (1.0)
#define BOWLGAMMA           (GAMMA*pow(log(rrg.num_vert)/(float)(rrg.num_vert), 1.0/(NUM_DIM)))
#define BOWLR               ( (BOWLGAMMA >= 0.1) ? BOWLGAMMA : 0.1)
#define PI                  (3.14156)

#define randf               (rand()/(RAND_MAX + 1.0))

#include "common.h"
#include "kdtree.h"
#include "figtree.h"

// limits of states
double max_states[NUM_DIM] = {1.0, 0.15, 0.15};
double min_states[NUM_DIM] = {0.0, -0.15, -0.15};

state sample()
{
    state s;
    
    for(int i=0; i< NUM_DIM; i++)
    {
        s.x[i] = min_states[i] + randf*( max_states[i] - min_states[i]);
    }
    return s;
}

kdtree *state_tree;
graph rrg;
vector<state> x, y, best_path;
vector<state> simx;
vector<state> xhat(max_states[0]/0.01 + 1);

#define HOW_FAR             (0.2)
double INIT_VAR[2] = {1e-4, 1e-4};
double OBS_VAR[2] = {1e-2, 1e-2};
double PRO_VAR[2] = {1e-4, 1e-4};
double OBS_VAR_KALMAN[2] = { (OBS_VAR[0]*2 + 2*HOW_FAR*HOW_FAR), (OBS_VAR[1]*2 + 2*HOW_FAR*HOW_FAR) };

state system(state s, double time, int is_clean)
{
    double var[2] = { PRO_VAR[0]/2*(exp(2*time) - 1), PRO_VAR[1]/2*(exp(2*time) - 1)};
    double tmp[2] = {0};
    double mean[2] = {0};
    multivar_normal(mean, var, tmp, 2);

    if(is_clean)
    {
        tmp[0] = 0;
        tmp[1] = 0;
    }
    state t;
    t.x[0] = s.x[0] + time;                         // add time
    t.x[1] = exp(-1*time)*s.x[1] + tmp[0];
    t.x[2] = exp(-1*time)*s.x[2] + tmp[1];
    return t;
}

state obs(state s, int is_clean)
{
    state t;
    double tmp[2] = {0};
    double mean[2] = {0};
    multivar_normal(mean, OBS_VAR, tmp, 2);
    
    if(is_clean)
    {
        tmp[0] = 0;
        tmp[1] = 0;
    }
    t.x[0] = s.x[0];                    // time is same
    t.x[1] = s.x[1] + tmp[0];             // add noise to obs
    t.x[2] = s.x[2] + tmp[1];             // add noise to obs
    return t;
}

void plot_rrg()
{
    DepthTag dtag;
    
    ofstream rrgout("rrg.dat");
    ofstream rrgpout("rrgp.dat");
    lout<<"write rrg"<<endl;

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
        rrgpout<<tstart->s.x[0]<<"\t"<<tstart->s.x[1]<<"\t"<<tstart->s.x[2]<<"\t"<<tstart->prob<<endl;
    }
    rrgout.close();
    rrgpout.close();
}

void plot_traj()
{
    DepthTag dtag;
    
    ofstream traj("traj.dat");

    traj<<"system"<<endl;
    lout<<"system"<<endl;
    for(int i=0; i< (int)x.size()-1; i++)
        traj<<x[i].x[0]<<"\t"<<x[i].x[1]<<"\t"<<x[i].x[2]<<endl;

    traj<<"observation"<<endl;
    lout<<"observation"<<endl;
    for(int i=0; i< (int)y.size()-1; i++)
        traj<<y[i].x[0]<<"\t"<<y[i].x[1]<<"\t"<<y[i].x[2]<<endl;

    traj<<"best_path"<<endl;
    lout<<"best_path"<<endl;
    for(unsigned int i=0; i< best_path.size(); i++)
        traj<< best_path[i].x[0]<<"\t"<< best_path[i].x[1]<<"\t"<< best_path[i].x[2]<<endl;
    
    traj<<"kf_path"<<endl;
    lout<<"kf_path"<<endl;
    for(unsigned int i=0; i< xhat.size(); i++)
        traj<< xhat[i].x[0]<<"\t"<< xhat[i].x[1]<<"\t"<< xhat[i].x[2]<<endl;

    traj.close();
}

vertex* nearest_vertex(state s)
{
    kdres *resn;
    resn = kd_nearest(state_tree, s.x);
    if(kd_res_end(resn)){
        lout<<"Error: no nearest"<<endl;
        exit(1);
    }
    vertex *v = (vertex*)kd_res_item_data(resn);
    kd_res_free(resn);

    return v;
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

// writes the best and normalizes
void normalize_edges(vertex *from)
{
    DepthTag dtag;
    
    double maxprob = 0;
    edge *besteout = NULL;
    double totprob = 0;
    for(unsigned int i=0; i< from->edgeout.size(); i++)
    {
        if(from->edgeout[i]->prob > maxprob)
        {
            maxprob = from->edgeout[i]->prob;
            besteout = from->edgeout[i];
        }
        totprob += from->edgeout[i]->prob;
    }
    // normalize
    for(unsigned int i=0; i< from->edgeout.size(); i++)
    {
        edge *etmp = from->edgeout[i];
        etmp->prob = etmp->prob / totprob;
    }
}


// change using particles
// start particles from node till obs_time, update prob using obs, resample
// propagate till next node, find prob of next node using figtree
#define NUM_PARTICLES       (50)
double get_edge_prob_particles(state xinit, double till_obs, double post_obs, state xend, state obs)
{
    DepthTag dtag;
    
    state parts[NUM_PARTICLES];
    double weights[NUM_PARTICLES] = {1.0/(double)NUM_PARTICLES};

    double totw = 0, totw2 = 0;
    //lout<< "xinit: "<< xinit.x[1] << endl;

    for(int i=0; i< NUM_PARTICLES; i++)
        parts[i] = (  system(xinit, till_obs, 0) );        // propagated particles directly
    
    for(int i=0; i< NUM_PARTICLES; i++)
    {
        //lout<< "px: "<< parts[i].x[1] << " obs: " << obs.x[1] << endl;
        double mean[NUM_DIM - 1] = { obs.x[1], obs.x[2]};
        double tocalci[NUM_DIM - 1] = { parts[i].x[1], parts[i].x[2]};

        weights[i] = normal_val( mean, OBS_VAR, tocalci, 2 ); 
        
        //lout<< "w: " << weights[i] << endl;

        totw += weights[i];
        totw2 += sq(weights[i]);
    }
       
    double mean1 = 0, mean2 = 0, var1 = 0, var2 = 0;
    double src[ NUM_PARTICLES];
    for(int i=0; i< NUM_PARTICLES; i++)
    {
        parts[i] = system( parts[i], post_obs, 0 );            // get new states
        
        mean1 += (parts[i].x[1]) * weights[i];
        mean2 += (parts[i].x[2]) * weights[i];
    }
    mean1 = mean1/totw;
    mean2 = mean2/totw;

    for(int i=0; i< NUM_PARTICLES; i++)
    {
        var1 += weights[i]*sq (parts[i].x[1] - mean1);
        var2 += weights[i]*sq (parts[i].x[2] - mean2);
    }
    var1 = totw*var1/(totw*totw - totw2);
    var2 = totw*var2/(totw*totw - totw2);
    //lout<<" mean: "<< mean <<" var: "<< var <<endl;
    
    double mean[NUM_DIM-1] = {mean1, mean2};
    double var[NUM_DIM-1] = {var1, var2};
    double tocalci[NUM_DIM-1] = {xend.x[1], xend.x[2]};
    
    double target_w = normal_val( mean, var, tocalci, 2);
    //cout<< "target_w: "<< target_w << endl;
    return target_w;
}

/*
double get_edge_prob_kalman(state xinit, double till_obs, double post_obs, state xend, state obs)
{
    DepthTag dtag;
    
    // propagate qdot = 2aq + bwb to till_obs
    double newvar = PRO_VAR/2*(exp(2*till_obs) - 1);
    //lout<<"old_var: "<<newvar<<endl;

    // reduce q by q bt rinv b q due to obs
    // get new xhat
    state xhat =  system(xinit, till_obs, 1);
    double L = newvar/(newvar + OBS_VAR_KALMAN);
    
    state obs_predict;
    obs_predict.x[0] = xhat.x[0];
    obs_predict.x[1] = xhat.x[1];
    
    xhat.x[1] = xhat.x[1] + L*(obs.x[1] - obs_predict.x[1]);
    newvar = (1 - L)*newvar;
    //lout<<"after obs: "<<newvar<<endl;

    // propagate post_obs using xhat and q
    // final gaussian is centered at final_state with var = newvar
    newvar = exp(-2*post_obs)*newvar + PRO_VAR/2*(exp(2*post_obs) -1);
    state final_state = system(xhat, post_obs, 1);
    //lout<<"after growing: "<<newvar<<" xhat: "<<final_state.x[1]<<" actual: "<<etmp->to->s.x[1]<<endl;

    // change probab of etmp then
    return normal_val(xhat.x[1], newvar, xend.x[1]);
}
*/

void write_next_prev(vertex *v)
{
    double max_prob = 0;
    for(unsigned int i=0; i< v->edgeout.size(); i++)
    {
        edge *etmp = v->edgeout[i];
        if( etmp->prob > max_prob)
        {
            max_prob = etmp->prob;
            v->next = etmp->to;
        }
    }

    max_prob = 0;
    for(unsigned int i=0; i< v->edgein.size(); i++)
    {
        edge *etmp = v->edgein[i];
        if( ( (etmp->prob) *(etmp->from->prob) ) > max_prob)
        {
            max_prob = ( (etmp->prob) *(etmp->from->prob) );
            v->prev = etmp->from;
            v->prob = max_prob;
        }
    }
}


// locate all vertices within bowlr and write prob of that vertex in weights
void update_obs_prob(state yt)
{
    DepthTag dtag;
    
    double *obs_state = yt.x;
    double pos[NUM_DIM] = {0};

    kdres *res;
    res = kd_nearest_range(state_tree, obs_state, BOWLR);
    while( !kd_res_end(res))
    {
        vertex *v = (vertex *)kd_res_item(res, pos);

        for(unsigned int i=0; i < v->edgeout.size(); i++)
        {
            edge *etmp = v->edgeout[i];
            double till_obs = obs_state[0] - v->s.x[0];

            // if edge time within observation
            if( (v->s.x[0] <= obs_state[0]) && (etmp->to->s.x[0] > obs_state[0]) )
            {
                double post_obs = etmp->to->s.x[0] - obs_state[0];
                //lout<<"till_obs: "<<till_obs<<" post_obs: "<<post_obs<<endl;
                
                // change by particles
                etmp->prob = get_edge_prob_particles( v->s, till_obs, post_obs, etmp->to->s, yt );
                
                // change by covariance propagation
                //etmp->prob = get_edge_prob_kalman( v->s, till_obs, post_obs, etmp->to->s, yt );
            }
        } 

        kd_res_next(res);

        //write_next_prev(v);
    }
    kd_res_free(res);
}

double update_edges(vertex *from, edge *e)
{
    DepthTag dtag;
    
    //lout<<"update edges: "<<edge_num<<" "<<from->s.x[0]<<" "<<from->s.x[1]<<endl;

    state newst = system(from->s, e->delt, 1);
    double newvar[NUM_DIM-1] = { PRO_VAR[0]/2*(exp(2*e->delt) - 1), PRO_VAR[1]/2*(exp(2*e->delt) - 1) };
    double xtmp[NUM_DIM - 1] = { newst.x[1], newst.x[2]};
    double tocalci[NUM_DIM-1] = {e->to->s.x[1], e->to->s.x[2]};
    e->prob = normal_val(xtmp, newvar, tocalci, 2);

    return e->prob;
}

void add_major_sample(vertex *v)
{
    double *toput = v->s.x;
    kd_insert(state_tree, toput, v);
    rrg.add_vertex(v);
    rrg.num_vert++;
}

void draw_edges(vertex *v)
{
    DepthTag dtag;
    
    double pos[NUM_DIM] = {0};
    double *toput = v->s.x;
    kdres *res;
    res = kd_nearest_range(state_tree, toput, BOWLR );
    //lout<<"got "<<kd_res_size(res)<<" states"<<endl;

    while( !kd_res_end(res))
    {
        vertex *v1 = (vertex *)kd_res_item(res, pos); 
        if( v != v1)
        {
            double e1t = v1->s.x[0] - v->s.x[0];
            double e2t = v->s.x[0] - v1->s.x[0];

            //lout<<"e1t: "<<e1t<<" e2t: "<<e2t<<endl;
            // make edges, no probab, update later
            // write edges
            if( e1t > 0)
            {
                edge *e1 = new edge(v, v1, e1t);
                v->edgeout.push_back(e1);
                v1->edgein.push_back(e1);
                update_edges(v, e1);
                
                //lout<<"wrote e: "<<v->s.x[0]<<" "<<v->s.x[1]<<" to "<<v1->s.x[0]<<" "<<v1->s.x[1]<<endl;  
            }
            else if( e2t > 0)
            {
                edge *e2 = new edge(v1, v, e2t);
                v->edgein.push_back(e2);
                v1->edgeout.push_back(e2);
                update_edges(v1, e2);
                
                //if( v1 == rrg.vlist[0] )
                //    lout<< " writing out edge from (0)" << endl;
                //lout<<"wrote e: "<<v1->s.x[0]<<" "<<v1->s.x[1]<<" to "<<v->s.x[0]<<" "<<v->s.x[1]<<endl;
            }
        }
        kd_res_next(res);
    }
    kd_res_free(res);

    //write_next_prev( v );
}

void get_best_path()
{
    DepthTag dtag;
    
    vertex *vcurr = NULL;

    vcurr = rrg.vlist[0];
    best_path.push_back(vcurr->s);
    //lout<< "vfirst: " << vcurr->s.x[0] << " " << vcurr->s.x[1] << endl; 
    while( fabs(vcurr->s.x[0] - max_states[0]) > 1e-2 )
    {
        
        if( vcurr->edgeout.size() != 0)
        {
            double max_prob = 0;
            vertex * next_tmp;
            //lout<< " num edge: " << vcurr->edgeout.size() << endl;
            for(unsigned int i=0; i< vcurr->edgeout.size(); i++)
            {
                edge *etmp = vcurr->edgeout[i];
                if( etmp-> prob > max_prob)
                {
                    next_tmp = etmp->to;
                    max_prob = etmp->prob;
                }
            }
            //lout<< " maxprob: " << max_prob << endl;
            vcurr = next_tmp;
            best_path.push_back(vcurr->s);
        }
        else
            break;
        
        lout<< vcurr->s.x[0]<<" "<< vcurr->s.x[1]<<endl;
    }
}

void get_kalman_path()
{
    DepthTag dtag;
    
    double tmp[2] = {0};
    double mean[2] = {0};
    multivar_normal( mean, INIT_VAR, tmp, 2);
    
    state start_state;
    start_state.x[0] = 0;
    start_state.x[1] = 0.1 + tmp[0];
    start_state.x[2] = 0.1 + tmp[1];
    
    xhat[0].x[0] = start_state.x[0];
    xhat[0].x[1] = start_state.x[1];
    xhat[0].x[2] = start_state.x[2];
    
    // create kalman filter output
    for(int j=1; j < NUM_DIM; j++)
    {
        double Q = INIT_VAR[j-1];
        for(int i= 0; i<= max_states[0]/0.01; i++)
        {
            // update xhat
            xhat[i].x[0] = x[i].x[0];
            state curr_obs = obs(xhat[i], 1);
            double S = y[i].x[j] - curr_obs.x[j];

            double L = Q/( Q + OBS_VAR[j-1] );
            xhat[i].x[j] += L*S;

            // update covar
            Q = (1 - L)*Q;

            // propagate
            xhat[i+1].x[j] = exp(-0.01)*(xhat[i].x[j]);
            Q = exp(-0.02)*Q + PRO_VAR[j-1]/2*(exp(0.02) - 1);
        }
    }
}

void get_hmmf_path()
{
    DepthTag dtag;
    
    double tmp[2] = {0};
    double mean[2] = {0};
    multivar_normal( mean, INIT_VAR, tmp, 2);
    state start_state;
    start_state.x[0] = 0;
    start_state.x[1] = 0.1 + tmp[0];
    start_state.x[2] = 0.1 + tmp[1];

    double start_time = get_msec();
    vertex *vfirst = new vertex(start_state);
    add_major_sample(vfirst);
    vfirst->prob = 1.0; vfirst->prev = NULL;

    for(unsigned int j = 0; j < y.size(); j++)
    {
        for(int i=0; i< dM; i++)
        {
            vertex *v = new vertex( sample() );
            
            //v->s.x[0] = y[j].x[0];                  // make the times equal
            add_major_sample(v);

            draw_edges(v);
        }
        update_obs_prob(y[j]);
        if( j % 20 == 0)
            lout<<"i: "<<j<<endl;
    }

    get_best_path();
    lout<<"exec time: "<< get_msec() - start_time <<" [ms]"<< endl;
}

int main()
{

    g_depth = 0;
    g_logDepth = 1;

    cout.precision(4);
    srand(time(0));
    state_tree = kd_create(NUM_DIM);

    state x0; x0.x[0] = 0;
    x0.x[1] = 0.1;
    x0.x[2] = 0.1;
    
    x.push_back(x0);
    y.push_back(obs(x0, 0));
    for(int i=1; i<= max_states[0]/0.01; i++)
    {
        // create new state from old
        state newstate = system(x.back(), 0.01, 0);

        // push_back
        x.push_back(newstate);
        y.push_back(obs(newstate, 0));
    }
    
    // get path from hmm filter
    get_hmmf_path();

    // get kalman filter output
    get_kalman_path();
    
    plot_rrg();
    plot_traj();

    //kd_free(state_tree);
    
    return 0;
}

