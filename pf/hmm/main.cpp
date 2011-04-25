// first dim is time, 1-d system
#define NUM_DIM         (2)

#define TMAX            (0.3)
#define TMIN            (0.0)

#define XMAX            (0.20)
#define XMIN            (0.00)

#define randf           (rand()/(RAND_MAX + 1.0))

#define EXTEND_DIST     (0.01)
#define dM              (100)
#define sinit           (1e-3)
#define sobs            (1e-5)
#define obs_max_uniform (sqrt(3*sobs))
#define spro            (1e-3)
#define BETA            (100)
#define GAMMA           (1.0)
#define BOWLGAMMA       (GAMMA*pow(log(rrg.num_vert)/(float)(rrg.num_vert), 1.0/(NUM_DIM)))
#define BOWLR           ( (BOWLGAMMA >= 0.01) ? BOWLGAMMA : 0.01)
//#define BOWLR           (0.1)
#define PI              (3.14156)

#include "common.h"
#include "kdtree.h"

// halton
int seed[NUM_DIM] = {0, 0};
int base[NUM_DIM] = {2, 3};
int step = 0;

kdtree *state_tree, *mini_tree;
graph rrg;
vector<state> x, y, best_path, xhat(TMAX/0.01 + 1);
vector<state> simx;

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
    return 1/sqrt(2*PI*var)*temp;
}

state system(state s, double time){

    state t;
    t.x[0] = s.x[0] + time;                 // add time
    t.x[1] = exp(-1*time)*s.x[1];           // xt+1 = xt + N(0, spro) (of continuous process)
    return t;
}

state obs(state s, int is_clean){
    state t;
    double tmp1=0, tmp2=0;
    randn(0, sobs, tmp1, tmp2);
    //tmp1 = randf*2*obs_max_uniform - obs_max_uniform;
    if(is_clean)
        tmp1 = 0;

    t.x[0] = s.x[0];                    // time is same
    t.x[1] = s.x[1] + tmp1;             // add noise to obs
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

    traj<<"observation"<<endl;
    cout<<"observation"<<endl;
    for(int i=0; i< (int)y.size()-1; i++)
        traj<<y[i].x[0]<<"\t"<<y[i].x[1]<<endl;

    traj<<"best_path"<<endl;
    cout<<"best_path"<<endl;
    for(unsigned int i=0; i< best_path.size(); i++)
        traj<< best_path[i].x[0]<<"\t"<< best_path[i].x[1]<<endl;
    
    traj<<"kf_path"<<endl;
    cout<<"kf_path"<<endl;
    for(unsigned int i=0; i< xhat.size(); i++)
        traj<< xhat[i].x[0]<<"\t"<< xhat[i].x[1]<<endl;

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
   */

// writes the best and normalizes
void normalize_edges(vertex *from)
{
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

    if(besteout != NULL)
        from->best_out = besteout;
    //else
    //    cout<<"best out is NULL"<<endl;

}

// locate all vertices within bowlr and write prob of that vertex in weights
void update_obs_prob(state yt)
{
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
            if( (v->s.x[0] <= obs_state[0]) && (etmp->to->s.x[0] >= obs_state[0]) )
            {
                double post_obs = etmp->to->s.x[0] - obs_state[0];
                //cout<<"till_obs: "<<till_obs<<" post_obs: "<<post_obs<<endl;

                // propagate qdot = 2aq + bwb to till_obs
                double newvar = spro/2*(exp(2*till_obs) - 1);
                //cout<<"old_var: "<<newvar<<endl;

                // reduce q by q bt rinv b q due to obs
                // get new xhat
                state xhat =  system(v->s, till_obs);
                double L = newvar/(newvar + sobs);
                state obs_predict = obs(xhat, 1);
                xhat.x[1] = xhat.x[1] + L*(obs_state[1] - obs_predict.x[1]);
                newvar = (1 - L)*newvar;
                //cout<<"after obs: "<<newvar<<endl;

                // propagate post_obs using xhat and q
                // final gaussian is centered at final_state with var = newvar
                newvar = exp(-2*post_obs)*newvar + spro/2*(exp(2*post_obs) -1);
                state final_state = system(xhat, post_obs);
                //cout<<"after growing: "<<newvar<<" xhat: "<<final_state.x[1]<<" actual: "<<etmp->to->s.x[1]<<endl;

                // change probab of etmp then
                etmp->prob = normal_val(xhat.x[1], newvar, etmp->to->s.x[1]);

                // rewire etmp->to
                vertex *vto = etmp->to;
                double max_prob = 0;
                for(unsigned int j=0; j< vto->edgein.size(); j++)
                {
                    edge *e_vto = vto->edgein[j];
                    if( (e_vto->prob)*(e_vto->from->prob) > max_prob)
                    {
                        max_prob = (e_vto->prob)*(e_vto->from->prob);
                        vto->prev = e_vto->from;
                    }
                }

            }
        } 

        // rewire v
        double max_prob = 0;
        for(unsigned int j=0; j< v->edgeout.size(); j++)
        {
            edge *e_vout = v->edgeout[j];
            if( (e_vout->prob) > max_prob)
            {
                max_prob = (e_vout->prob);
                v->next = e_vout->to;
            }
        }
        
        kd_res_next(res);
    }
    kd_res_free(res);
}

double update_edges(vertex *from, edge *e)
{
    //cout<<"update edges: "<<edge_num<<" "<<from->s.x[0]<<" "<<from->s.x[1]<<endl;

    state newst = system(from->s, e->delt);
    double newvar = spro/2*(exp(2*e->delt) - 1);
    double xtmp = newst.x[1];
    e->prob = normal_val(xtmp, newvar, e->to->s.x[1]);

    return e->prob;
}

void add_major_sample(vertex *v)
{
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
}

void draw_edges(vertex *v)
{
    double pos[NUM_DIM] = {0};
    double *toput = v->s.x;
    kdres *res;
    res = kd_nearest_range(state_tree, toput, BOWLR);
    //cout<<"got "<<kd_res_size(res)<<" states"<<endl;

    while( !kd_res_end(res))
    {
        vertex *v1 = (vertex *)kd_res_item(res, pos); 
        if( v != v1)
        {
            double e1t = v1->s.x[0] - v->s.x[0];
            double e2t = v->s.x[0] - v1->s.x[0];

            //cout<<"e1t: "<<e1t<<" e2t: "<<e2t<<endl;
            // make edges, no probab, update later
            // write edges
            if( e1t > 0)
            {
                edge *e1 = new edge(v, v1, e1t);
                v->edgeout.push_back(e1);
                v1->edgein.push_back(e1);
                double e1prob = update_edges(v, e1);

                // outgoing edge from v, write next
                if( v->next != NULL)
                {
                    if( e1prob > (v->best_out->prob))
                    {
                        v->next = v1;
                        v->best_out = e1;
                    }
                }
                else
                {
                    if( v1 == rrg.vlist[0])
                        cout<<"writing out edge to first vert here"<<endl;
                    v->next = v1;
                    v->best_out = e1;
                }

                //cout<<"wrote e: "<<v->s.x[0]<<" "<<v->s.x[1]<<" to "<<v1->s.x[0]<<" "<<v1->s.x[1]<<endl;  
            }
            else if( e2t > 0)
            {
                edge *e2 = new edge(v1, v, e2t);
                v->edgein.push_back(e2);
                v1->edgeout.push_back(e2);
                double e2prob = update_edges(v1, e2);

                // incoming edge into v, write prev
                if( v->prev != NULL)
                {
                    if( (v->prev->prob)*(v->best_in->prob) < (v1->prob)*(e2prob))
                    {
                        v->prev = v1;
                        v->best_in = e2;
                        v->prob = (v1->prob)*(e2prob);
                    }
                }
                else
                {
                    v->prev = v1;
                    v->best_in = e2;
                    v->prob = (v1->prob)*(e2prob);
                }

                //cout<<"wrote e: "<<v1->s.x[0]<<" "<<v1->s.x[1]<<" to "<<v->s.x[0]<<" "<<v->s.x[1]<<endl;
            }

        }
        kd_res_next(res);
    }
    kd_res_free(res);
}

void get_best_path(int from_first)
{
    if(from_first)
    {
        state start_state = rrg.vlist[0]->s;
        vertex *vcurr = NULL;
        // find best node near zero
        kdres *res;
        res = kd_nearest_range(state_tree, start_state.x, (XMAX - XMIN)/2);
        double pos[2];
        double max_prob = 0;

        while( !kd_res_end(res) )
        {
            vertex *vtmp = (vertex *)kd_res_item(res, pos);
            if( vtmp->s.x[0] < 0.01)
            {
                if(vtmp->prob > max_prob)
                {
                    max_prob = vtmp->prob;
                    vcurr = vtmp;
                }
            }
            kd_res_next(res);
        }
        kd_res_free(res);
        cout<<"got first as: "<< vcurr->s.x[0]<<" "<< vcurr->s.x[1]<<endl;

        vcurr = rrg.vlist[0];
        while( (vcurr->s.x[0] < TMAX) && (vcurr->next != NULL) )
        {
            cout<< vcurr->s.x[0]<<" "<< vcurr->s.x[1]<<endl;
            best_path.push_back(vcurr->s);
            vcurr = vcurr->next;
        }
    }
    else
    {
        vertex *vlast = NULL;
        // find best node at TMAX
        state last; last.x[0] = TMAX; last.x[1] = (XMAX + XMIN)/2;
        kdres *res;
        res = kd_nearest_range(state_tree, last.x, (XMAX - XMIN)/2);
        double pos[2];
        double max_prob = 0;

        while( !kd_res_end(res) )
        {
            vertex *vtmp = (vertex *)kd_res_item(res, pos);
            double tmp = fabs(vtmp->s.x[0] - TMAX); 
            if( tmp < 0.1)
            {
                if(vtmp->prob > max_prob)
                {
                    max_prob = vtmp->prob;
                    vlast = vtmp;
                }
            }
            kd_res_next(res);
        }
        kd_res_free(res);
        cout<<"got last as: "<< vlast->s.x[0]<<" "<< vlast->s.x[1]<<endl;

        vertex *vcurr = vlast;
        while( (vcurr->s.x[0] != 0) && (vcurr->prev != NULL) )
        {
            cout<< vcurr->s.x[0]<<" "<< vcurr->s.x[1]<<endl;
            best_path.push_back(vcurr->s);
            vcurr = vcurr->prev;

        }

    }
}

int main()
{
    cout.precision(4);
    srand(time(0));
    state_tree = kd_create(NUM_DIM);
    mini_tree = kd_create(NUM_DIM);

    state x0; x0.x[0] = 0; x0.x[1] = 0.1;
    double tmp1, tmp2;
    randn(0, sinit, tmp1, tmp2);
    state start_state;
    start_state.x[0] = 0; start_state.x[1] = x0.x[1] + sqrt(sinit)*tmp1;
    
    x.push_back(x0);
    y.push_back(obs(x0, 0));
    for(int i=1; i<= TMAX/0.01; i++)
    {
        // create new state from old
        double tmp1=0, tmp2=0;
        double var = spro/2*(exp(2*0.01) - 1);
        randn(0, var, tmp1, tmp2);
        state newstate = system(x.back(), 0.01);
        newstate.x[1] += tmp1;

        // push_back
        x.push_back(newstate);
        y.push_back(obs(newstate, 0));
    }
    
    /*
    xhat[0].x[1] = start_state.x[1];
    xhat[0].x[0] = start_state.x[0];

    // create kalman filter output
    double Q = sinit;
    for(int i= 0; i<= TMAX/0.01; i++)
    {
        // update xhat
        xhat[i].x[0] = x[i].x[0];
        state curr_obs = obs(xhat[i], 1);
        double S = y[0].x[1] - curr_obs.x[1];
        double L = Q/(Q + sobs);
        xhat[i].x[1] += L*S;
        
        // update covar
        Q = (1 - L)*Q;

        // propagate
        xhat[i+1].x[1] = exp(-0.01)*xhat[i].x[1];
        Q = exp(-0.02)*Q + spro/2*(exp(0.02) - 1);
    }
    */

    double start_time = get_msec();

    vertex *vfirst = new vertex(x0);
    add_major_sample(vfirst);
    vfirst->prob = 1.0; vfirst->prev = NULL;

    for(unsigned int j = 0; j < y.size(); j++)
    {
        for(int i=0; i< 400; i++)
        {
            vertex *v = new vertex(sample());
            //v->s.x[0] = y[j].x[0];
            add_major_sample(v);

            draw_edges(v);
        }
        update_obs_prob(y[j]);
        cout<<"i: "<<j<<endl;
    }
    
    get_best_path(1);

    cout<<"exec time: "<< get_msec() - start_time <<endl;

    plot_rrg();
    plot_traj();

    kd_free(state_tree);
    kd_free(mini_tree); 
    return 0;
}

