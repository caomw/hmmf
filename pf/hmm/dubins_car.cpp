// first dim is time, 3-d system
#define NUM_DIM             (3)
#define PI                  (3.14156)

#define NUM_PARTICLES       (100)
#define dM                  (50)
#define GAMMA               (1.0)
#define BOWLR               (GAMMA*pow( (max_states[0]/DT)*(log(rrg.num_vert) - log(max_states[0]/DT))/(float)(rrg.num_vert), 1.0/(NUM_DIM-1.0)))
#define EXTEND_DIST         (1.0)
#define DT                  (0.01)
#define TIME_SCALE          (1)

#define randf               (rand()/(RAND_MAX + 1.0))

#include "common.h"
#include "kdtree.h"
#include "figtree.h"

// limits of states
double max_states[NUM_DIM] = {1.0, 1.0, 1.0};
double min_states[NUM_DIM] = {0.0, 0, 0};

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
vector<state> xhat( (max_states[0]/DT) + 1);

double INIT_VAR[3] = {1e-4, 1e-4, 1e-4};
double OBS_VAR[3] = {1e-4, 1e-4, 1e-4};
double PRO_VAR[3] = {1e-4, 1e-4, 1e-4};                       // noise of the continuous process

state system(state s, double time, int is_clean)
{
    /*
    double dt = DT;
    double var[2] = { PRO_VAR[0]*dt, PRO_VAR[1]*dt};
    double tmp[2] = {0};
    double mean[2] = {0};
    
    state t;  t.x[0] = s.x[0]; t.x[1] = s.x[1]; t.x[2] = s.x[2]; t.x[3] = s.x[3];

    while(t.x[0] < (s.x[0] + time) )
    {
        multivar_normal(mean, var, tmp, 2);

        if(is_clean)
        {
            tmp[0] = 0;
            tmp[1] = 0;
        }
        double u1 = 1.0 + tmp[0];
        double u2 = 1.0 + tmp[1];

        double th_old = t.x[3];

        t.x[3] += u2*dt;
        if( t.x[3] > 2*PI)
            t.x[3] = t.x[3] - 2*PI;
        else if( t.x[3] < 0)
            t.x[3] = t.x[3] + 2*PI;

        t.x[2] += u1*sin(th_old)*dt;
        t.x[1] += u1*cos(th_old)*dt;

        t.x[0] += dt;

    }
    */   
    
    // single integrator
    state t;  t.x[0] = s.x[0]; t.x[1] = s.x[1]; t.x[2] = s.x[2];
    double var[2] = { PRO_VAR[0]/2*( exp(2*time) - 1), PRO_VAR[1]/2*( exp(2*time) - 1)};
    double tmp[2] = {0};
    double mean[2] = {0};
    multivar_normal(mean, var, tmp, 2);
   
    if ( is_clean)
    {
        tmp[0] = 0; tmp[1] = 0; tmp[2] = 0;
    }

    t.x[0] = t.x[0] + time;
    t.x[1] = exp(-1*time)*t.x[1] + tmp[0];
    t.x[2] = exp(-1*time)*t.x[2] + tmp[1];
     
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
    t.x[0] = s.x[0];                        // time is same
    t.x[1] = s.x[1] + tmp[0];               // add noise to obs
    t.x[2] = s.x[2] + tmp[1];             
    
    return t;
}

void plot_rrg()
{
    DepthTag dtag;
    
    ofstream rrgout("rrg.dat");
    ofstream rrgpout("rrgp.dat");
    lout<<"write rrg"<<endl;
    cout<<"rrg size: "<< rrg.vlist.size() << endl;
    for(vector<vertex*>::iterator i = rrg.vlist.begin(); i != rrg.vlist.end(); i++)
    {
        vertex *tstart = (*i);
        for(vector<edge*>::iterator eo = tstart->edgeout.begin(); eo != tstart->edgeout.end(); eo++)
        {
            vertex *tend = (*eo)->to;
            edge *etmp = (*eo);

            //draw the edge
            rrgout<<tstart->s.x[0]<<"\t"<<tstart->s.x[1]<<"\t"<<tend->s.x[0]<<"\t"<<tend->s.x[1]<<"\t"<<etmp->prob<<"\t"<<etmp->delt<<endl;
        }
        rrgpout<<tstart->s.x[0]<<"\t"<<tstart->s.x[1]<<"\t"<<tstart->s.x[2]<<endl;
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
    for(int i=0; i< (int)x.size(); i++)
    {
        for(int j=0; j< NUM_DIM; j++)
        {
            traj<<x[i].x[j]<<"\t";
        }
        traj<<endl;
    }
    traj<<"observation"<<endl;
    lout<<"observation"<<endl;
    for(int i=0; i< (int)y.size(); i++)
    {
        for(int j=0; j< NUM_DIM; j++)
        {
            traj<<y[i].x[j]<<"\t";
        }
        traj<<endl;
    }

    traj<<"best_path"<<endl;
    lout<<"best_path"<<endl;
    for(int i=0; i< (int)best_path.size(); i++)
    {
        for(int j=0; j< NUM_DIM; j++)
        {
            traj << best_path[i].x[j]<<"\t";
        }
        traj<<endl;
    }
    
    lout<<"plotting traj done" << endl;
    traj.close();
}

vertex* nearest_vertex(state s)
{
    kdres *resn;
    resn = kd_nearest(state_tree, s.x);
    if(kd_res_end(resn))
    {
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

double get_edge_prob_particles(state xinit, double till_obs, double post_obs, state xend, state obs)
{
    DepthTag dtag;
    
    state parts[NUM_PARTICLES];
    double weights[NUM_PARTICLES] = {0};

    double totw = 0, totw2 = 0;
    //lout<< "xinit: "<< xinit.x[0] << endl;

    for(int i=0; i< NUM_PARTICLES; i++)
    {
        parts[i] = (  system(xinit, till_obs, 0) );        // propagated particles directly
        weights[i] = 1.0/(double)NUM_PARTICLES ; 
    }
    for(int i=0; i < NUM_PARTICLES; i++)
    {
        if ( post_obs >= 0)
        {
            //lout<< "px: "<< parts[i].x[1] << " obs: " << obs.x[1] << endl;
            double mean[NUM_DIM - 1] = { obs.x[1], obs.x[2] };
            double tocalci[NUM_DIM - 1] = { parts[i].x[1], parts[i].x[2] };

            weights[i] = normal_val( mean, OBS_VAR, tocalci, NUM_DIM-1 ); 
        }
        //cout<< "w: " << weights[i] << endl;

        totw += weights[i];
        totw2 += sq(weights[i]);
    }
    //cout<< "post_obs: "<< post_obs<< " tot: " << totw << " " << totw2 << endl;

    double mean[NUM_DIM-1] = {0};
    double var[NUM_DIM-1] = {0};

    for(int i=0; i< NUM_PARTICLES; i++)
    {
        if ( post_obs >= 0)
            parts[i] = system( parts[i], post_obs, 0 );            // get new states

        for(int j=1; j< NUM_DIM; j++)
        {
            mean[j-1] += (parts[i].x[j]) * weights[i];
        }
    }
    for(int j=1; j< NUM_DIM; j++)
    {
        mean[j-1] = mean[j-1]/totw;
        //cout<< "mean: "<< j << " " << mean[j-1] << endl;
    }
    for(int i=0; i< NUM_PARTICLES; i++)
    {
        for(int j=1; j< NUM_DIM; j++)
        {
            var[j-1] += weights[i]*sq (parts[i].x[j] - mean[j]);
        }
    }
    for(int j=1; j< NUM_DIM; j++)
    {
        var[j-1] = totw*var[j-1]/(totw*totw - totw2);
        //cout<< "var: "<< j << " " << var[j-1] << endl;
    }
    
    //lout<<" mean: "<< mean <<" var: "<< var <<endl;
    
    double tocalci[NUM_DIM-1] = { xend.x[1], xend.x[2] } ;
    double target_w = normal_val( mean, var, tocalci, NUM_DIM-1);
   

    //if( post_obs < 0)
    //    cout<< "target_w: "<< target_w << endl;
    
    return target_w;
}

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
                
                if( etmp->prob < 1e-6)
                {
                    //cout<< "obs made edge prob low: deleting" << endl;
                    delete etmp;
                }
            }
        } 

        kd_res_next(res);

        //write_next_prev(v);
    }
    kd_res_free(res);
}

void draw_edges(vertex *v)
{
    double edge_prob_threshold = 1e-3;
    DepthTag dtag;
    
    double pos[NUM_DIM] = {0};
    double tolook[NUM_DIM] = {0};
    for(int i=0; i< NUM_DIM; i++)
        tolook[i] = v->s.x[i];
    tolook[0] = tolook[0]/TIME_SCALE;

    kdres *res;
    res = kd_nearest_range(state_tree, tolook, BOWLR );
    //cout<<"got "<<kd_res_size(res)<<" states"<<endl;
    state dummy;

    while( !kd_res_end(res))
    {
        vertex *v1 = (vertex *)kd_res_item(res, pos); 
        if( v != v1)
        {
            double e1t = v1->s.x[0] - v->s.x[0];
            double e2t = v->s.x[0] - v1->s.x[0];

            // cout<<"e1t: "<<e1t<<" e2t: "<<e2t<<endl;
            // make edges, no probab, update later
            // write edges
            if( e1t > 0)
            {
                edge *e1 = new edge(v, v1, e1t);
                v->edgeout.push_back(e1);
                v1->edgein.push_back(e1);
                double prob = get_edge_prob_particles(v->s, e1->delt, -1, v1->s, dummy);
                e1->prob = prob;
                //cout<< "e->prob: " << prob << endl;       
                
                if( prob < edge_prob_threshold)
                {
                    //cout<< "deleting edge" << endl;
                    delete e1;
                }
                //cout<<"wrote e: "<<v->s.x[0]<<" "<<v->s.x[1]<<" to "<<v1->s.x[0]<<" "<<v1->s.x[1]<<endl;  
            }
            else if( e2t > 0)
            {
                edge *e2 = new edge(v1, v, e2t);
                v->edgein.push_back(e2);
                v1->edgeout.push_back(e2);
                double prob = get_edge_prob_particles(v1->s, e2->delt, -1, v->s, dummy);
                e2->prob = prob;
                //cout<< "e->prob: " << prob << endl;               
                
                if( prob < edge_prob_threshold)
                {
                    //cout<< "deleting edge" << endl;
                    delete e2;
                }
                //if( v1 == rrg.vlist[0] )
                //    lout<< " writing out edge from (0)" << endl;
                //cout<<"wrote e: "<<v1->s.x[0]<<" "<<v1->s.x[1]<<" to "<<v->s.x[0]<<" "<<v->s.x[1]<<endl;
            }
        }
        kd_res_next(res);
    }
    kd_res_free(res);
    
    //cout<<"getchar: "; getchar();
    
    //write_next_prev( v );
}

void add_sample(vertex *v)
{
    /*
    if( rrg.vlist.size() != 0)
    {
        vertex *near = nearest_vertex( v->s );

        double d = dist( v->s, near->s);
        if ( d > EXTEND_DIST)
        {
            for(int i=0; i< NUM_DIM; i++)
                v->s.x[i] = near->s.x[i] + (v->s.x[i] - near->s.x[i])*EXTEND_DIST/d;
        }
    }
    */
    draw_edges (v);
    //cout<<"i: "<< rrg.vlist.size()-1 << " edges- in: " << v->edgein.size() << " out: "<< v->edgeout.size() << endl;

    double toput[NUM_DIM] = {0};
    for(int i=0; i< NUM_DIM; i++)
        toput[i] = v->s.x[i];
    toput[0] = toput[0]/TIME_SCALE;

    if ( 0 )
    {
        if( (v->edgeout.size() == 0) || (v->edgein.size() == 0) )
        {
            //cout<< " didnot add to graph " <<endl;
        }
        else
        {
            rrg.add_vertex(v);
            rrg.num_vert++;
            kd_insert(state_tree, toput, v);
        }
    }
    else
    {
        rrg.add_vertex(v);
        rrg.num_vert++;
        kd_insert(state_tree, toput, v);
    }

    //cout<<"getchar: "; getchar();
}

void get_best_path()
{
    DepthTag dtag;
    
    vertex *vcurr = rrg.vlist[0];
    best_path.push_back(vcurr->s);
    cout<< "vfirst: " << vcurr->s.x[0] << " " << vcurr->s.x[1] << endl; 
    while( fabs(vcurr->s.x[0] - max_states[0]) > 1e-2 )
    {
        if( vcurr->edgeout.size() != 0)
        {
            double max_prob = 0;
            vertex * next_tmp;
            //cout<< "num edge: " << vcurr->edgeout.size() << endl;
            for(unsigned int i=0; i< vcurr->edgeout.size(); i++)
            {
                edge *etmp = vcurr->edgeout[i];
                //cout<<"etmp_prob: "<< etmp->prob << endl;
                if( etmp-> prob > max_prob)
                {
                    next_tmp = etmp->to;
                    max_prob = etmp->prob;
                }
            }
            //cout<< "maxprob: " << max_prob << endl;
            vcurr = next_tmp;
            best_path.push_back(vcurr->s);
            //cout<<"getchar: "; getchar(); cout<< endl;
        }
        else
        {
            cout << "Couldn't find any outgoing edge, quitting" << endl;
            break;
        }
        lout<< vcurr->s.x[0]<<" "<< vcurr->s.x[1]<<endl;
    }
    
    /*
    vertex *vcurr = NULL;
    double pos[NUM_DIM] = {0};
    kdres *res;
    res = kd_nearest_range(state_tree, vcurr->s, BOWLR);
    while( !kd_res_end(res))
    {
        vertex *vtmp = (vertex *)kd_res_item(res, pos);
        if( fabs(vtmp->s.x[0] - min_states[0]) < 1e-2)
        {

        }

        kd_res_nex();
    }
    kd_res_free(res);
    */
}

void print_rrg()
{
    for(unsigned int i=0; i< rrg.vlist.size(); i++)
    {
        vertex *v = rrg.vlist[i];
        cout<<"node: " << i << endl << "ei: " << endl;
        for(int j=0; j < v->edgein.size(); j++)
        {
            cout<<"\t "<< j << " " << v->edgein[j]->prob << endl;
        }
        cout<<"eo: " << endl;
        for(int j=0; j < v->edgeout.size(); j++)
        {
            cout<<"\t "<< j << " " << v->edgeout[j]->prob << endl;
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
    for(int i=1; i<NUM_DIM; i++)
        start_state.x[i] = 1.0 + tmp[i-1];

    double start_time = get_msec();
    vertex *vfirst = new vertex(start_state);
    add_sample(vfirst);
    vfirst->prob = 1.0; vfirst->prev = NULL;

    for(unsigned int j = 0; j < y.size(); j++)
    {
        for(int i=0; i< dM; i++)
        {
            vertex *v = new vertex( sample() );
            
            //v->s.x[0] = y[j].x[0];                  // make the times equal
            add_sample(v);
        }
        update_obs_prob(y[j]);

        //print_rrg();
        //cout<<"getchar: "; getchar();
        
        if( j % 10 == 0)
            lout<<"i: "<<j<<" bowl: "<< BOWLR << endl;
    }
    /* 
    for(int i=0; i< rrg.vlist.size(); i++)
    {
        vertex *v = rrg.vlist[i];
        cout<<"v: "<< i <<" ei: "<< v->edgein.size() << " eo: " << v->edgeout.size() << endl;
    }
    */
    lout<<"\nfinished exec, getting best path" << endl;
    get_best_path();
    lout<<"exec time: "<< (get_msec() - start_time)/1000.0 <<" [s]"<< endl;
}

int main()
{

    g_depth = 0;
    g_logDepth = 1;

    srand(time(0));
    state_tree = kd_create(NUM_DIM);

    state x0; x0.x[0] = 0;
    for(int i=1; i<NUM_DIM; i++)
        x0.x[i] = 1.0;
    
    x.push_back(x0);
    y.push_back(obs(x0, 0));
    for(int i=1; i<= max_states[0]/DT; i++)
    {
        // create new state from old
        state newstate = system(x.back(), DT, 0);

        // push_back
        x.push_back(newstate);
        y.push_back(obs(newstate, 0));
    }
    
    // get path from hmm filter
    get_hmmf_path();

    plot_rrg();
    plot_traj();

    kd_free(state_tree);
    
    return 0;
}

