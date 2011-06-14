// first dim is time
#define NUM_DIM             (2)
#define PI                  (3.14156)

#define NUM_PARTICLES       (100)
#define dM                  (100)
#define GAMMA               (0.5)
//#define BOWLR               (GAMMA*pow( log( rrg.num_vert)/ rrg.num_vert, 1.0/NUM_DIM) )
#define BOWLR               (0.02)
#define DT                  (0.01)
#define EXTEND_DIST         (max_states[0])
#define TIME_SCALE          (1.0)
#define BOWLR_OBS           (0.5)

#define randf               (rand()/(RAND_MAX + 1.0))

#define MIN_PROB_LOG        (0)

#include "common.h"
#include "kdtree.h"
#include "figtree.h"

// limits of states
double max_states[NUM_DIM] = {0.2, 0.6};
double min_states[NUM_DIM] = {0.0, 0.2};

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
vector<vertex *> open_vertices;


double INIT_VAR[3] = {1e-4, 1e-4, 1e-4};
double OBS_VAR[3] = {1e-4, 1e-4, 1e-4};
double PRO_VAR[3] = {1e-4, 1e-4, 1e-4};                       // noise of the continuous process

state system(state s, double time, int is_clean)
{
    state t;  t.x[0] = s.x[0]; t.x[1] = s.x[1];
    double var = PRO_VAR[0]/2*( exp(2*time) - 1);
    double tmp = 0;
    double mean = 0;
    multivar_normal( &mean, &var, &tmp, 1);
   
    if ( is_clean)
    {
        tmp = 0;
    }

    t.x[0] = t.x[0] + time;
    t.x[1] = exp(-1*time)*t.x[1] + tmp;
     
    return t;
}

state obs(state s, int is_clean)
{
    state t;
    double tmp = 0;;
    double mean = 0;
    multivar_normal( &mean, OBS_VAR, &tmp, 1);
    
    if(is_clean)
    {
        tmp = 0;
    }
    t.x[0] = s.x[0];                        // time is same
    t.x[1] = s.x[1] + tmp;

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
        rrgpout<<tstart->s.x[0]<<"\t"<<tstart->s.x[1]<<endl;
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

int edge_index_connected_vertex(vertex *end, vertex *tocheck)
{
    int index = -1;
    for(int i=0; i< (int)end->edgein.size(); i++)
    {
        if( (end->edgein[i])->from == tocheck)
            index = i;
    }
    return index;
}

// writes the best and normalizes
void normalize_edges(vertex *from)
{
    DepthTag dtag;
    double totprob = 0;
    for(unsigned int i=0; i< from->edgeout.size(); i++)
    {
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
            double mean = obs.x[1];
            double tocalci = parts[i].x[1];

            weights[i] = normal_val( &mean, OBS_VAR, &tocalci, NUM_DIM-1 ); 
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
    
    double tocalci[NUM_DIM-1] = { xend.x[1] } ;
    double target_w = normal_val( mean, var, tocalci, NUM_DIM-1) * (PI*BOWLR*BOWLR)/log(rrg.num_vert);
   

    //if( post_obs < 0)
    //    cout<< "target_w: "<< target_w << endl;
    
    return target_w;
}

void write_next_prev(vertex *v)
{
    double max_prob = MIN_PROB_LOG;
    for(unsigned int i=0; i< v->edgeout.size(); i++)
    {
        edge *etmp = v->edgeout[i];
        if( etmp->prob > max_prob)
        {
            max_prob = etmp->prob;
            v->next = etmp->to;
        }
    }

    max_prob = MIN_PROB_LOG;
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

void write_viterbi( vertex *v )
{
    double max_prob = MIN_PROB_LOG;
    for(unsigned int i=0; i< v->edgein.size(); i++)
    {
        edge *etmp = v->edgein[i];
        if( ( (etmp->prob) * (etmp->from->prob) ) > max_prob)
        {
            max_prob = ( (etmp->prob) * (etmp->from->prob) );
            v->prev = etmp->from;
            v->prob = max_prob;
        }
    }
}

// locate all vertices within bowlr and write prob of that vertex in weights
void update_obs_prob(state yt)
{
    DepthTag dtag;
    double edge_prob_threshold = -1e6;

    double *obs_state = yt.x;
    double pos[NUM_DIM] = {0};

    kdres *res;
    res = kd_nearest_range(state_tree, obs_state, BOWLR_OBS);
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
                etmp->prob = (etmp->prob) + ( get_edge_prob_particles( v->s, till_obs, post_obs, etmp->to->s, yt ) );
                
                if( etmp->prob < edge_prob_threshold)
                {
                    //cout<< "obs made edge prob low: deleting" << endl;
                    delete etmp;
                }
            }
        } 
        

        kd_res_next(res);

        //write_next_prev(v);
    }
    
    /*
    kd_res_rewind( res);
    while( !kd_res_end(res))
    {
        vertex *v = (vertex *)kd_res_item(res, pos);
        write_viterbi(v);
        kd_res_next(res);
    }
    */
    kd_res_free(res);
}

void draw_edges(vertex *v)
{
    double edge_prob_threshold = -1e6;
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
    
    //write_viterbi(v);
    //cout<<"getchar: "; getchar();
    
    //write_next_prev( v );
}

void add_sample(vertex *v)
{
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
    
    draw_edges (v);
    //cout<<"i: "<< rrg.vlist.size()-1 << " edges- in: " << v->edgein.size() << " out: "<< v->edgeout.size() << endl;

    double toput[NUM_DIM] = {0};
    for(int i=0; i< NUM_DIM; i++)
        toput[i] = v->s.x[i];
    toput[0] = toput[0]/TIME_SCALE;
   
    
    bool add_vertex = 1;
    /*
    if( rrg.num_vert > 2*dM)
    {
        if ( v->edgein.size() != 0)
            add_vertex = 1;
    }
    else
    {
        add_vertex = 1;
    }
    */
    if ( add_vertex)
    {
        rrg.add_vertex(v);
        rrg.num_vert++;
        kd_insert(state_tree, toput, v);
    }

    //cout<<"getchar: "; getchar();
}

void get_best_path()
{
    double pos[NUM_DIM] = {0};
    kdres *res;
    res = kd_nearest_range(state_tree, (x.back()).x, BOWLR_OBS);
    double max_prob = MIN_PROB_LOG;
    vertex *vcurr = NULL;
    cout<< "did kdtree query with: "<< 1.0 << " size: " << kd_res_size(res) << endl;

    while( !kd_res_end(res))
    {
        vertex *vtmp = (vertex *)kd_res_item(res, pos);
        if( fabs(vtmp->s.x[0] - max_states[0]) < DT)
        {
            if( vtmp->prob > max_prob)
            {
                vcurr = vtmp;
                max_prob = vtmp->prob;
            }
            //cout<< vtmp->prob << endl;
        }

        kd_res_next(res);
    }
    kd_res_free(res);
    
    if( vcurr != NULL)
        cout<<"Found last: "<< vcurr->s.x[0]<<" "<< vcurr->s.x[1] << " " << vcurr->prob<< endl;
    else
    {
        cout<<"Couldn't find vcurr, quitting" << endl;
        return;
    }
    best_path.clear();
    while( vcurr != NULL)
    {
        best_path.push_back( vcurr->s );
        
        /*
        cout<<" all edges are" << endl;
        for(unsigned int i=0; i< vcurr->edgein.size(); i++)
        {
            cout<<"edge i: "<< i << " : " << (vcurr->edgein[i]->prob) * (vcurr->edgein[i]->from->prob) << endl;
        }
        cout<<"chose: "<< edge_index_connected_vertex(vcurr, vcurr->prev) << endl;
        cout<<"getchar: "; getchar();
        */

        vcurr = vcurr->prev;
         
    }
    cout<< "wrote best_path" << endl;
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


void create_rrg()
{
    DepthTag dtag;
    
    state start_state;
    start_state.x[0] = 0;
    for(int i=1; i<NUM_DIM; i++)
        start_state.x[i] = 0.5;

    double start_time = get_msec();
    double mean = start_state.x[1];
    
    // add vertices at the beginning
    for(int i=0; i< dM; i++)
    {
        vertex *v = new vertex( sample() );
        v->s.x[0] = y[0].x[0];
        double tocalci = v->s.x[1];
        v->prob = normal_val( &mean, INIT_VAR, &tocalci, NUM_DIM-1);
        //cout<< "vinit->prob: "<< v->s.x[1]<<" " << v->prob << endl;
        v->prev = NULL;
        
        add_sample(v);
        open_vertices.push_back( v );
    }
    // start filter here
    for(unsigned int j = 0; j < y.size(); j++)
    {
        for(int i=0; i< dM; i++)
        {
            vertex *v = new vertex( sample() );
            //v->s.x[0] = y[j].x[0];
            add_sample(v);
        }
        //print_rrg();
        //cout<<"getchar: "; getchar();
        
        if( j % 5 == 0)
            lout<<"i: "<<j<<" bowl: "<< BOWLR << endl;
    }
    
    for(unsigned int j = 0; j < y.size(); j++)
    {
        update_obs_prob(y[j]);
    }
    
    lout<<"exec time: "<< (get_msec() - start_time)/1000.0 <<" [s]"<< endl;
}

void normalize_rrg()
{
    cout<<"Normalizing" << endl;
    for(unsigned int i=0; i< rrg.vlist.size(); i++)
    {
        vertex *v = rrg.vlist[i];
        normalize_edges( v );
    }
    
    cout<<"Doing viterbi: "<< rrg.num_vert<<" vertices"<<endl; 
    
    // do Dijkstra's
    unsigned int tot_is_open = open_vertices.size();
    int iter = 0;
    cout<<"iter: "<<iter<<" tot_open: "<< tot_is_open<<endl;
    while( tot_is_open != 0)
    {
        vertex *v = open_vertices.front();
        open_vertices.erase( open_vertices.begin() );
        write_viterbi( v);

        for(int i=0; i< v->edgeout.size(); i++)
        {
            if( v->edgeout[i]->to != open_vertices.back() )
            {
                vector<vertex *>::iterator opv_iter;
                opv_iter = find( open_vertices.begin(), open_vertices.end(), v->edgeout[i]->to);
                
                if(opv_iter == open_vertices.end())
                    open_vertices.push_back( v->edgeout[i]->to );
            }
        }
        tot_is_open = open_vertices.size();
        iter++;
        /* 
        if( iter %100 == 0)
            cout<<"iter: "<<iter<<" tot_open: "<< tot_is_open<<endl;
        //cout<<"getchar: "; getchar();
        */
    }
}


int main()
{
    g_depth = 0;
    g_logDepth = 1;

    srand( 0 );

    state_tree = kd_create(NUM_DIM);

    state x0; x0.x[0] = 0;
    for(int i=1; i<NUM_DIM; i++)
        x0.x[i] = 0.5;
    
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
    create_rrg();

    normalize_rrg();

    lout<<"\nfinished exec, getting best path" << endl;
    get_best_path();
    
    plot_rrg();
    plot_traj();

    kd_free(state_tree);
    
    return 0;
}

