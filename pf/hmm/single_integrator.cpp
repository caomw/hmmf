// first dim is time
#define NUM_DIM             (2)
#define PI                  (3.14156)

#define dM                  (50)
#define NUM_PARTICLES       (50)
#define GAMMA               ((max_states[1] - min_states[1]))
#define BOWLR               (GAMMA*pow( log( rrg.num_vert + 1)/ (rrg.num_vert + 1.0), 1.0/NUM_DIM) )
#define DT                  (0.01)
#define TIME_SCALE          (1.0)
#define EXTEND_DIST         (0.1)

#define randf               (rand()/(RAND_MAX + 1.0))

#define MIN_PROB            (0)
#define EDGE_PROB_THRESH    (1e-50)
#define NODE_PROB_THRESH    (1e-50)

#include "common.h"
#include "kdtree.h"

// limits of states
double max_states[4] = {0.5, 1.0, 1.0, 1.0};
double min_states[4] = {0.0, -0.2, -0.2, 0.0};

graph rrg;
vector<state> x, y, best_path;
vector<state> xkf( (max_states[0]/DT) + 1);
vector<vertex *> open_vertices;
double largest_edge = 0;

double INIT_VAR[3] = {1e-4, 1e-4, 1e-4};
double OBS_VAR[3] = {1e-2, 1e-2, 1e-2};
double PRO_VAR[3] = {1e-2, 1e-2, 1e-2};                       // noise of the continuous process


int is_free(state s)
{
    //return 1;

    bool retflag = 0;

    // obs 1
    if( (s.x[0] >= 0.127) && (s.x[0] <= 0.26) )
    {
        if( (s.x[1] >= 0) && (s.x[1] <= .217) )
            retflag = 0;
        else
            retflag = 1;
    }
    else
        retflag = 1;

    if (retflag == 0)
        return 0;

    // obs 2
    if( (s.x[0] >= 0.1) && (s.x[0] <= 0.2) )
    {
        if( (s.x[1] >= .32) && (s.x[1] <= .5) )
            retflag = 0;
        else
            retflag = 1;
    }
    else
        retflag = 1;

    return retflag;
}


state system(state s, double time, int is_clean)
{
    state t;

    double *var = new double[NUM_DIM-1];
    double *mean = new double[NUM_DIM-1];
    double *tmp = new double[NUM_DIM-1];

    for(int i=0; i<NUM_DIM; i++)
    {
        t.x[i] = s.x[i];
    }

    for(int i=0; i<NUM_DIM-1; i++)
    {
        tmp[i] = 0;
        mean[i] = 0;
        var[i] = PRO_VAR[i]/6*( exp(6*time) - 1);
    }
    if( !is_clean)  
        multivar_normal( mean, var, tmp, NUM_DIM-1);

    t.x[0] = t.x[0] + time;
    for(int i=1; i<NUM_DIM; i++)
        t.x[i] = exp(-3*time)*t.x[i] + tmp[i-1];

    delete[] var;
    delete[] mean;
    delete[] tmp;

    return t;
}

state observation(state s, int is_clean)
{
    state t;

    double *tmp = new double[NUM_DIM-1];
    double *mean = new double[NUM_DIM-1];

    if( !is_clean)  
        multivar_normal( mean, OBS_VAR, tmp, NUM_DIM-1);

    if(is_clean)
    {
        for(int i=0; i<NUM_DIM-1; i++)
            tmp[i] = 0;
    }

    t.x[0] = s.x[0];                        // time is same
    for(int i=1; i<NUM_DIM; i++)
        t.x[i] = s.x[i] + tmp[i-1];

    delete[] mean;
    delete[] tmp;

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
        for(int i=0; i<NUM_DIM; i++)
            rrgpout<<tstart->s.x[i]<<"\t";

        rrgpout<< tstart->prob << "\t";;
        rrgpout<<endl;
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

    traj<<"kf_path"<<endl;
    lout<<"kf_path"<<endl;
    for(int i=0; i< (int)xkf.size()-1; i++)
    {
        for(int j=0; j< NUM_DIM; j++)
        {
            traj << xkf[i].x[j]<<"\t";
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

bool is_normalize_edges(vertex *from)
{ 
    double totprob = 0;
    for(unsigned int i=0; i< from->edgeout.size(); i++)
    {
        totprob += from->edgeout[i]->prob;
    }
    if( fabs(totprob - 1) > 1e-6)
        return 0;
    else
        return 1;
}

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

double get_edge_prob_obs(state xinit, double till_obs, state obs)
{
    DepthTag dtag;

    state parts[NUM_PARTICLES];

    //lout<< "xinit: "<< xinit.x[0] << endl;

    for(int i=0; i< NUM_PARTICLES; i++)
    {
        parts[i] = (  system(xinit, till_obs, 0) );        // propagated particles directly
    }

    double *mean = new double[NUM_DIM-1];
    double *var = new double[NUM_DIM-1];
    double *tocalci = new double[NUM_DIM-1];
    double *tocalci_obs = new double[NUM_DIM-1];
    for(int j=1; j<NUM_DIM; j++)
    {
        tocalci_obs[j-1] = obs.x[j];
    }

    for(int i=0; i < NUM_PARTICLES; i++)
    {
        //lout<< "px: "<< parts[i].x[1] << " obs: " << obs.x[1] << endl;
        for(int j=1; j<NUM_DIM; j++)
        {
            mean[j-1] += parts[i].x[j];
        }
    }

    for(int i=0; i < NUM_PARTICLES; i++)
    {
        for(int j=1; j<NUM_DIM; j++)
        {
            var[j-1] += sq(parts[i].x[j] - mean[j-1]);
        } 
    } 
    for(int j=1; j<NUM_DIM; j++)
    {
        mean[j-1] = mean[j-1]/NUM_PARTICLES;
        var[j-1] = var[j-1] / (NUM_PARTICLES-1);
    }

    double retval = 0;
    for(int i=0; i < 1000; i++)
    {
        for(int j=1; j<NUM_DIM; j++)
        {
            tocalci[j-1] = randf*(6*var[j-1]) - 3*var[j-1] + mean[j-1];
        } 
        //retval += normal_val( mean, var, tocalci, NUM_DIM-1) * ( normal_val(tocalci, OBS_VAR, tocalci_obs, NUM_DIM-1) );
        retval += normal_val(tocalci, OBS_VAR, tocalci_obs, NUM_DIM-1);
    }

    double volume_parts = 0;
    if (NUM_DIM == 2)
        volume_parts = 2*3*sqrt(var[0]);
    else if( NUM_DIM == 3)
        volume_parts = 4*PI*9*sqrt((var[0])*(var[1]));
    else if( NUM_DIM == 4)
        volume_parts = 8*4*PI/3*27*sqrt((var[0])*(var[1])*(var[2]));

    retval = volume_parts*retval/(1000.0);
    //cout << "retval: " << retval << endl;

    delete[] mean;
    delete[] var;
    delete[] tocalci;
    delete[] tocalci_obs;

    return retval;
}

double get_edge_prob_system(state xinit, double till_end, state xend)
{
    DepthTag dtag;

    state parts[NUM_PARTICLES];

    //lout<< "xinit: "<< xinit.x[0] << endl;

    for(int i=0; i< NUM_PARTICLES; i++)
    {
        parts[i] = (  system(xinit, till_end, 0) );        // propagated particles directly
    }

    double *mean = new double[NUM_DIM-1];
    double *var = new double[NUM_DIM-1];
    double *tocalci = new double[NUM_DIM-1];

    //cout<< "post_obs: "<< post_obs<< " tot: " << totw << " " << totw2 << endl;

    for(int i=0; i<NUM_DIM-1; i++)
    {
        mean[i] = 0;
        var[i] = 0;
        tocalci[i] = 0;
    }

    for(int i=0; i< NUM_PARTICLES; i++)
    {
        for(int j=1; j< NUM_DIM; j++)
        {
            mean[j-1] += (parts[i].x[j]);
        }
    }
    for(int j=1; j< NUM_DIM; j++)
    {
        mean[j-1] = mean[j-1]/(NUM_PARTICLES);
        //cout<< "mean: "<< j << " " << mean[j-1] << endl;
    }
    for(int i=0; i< NUM_PARTICLES; i++)
    {
        for(int j=1; j< NUM_DIM; j++)
        {
            var[j-1] += sq (parts[i].x[j] - mean[j-1]);
        }
    }
    for(int j=1; j< NUM_DIM; j++)
    {
        var[j-1] = var[j-1]/(NUM_PARTICLES - 1);
        //cout<< "var: "<< j << " " << var[j-1] << endl;
    }

    //lout<<" mean: "<< mean <<" var: "<< var <<endl;

    for(int j=1; j< NUM_DIM; j++)
        tocalci[j-1] = xend.x[j];

    double target_w = normal_val( mean, var, tocalci, NUM_DIM-1);

    //if( post_obs < 0)
    //    cout<< "target_w: "<< target_w << endl;

    delete[] mean;
    delete[] var;
    delete[] tocalci;

    return target_w;
}

void write_next_prev(vertex *v)
{
    double max_prob = MIN_PROB;
    for(unsigned int i=0; i< v->edgeout.size(); i++)
    {
        edge *etmp = v->edgeout[i];
        if( etmp->prob > max_prob)
        {
            max_prob = etmp->prob;
            v->next = etmp->to;
        }
    }

    max_prob = MIN_PROB;
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
    double max_prob = MIN_PROB;
    for(unsigned int i=0; i< v->edgein.size(); i++)
    {
        edge *etmp = v->edgein[i];
        if( ( (etmp->prob) * (etmp->from->prob) ) > max_prob)
        {
            max_prob = ( (etmp->prob) * (etmp->from->prob) );
            v->prev = etmp->from;
            v->prob = max_prob;

            /*
               if(v->prob > 1)
               {
               cout<<"etmp->prob: " << etmp->prob << endl;
               cout<<"from->prob: "<< etmp->from->prob << endl;
               cout<<"rrg.num_vert: "<< rrg.num_vert << endl;
               getchar();
               }
               */
        }
    }
    /*
       if( (v->prob < NODE_PROB_THRESH) && (rrg.num_vert > 10*dM) )
       rrg.remove_vertex(v);
       */
}

// locate all vertices within bowlr and write prob of that vertex in weights
void update_obs_prob(state yt)
{
    DepthTag dtag;

    double *obs_state = yt.x;
    double obs_time_tolook = yt.x[0] - 2*BOWLR;
    double pos[NUM_DIM] = {0};

    kdres *res;
    res = kd_nearest_range(time_tree, &(obs_time_tolook), 2*BOWLR + 1e-3);
    while( !kd_res_end(res))
    {
        bool updated_edges = 0;
        vertex *v = (vertex *)kd_res_item(res, pos);

        for(unsigned int i=0; i < v->edgeout.size(); i++)
        {
            edge *etmp = v->edgeout[i];
            double till_obs = obs_state[0] - v->s.x[0];

            // if edge time within observation
            if( (v->s.x[0] <= obs_state[0]) && (etmp->to->s.x[0] > obs_state[0]) )
            {
                updated_edges = 1;

                //lout<<"till_obs: "<<till_obs<<" post_obs: "<<post_obs<<endl;

                // change by particles
                etmp->prob = (etmp->prob) * ( get_edge_prob_obs( v->s, till_obs, yt ) );

                if( etmp->prob < EDGE_PROB_THRESH)
                {
                    //cout<< "obs made edge prob low: deleting" << endl;
                    delete etmp;
                }
            }
        } 

        /*
           if(updated_edges)
           normalize_edges(v);
           */

        kd_res_next(res);
    }

    kd_res_free(res);
}

int is_edge_free( edge *etmp)
{
    //return 1;

    state init = etmp->from->s;
    state end = etmp->to->s;

    for(int i=0; i< 11; i++)
    {
        state stmp;
        for(int j=0; j< NUM_DIM; j++)
            stmp.x[j] = init.x[j] + (end.x[j] - init.x[j])/10*i;

        if( is_free( stmp) == 0)
            return 0;
    }
    return 1;
}

void draw_edges(vertex *v)
{
    DepthTag dtag;

    double pos[NUM_DIM] = {0};
    double tolook[NUM_DIM] = {0};
    for(int i=0; i< NUM_DIM; i++)
        tolook[i] = v->s.x[i];
    tolook[0] = tolook[0]/TIME_SCALE;

    kdres *res;
    res = kd_nearest_range(state_tree, tolook, BOWLR );
    //cout<<"got "<<kd_res_size(res)<<" states"<<endl;

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
                double prob = get_edge_prob_system(v->s, e1->delt, v1->s);
                e1->prob = prob;
                //cout<< "e->prob: " << prob << endl;       

                if( ( prob < EDGE_PROB_THRESH) || (is_edge_free(e1) == 0) )
                {
                    //cout<< "deleting edge" << endl;
                    delete e1;
                }
                else
                {
                    if( e1t > largest_edge)
                        largest_edge = e1t;
                }
                //cout<<"wrote e: "<<v->s.x[0]<<" "<<v->s.x[1]<<" to "<<v1->s.x[0]<<" "<<v1->s.x[1]<<endl;  
            }
            else if( e2t > 0)
            {
                edge *e2 = new edge(v1, v, e2t);
                v->edgein.push_back(e2);
                v1->edgeout.push_back(e2);
                double prob = get_edge_prob_system(v1->s, e2->delt, v->s);
                e2->prob = prob;
                //cout<< "e->prob: " << prob << endl;               

                if( ( prob < EDGE_PROB_THRESH) || (is_edge_free(e2) == 0) )
                {
                    //cout<< "deleting edge" << endl;
                    delete e2;
                }
                else
                {
                    if( e2t > largest_edge)
                        largest_edge = e2t;
                }

                //if( v1 == rrg.vlist[0] )
                //    lout<< " writing out edge from (0)" << endl;

                //normalize_edges(v1);
                //cout<<"wrote e: "<<v1->s.x[0]<<" "<<v1->s.x[1]<<" to "<<v->s.x[0]<<" "<<v->s.x[1]<<endl;
            }
        }
        kd_res_next(res);
    }
    kd_res_free(res);

    //normalize_edges(v);

    //cout<<"getchar: "; getchar();

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

    /*
    // concentrate around best trajectory
    if(best_path.size() != 0)
    {
    double stime = v->s.x[0];
    for(unsigned int i=0; i< best_path.size()-1; i++)
    {
    if( (best_path[i].x[0] >= stime) && ( best_path[i+1].x[0] <= stime) )
    {
    state start = best_path[i+1];
    state end = best_path[i];
    double d = dist(start, end);
    double extend = stime - start.x[0];

    //cout<<"concentrating: " << start.x[0] << " "<< end.x[0] << endl;

    for(int i=0; i< NUM_DIM; i++)
    v->s.x[i] = start.x[i] + (end.x[i] - start.x[i])*extend/d;

    state interpolate = v->s;
    bool sample_free = 0;
    int iter = 0;
    while ( sample_free == 0)
    {
    //cout<< iter++ << endl;

    double tmp[NUM_DIM-1] = {0};

    for(int i=1; i<NUM_DIM; i++)
    {
    tmp[i-1] = randf*0.5 - 0.25;
    v->s.x[i] = interpolate.x[i] + tmp[i-1];
    }

    if(is_free(v->s))
    sample_free = 1;
    }
    break;
    }
    }
    }
    */

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
        kd_insert(time_tree, &(toput[0]), v);

    }

    //cout<<"getchar: "; getchar();
}

void get_best_path()
{
    double pos[NUM_DIM] = {0};
    kdres *res;
    res = kd_nearest_range(state_tree, (x.back()).x, 0.5);
    double max_prob = MIN_PROB;
    vertex *vcurr = NULL;
    //cout<< "did kdtree query with: "<< 1.0 << " size: " << kd_res_size(res) << endl;

    while( !kd_res_end(res))
    {
        vertex *vtmp = (vertex *)kd_res_item(res, pos);
        if( fabs(vtmp->s.x[0] - max_states[0]) <= DT)
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
    {
        cout<<"Found last: "<< vcurr->s.x[0]<<" "<< vcurr->s.x[1] << " " << vcurr->prob<< endl;
    }
    else
    {
        cout<<"Couldn't find vcurr" << endl;
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
    //cout<< "wrote best_path" << endl;
}

void print_rrg()
{
    for(unsigned int i=0; i< rrg.vlist.size(); i++)
    {
        vertex *v = rrg.vlist[i];
        cout<<"node: " << i << " state: " << v->s.x[0] << " " << v->s.x[1] << " " << v-> prob << endl;

        cout << "ei: " << endl;
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

    double *mean = new double[NUM_DIM-1];
    double *tocalci = new double[NUM_DIM-1];

    for(int i=1; i<NUM_DIM; i++)
        mean[i-1] = start_state.x[i];

    double start_time = get_msec();

    // add vertices at the beginning
    for(int i=0; i< 100; i++)
    {
        vertex *v = new vertex( sample() );
        v->s.x[0] = y[0].x[0];

        for(int i=1; i<NUM_DIM; i++)
            tocalci[i-1] = v->s.x[i];

        v->prob = normal_val( mean, INIT_VAR, tocalci, NUM_DIM-1)* (max_states[1] - min_states[1])/100;
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
    }

    for(unsigned int j = 0; j < y.size(); j++)
    {
        update_obs_prob(y[j]);
    }

    lout<<"exec time: "<< (get_msec() - start_time)/1000.0 <<" [s]"<< endl;
}

void create_rrg_inc_obs_sampling()
{
    DepthTag dtag;
    vector<vertex *> last_added_vertices;
    vector<vertex *> curr_added_vertices;

    state start_state;
    start_state.x[0] = 0;
    for(int i=1; i<NUM_DIM; i++)
        start_state.x[i] = 0.5;

    double *mean = new double[NUM_DIM-1];
    double *tocalci = new double[NUM_DIM-1];

    for(int i=1; i<NUM_DIM; i++)
        mean[i-1] = start_state.x[i];

    double start_time = get_msec();

    // add vertices at the beginning
    for(int i=0; i< 100; i++)
    {
        vertex *v = new vertex( sample() );
        v->s.x[0] = y[0].x[0];

        for(int i=1; i<NUM_DIM; i++)
            tocalci[i-1] = v->s.x[i];

        v->prob = normal_val( mean, INIT_VAR, tocalci, NUM_DIM-1)* (max_states[1] - min_states[1])/100;
        //cout<< "vinit->prob: "<< v->s.x[1]<<" " << v->prob << endl;
        v->prev = NULL;

        add_sample(v);
        open_vertices.push_back( v );
        last_added_vertices.push_back(v);
    }
    update_obs_prob(y[0]);

    // start filter here
    for(unsigned int j = 1; j < y.size(); j++)
    {
        curr_added_vertices.clear();
        for(int i=0; i< dM; i++)
        {
            vertex *v = new vertex( sample() );
            v->s.x[0] = y[j].x[0];
            add_sample(v);

            curr_added_vertices.push_back(v);
        }
        update_obs_prob(y[j]);

        for(unsigned int i=0; i< last_added_vertices.size(); i++)
        {
            vertex *v = last_added_vertices[i];
            normalize_edges (v);
            //cout<<"j: "<< j << " prob was: "<< v->prob;
            write_viterbi(v);
            //cout<<" is: " << v->prob << endl;
        }
        //cout<< "size of lists: " << last_added_vertices.size() << " " << curr_added_vertices.size() << endl;

        last_added_vertices.clear();
        last_added_vertices = curr_added_vertices;

        //print_rrg();
        //cout<<"getchar: "; getchar();
    }

    // one last update
    for(unsigned int i=0; i< last_added_vertices.size(); i++)
    {
        vertex *v = last_added_vertices[i];
        normalize_edges (v);
        //cout<<"j: "<< j << " prob was: "<< v->prob;
        write_viterbi(v);
        //cout<<" is: " << v->prob << endl;
    }

    //print_rrg();
    get_best_path();
    lout<<"exec time: "<< (get_msec() - start_time)/1000.0 <<" [s]"<< endl;
}

void create_rrg_inc_uniform_sampling()
{
    DepthTag dtag;

    state start_state;
    start_state.x[0] = 0;
    for(int i=1; i<NUM_DIM; i++)
        start_state.x[i] = 0.5;

    double start_time = get_msec();

    double *mean = new double[NUM_DIM-1];
    double *tocalci = new double[NUM_DIM-1];

    for(int i=1; i<NUM_DIM; i++)
        mean[i-1] = start_state.x[i];

    // add vertices at the beginning
    for(int i=0; i< 100; i++)
    {
        vertex *v = new vertex( sample() );
        v->s.x[0] = y[0].x[0];

        for(int i=1; i<NUM_DIM; i++)
            tocalci[i-1] = v->s.x[i];

        v->prob = normal_val( mean, INIT_VAR, tocalci, NUM_DIM-1)* (max_states[1] - min_states[1])/100;
        //cout<< "vinit->prob: "<< v->s.x[1]<<" " << v->prob << endl;
        v->prev = NULL;
        add_sample(v);
    }
    update_obs_prob(y[0]);

    // start filter here
    for(unsigned int j = 0; j < y.size(); j++)
    {
        for(int i=0; i< dM; i++)
        {
            vertex *v = new vertex( sample() );
            add_sample(v);
        }
        update_obs_prob(y[j]);

        int num_updates = 0;
        while( num_updates < min(100*log(rrg.num_vert), (double)rrg.num_vert) )
        {
            int which = randf*(rrg.num_vert);
            vertex *v = rrg.vlist[which];

            normalize_edges(v);
            write_viterbi(v);

            num_updates++;
        }

        //cout<< "size of lists: " << last_added_vertices.size() << " " << curr_added_vertices.size() << endl;

        if( j%10 == 0)
        {
            get_best_path();
            cout<<"i: "<< j<< " bowl: "<< BOWLR << endl;
        }
        //print_rrg();
        //cout<<"getchar: "; getchar();
    }
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

        bool all_prev_closed = 1;
        for(int i=0; i< v->edgein.size(); i++)
        {
            if (v->edgein[i]->from->is_closed == 0)
                all_prev_closed = 0;
        }
        if(all_prev_closed)
        {
            write_viterbi( v);
            v->is_closed = 1;
        }
        for(int i=0; i< v->edgeout.size(); i++)
        {
            if( v->edgeout[i]->to != open_vertices.back() )
            {
                vector<vertex *>::iterator opv_iter;
                opv_iter = find( open_vertices.begin(), open_vertices.end(), v->edgeout[i]->to);

                if( (opv_iter == open_vertices.end()) && ( v->edgeout[i]->to->is_closed == 0) )
                    open_vertices.push_back( v->edgeout[i]->to );
            }
        }
        tot_is_open = open_vertices.size();
        iter++;

        if( iter % 5000 == 0)
            cout<<"iter: "<<iter<<" tot_open: "<< tot_is_open<<endl;

        //cout<<"getchar: "; getchar();
    }
    get_best_path();
}

void get_kalman_path()
{
    state start_state;
    start_state.x[0] = 0;
    xkf[0].x[0] = start_state.x[0];
    for(int i= 1; i< NUM_DIM; i++)
    {
        start_state.x[i] = 0.5;
        xkf[0].x[i] = start_state.x[i];
    }

    for(int dim=1; dim< NUM_DIM; dim++)
    {
        // create kalman filter output
        double Q = INIT_VAR[dim];
        for(int i= 0; i< max_states[0]/DT; i++)
        {
            // update xkf
            xkf[i].x[0] = x[i].x[0];
            state curr_obs = observation(xkf[i], 1);
            double S = y[i].x[dim] - curr_obs.x[dim];
            double L = Q/(Q + OBS_VAR[dim]);
            xkf[i].x[dim] += L*S;

            // update covar
            Q = (1 - L)*Q;

            // propagate
            xkf[i+1].x[dim] = exp(-3*DT)*(xkf[i].x[dim]);
            Q = exp(-6*DT)*Q + PRO_VAR[dim]/6*(exp(6*DT) - 1);
        }
    }
}

int main()
{
    g_depth = 0;
    g_logDepth = 1;

    srand( 0 );

    state_tree = kd_create(NUM_DIM);
    time_tree = kd_create(1);

    state x0; x0.x[0] = 0;
    for(int i=1; i<NUM_DIM; i++)
        x0.x[i] = 0.5;

    x.push_back(x0);
    y.push_back(observation(x0, 0));
    for(int i=1; i<= max_states[0]/DT; i++)
    {
        // create new state from old
        state newstate = system(x.back(), DT, 0);

        // push_back
        x.push_back(newstate);
        y.push_back(observation(newstate, 0));
    }

    double start_time = get_msec();

    // do batch
    // create_rrg();
    // normalize_rrg();

    // do incrementally
    // create_rrg_inc_obs_sampling();
    create_rrg_inc_uniform_sampling();

    lout<<"\nfinished exec, getting best path" << endl;
    lout<<"exec time: "<< (get_msec() - start_time)/1000.0 <<" [s]"<< endl;
    //lout<<"largest_edge: " << largest_edge << endl;    
    // kalman filter output
    get_kalman_path();

    plot_rrg();
    plot_traj();

    kd_free(state_tree);
    kd_free(time_tree);

    return 0;
}

