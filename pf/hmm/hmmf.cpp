#include "hmmf.h"


Vertex::Vertex(State& st)
{
    s = st;

    prev = NULL;
    next = NULL;
    prob_best_path = 1e-50;
    is_open_queue = 0;
    is_closed_queue = 0;

    edges_in.clear();
    edges_out.clear();
}

Edge::Edge(Vertex *f, Vertex *t, double prob){
    from = f;
    to = t;
    transition_prob = prob;
    transition_time = t->s[0] - f->s[0];
}

Graph::Graph(System& sys) {
    
    system = &sys;

    vlist.clear();
    num_vert = 0;
    num_particles = 25;
    samples_per_obs = 5*(NUM_DIM);
    num_observations = 0;

    state_tree = kd_create(NUM_DIM);
    time_tree = kd_create(1);
    
    gamma = 2.0;
};

Graph::~Graph()
{
    for(vector<Vertex*>::reverse_iterator i = vlist.rbegin(); i != vlist.rend(); i++)
    {
        delete *i;
        num_vert--;
    }
    for(list<Edge*>::reverse_iterator i = elist.rbegin(); i != elist.rend(); i++)
    {
        delete *i;
    }

    vlist.clear();
    truth.clear();
    obs.clear();
    kalman_path.clear();
    best_path.clear();

    kd_free(state_tree);
    kd_free(time_tree);
}

void Graph::remove_vertex(Vertex *v)
{
    // remove from vlist
    //vlist.remove(v); 
    //num_vert--;
    //cout<<"-------------------\nremoving edges: "<< v->edges_out.size() << " " << v->edges_in.size() << endl;
    
    /*
    cout<<"**removing incoming" << endl;
    // delete all edges from v and its neighbors
    for(list<Edge *>::reverse_iterator i = v->edges_in.rbegin(); i != v->edges_in.rend(); i++)
        remove_edge(*i);
    
    cout<<"**removing outgoing" << endl;
    for(list<Edge *>::reverse_iterator i = v->edges_out.rbegin(); i != v->edges_out.rend(); i++)
        remove_edge(*i);
    */

    // finally delete v
    delete v;
    //cout<<"deleted: num_vert - " << num_vert << endl;
}

void Graph::remove_edge(Edge *e)
{
    cout<<"--- inside remove_edge" << endl;
    // remove from, from + to lists and then call destructor
    e->from->edges_out.remove(e);
    cout<<"removed from, do to" << endl;
    if( e->to )
        cout<<"e->to exists " << e->to->edges_in.size() << endl;
    else
        cout<<" e->to does not exist: panic" << endl;
    e->to->edges_in.remove(e);
    cout<<"removed to, delete" << endl;

    delete e;
}

void Graph::plot_graph()
{
    ofstream rrgout("data/rrg.dat");
    ofstream rrgpout("data/rrgp.dat");
    //cout<<"writing rrg"<<endl;
    //cout<<"rrg size: "<< vlist.size() << endl;
    for(vector<Vertex*>::iterator i = vlist.begin(); i != vlist.end(); i++)
    {
        Vertex *tstart = (*i);
        for(list<Edge*>::iterator eo = tstart->edges_out.begin(); eo != tstart->edges_out.end(); eo++)
        {
            Vertex *tend = (*eo)->to;
            Edge *etmp = (*eo);

            //draw the edge
            rrgout<<tstart->s.x[0]<<"\t"<<tstart->s.x[1]<<"\t"<<tend->s.x[0]<<"\t"<<tend->s.x[1]<<"\t"<<etmp->transition_prob<<"\t"<<etmp->transition_time<<endl;
        }
        for(int i=0; i<NUM_DIM; i++)
            rrgpout<< tstart->s.x[i]<<"\t";

        rrgpout<< tstart->prob_best_path << "\t";
        rrgpout<<endl;
    }
    rrgout.close();
    rrgpout.close();
}

void Graph::plot_trajectory()
{
    ofstream traj("data/traj.dat");

    traj<<"system"<<endl;
    for(list<State>::iterator i= truth.begin(); i != truth.end(); i++)
    {
        State& curr = *i;
        for(int j=0; j< NUM_DIM; j++)
        {
            traj<< curr.x[j]<<"\t";
        }
        traj<<endl;
    }
    
    traj<<"observation"<<endl;
    for(list<State>::iterator i= obs.begin(); i != obs.end(); i++)
    {
        State& curr = *i;
        for(int j=0; j< NUM_DIM; j++)
        {
            traj<< curr.x[j]<<"\t";
        }
        traj<<endl;
    }

    traj<<"best_path"<<endl;
    for(list<State>::iterator i= best_path.begin(); i != best_path.end(); i++)
    {
        State& curr = *i;
        for(int j=0; j< NUM_DIM; j++)
        {
            traj<< curr.x[j]<<"\t";
        }
        traj<<endl;
    }
    
    traj<<"kf_path"<<endl;
    for(list<State>::iterator i= kalman_path.begin(); i != kalman_path.end(); i++)
    {
        State& curr = *i;
        for(int j=0; j< NUM_DIM; j++)
        {
            traj<< curr.x[j]<<"\t";
        }
        traj<<endl;
    }
    cout<<"trajectory plotting done" << endl;
    traj.close();
}

        
int Graph::insert_into_kdtree(Vertex *v)
{
    double key[NUM_DIM];
    system->get_key(v->s, key);
        
    kd_insert(state_tree, key, v);
    kd_insert(time_tree, &(key[0]), v);

    return 0;
}

Vertex* Graph::nearest_vertex(State s)
{
    assert(num_vert > 0);
    double key[NUM_DIM];
    system->get_key(s, key);

    kdres *res;
    res = kd_nearest(state_tree, key);
    if(kd_res_end(res))
    {
        cout<<"Error: no nearest"<<endl;
        exit(1);
    }
    Vertex *v = (Vertex*)kd_res_item_data(res);
    kd_res_free(res);

    return v;
}

void Graph::normalize_edges(Vertex *from)
{
    int nedges = from->edges_out.size();
    if( nedges == 0)
        return;
    else if ( nedges == 1)
    {
        Edge *etmp = from->edges_out.front();
        etmp->transition_prob = 1.0; 
    }
    else
    {
        double totprob = 0;
        for(list<Edge *>::iterator i = from->edges_out.begin(); i != from->edges_out.end(); i++)
        {
            totprob += (*i)->transition_prob;
        }
        
        if(totprob > 1e-100)
        {
            for(list<Edge *>::iterator i = from->edges_out.begin(); i != from->edges_out.end(); i++)
            {
                Edge *etmp = *i;
                etmp->transition_prob = etmp->transition_prob / totprob;
                //cout<<"wrote edge prob: "<< etmp->transition_prob << endl;
                if(etmp->transition_prob != etmp->transition_prob)
                {
                    cout<<"found a nan: "<< totprob << " nedges: "<< nedges << endl;
                    getchar();
                }
            }
        }
    }
    /*
    cout<<"getchar(): ";
    getchar();
    cout<<endl;
    */
}

int Graph::write_transition_prob(Edge *e)
{
    State parts[num_particles];

    for(int i=0; i< num_particles; i++)
    {
        parts[i] = system->integrate(e->from->s, e->transition_time, false);        // propagate noisy particles directly
    }

    double *mean = new double[NUM_DIM-1];
    double *var = new double[NUM_DIM-1];
    double *tocalci = new double[NUM_DIM-1];

    for(int i=0; i<NUM_DIM-1; i++)
    {
        mean[i] = 0;
        var[i] = 0;
        tocalci[i] = 0;
    }

    for(int i=0; i< num_particles; i++)
    {
        for(int j=1; j< NUM_DIM; j++)
        {
            mean[j-1] += (parts[i].x[j]);
        }
    }
    for(int j=1; j< NUM_DIM; j++)
    {
        mean[j-1] = mean[j-1]/(num_particles);
        //cout<< "mean: "<< j << " " << mean[j-1] << endl;
    }
    for(int i=0; i< num_particles; i++)
    {
        for(int j=1; j< NUM_DIM; j++)
        {
            var[j-1] += sq (parts[i].x[j] - mean[j-1]);
        }
    }
    for(int j=1; j< NUM_DIM; j++)
    {
        var[j-1] = var[j-1]/(num_particles - 1);
        //cout<< "var: "<< j << " " << var[j-1] << endl;
    }

    for(int j=1; j< NUM_DIM; j++)
        tocalci[j-1] = e->to->s.x[j];

    e->transition_prob = normal_val( mean, var, tocalci, NUM_DIM-1);
    
    //if( post_obs < 0)
    //    cout<< "target_w: "<< target_w << endl;

    delete[] mean;
    delete[] var;
    delete[] tocalci;

    return 0;
}

int Graph::write_observation_prob(Edge *e, State& obs)
{
    State parts[num_particles];

    for(int i=0; i< num_particles; i++)
    {
        parts[i] = system->integrate(e->from->s, (obs.x[0] - e->from->s.x[0]), false);        // propagate particles directly
    }

    double *mean = new double[NUM_DIM-1];
    double *var = new double[NUM_DIM-1];
    double *tocalci = new double[NUM_DIM-1];
    double *tocalci_obs = new double[NUM_DIM-1];
    for(int j=1; j<NUM_DIM; j++)
    {
        tocalci_obs[j-1] = obs.x[j];
    }

    for(int i=0; i < num_particles; i++)
    {
        //cout<< "px: "<< parts[i].x[1] << " obs: " << obs.x[1] << endl;
        for(int j=1; j<NUM_DIM; j++)
        {
            mean[j-1] += parts[i].x[j];
        }
    }

    for(int i=0; i < num_particles; i++)
    {
        for(int j=1; j<NUM_DIM; j++)
        {
            var[j-1] += sq(parts[i].x[j] - mean[j-1]);
        } 
    } 
    for(int j=1; j<NUM_DIM; j++)
    {
        mean[j-1] = mean[j-1]/(double)num_particles;
        var[j-1] = var[j-1] /(double)(num_particles-1);
    }
    
    // monte-carlo for getting observation probability
    double retval = 0;
    for(int i=0; i < 1000; i++)
    {
        for(int j=1; j<NUM_DIM; j++)
        {
            tocalci[j-1] = RANDF*(6*var[j-1]) - 3*var[j-1] + mean[j-1];
        } 
        //retval += normal_val( mean, var, tocalci, NUM_DIM-1) * ( normal_val(tocalci, OBS_VAR, tocalci_obs, NUM_DIM-1) );
        retval += normal_val(tocalci, system->obs_noise, tocalci_obs, NUM_DIM-1);
    }

    double volume_parts = 0;
    if (NUM_DIM == 2)
        volume_parts = 2*3*sqrt(var[0]);
    else if( NUM_DIM == 3)
        volume_parts = 4*M_PI*9*sqrt((var[0])*(var[1]));
    else if( NUM_DIM == 4)
        volume_parts = 8*4*M_PI/3*27*sqrt((var[0])*(var[1])*(var[2]));

    retval = volume_parts*retval/(1000.0);
    
    //cout<<"obs: prob was: " << e->transition_prob;
    e->transition_prob = e->transition_prob * retval;
    //cout<<" became: " << e->transition_prob << endl;

    delete[] mean;
    delete[] var;
    delete[] tocalci;
    delete[] tocalci_obs;

    return 0;
}


void Graph::update_viterbi( Vertex *v )
{
    double max_prob = -DBL_MAX;
    for(list<Edge*>::iterator i= v->edges_in.begin(); i != v->edges_in.end(); i++)
    {
        Edge *etmp = *i;;
        if( ( (etmp->transition_prob) * (etmp->from->prob_best_path) ) > max_prob)
        {
            max_prob = ( (etmp->transition_prob) * (etmp->from->prob_best_path) );
            v->prev = etmp->from;
            v->prob_best_path = max_prob;
            
            //assert(v->prob_best_path < 1.0);
            
            /*
            if(v->prob_best_path > 1)
            {
                cout<<"etmp->prob: " << etmp->transition_prob << endl;
                cout<<"from->prob: "<< etmp->from->prob_best_path << endl;
                cout<<"rrg.num_vert: "<< num_vert << endl;
                getchar();
            }
            */
        }
    }
}

void Graph::update_observation_prob(State& yt)
{
    /*
    double bowlr = gamma * pow( log(num_vert)/(num_vert), 1/(double)NUM_DIM);
    double *obs_state = yt.x;
    double obs_time_tolook = yt.x[0] - 2*bowlr;
    double pos[NUM_DIM] = {0};

    kdres *res;
    res = kd_nearest_range(time_tree, &(obs_time_tolook), 2*bowlr + 1e-3);
    //cout<<"yt[0]: " << yt.x[0]<<" inside bowl: "<< kd_res_size(res) << endl;
    
    while( !kd_res_end(res))
    {
        bool updated_edges = 0;
        Vertex *v = (Vertex *)kd_res_item(res, pos);

        for(list<Edge*>::iterator i= v->edges_out.begin(); i != v->edges_out.end(); i++)
        {
            Edge *etmp = *i;
            double till_obs = obs_state[0] - v->s.x[0];

            // if edge time within observation
            if( (v->s.x[0] <= obs_state[0]) && (etmp->to->s.x[0] > obs_state[0]) )
            {
                updated_edges = 1;

                //cout<<"till_obs: "<<till_obs<<" post_obs: "<<post_obs<<endl;

                // change by particles
                write_observation_prob(etmp, yt);
            }
        } 

        if(updated_edges)
            normalize_edges(v);

        kd_res_next(res);
    }

    kd_res_free(res);
    */

    double *obs_state = yt.x;
    for(vector<Vertex*>::iterator i=vlist.begin(); i!= vlist.end(); i++)
    {
        bool updated_edges = 0;
        Vertex *v = *i;

        for(list<Edge*>::iterator i= v->edges_out.begin(); i != v->edges_out.end(); i++)
        {
            Edge *etmp = *i;
            double till_obs = obs_state[0] - v->s.x[0];

            // if edge time within observation
            if( (v->s.x[0] <= obs_state[0]) && (etmp->to->s.x[0] > obs_state[0]) )
            {
                updated_edges = 1;

                //cout<<"till_obs: "<<till_obs<<" post_obs: "<<post_obs<<endl;

                // change by particles
                write_observation_prob(etmp, yt);
            }
        } 

        if(updated_edges)
            normalize_edges(v);
    }

}

bool Graph::is_edge_free( Edge *etmp)
{
    //return 1;
    
    int num_discrete = 10;
    State& init = etmp->from->s;
    State& end = etmp->to->s;

    for(int i=0; i< num_discrete+1; i++)
    {
        State stmp;
        for(int j=0; j< NUM_DIM; j++)
            stmp.x[j] = init.x[j] + (end.x[j] - init.x[j])/num_discrete*i;

        if( system->is_free(stmp) == 0)
            return false;
    }
    return true;
}

void Graph::add_sample()
{
    State stmp = system->sample();
    Vertex *v = new Vertex(stmp);

    /*
       if( rrg.vlist.size() != 0)
       {
       Vertex *near = nearest_vertex( v->s );

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

    if(v->s.is_free())
    sample_free = 1;
    }
    break;
    }
    }
    }
    */

    //cout<<"i: "<< rrg.vlist.size()-1 << " edges- in: " << v->edgein.size() << " out: "<< v->edgeout.size() << endl;

    vlist.push_back(v);
    num_vert++;
    insert_into_kdtree(v);

    connect_edges(v);
    update_viterbi(v);
}

int Graph::connect_edges(Vertex *v)
{
    double key[NUM_DIM] ={0};
    system->get_key(v->s, key);

    double bowlr = gamma * pow( log(num_vert)/(num_vert), 1/(double)NUM_DIM);
    
    kdres *res;
    res = kd_nearest_range(state_tree, key, bowlr );
    //cout<<"got "<<kd_res_size(res)<<" states in bowlr= "<< bowlr << endl;

    double pos[NUM_DIM] = {0};
    while( !kd_res_end(res))
    {
        Vertex *v1 = (Vertex *)kd_res_item(res, pos); 
        if( v != v1)
        {
            // cout<<"e1t: "<<e1t<<" e2t: "<<e2t<<endl;
            // make edges, no probab, update later
            // write edges
            if( (v1->s.x[0] - v->s.x[0]) > 0)
            {
                Edge *e1 = new Edge(v, v1);
                write_transition_prob(e1);
                
                if( is_edge_free(e1) )
                {
                    elist.push_back(e1);

                    v->edges_out.push_back(e1);
                    v1->edges_in.push_back(e1);
                }
                else
                    delete e1;
                //cout<<"wrote e: "<<v->s.x[0]<<" "<<v->s.x[1]<<" to "<<v1->s.x[0]<<" "<<v1->s.x[1]<<endl;  
            }
            else if( (v->s.x[0] - v1->s.x[0]) > 0)
            {
                Edge *e2 = new Edge(v1, v);
                write_transition_prob(e2);

                if(  is_edge_free(e2) )
                {
                    elist.push_back(e2);
                    
                    v->edges_in.push_back(e2);
                    v1->edges_out.push_back(e2);

                    normalize_edges(v1);
                }
                else
                    delete e2;
            }
        }
        kd_res_next(res);
    }
    kd_res_free(res);

    normalize_edges(v);
    
    //cout<<"getchar: "; getchar();
    return 0;
}


void Graph::get_best_path()
{
    double pos[NUM_DIM] = {0};
    double to_query[NUM_DIM];
    cout<<"truth: "<< truth.back().x[0] <<" " << truth.back().x[1] << endl;
    system->get_key(truth.back(), to_query);

    kdres *res;
    res = kd_nearest_range(time_tree, &(to_query[0]), 0.5);
    double max_prob = -DBL_MAX;
    Vertex *vcurr = NULL;
    //cout<< "did kdtree query with: "<< truth.back().x[1] << " size: " << kd_res_size(res) << endl;

    while( !kd_res_end(res))
    {
        Vertex *vtmp = (Vertex *)kd_res_item(res, pos);
        if( fabs(vtmp->s[0] - system->max_states[0]) <= system->sim_time_delta)
        {
            if( vtmp->prob_best_path > max_prob)
            {
                vcurr = vtmp;
                max_prob = vtmp->prob_best_path;
            }
            //cout<< vtmp->prob_best_path << endl;
        }
        kd_res_next(res);
    }
    kd_res_free(res);

    if( vcurr != NULL)
    {
        cout<<"Found last: "<< vcurr->s.x[0]<<" "<< vcurr->s.x[1] << " " << vcurr->prob_best_path<< endl;
    }
    else
    {
        cout<<"Couldn't find vcurr" << endl;
        return;
    }

    // populate best_path
    best_path.clear();
    while( vcurr != NULL)
    {
        best_path.push_back( vcurr->s );
        
        /*
        //cout<<" all edges are" << endl;
        for(unsigned int i=0; i< vcurr->edgein.size(); i++)
        {
            cout<<"edge i: "<< i << " : " << (vcurr->edgein[i]->prob) * (vcurr->edgein[i]->from->prob) << endl;
        }
        //cout<<"chose: "<< edge_index_connected_vertex(vcurr, vcurr->prev) << endl;
        //cout<<"getchar: "; getchar();
        */

        vcurr = vcurr->prev;

    }
    cout<< "wrote best_path" << endl;
}

void Graph::print_rrg()
{
    int counte = 0;
    int countv = 0;
    for(vector<Vertex*>::iterator i = vlist.begin(); i != vlist.end(); i++)
    {
        Vertex *v = *i;
        cout<<"node: " << countv++ << " state: " << v->s.x[0] << " " << v->s.x[1] << " " << v->prob_best_path << endl;

        counte = 0;
        cout << "ei: " << endl;
        for(list<Edge *>::iterator j = v->edges_in.begin(); j != v->edges_in.end(); j++)
        {
            cout<<"\t "<< counte++ << " " << (*j)->transition_prob << endl;
        }
        
        counte = 0;
        cout<<"eo: " << endl;
        for(list<Edge *>::iterator j = v->edges_out.begin(); j != v->edges_out.end(); j++)
        {
            cout<<"\t "<< counte++ << " " << (*j)->transition_prob << endl;
        }
    }
}

void Graph::propagate_system()
{
    num_observations = ( system->max_states[0] - system->min_states[0])/system->sim_time_delta;
    truth.push_back(system->init_state);

    for(int i =0; i< num_observations; i++)
    {
        State snext = system->integrate( truth.back(), system->sim_time_delta, false);
        truth.push_back( snext);
        State next_obs = system->observation( snext, false);
        obs.push_back(next_obs);
    }
}

void Graph::put_init_samples()
{
    double totprob = 0;
    for(int i=0; i < 100; i++)
    {
        State stmp = system->sample();
        Vertex *v = new Vertex(stmp);
        
        v->s.x[0] = system->min_states[0];
        v->prob_best_path = normal_val( &(system->init_state.x[1]), system->init_var, &(v->s.x[1]), NUM_DIM-1);      
        totprob += v->prob_best_path;

        vlist.push_back(v);
        num_vert++;
        insert_into_kdtree(v);
    }

    // normalize the init_samples
    for(vector<Vertex*>::iterator i = vlist.begin(); i != vlist.end(); i++)
    {
        Vertex* v = *i;
        v->prob_best_path = v->prob_best_path/totprob;
    }
}

bool Graph::is_everything_normalized()
{
    for(vector<Vertex*>::iterator i = vlist.begin(); i != vlist.end(); i++)
    {
        Vertex* v = *i;
        if(v->edges_out.size() > 0)
        {
            double totprob = 0;
            for(list<Edge *>::iterator j = v->edges_out.begin(); j != v->edges_out.end(); j++)
            {
                totprob += (*j)->transition_prob;
            }

            if( fabs(totprob - 1.0) > 0.1)
            {
                cout<<"offending vertex prob: "<< totprob << " nedges: "<< v->edges_out.size() << endl;
                return false;        
            }
        }
    }
    return true;
}

// returns 
int Graph::simulate_trajectory()
{
    list<State> seq;

    State stmp = system->init_state;
    stmp.x[0] = system->min_states[0];
    double init_prob = normal_val( &(system->init_state.x[1]), system->init_var, &(stmp.x[1]), NUM_DIM-1);
    
    bool can_go_forward = false;
    
    Vertex *vcurr = nearest_vertex(stmp);
    if(vcurr->edges_out.size() > 0)
        can_go_forward = true;

    seq.push_back(vcurr->s);
    double traj_prob = vcurr->prob_best_path;
    
    while(can_go_forward)
    {
        double rand_tmp = RANDF;
        double runner = 0;
        Edge *which_edge = NULL;
        for(list<Edge*>::iterator eo = vcurr->edges_out.begin(); eo != vcurr->edges_out.end(); eo++)
        {
            runner += (*eo)->transition_prob;
            if(runner > rand_tmp)
            {
                vcurr = (*eo)->to;
                which_edge = *eo;
                break;
            }
        }
        if( (vcurr->edges_out.size() > 0) && (which_edge == NULL))
        {
            cout<<"error- rand: "<< rand_tmp<<" runner: "<< runner <<" num_edges: "<< vcurr->edges_out.size() << endl;
            getchar();
        }

        if(vcurr->edges_out.size() == 0)
        {
            can_go_forward = false;
            break;
        }

        seq.push_back(vcurr->s);
        traj_prob *= which_edge->transition_prob;
    
    }

    if( fabs(vcurr->s.x[0] - system->max_states[0]) < system->sim_time_delta)
    {
        cout<<"traj_prob: "<<traj_prob << endl; 
        sanity_trajectories.push_back( seq);
        sanity_probabilities.push_back(traj_prob);
        return 0;
    }
    return 1;
}


