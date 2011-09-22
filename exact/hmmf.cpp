#include "hmmf.h"


Vertex::Vertex(State& st)
{
    s = st;

    prev = NULL;
    next = NULL;
    prob_best_path = -1;
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
    num_particles = 50;
    obs_interval = 5;

    state_tree = kd_create(NUM_DIM);
    time_tree = kd_create(1);

    goal_width = 10*system->sim_time_delta;
    gamma = 2.2;

    curr_best_cost = DBL_MAX;
    curr_best_vertex = NULL;
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
    vector<Vertex*>::iterator which = find(vlist.begin(), vlist.end(), v);
    if( which != vlist.end())
    {
        vlist.erase(which); 
        num_vert--;
    }
    //cout<<"-------------------\nremoving edges: "<< v->edges_out.size() << " " << v->edges_in.size() << endl;
    
    //cout<<"**removing incoming" << endl;
    // delete all edges from v and its neighbors
    for(list<Edge *>::reverse_iterator i = v->edges_in.rbegin(); i != v->edges_in.rend(); i++)
        remove_edge(*i);
    
    //cout<<"**removing outgoing" << endl;
    for(list<Edge *>::reverse_iterator i = v->edges_out.rbegin(); i != v->edges_out.rend(); i++)
        remove_edge(*i);

    // finally delete v
    delete v;
    //cout<<"deleted: num_vert - " << num_vert << endl;
}

void Graph::remove_edge(Edge *e)
{
    elist.remove(e);

    //cout<<"--- inside remove_edge" << endl;
    // remove from, from + to lists and then call destructor
    e->from->edges_out.remove(e);
    /*
    cout<<"removed from, do to" << endl;
    if( e->to )
        cout<<"e->to exists " << e->to->edges_in.size() << endl;
    else
        cout<<" e->to does not exist: panic" << endl;
    */
    e->to->edges_in.remove(e);
    //cout<<"removed to, delete" << endl;

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
            //rrgout<<tstart->s.x[0]<<"\t"<<tstart->s.x[1]<<"\t"<<tend->s.x[0]<<"\t"<<tend->s.x[1]<<"\t"<<etmp->transition_prob<<"\t"<<etmp->transition_time<<endl;
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
            Edge *etmp = *i;
            if(etmp->transition_prob != etmp->transition_prob)
            {
                cout<<"found a nan: "<< totprob << " nedges: "<< nedges << endl;
                getchar();
            }
            totprob += (*i)->transition_prob;
        }
        
        if(totprob > 1.0/DBL_MAX)
        {
            for(list<Edge *>::iterator i = from->edges_out.begin(); i != from->edges_out.end(); i++)
            {
                Edge *etmp = *i;
                etmp->transition_prob = etmp->transition_prob / totprob;
                //cout<<"wrote edge prob: "<< etmp->transition_prob << endl;
            }
        }
        else
        {
            cout<<"totprob is: "<< totprob << " [DITCH]" << endl;
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

// importance resampling
int pfilter_resample(State* parts, double *weights, int num_particles)
{
    double totweight = 0;
    for(int i=0; i<num_particles; i++)
        totweight += weights[i];

    double* cum = new double[num_particles];
    double curr_tot = 0;
    for(int i=0; i<num_particles; i++)
    {
        weights[i] = weights[i]/totweight;
        curr_tot += weights[i];
        cum[i] = curr_tot;

        //reset to equal weights
        weights[i] = 1.0/num_particles;
    }

    State* newparts = new State[num_particles];
    double u = RANDF/(double)num_particles;
    int i = 0;
    for(int j=0; j<num_particles; j++)
    {
        double tocheck = u + j/(double)num_particles;
        while(tocheck > cum[i])
            i++;

        newparts[j] = parts[i];
    }
    for(int j=0; j<num_particles; j++)
        parts[j] = newparts[j]; 


    delete[] cum;
    delete[] newparts;
    
    // weights are already set above
    return 0;
}


// particle filter step incorporating all observations within e->transition_time
int Graph::write_transition_prob_with_obs(Edge *e)
{
    vector<State> obs_states;
    int obs_size = 0;

    double stime = e->from->s.x[0];
    double etime = e->to->s.x[0];
    for(list<State>::iterator i = obs.begin(); i != obs.end(); i++)
    {
        State& curr_state = *i;
        if( (curr_state.x[0] < etime) && (curr_state.x[0] > stime) )
        {
            obs_states.push_back(curr_state);
            obs_size++;
        }
    }
    
    /*
    cout<<"obs: ";
    for(int i=0;i < obs_size; i++)
    {
        cout<< obs_states[i].x[0]<<" ";
    }
    cout<<endl;
    */

    State parts[num_particles];
    double weights[num_particles];
    for(int i=0; i< num_particles; i++)
    {
        parts[i] = e->from->s;
        weights[i] = 1/(double)num_particles;
    }

    //cout<<"stime: "<< stime<<" transition_time: "<< e->transition_time << " etime: "<< etime << endl;
    
    for(int obs_iter=0; obs_iter<obs_size; obs_iter++)
    {
        double delta_t = obs_states[obs_iter].x[0] - parts[0].x[0];           // for how much time to propagate
        //cout<<"parts: "<< parts[0].x[0]<<" delta_t: "<< delta_t << " obs: "<< obs_states[obs_iter].x[0]<<endl;

        for(int i=0; i< num_particles; i++)
        {
            parts[i] = system->integrate(parts[i], delta_t, false);
        }

        for(int i=0; i< num_particles; i++)
        {
            weights[i] = weights[i] * normal_val( &(obs_states[obs_iter].x[1]), system->obs_noise, &(parts[i].x[1]), NUM_DIM-1);
        }

        pfilter_resample(parts, weights, num_particles);

        /*
        for(int i=0; i< num_particles; i++)
        {
            cout<<"weight: "<< weights[i]<< " state: ";
            for(int j=0; j<NUM_DIM; j++)
            {
                cout<<parts[i].x[j]<<" ";
            }
            cout<<endl;
        }
        cout<<"press key: "; getchar(); cout<<endl;
        */
    }
    
    //cout<<"press key: "; getchar(); cout<<endl;
    double delta_t = e->to->s.x[0] - parts[0].x[0];
    if (delta_t > 0)
    {
        for(int i=0; i< num_particles; i++)
        {
            parts[i] = system->integrate(parts[i], delta_t, false);
        }
    }
    else
        cout<<"pfilter delta_t was < 0"<<endl;

    double *mean = new double[NUM_DIM-1];
    double *var = new double[NUM_DIM-1];
    double *tocalci = new double[NUM_DIM-1];
    double totw = 0, totw2 = 0;

    for(int i=0; i<NUM_DIM-1; i++)
    {
        mean[i] = 0;
        var[i] = 0;
        tocalci[i] = 0;
    }

    /*
    for(int i=0; i< num_particles; i++)
    {
        cout<<"weight: "<< weights[i]<< " state: ";
        for(int j=0; j<NUM_DIM; j++)
        {
            cout<<parts[i].x[j]<<" ";
        }
        cout<<endl;
    }
    //cout<<"press key: "; getchar(); cout<<endl;
    */

    for(int i=0; i< num_particles; i++)
    {
        //cout<<"weights[i]: "<< weights[i] << endl;
        for(int j=1; j< NUM_DIM; j++)
        {
            mean[j-1] = mean[j-1] + (parts[i].x[j]) * weights[i];
        }
        totw += weights[i];
        totw2 += sq(weights[i]);
    }

    for(int i=0; i< num_particles; i++)
    {
        for(int j=1; j< NUM_DIM; j++)
        {
            var[j-1] = var[j-1] + sq(parts[i].x[j] - mean[j-1])*weights[i];
        }
    }

    //cout<<"totw: " << totw <<" totw2: "<< totw2 << endl;
    for(int j=1; j< NUM_DIM; j++)
    {
        mean[j-1] = mean[j-1]/totw;
        var[j-1] = totw*var[j-1]/(totw*totw - totw2);
        
        if(var[j-1] < 1e-100)
            var[j-1] = 1e-100;

        //cout<< "mean: "<< j << " " << mean[j-1] << endl;
        //cout<< "var: "<< j << " " << var[j-1] << endl;
    }

    for(int j=1; j< NUM_DIM; j++)
        tocalci[j-1] = e->to->s.x[j];

    e->transition_prob = normal_val( mean, var, tocalci, NUM_DIM-1);
    //cout<<"e->tprob: "<< e->transition_prob << endl;

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

void Graph::propagate_viterbi(Vertex* v)
{
    int max_size = -1;
    list<Vertex*> myq;
    int myq_size = 0;
    
    myq.push_back(v);
    myq_size++;

    for(list<Edge*>::iterator i = v->edges_out.begin(); i != v->edges_out.end(); i++)
    {
        myq.push_back((*i)->to);
        myq_size++;
    }

    while(myq_size)
    {
        Vertex* vtmp = myq.front();
        update_viterbi(vtmp);
        
        myq.pop_front();
        myq_size--;

        for(list<Edge*>::iterator i = vtmp->edges_out.begin(); i != vtmp->edges_out.end(); i++)
        {
            Edge* etmp = *i;
            if( (vtmp->prob_best_path)*(etmp->transition_prob) > (etmp->to->prob_best_path) )
            {
                myq.push_back(etmp->to);
                myq_size++;
            }
        }
        if( myq_size > max_size)
            max_size = myq_size;
    }
    //cout<<"max_size: "<< max_size << endl;
}

void Graph::update_viterbi( Vertex *v )
{
    double max_prob = -1;
    for(list<Edge*>::iterator i= v->edges_in.begin(); i != v->edges_in.end(); i++)
    {
        Edge *etmp = *i;
        double tocalci = (etmp->transition_prob)*(etmp->from->prob_best_path);
        if( tocalci > max_prob )
        {
            max_prob = tocalci;
            v->prev = etmp->from;
            v->prob_best_path = tocalci;
        }
    }
}

void Graph::update_goal_viterbi()
{
    State last_state;
    last_state.x[0] = system->max_states[0];

    double to_query[NUM_DIM];
    double pos[NUM_DIM] = {0};
    system->get_key( last_state, to_query);

    kdres *res;
    res = kd_nearest_range(time_tree, &(to_query[0]), goal_width);

    while(!kd_res_end(res))
    {
        Vertex *v = (Vertex *)kd_res_item(res, pos);
        
        update_viterbi(v);

        kd_res_next(res);
    }
    kd_res_free(res);
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

    double extend_dist = 0.1;
    if( vlist.size() != 0)
    {
        Vertex *near = nearest_vertex( v->s );

        double d = dist( v->s, near->s);
        if ( d > extend_dist)
        {
            for(int i=0; i< NUM_DIM; i++)
                v->s.x[i] = near->s.x[i] + (v->s.x[i] - near->s.x[i])*extend_dist/d;
        }
    }
    
    connect_edges(v);
    
    vlist.push_back(v);
    num_vert++;
    insert_into_kdtree(v);
}

int Graph::connect_edges(Vertex *v)
{
    double key[NUM_DIM] ={0};
    system->get_key(v->s, key);

    double bowlr = gamma * pow( log(num_vert)/(num_vert), 1/(double)NUM_DIM);

    kdres *res;
    res = kd_nearest_range(state_tree, key, bowlr );
    //cout<<"got "<<kd_res_size(res)<<" states in bowlr= "<< bowlr << endl;

    bool added_outgoing = false;
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
                write_transition_prob_with_obs(e1);
                //cout<<"e1 prob: " << e1->transition_prob << endl;

                if( is_edge_free(e1))
                {
                    added_outgoing = true;
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
                write_transition_prob_with_obs(e2);
                //cout<<"e2 prob: " << e2->transition_prob << endl;

                if(  is_edge_free(e2) )
                {
                    elist.push_back(e2);

                    v->edges_in.push_back(e2);
                    v1->edges_out.push_back(e2);

                    normalize_edges(v1);
                    propagate_viterbi(v1);
                }
                else
                    delete e2;
            }
        }
        kd_res_next(res);
    }
    kd_res_free(res);

    normalize_edges(v);

    if(added_outgoing)
        propagate_viterbi(v);

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
    res = kd_nearest_range(time_tree, &(to_query[0]), goal_width);

    double max_prob = -1;
    Vertex *vcurr = NULL;
    cout<< "did kdtree query with: "<< truth.back().x[0] << " size: " << kd_res_size(res) << endl;

    while( !kd_res_end(res))
    {
        Vertex *vtmp = (Vertex *)kd_res_item(res, pos);
        //cout<< vtmp->prob_best_path<<" ";

        if( vtmp->prob_best_path > max_prob)
        {
            vcurr = vtmp;
            max_prob = vtmp->prob_best_path;
        }
        kd_res_next(res);
    }
    kd_res_free(res);

    //Vertex *vcurr = curr_best_vertex;
    if( vcurr != NULL)
    {
        cout<<"found last: ";
        for(int i=0; i< NUM_DIM; i++)
            cout<< vcurr->s.x[i]<<" ";
        cout<<"prob: "<< vcurr->prob_best_path << endl;
    }
    else
    {
        cout<<"Couldn't find vcurr" << endl;
        return;
    }

    // populate best_path
    cout<<"writing best path"<<endl;
    best_path.clear();
    while( vcurr != NULL)
    {
        best_path.push_back( vcurr->s );
        cout<< vcurr->prob_best_path << endl;
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


void Graph::get_kalman_path()
{
    system->get_kalman_path(obs, kalman_path);
    return;
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
        double totprob = 0;
        cout<<"eo: " << endl;
        for(list<Edge *>::iterator j = v->edges_out.begin(); j != v->edges_out.end(); j++)
        {
            cout<<"\t "<< counte++ << " " << (*j)->transition_prob << endl;
            totprob += (*j)->transition_prob;
        }
        cout<<"totprob: "<< totprob << endl;
    }
}

void Graph::propagate_system()
{
    int num_system = ( system->max_states[0] - system->min_states[0])/system->sim_time_delta;
    int count = obs_interval;

    truth.push_back(system->init_state);

    for(int i =0; i< num_system; i++)
    {
        State snext = system->integrate( truth.back(), system->sim_time_delta, false);
        truth.push_back( snext);

        count--;
        if(count == 0)
        {
            count = obs_interval;
            State next_obs = system->observation( snext, false);
            obs.push_back(next_obs);
        }
    }
}

void Graph::put_init_samples()
{
    double totprob = 0;
    for(int i=0; i < 50; i++)
    {
        State stmp = system->sample();
        Vertex *v = new Vertex(stmp);


        v->s.x[0] = system->min_states[0];
        multivar_normal( &(system->init_state.x[1]), system->init_var, &(v->s.x[1]), NUM_DIM-1);

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
        double tmp = v->prob_best_path/totprob;
        v->prob_best_path = tmp;
    }

    /*
    // just add one vertex at the mean
    State stmp = system->init_state;
    Vertex* v = new Vertex(stmp);
    v->cost_best_path = 0;
    vlist.push_back(v);
    num_vert++;
    insert_into_kdtree(v);
    */

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
    // keep in multiplicative form for simulation
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
        /*
           if( (vcurr->edges_out.size() > 0) && (which_edge == NULL))
           {
           cout<<"error- rand: "<< rand_tmp<<" runner: "<< runner <<" num_edges: "<< vcurr->edges_out.size() << endl;
           getchar();
           }
           */
        if(vcurr->edges_out.size() == 0)
        {
            can_go_forward = false;
            break;
        }

        if(which_edge != NULL)
        {
            seq.push_back(vcurr->s);
            traj_prob *= which_edge->transition_prob;
        }
    }

    //if( fabs(vcurr->s.x[0] - system->max_states[0]) < system->sim_time_delta)
    if(1)
    {
        //cout<<"traj_prob: "<<traj_prob << endl; 
        monte_carlo_trajectories.push_back( seq);
        monte_carlo_probabilities.push_back(traj_prob);
        return 0;
    }
    return 1;
}


int Graph::prune_graph()
{
    int npruned = 0;
    for(vector<Vertex*>::reverse_iterator i = vlist.rbegin(); i!= vlist.rend(); i++)
    {
        Vertex* v = *i;
        if( (v->edges_in.size() == 0) && ( fabs(v->s.x[0] - system->min_states[0]) > system->sim_time_delta) )
        {
            remove_vertex(v);
            npruned++;
        }
    }
    return npruned;
}

void Graph::plot_monte_carlo_trajectories()
{
    ofstream mout("data/monte_carlo.dat");
    int count = 0;

    double totprob = 0;
    for(list<double>::iterator i = monte_carlo_probabilities.begin(); i!= monte_carlo_probabilities.end(); \
            i++)
    {
        totprob += (*i);
    }

    list<double>::iterator prob_iter = monte_carlo_probabilities.begin();
    for(list< list<State> >::iterator i= monte_carlo_trajectories.begin(); i != monte_carlo_trajectories.end(); \
            i++)
    {
        list<State> curr_traj = (*i);
        double curr_prob = (*prob_iter)/totprob;
        mout<<count<<"\t"<< curr_prob <<"\t"<<endl;
        for(list<State>::iterator j = curr_traj.begin(); j != curr_traj.end(); j++)
        {
            State& curr_state = *j;
            for(int k=0; k< NUM_DIM; k++)
            {
                mout<< curr_state.x[k]<<"\t";
            }
            for(int k=NUM_DIM; k< 4; k++)
            {
                mout<< 0 <<"\t";
            }

            mout<<endl;
        }

        prob_iter++;
        count++;
    }

    mout.close();
}   

