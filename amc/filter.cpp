#include "filter.h"

Vertex::Vertex(State& st)
{
    s = st;

    prev = NULL;
    prob_best_path = -1;
    prob_best_path_buffer = -1;
    obs_update_time = 0;

    self_transition_prob = 0;

    edges_in.clear();
    edges_out.clear();
}

Edge::Edge(Vertex *f, Vertex *t, double prob, double trans_time){
    from = f;
    to = t;
    transition_prob = prob;
    transition_time = trans_time;
}

Graph::Graph(System& sys) {

    system = &sys;

    vlist.clear();
    num_vert = 0;
    obs_interval = 1;
    obs_curr_index = 0;
    delta = system->sim_time_delta;
    max_obs_time = 10;

    min_holding_time = delta;
    seeding_finished = false;

    state_tree = kd_create(NUM_DIM);

    gamma = 2.1*pow(1+1/(double)NUM_DIM, 1/(double)NUM_DIM);
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
    obs_times.clear();
    kalman_path.clear();
    pf_path.clear();
    best_path.clear();

    monte_carlo_trajectories.clear();
    monte_carlo_probabilities.clear();
    monte_carlo_times.clear();
    actual_trajectories.clear();

    kd_free(state_tree);
}

int Graph::vertex_delete_edges(Vertex* v)
{
    for(list<Edge *>::iterator i = v->edges_out.begin(); i != v->edges_out.end(); i++)
    {
        Edge* etmp = (*i);
        elist.erase(etmp->elist_iter);
        etmp->to->edges_in.erase(etmp->to_iter);
        delete etmp;
    }
    v->edges_out.clear();
    return 0;
}

void Graph::remove_edge(Edge *e)
{
    elist.erase(e->elist_iter);

    //cout<<"--- inside remove_edge" << endl;
    // remove from, from + to lists and then call destructor
    e->from->edges_out.erase(e->from_iter);

    e->to->edges_in.erase(e->to_iter);
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
#if 0
        for(list<Edge*>::iterator eo = tstart->edges_out.begin(); eo != tstart->edges_out.end(); eo++)
        {
            Vertex *tend = (*eo)->to;
            Edge *etmp = (*eo);

            //draw the edge
            rrgout<<tstart->s.x[0]<<"\t"<<tstart->s.x[1]<<"\t"<<tend->s.x[0]<<"\t"<<tend->s.x[1]<<"\t"<<etmp->transition_prob<<"\t"<<etmp->transition_time<<endl;
        }
#endif
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

    int count = 0;
    traj<<"system"<<endl;
    double curr_time =0;
    for(list<State>::iterator i= truth.begin(); i != truth.end(); i++)
    {
        traj<< curr_time<<"\t";
        State& curr = *i;
        for(int j=0; j< NUM_DIM; j++)
        {
            traj<< curr.x[j]<<"\t";
        }
        traj<<endl;

        curr_time += system->sim_time_delta;
    }

    count = 0;
    traj<<"observation"<<endl;
    for(vector<State>::iterator i= obs.begin(); i != obs.end(); i++)
    {
        traj<< obs_times[count]<<"\t";
        State& curr = *i;
        for(int j=0; j< NUM_DIM; j++)
        {
            traj<< curr.x[j]<<"\t";
        }
        traj<<endl;
        count++;
    }

    curr_time =0;
    traj<<"best_path"<<endl;
    for(list<State>::iterator i= best_path.begin(); i != best_path.end(); i++)
    {
        traj<< curr_time<<"\t";
        State& curr = *i;
        for(int j=0; j< NUM_DIM; j++)
        {
            traj<< curr.x[j]<<"\t";
        }
        traj<<endl;
        curr_time += system->sim_time_delta;
    }

    curr_time =0;
    traj<<"kf_path"<<endl;
    for(list<State>::iterator i= kalman_path.begin(); i != kalman_path.end(); i++)
    {
        traj<< curr_time<<"\t";
        State& curr = *i;
        for(int j=0; j< NUM_DIM; j++)
        {
            traj<< curr.x[j]<<"\t";
        }
        traj<<endl;
        curr_time += system->sim_time_delta;
    }
    curr_time =0;
    traj<<"pf_path"<<endl;
    for(list<State>::iterator i= pf_path.begin(); i != pf_path.end(); i++)
    {
        traj<< curr_time<<"\t";
        State& curr = *i;
        for(int j=0; j< NUM_DIM; j++)
        {
            traj<< curr.x[j]<<"\t";
        }
        traj<<endl;
        curr_time += system->sim_time_delta;
    }

    curr_time =0;
    traj<<"kf_covar"<<endl;
    for(list<State>::iterator i= kalman_covar.begin(); i != kalman_covar.end(); i++)
    {
        traj<< curr_time<<"\t";
        State& curr = *i;
        for(int j=0; j< NUM_DIM; j++)
        {
            traj<< curr.x[j]<<"\t";
        }
        traj<<endl;
        curr_time += system->sim_time_delta;
    }

    cout<<"trajectory plotting done" << endl;
    traj.close();
}


int Graph::insert_into_kdtree(Vertex *v)
{
    double key[NUM_DIM];
    system->get_key(v->s, key);

    kd_insert(state_tree, key, v);

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
    else if(nedges > 0)
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

void Graph::propagate_density(Vertex* v)
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
        update_density_implicit(vtmp);

        myq.pop_front();
        myq_size--;

        for(list<Edge*>::iterator i = vtmp->edges_out.begin(); i != vtmp->edges_out.end(); i++)
        {
            Edge* etmp = *i;
            // add to queue if not updated with latest observation
            if( etmp->to->obs_update_time < obs_times[obs_curr_index] )
            {
                myq.push_back(etmp->to);
                myq_size++;
            }
            //else
            //cout<<"didn't add: " << etmp->to->obs_update_time << " " << obs_times[obs_curr_index] << endl;
        }
        if( myq_size > max_size)
            max_size = myq_size;
    }
    //cout<<"max_size: "<< max_size << endl;
}

bool Graph::update_density_implicit(Vertex* v)
{
    State gx = system->observation(v->s, true);
    double *obs_var = new double[NUM_DIM_OBS];
    system->get_obs_variance(v->s, obs_var);

    double sum = 0;
    for(list<Edge*>::iterator i = v->edges_in.begin(); i!= v->edges_in.end(); i++)
    {
        Edge* etmp = *i;
        sum = sum + etmp->transition_prob_delta * (etmp->from->prob_best_path);
    }
    sum = sum + v->self_transition_prob *(v->prob_best_path); 

    double sumt = sum;
    double tomult = normal_val(obs[obs_curr_index].x, obs_var, gx.x, NUM_DIM_OBS);

    v->prob_best_path_buffer = sum*tomult;
    
    if(v->prob_best_path_buffer != v->prob_best_path_buffer)
    {
        cout<<"implicit update nan prob_best_path - diagnose - sum: " << sumt <<" vprob_bp: "<< v->prob_best_path <<endl;
        getchar();
    }
    
    v->obs_update_time = obs_times[obs_curr_index];
    
    delete[] obs_var;
    
    if(v->prob_best_path_buffer > 100)
        return true;
    
    return false;
}

void Graph::average_density(Vertex* v)
{
    double totprob = 0;
    int count =0;
    for(list<Edge*>::iterator i = v->edges_in.begin(); i!= v->edges_in.end(); i++)
    {
        Edge* etmp = *i;
        totprob += (etmp->to->prob_best_path);
        count++;
    }
    v->prob_best_path = totprob/(double)count;
}

int Graph::calculate_delta()
{
    double min_time = 1000;
    for(unsigned int i=0; i< num_vert; i++)
    {
        Vertex* v = vlist[i];
        if( v->holding_time < min_time)
            min_time = v->holding_time;
    }
    delta = min(0.01, 0.99*min_time);
    //cout<<"delta: "<< delta << endl;
    return 0;
}

int Graph::calculate_probabilities_delta_all()
{
    for(unsigned int i=0; i< num_vert; i++)
    {
        Vertex* v = vlist[i];
        calculate_probabilities_delta(v);
    }
    return 0;
}

int Graph::calculate_probabilities_delta(Vertex* v)
{
    double holding_time = v->holding_time;
    double pself = 1*(1 - delta/holding_time);
    /*
    if(holding_time > delta)
    {
        cout<<"add more nodes to reduce holding time" << endl;
        getchar();
    }
    double pself = holding_time/delta;
    */
    // double pself = 0;
    for(list<Edge*>::iterator j = v->edges_out.begin(); j != v->edges_out.end(); j++)
    {
        Edge* etmp = (*j);
        etmp->transition_prob_delta = (etmp->transition_prob)*(1-pself);
    }

    v->self_transition_prob = pself;
    v->holding_time_delta = holding_time;
    //cout<<"ht_del: "<< v->holding_time<<endl;
    return 0;
}

void Graph::normalize_density()
{
    double totprob = 0;
    for(unsigned int j=0; j < num_vert; j++)
    {
        Vertex* v = vlist[j];
        totprob += (v->prob_best_path);
    }
    for(unsigned int j=0; j < num_vert; j++)
    {
        Vertex* v = vlist[j];
        v->prob_best_path = v->prob_best_path/totprob;
    }
}

void Graph::update_density_implicit_no_obs(Vertex* v)
{
    double sum = 0;
    for(list<Edge*>::iterator i = v->edges_in.begin(); i!= v->edges_in.end(); i++)
    {
        Edge* etmp = *i;
        sum = sum + etmp->transition_prob_delta * (etmp->from->prob_best_path);
    }
    sum = v->self_transition_prob * (sum + v->prob_best_path); 
    
    v->prob_best_path_buffer = sum;
}

void Graph::update_density_implicit_no_obs_all()
{
    for(unsigned int j = 0; j< num_vert; j++)
    {
        Vertex* v = vlist[j];
        update_density_implicit_no_obs(v);
    }
    buffer_prob_copy();
}

void Graph::update_density_implicit_all(double covar)
{
    //double key_covar = 0.1;
    double key_covar = 3*sqrt(covar + system->process_noise[0]*delta)/fabs(system->max_states[0] - system->min_states[0]);
    key_covar = max(key_covar, 0.1);
    //cout<<key_covar<<endl;
    bool to_normalize = false;
    State last_mean = best_path.back();
    State next_mean = system->integrate(last_mean, delta, true);
    
    // zero whole buffer
    for(unsigned int j = 0; j< num_vert; j++)
        vlist[j]->prob_best_path_buffer = 0;
#if 1
    double key[NUM_DIM] ={0};
    system->get_key(next_mean, key);
    kdres *res;
    res = kd_nearest_range(state_tree, key, key_covar);
    double pos[NUM_DIM] = {0};
    while( !kd_res_end(res) )
    {
        Vertex* v1 = (Vertex* ) kd_res_item(res, pos);
        bool retval = update_density_implicit(v1);
        to_normalize |= retval;
        kd_res_next(res);
    }
    kd_res_free(res);
    buffer_prob_copy();
#else
    to_normalize = true;
    for(unsigned int i=0; i< num_vert; i++)
    {
        Vertex* vtmp = vlist[i];
        update_density_implicit(vtmp);
    }
    buffer_prob_copy();
#endif
    if(to_normalize)
        normalize_density();
}

void Graph::buffer_prob_copy()
{
    for(unsigned int j = 0; j< num_vert; j++)
    {
        Vertex* v = vlist[j];
        v->prob_best_path = v->prob_best_path_buffer;
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

Vertex* Graph::add_sample(bool is_seed)
{
    State stmp;
    
    if(is_seed)
        multivar_normal(system->init_state.x, system->init_var, stmp.x, NUM_DIM);
    else
        stmp = system->sample();
    Vertex *v = new Vertex(stmp);

    /*
       if(num_vert != 0)
       {
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
       }
       */

    vlist.push_back(v);
    num_vert++;
    insert_into_kdtree(v);

    return v;
}

void Graph::approximate_density(Vertex* v)
{
#if 1
    v->prob_best_path = 0;

#endif
#if 0
    double tprob = 0;

    double key[NUM_DIM] ={0};
    system->get_key(v->s, key);

    double bowlr = gamma/10 * pow( log(num_vert)/(num_vert), 1.0/(double)NUM_DIM);

    kdres *res;
    res = kd_nearest_range(state_tree, key, bowlr );
    //cout<<"approx size: "<< kd_res_size(res) << " bowlr: " << bowlr << endl;

    double pos[NUM_DIM] = {0};
    while( !kd_res_end(res) )
    {
        Vertex* v1 = (Vertex* ) kd_res_item(res, pos);

        if(v1 != v)
        {
            tprob = tprob + (v1->prob_best_path);
            //cout<<"v1_pbp: "<< v1->prob_best_path << endl;
        }
        kd_res_next(res);
    }

    v->prob_best_path = tprob/((double)kd_res_size(res));
    //cout<<"v_pbp: "<< v->prob_best_path << endl; getchar();

    kd_res_free(res);
    
    double factor = 1 - v->prob_best_path;
    for(unsigned int i=0; i< num_vert; i++)
    {
        Vertex* vtmp = vlist[i];
        if(vtmp != v)
        {
            vtmp->prob_best_path = vtmp->prob_best_path*factor;
        }
    }
#endif

#if 0
    State gx = system->observation(v->s, true);

    double sum = 0;
    //cout<<"Ry: " << Ry << endl;
    for(list<Edge*>::iterator i = v->edges_in.begin(); i!= v->edges_in.end(); i++)
    {
        Edge* etmp = *i;
        sum = sum + etmp->transition_prob_delta * (etmp->from->prob_best_path);
    }
    
    v->prob_best_path = sum*normal_val(obs[obs_curr_index].x,
            system->obs_noise, gx.x, NUM_DIM_OBS);
#endif
}

int Graph::reconnect_edges_neighbors(Vertex* v)
{
#if 1
    double key[NUM_DIM] ={0};
    system->get_key(v->s, key);

    double bowlr = gamma/5 * pow( log(num_vert)/(num_vert), 1.0/(double)NUM_DIM);

    kdres *res;
    res = kd_nearest_range(state_tree, key, bowlr );
    //int pr = kd_res_size(res);
    //cout<<"reconnecting "<<kd_res_size(res)<<" nodes, radius: " << bowlr << endl;

    double pos[NUM_DIM] = {0};
    while( !kd_res_end(res) )
    {
        Vertex* v1 = (Vertex* ) kd_res_item(res, pos);

        if(v1 != v)
        {
            // remove old edges
            
            vertex_delete_edges(v1);
            connect_edges_approx(v1);
        }

        kd_res_next(res);
    }
    kd_res_free(res);
    
#endif

#if 0
    Vertex* v1 = nearest_vertex(v->s);
    vertex_delete_edges(v1);
    connect_edges_approx(v1);
#endif

    return 0;
}

int Graph::connect_edges_approx(Vertex* v)
{
    double key[NUM_DIM] ={0};
    system->get_key(v->s, key);

    double bowlr = gamma * pow( log(num_vert)/(double)(num_vert), 1.0/(double)NUM_DIM);
    //cout<<"bowlr: " << bowlr << endl;

    double holding_time = system->get_holding_time(v->s, gamma, num_vert);
    //cout<< holding_time << endl;
    v->holding_time = holding_time;

    kdres *res;
    res = kd_nearest_range(state_tree, key, bowlr );
    //cout<<"got "<<kd_res_size(res)<<" states in bowlr= "<< bowlr << endl;
    //int pr = kd_res_size(res);

    double *sys_var = new double[NUM_DIM];
    State stmp = system->integrate(v->s, holding_time, true);
    system->get_variance(v->s, holding_time, sys_var);

    double sum_prob = 0;
    double pos[NUM_DIM] = {0};
    while( !kd_res_end(res) )
    {
        Vertex* v1 = (Vertex* ) kd_res_item(res, pos);

        if(v1 != v)
        {
            double prob_tmp = normal_val(stmp.x, sys_var, v1->s.x, NUM_DIM);
            
            if(prob_tmp > 0)
            {
                Edge *e1 = new Edge(v, v1, prob_tmp, holding_time);

                if(1)
                {
                    sum_prob += prob_tmp;
                    
                    elist.push_back(e1);
                    v->edges_out.push_back(e1);
                    v1->edges_in.push_back(e1);

                    e1->elist_iter = elist.end();   e1->elist_iter--;
                    e1->from_iter = v->edges_out.end(); e1->from_iter--;
                    e1->to_iter = v1->edges_in.end();   e1->to_iter--;
                }
                else
                    delete e1;
            }
        }
        kd_res_next(res);
    }
    kd_res_free(res);

    normalize_edges(v);

    delete[] sys_var;

    return 0;
}

void Graph::get_kalman_path()
{
    system->get_kalman_path(obs, obs_times, kalman_path, kalman_covar);
    return;
}

void Graph::get_pf_path(int nparticles)
{
    system->get_pf_path(obs, obs_times, pf_path, nparticles);
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
    truth.clear();
    obs.clear();
    obs_times.clear();
    kalman_path.clear();
    pf_path.clear();
    best_path.clear();

    double curr_time = 0;
    double max_time = max_obs_time;

    State start_state = system->init_state_real;
    truth.push_back(start_state);

    while(curr_time < max_time)
    {
        State snext = system->integrate( truth.back(), system->sim_time_delta, false, true);
        truth.push_back( snext);
        curr_time += system->sim_time_delta;

        State next_obs = system->observation( snext, false);
        obs.push_back(next_obs);
        obs_times.push_back(curr_time);
    }
}

void Graph::put_init_samples(int howmany)
{
    for(int i=0; i < howmany; i++)
    {
        add_sample();
    }

    for(vector<Vertex*>::iterator i = vlist.begin(); i != vlist.end(); i++)
    {
        Vertex* v = *i;
        v->prob_best_path = normal_val( system->init_state.x, system->init_var, v->s.x, NUM_DIM);      
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

int Graph::simulate_trajectory_implicit()
{
    list<State> seq;
    list<State> actual_seq;
    list<double> seq_times;

    State stmp = system->init_state;

    bool can_go_forward = false;

    Vertex *vcurr = nearest_vertex(stmp);

    if(vcurr->edges_out.size() > 0)
        can_go_forward = true;
    else
        cout<<"vcurr has zero edges" << endl;

    double curr_time = 0;
    double max_time = max_obs_time;
    
    double traj_prob = 1.0;
    seq.push_back(vcurr->s);
    actual_seq.push_back(vcurr->s);

    seq_times.push_back(curr_time);

    while(curr_time < max_time)
    {
        State next_actual_state = system->integrate(vcurr->s, vcurr->holding_time_delta, false, true);
        if(next_actual_state.x[0] > system->max_states[0])
            next_actual_state.x[0] = system->max_states[0];
        else if(next_actual_state.x[0] < system->min_states[0])
            next_actual_state.x[0] = system->min_states[0];

        curr_time += (vcurr->holding_time_delta);
        double rand_tmp = RANDF - vcurr->self_transition_prob;
        if(rand_tmp > 0)
        { 
            double runner = 0;
            for(list<Edge*>::iterator eo = vcurr->edges_out.begin(); \
                    eo != vcurr->edges_out.end(); eo++)
            {
                runner += (*eo)->transition_prob_delta;
                if(runner > rand_tmp)
                {
                    vcurr = (*eo)->to;
                    break;
                }
            }
        }
        //cout<<vcurr->holding_time_delta<<" ";

        seq.push_back(vcurr->s);
        actual_seq.push_back(next_actual_state);
        seq_times.push_back(curr_time);
        if(curr_time > max_time)
        {
            //cout<<"finished one sim" << endl;
            break;
        }
    }

    if(1)
    {
        //cout<<"traj_prob: "<<traj_prob << endl;
        monte_carlo_times.push_back(seq_times);
        monte_carlo_trajectories.push_back( seq);
        actual_trajectories.push_back(actual_seq);
        monte_carlo_probabilities.push_back(traj_prob);
        return 0;
    }
    return 1;
}

void Graph::analyse_monte_carlo_trajectories(double& err_m1, double& err_m2)
{
    list< list<State> >::iterator traj_iter = monte_carlo_trajectories.begin();
    list< list<State> >::iterator actual_traj_iter = actual_trajectories.begin();
    double totprob = 0;
    int num_trajs = 0;
    double traj_m1_1 = 0, traj_m1_2 = 0;
    double actual_m1_1 = 0, actual_m1_2 = 0;
    double traj_m2_1 = 0, traj_m2_2 = 0;
    double actual_m2_1 = 0, actual_m2_2 = 0;
    for(list<double>::iterator i = monte_carlo_probabilities.begin(); \
            i!= monte_carlo_probabilities.end(); i++)
    {
        totprob += (*i);
        num_trajs += 1;
        
        State& last_state = (*traj_iter).back();
        State& actual_last_state = (*actual_traj_iter).back();
        traj_m1_1 = traj_m1_1 + last_state.x[0]; // *(*i);
        traj_m1_2 = traj_m1_2 + last_state.x[1]; // *(*i);
        actual_m1_1 = actual_m1_1 + actual_last_state.x[0];
        actual_m1_2 = actual_m1_2 + actual_last_state.x[1];
        
        traj_m2_1 = traj_m2_1 + sq(last_state.x[0]); //*(*i);
        traj_m2_2 = traj_m2_2 + sq(last_state.x[1]); //*(*i);
        actual_m2_1 = actual_m2_1 + sq(actual_last_state.x[0]); //*(*i);
        actual_m2_2 = actual_m2_2 + sq(actual_last_state.x[1]); //*(*i);

        traj_iter++;
        actual_traj_iter++;
    }
    traj_m1_1 = traj_m1_1/(float)num_trajs;
    actual_m1_1 = actual_m1_1/(float)num_trajs;
    traj_m2_1 = traj_m2_1/(float)num_trajs;
    actual_m2_1 = actual_m2_1/(float)num_trajs;
    traj_m1_2 = traj_m1_2/(float)num_trajs;
    actual_m1_2 = actual_m1_2/(float)num_trajs;
    traj_m2_2 = traj_m2_2/(float)num_trajs;
    actual_m2_2 = actual_m2_2/(float)num_trajs;
   
    err_m1 = sqrt(sq(traj_m1_1 - actual_m1_1) + sq(traj_m1_2 - actual_m1_2));
    err_m2 = sqrt(sq(traj_m2_1 - actual_m2_1) + sq(traj_m2_2 - actual_m2_2));
    cout<<num_vert<<" "<< sqrt(sq(traj_m1_1 - actual_m1_1) + sq(traj_m1_2 - actual_m1_2)) << " " << sqrt(sq(traj_m2_1 - actual_m2_1) + sq(traj_m2_2 - actual_m2_2)) << endl;
}

void Graph::clear_monte_carlo_trajectories()
{
    monte_carlo_trajectories.clear();
    monte_carlo_probabilities.clear();
    monte_carlo_times.clear();
    actual_trajectories.clear();
}

void Graph::plot_monte_carlo_trajectories()
{
    ofstream mout("data/monte_carlo.dat");
    int count = 0;

    double totprob = 0;
    for(list<double>::iterator i = monte_carlo_probabilities.begin(); \
            i!= monte_carlo_probabilities.end(); i++)
    {
        totprob += (*i);
    }

    list<double>::iterator prob_iter = monte_carlo_probabilities.begin();
    list< list<double> >::iterator times_iter = monte_carlo_times.begin();

    for(list< list<State> >::iterator i= monte_carlo_trajectories.begin(); \
            i != monte_carlo_trajectories.end(); i++)
    {
        list<State> curr_traj = (*i);
        double curr_prob = (*prob_iter);

        list<double> time_seq = (*times_iter);
        list<double>::iterator time_seq_iter = time_seq.begin();
        for(list<State>::iterator j = curr_traj.begin(); j != curr_traj.end(); j++)
        {
            mout<< (*time_seq_iter) <<"\t";
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
            time_seq_iter++;
        }
        mout<<count<<"\t"<< curr_prob <<"\t"<<endl;

        prob_iter++;
        times_iter++;
        count++;
    }
    mout.close();

    ofstream aout("data/actual_traj.dat");
    count = 0;

    times_iter = monte_carlo_times.begin();

    for(list< list<State> >::iterator i= actual_trajectories.begin(); \
            i != actual_trajectories.end(); i++)
    {
        list<State> curr_traj = (*i);

        list<double> time_seq = (*times_iter);
        list<double>::iterator time_seq_iter = time_seq.begin();
        for(list<State>::iterator j = curr_traj.begin(); j != curr_traj.end(); j++)
        {
            aout<< (*time_seq_iter) <<"\t";
            State& curr_state = *j;
            for(int k=0; k< NUM_DIM; k++)
            {
                aout<< curr_state.x[k]<<"\t";
            }
            for(int k=NUM_DIM; k< 4; k++)
            {
                aout<< 0 <<"\t";
            }

            aout<<endl;
            time_seq_iter++;
        }
        aout<<count<<"\t"<< 1.0 <<"\t"<<endl;

        times_iter++;
        count++;
    }
    aout.close();

} 

void Graph::plot_monte_carlo_density(char* filename)
{
    ofstream dout(filename);
    int count = 0;

    double totprob = 0;
    for(list<double>::iterator i = monte_carlo_probabilities.begin(); \
            i!= monte_carlo_probabilities.end(); i++)
    {
        totprob += (*i);
    }
    //cout<<"totprob: "<< totprob << endl;

    list<double>::iterator prob_iter = monte_carlo_probabilities.begin();
    for(list< list<State> >::iterator i= monte_carlo_trajectories.begin(); \
            i != monte_carlo_trajectories.end(); i++)
    {
        list<State> curr_traj = (*i);
        double curr_prob = (*prob_iter);
        curr_prob = curr_prob/totprob; 
        //cout<<"curr_prob: "<< curr_prob << endl;

        dout << curr_prob <<"\t";
        State& curr_state = curr_traj.back();
        for(int k=0; k< NUM_DIM; k++)
        {
            dout<< curr_state.x[k]<<"\t";
        }
        for(int k=NUM_DIM; k< 4; k++)
        {
            dout<< 0 <<"\t";
        }

        dout<<endl;

        prob_iter++;
        count++;
    }

    dout.close();
}

