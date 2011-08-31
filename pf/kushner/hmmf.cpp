#include "hmmf.h"


Vertex::Vertex(State& st)
{
    s = st;

    prev = NULL;
    prob_best_path = -1;
    prob_best_path_buffer = -1;
    obs_update_time = 0;
    
    self_transition_prob = -1;

    edges_in.clear();
    edges_out.clear();
}


Edge::Edge(Vertex *f, Vertex *t, float prob, float trans_time){
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
    max_obs_time = 0.3;
    delta = system->sim_time_delta;

    state_tree = kd_create(NUM_DIM);

    gamma = 2.5;
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
    best_path.clear();

    kd_free(state_tree);
}

int Graph::vertex_delete_edges(Vertex* v)
{
    for(list<Edge *>::reverse_iterator i = v->edges_out.rbegin(); i != v->edges_out.rend(); i++)
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
    float curr_time =0;
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
    
    count = 0;
    traj<<"kf_path"<<endl;
    for(list<State>::iterator i= kalman_path.begin(); i != kalman_path.end(); i++)
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
        update_density_explicit(vtmp);
        
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

void Graph::update_density_explicit(Vertex* v)
{
    State gx = system->observation(v->s, true);

    float sum = 0;
    //cout<<"Ry: " << Ry << endl;
    for(list<Edge*>::iterator i = v->edges_in.begin(); i!= v->edges_in.end(); i++)
    {
        Edge* etmp = *i;
        sum = sum + normal_val(obs[obs_curr_index].x, system->obs_noise, gx.x, NUM_DIM) * \
              etmp->transition_prob * (etmp->from->prob_best_path); 
    }
    v->prob_best_path_buffer = sum;

    v->obs_update_time = obs_times[obs_curr_index];
}

void Graph::update_density_implicit(Vertex* v)
{
    State gx = system->observation(v->s, true);

    float sum = 0;
    //cout<<"Ry: " << Ry << endl;
    for(list<Edge*>::iterator i = v->edges_in.begin(); i!= v->edges_in.end(); i++)
    {
        Edge* etmp = *i;
        sum = sum + normal_val(obs[obs_curr_index].x, system->obs_noise, gx.x, NUM_DIM) * \
              etmp->transition_prob_delta * (etmp->from->prob_best_path); 
    }
    v->prob_best_path_buffer = sum;

    v->obs_update_time = obs_times[obs_curr_index];
}

int Graph::calculate_delta()
{
    float max_time = -1000;
    for(unsigned int i=0; i< num_vert; i++)
    {
        Vertex* v = vlist[i];
        if( v->holding_time > max_time)
            max_time = v->holding_time;
    }
    delta = max_time;
    
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
    float holding_time = v->holding_time;
    for(list<Edge*>::iterator j = v->edges_out.begin(); j != v->edges_out.end(); j++)
    {
        Edge* etmp = (*j);
        etmp->transition_prob_delta = (etmp->transition_prob)*(1 - holding_time/(holding_time + delta));
    }

    v->self_transition_prob = holding_time/(holding_time + delta);
    v->holding_time_delta = delta*holding_time/(holding_time + delta);
    
    return 0;
}

float Graph::make_holding_time_constant()
{
    float min_holding_time = system->get_min_holding_time(gamma, num_vert);
        
    cout<<"min_hold: "<< min_holding_time << endl;

    for(unsigned int i=0; i< num_vert; i++)
    {
        Vertex* v = vlist[i];
        float pself = 1 - min_holding_time/(v->holding_time);
        v->holding_time = min_holding_time;

        for(list<Edge*>::iterator j = v->edges_out.begin(); j != v->edges_out.end(); j++)
        {
            Edge* etmp = *j;
            etmp->transition_prob = (etmp->transition_prob)*(1 - pself);
            etmp->transition_time = min_holding_time;
        }
        
        //add new edge to itself
        Edge* new_edge = new Edge(v, v, pself, min_holding_time);
        
        elist.push_back(new_edge);
        v->edges_out.push_back(new_edge);
        v->edges_in.push_back(new_edge);
        
        new_edge->elist_iter = elist.end();         new_edge->elist_iter--;
        new_edge->from_iter = v->edges_out.end();   new_edge->from_iter--;
        new_edge->to_iter = v->edges_in.end();      new_edge->to_iter--;
        
        // normalize
        //normalize_edges(v);
    }
    return min_holding_time;
}

void Graph::update_viterbi_neighbors(Vertex* v)
{
    double key[NUM_DIM] ={0};
    system->get_key(v->s, key);

    float bowlr = gamma * pow( log(num_vert)/(num_vert), 1.0/(float)NUM_DIM);

    kdres *res;
    res = kd_nearest_range(state_tree, key, bowlr );
    //cout<<"reconnecting "<<kd_res_size(res)<<" nodes, radius: " << bowlr << endl;
    
    float average = 0;
    double pos[NUM_DIM] = {0};
    while( !kd_res_end(res) )
    {
        Vertex* v1 = (Vertex* ) kd_res_item(res, pos);
       
        if(v1 != v)
        {
            // remove old edges
            average += (v1->prob_best_path);
        }

        kd_res_next(res);
    }

    v->prob_best_path = average/(float)kd_res_size(res);

    kd_res_free(res);
}

void Graph::update_viterbi(Vertex* v)
{
    float max_prob = -1;
    float tmp;
    State gx = system->observation(v->s, true);
    float obs_prob = normal_val(obs[obs_curr_index].x, system->obs_noise, gx.x, NUM_DIM);
    for(list<Edge*>::iterator i = v->edges_in.begin(); i!= v->edges_in.end(); i++)
    {
        Edge* etmp = *i;

        tmp = (etmp->from->prob_best_path)*(etmp->transition_prob)*obs_prob;

        if(tmp > max_prob)
        {
            max_prob = tmp;
            v->prev = etmp->from;
            v->prob_best_path_buffer = max_prob;
        }
    }
    v->obs_update_time = obs_times[obs_curr_index];
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
            
            State gx = system->observation(etmp->to->s, true);
            //float tmp = (vtmp->prob_best_path)*(etmp->transition_prob)* normal_val(obs[obs_curr_index].x, system->obs_noise, gx.x, NUM_DIM);
            // add to queue if not updated with latest observation
            if( etmp->to->obs_update_time < (obs_times[obs_curr_index] - 1e-10) )
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


void Graph::normalize_density()
{
    float totprob = 0;
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

Vertex* Graph::add_sample()
{
    State stmp = system->sample();
    Vertex *v = new Vertex(stmp);
    
    /*
    float extend_dist = 0.1;
    if( vlist.size() != 0)
    {
        Vertex *near = nearest_vertex( v->s );

        float d = dist( v->s, near->s);
        if ( d > extend_dist)
        {
            for(int i=0; i< NUM_DIM; i++)
                v->s.x[i] = near->s.x[i] + (v->s.x[i] - near->s.x[i])*extend_dist/d;
        }
    }
    */

    vlist.push_back(v);
    num_vert++;
    insert_into_kdtree(v);

    return v;
}

int Graph::reconnect_edges_neighbors(Vertex* v)
{
#if 1
    double key[NUM_DIM] ={0};
    system->get_key(v->s, key);

    float bowlr = gamma/2 * pow( log(num_vert)/(num_vert), 1.0/(float)NUM_DIM);

    kdres *res;
    res = kd_nearest_range(state_tree, key, bowlr );
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
    connect_edges_approx(v1);
#endif

    return 0;
}


int Graph::connect_edges_approx(Vertex* v)
{
    double key[NUM_DIM] ={0};
    system->get_key(v->s, key);

    float bowlr = gamma * pow( log(num_vert)/(num_vert), 1.0/(float)NUM_DIM);
    //cout<<"bowlr: " << bowlr << endl;
    
    float holding_time = system->get_holding_time(v->s, gamma, num_vert);
    //cout<< holding_time << endl;
    v->holding_time = holding_time;
     
    kdres *res;
    res = kd_nearest_range(state_tree, key, bowlr );
    //cout<<"got "<<kd_res_size(res)<<" states in bowlr= "<< bowlr << endl;
    //int pr = kd_res_size(res);
    
    float *var = new float[NUM_DIM];
    State stmp = system->integrate(v->s, holding_time, true);
    system->get_variance(v->s, holding_time, var);
        
    float sum_prob = 0;
    double pos[NUM_DIM] = {0};
    while( !kd_res_end(res) )
    {
        Vertex* v1 = (Vertex* ) kd_res_item(res, pos);
       
        if(v1 != v)
        {
            float prob_tmp = normal_val(stmp.x, var, v1->s.x, NUM_DIM);
            if(prob_tmp > 0)
            {
                Edge *e1 = new Edge(v, v1, prob_tmp, holding_time);
                sum_prob += prob_tmp;

                if( 1)
                {
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
    
    calculate_probabilities_delta(v);

    delete[] var;

    return 0;
}

void Graph::get_kalman_path()
{
    system->get_kalman_path(obs, obs_times, kalman_path);
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
        float totprob = 0;
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

    float curr_time = 0;
    float max_time = max_obs_time;
    int count = obs_interval;

    truth.push_back(system->init_state);

    while(curr_time < max_time)
    {
        State snext = system->integrate( truth.back(), system->sim_time_delta, false);
        truth.push_back( snext);

        count--;
        if(count == 0)
        {
            count = obs_interval;
            State next_obs = system->observation( snext, false);
            obs.push_back(next_obs);
            obs_times.push_back(curr_time);
        }
        curr_time += system->sim_time_delta;
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
            float totprob = 0;
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
    list<float> seq_times;

    
    State stmp = system->init_state;
    multivar_normal(system->init_state.x, system->init_var, stmp.x, NUM_DIM);
    
    // keep in multiplicative form for simulation
    float traj_prob = normal_val( system->init_state.x, system->init_var, stmp.x, NUM_DIM);
    float curr_time = 0;
    float max_time = max_obs_time;
    
    bool can_go_forward = false;

    Vertex *vcurr = nearest_vertex(stmp);
    
    if(vcurr->edges_out.size() > 0)
        can_go_forward = true;
    else
        cout<<"vcurr has zero edges" << endl;

    seq.push_back(vcurr->s);
    seq_times.push_back(curr_time);

    while(curr_time < max_time)
    {
        while(1)
        {
            float rand_tmp = RANDF - vcurr->self_transition_prob;
            if(rand_tmp < 0)
            {
                // jump in time
                
                traj_prob = traj_prob * vcurr->self_transition_prob;
                curr_time += (vcurr->holding_time_delta);
                break;
            }
            else
            { 
                float runner = 0;
                for(list<Edge*>::iterator eo = vcurr->edges_out.begin(); \
                        eo != vcurr->edges_out.end(); eo++)
                {
                    runner += (*eo)->transition_prob_delta;
                    if(runner > rand_tmp)
                    {
                        traj_prob = traj_prob * ((*eo)->transition_prob_delta);
                        vcurr = (*eo)->to;
                        break;
                    }
                }
            }
        }
        
        seq.push_back(vcurr->s);
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
        monte_carlo_probabilities.push_back(traj_prob);
        return 0;
    }
    return 1;
}


// returns 
int Graph::simulate_trajectory_explicit()
{
    list<State> seq;

    State stmp = system->init_state;
    multivar_normal(system->init_state.x, system->init_var, stmp.x, NUM_DIM);
    
    bool can_go_forward = false;

    Vertex *vcurr = nearest_vertex(stmp);
    
    if(vcurr->edges_out.size() > 0)
        can_go_forward = true;
    else
        cout<<"vcurr has zero edges" << endl;

    seq.push_back(vcurr->s);
    
    // keep in multiplicative form for simulation
    float traj_prob = normal_val( system->init_state.x, system->init_var, stmp.x, NUM_DIM);
    float curr_time = 0;
    float max_time = max_obs_time;

    while(can_go_forward)
    {
        float rand_tmp = RANDF;
        float runner = 0;
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
        
        if(vcurr->edges_out.size() == 0)
        {
            can_go_forward = false;
            //cout<<"break line 684 size: " << seq.size() << endl;
            break;
        }

        if(which_edge != NULL)
        {
            seq.push_back(vcurr->s);
            traj_prob *= which_edge->transition_prob;
            
            curr_time += which_edge->transition_time;
            if(curr_time > max_time)
            {
                //cout<<"finished one sim" << endl;
                break;
            }
        }
        //cout<<curr_time << " ";
    }

    if(1)
    {
        //cout<<"traj_prob: "<<traj_prob << endl; 
        monte_carlo_trajectories.push_back( seq);
        monte_carlo_probabilities.push_back(traj_prob);
        return 0;
    }
    return 1;
}


void Graph::plot_monte_carlo_trajectories()
{
    ofstream mout("data/monte_carlo.dat");
    int count = 0;

    float totprob = 0;
    for(list<float>::iterator i = monte_carlo_probabilities.begin(); \
            i!= monte_carlo_probabilities.end(); i++)
    {
        totprob += (*i);
    }

    list<float>::iterator prob_iter = monte_carlo_probabilities.begin();
    list< list<float> >::iterator times_iter = monte_carlo_times.begin();

    for(list< list<State> >::iterator i= monte_carlo_trajectories.begin(); \
            i != monte_carlo_trajectories.end(); i++)
    {
        list<State> curr_traj = (*i);
        float curr_prob = (*prob_iter)/totprob;

        list<float> time_seq = (*times_iter);
        list<float>::iterator time_seq_iter = time_seq.begin();
        mout<<count<<"\t"<< curr_prob <<"\t"<<endl;
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

        prob_iter++;
        times_iter++;
        count++;
    }

    mout.close();
}   

