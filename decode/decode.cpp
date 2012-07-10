#include "decode.h"

Vertex::Vertex(State& st)
{
    s = st;

    prev = NULL;
    prob_best_path = -100;

    edges_in.clear();
    edges_out.clear();
}

Edge::Edge(Vertex *f, Vertex *t, double prob, double trans_time){
    from = f;
    to = t;
    transition_prob = prob;
    transition_time = trans_time;
}

Graph::Graph(System& sys) 
{
    system = &sys;

    vlist.clear();
    num_vert = 0;
    obs_curr_index = 0;
    expt_time = system->max_states[0];
    
    sim_delta = 1e-3;
    delta = sim_delta;

    min_holding_time = delta;
    seeding_finished = false;

    state_tree = kd_create(NUM_DIM);
    time_tree = kd_create(1);

    double factor = 1;
    if(NUM_DIM == 2)
        factor = M_PI;
    else if(NUM_DIM == 3)
        factor = 4/3*M_PI;
    else if(NUM_DIM == 4)
        factor = 0.5*M_PI*M_PI;
    else if(NUM_DIM == 5)
        factor = 8*M_PI*M_PI/15;
    
    factor = 1;
    gamma = 2.2*pow( (1+1/(double)NUM_DIM), 1/(double)NUM_DIM) *pow(factor, -1/(double)NUM_DIM);
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

    kd_free(state_tree);
    kd_free(time_tree);
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
    
#if 0
    for(list<Edge *>::reverse_iterator i = v->edges_in.rbegin(); i != v->edges_in.rend(); i++)
    {
        Edge* etmp = (*i);
        elist.erase(etmp->elist_iter);
        etmp->from->edges_out.erase(etmp->from_iter);
        delete etmp;
    }
    v->edges_in.clear();
#endif

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

        curr_time += sim_delta;
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

    traj<<"best_path"<<endl;
    for(list<State>::iterator i= best_path.begin(); i != best_path.end(); i++)
    {
        State& curr = *i;
        traj<< curr.x[0] <<"\t";
        for(int j=0; j< NUM_DIM; j++)
        {
            traj<< curr.x[j]<<"\t";
        }
        traj<<endl;
    }

#if 0
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
#endif

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

double Graph::get_obs_prob_vertex(Vertex* v)
{
    State gx = system->observation(v->s, true);
    State closest_obs = obs[(int)((v->s.x[0] - system->min_states[0])/sim_delta)];
    double toret = normal_val( &(closest_obs.x[1]), system->obs_noise, &(gx.x[1]), NUM_DIM_OBS-1);
    //cout<<toret<<endl;
    return toret;
}

double Graph::get_obs_prob_edge(Edge* etmp)
{
    return 0;
}

void Graph::update_viterbi(Vertex* v)
{
#if 1
    double max_prob = -1;
    double tmp;
    
    for(list<Edge*>::iterator i = v->edges_in.begin(); i!= v->edges_in.end(); i++)
    {
        Edge* etmp = *i;

        tmp = (etmp->from->prob_best_path)*(etmp->transition_prob)*get_obs_prob_vertex(v);

        if( (tmp > max_prob) && (tmp > 0))
        {
            max_prob = tmp;
            v->prev = etmp->from;
            v->prob_best_path = max_prob;
        }
    }
    //cout<<"update_viterbi- max_prob: "<< max_prob << endl;
#endif
}

void Graph::propagate_viterbi(Vertex* v)
{
#if 1
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
            Vertex* vto = etmp->to;
            
            double tmp = (vtmp->prob_best_path)*(etmp->transition_prob)*get_obs_prob_vertex(vto);
#if 1 
            // add to queue if its probability is less than the current one
            if( vto->prob_best_path < tmp )
            {
                myq.push_back(etmp->to);
                myq_size++;
            }
#endif
        }
        if( myq_size > max_size)
            max_size = myq_size;
    }
    //cout<<"max_size: "<< max_size << endl;
#endif
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
    State stmp;
    
    float p = RANDF;

    if(p < 0.05)
    {
        multivar_normal(&(system->init_state.x[1]), system->init_var, &(stmp.x[1]), NUM_DIM-1);
        stmp.x[0] = system->min_states[0];
    }
    else if(p < 0.15)
    {
        stmp = system->sample();
        stmp.x[0] = system->max_states[0];
    }
    else
        stmp = system->sample();
    
    Vertex *v = new Vertex(stmp);

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

    double bowlr = gamma/2 * pow( log(num_vert)/(num_vert), 1.0/(double)NUM_DIM);

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

    kdres *res;
    res = kd_nearest_range(state_tree, key, bowlr );
    //cout<<"got "<<kd_res_size(res)<<" states in bowlr= "<< bowlr << endl;
    //int pr = kd_res_size(res);
   
    double *next_state = new double[NUM_DIM];
    double *sys_var = new double[NUM_DIM];

    double pos[NUM_DIM] = {0};
    while( !kd_res_end(res) )
    {
        Vertex* v1 = (Vertex* ) kd_res_item(res, pos);

        if(v1 != v)
        {
            double delta_t = v1->s.x[0] - v->s.x[0];
            if(delta_t > 1e-5)
            {
                system->get_fdt(v->s, delta_t, next_state);
                system->get_variance(v->s, delta_t, sys_var);

                double prob_tmp = normal_val(next_state, sys_var, v1->s.x, NUM_DIM);

                if(prob_tmp > 0)
                {
                    Edge *e1 = new Edge(v, v1, prob_tmp, delta_t);
                    
                    if(1)
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
        }
        kd_res_next(res);
    }
    kd_res_free(res);

    normalize_edges(v);

    delete[] next_state;
    delete[] sys_var;

    return 0;
}


bool compare_vlist_times_pair(pair<double, Vertex*> p1, pair<double, Vertex*> p2)
{
    if(p1.first < p2.first)
        return true;
    else
        return false;
}

void Graph::get_best_path()
{
    list< pair<double, Vertex*> > vlist_times;
    for(unsigned int i=0; i< vlist.size(); i++)
    {
        Vertex* vcurr = vlist[i];
        vlist_times.push_back( make_pair(vcurr->s.x[0],vcurr));
    }
    vlist_times.sort(compare_vlist_times_pair);
     
    list< pair<double, Vertex*> >::iterator runner = vlist_times.begin();
    Vertex* best_end_vertex = NULL;
    double max_prob = -1;

    while(runner != vlist_times.end())
    {
        Vertex* vcurr = (*runner).second;
    
        if( vcurr->prob_best_path < 0 )
        {
            update_viterbi(vcurr);
        }
        
        if( (system->max_states[0] - vcurr->s.x[0]) < 1e-5)
        {
            if(vcurr->prob_best_path > max_prob)
            {
                max_prob = vcurr->prob_best_path;
                best_end_vertex = vcurr;
            }
        }
        //cout<< vcurr->s.x[0] << " " << vcurr->prob_best_path << endl; 
        runner++;
    }
    //cout<<"Found best_vertex: "<< best_end_vertex->prob_best_path << endl;

    best_path.clear();
    while(best_end_vertex != NULL)
    {
        best_path.push_back(best_end_vertex->s);
        best_end_vertex = best_end_vertex->prev;
    }
}

void Graph::get_kalman_path()
{
    system->get_kalman_path(obs, obs_times, kalman_path, kalman_covar);
    return;
}

void Graph::get_pf_path()
{
    system->get_pf_path(obs, obs_times, pf_path);
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
    
    double delta_t = sim_delta;
    double curr_time = 0;

    truth.push_back(system->init_state);

    while(curr_time < expt_time)
    {
        State snext = system->integrate( truth.back(), delta_t, false);
        truth.push_back( snext);

        State next_obs = system->observation( snext, false);
        obs.push_back(next_obs);
        obs_times.push_back(curr_time);
        curr_time += delta_t;
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

// returns 
int Graph::simulate_trajectory()
{
    list<State> seq;
    list<double> seq_times;

    State stmp = system->init_state;
    //multivar_normal(&(system->init_state.x[1]), system->init_var, &(stmp.x[1]), NUM_DIM-1);

    bool can_go_forward = false;

    Vertex *vcurr = nearest_vertex(stmp);

    if(vcurr->edges_out.size() > 0)
        can_go_forward = true;
    else
        cout<<"vcurr has zero edges" << endl;

    seq.push_back(vcurr->s);
    seq_times.push_back(0);

    // keep in multiplicative form for simulation
    double traj_prob = vcurr->prob_best_path;
    double curr_time = 0;
    double max_time = expt_time;

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

        if(vcurr->edges_out.size() == 0)
        {
            can_go_forward = false;
            cout<<"break line 684 size: " << seq.size() << endl;
            break;
        }

        if(which_edge != NULL)
        {
            seq.push_back(vcurr->s);
            seq_times.push_back(curr_time);
            traj_prob = traj_prob * which_edge->transition_prob;

            curr_time += which_edge->transition_time;
            if(curr_time > max_time)
            {
                cout<<"finished one sim" << endl;
                break;
            }
        }
        //cout<<curr_time << " ";
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
        double curr_prob = (*prob_iter)/totprob;

        list<double> time_seq = (*times_iter);
        list<double>::iterator time_seq_iter = time_seq.begin();
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

