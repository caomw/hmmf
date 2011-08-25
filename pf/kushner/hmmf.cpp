#include "hmmf.h"


Vertex::Vertex(State& st)
{
    s = st;

    prev = NULL;
    next = NULL;
    prob_best_path = -1;

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

    state_tree = kd_create(NUM_DIM);

    gamma = 2.2;
    gamma_t = gamma/25/3/(system->max_states[0]);
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
        float totprob = 0;
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
    float max_prob = -1;
    for(list<Edge*>::iterator i= v->edges_in.begin(); i != v->edges_in.end(); i++)
    {
        Edge *etmp = *i;
        float tocalci = (etmp->transition_prob)*(etmp->from->prob_best_path);
        if( tocalci > max_prob )
        {
            max_prob = tocalci;
            v->prev = etmp->from;
            v->prob_best_path = tocalci;
        }
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
    
    connect_edges_approx(v);

    vlist.push_back(v);
    num_vert++;
    insert_into_kdtree(v);
}

// sys: xdot = -3x + w
// obs: y = x + n
int Graph::connect_edges_approx(Vertex* v)
{
    double key[NUM_DIM] ={0};
    system->get_key(v->s, key);

    float bowlr = gamma * pow( log(num_vert)/(num_vert), 1/(float)NUM_DIM);
    float holding_time = gamma_t * pow( log(num_vert)/(num_vert), 1/(float)NUM_DIM);
    //cout<<"holding time: "<< holding_time << endl;
    
    kdres *res;
    res = kd_nearest_range(state_tree, key, bowlr );
    //cout<<"got "<<kd_res_size(res)<<" states in bowlr= "<< bowlr << endl;

    int pr = kd_res_size(res);
        
    MatrixXf Z1(NUM_DIM, pr);
    MatrixXf Z2( (NUM_DIM)*(NUM_DIM), pr);
    VectorXf q0(pr, 1);
    VectorXf q1(pr, 1);
    VectorXf fdt( (NUM_DIM), 1);
    VectorXf FFdt( (NUM_DIM)*(NUM_DIM), 1 );
    
    for(int i=0; i <NUM_DIM; i++)
        fdt[i] = -3*(v->s.x[i])*holding_time;
    
    for(int i=0; i< (NUM_DIM)*(NUM_DIM); i++)
    {
        if (i % (NUM_DIM) == 0)
            FFdt(i) = (system->process_noise[i/(NUM_DIM)])*holding_time*holding_time;       // already covariance
        else
            FFdt(i) = 0;
    }

    double pos[NUM_DIM] = {0};
    int which_num = 0;
    while( !kd_res_end(res) )
    {
        Vertex* v1 = (Vertex* )kd_res_item(res, pos);

        for(int i=0; i< NUM_DIM; i++)
        {
            Z1(i, which_num) = (v1->s.x[i] - v->s.x[i]);

            for(int j=0; j< NUM_DIM; j++)
            {
                Z2(i*(NUM_DIM) + j, which_num) = (v1->s.x[i] - v->s.x[i])*(v1->s.x[j] - v->s.x[j]); 
            }
        }

        which_num++;
        kd_res_next(res);
    }
    
    q0 = Z1.transpose() * (Z1 * Z1.transpose()).inverse() * fdt;        // under-determined


    MatrixXf Z3( (NUM_DIM) + (NUM_DIM)*(NUM_DIM) , pr);
    Z3.topLeftCorner( NUM_DIM, pr) = Z1;
    Z3.bottomLeftCorner( (NUM_DIM)*(NUM_DIM), pr) = Z2;
    
    VectorXf b3( NUM_DIM + (NUM_DIM)*(NUM_DIM) );
    
    for(int i=0; i< NUM_DIM; i++)
        b3[i] = 0;
    
    for(int i= NUM_DIM; i< NUM_DIM*NUM_DIM + NUM_DIM; i++)
    {
        b3[i] = FFdt(i - NUM_DIM);
    }

    // Z3*q1 = b3 is overdetermined now, add a constraint q1 > 0 to it
    
    Vector = (Z3.transpose() * Z3).inverse() * Z3.transpose() * b3;

    float sum = 0;
    for(int i=0; i< pr; i++)
    {
        sum += q0[i];
        sum += q1[i];
    }
    
    float sum_prob = 0;
    which_num = 0;
    kd_res_rewind(res);
    while( !kd_res_end(res) )
    {
        Vertex* v1 = (Vertex* ) kd_res_item(res, pos);

        Edge *e1 = new Edge(v, v1, (q0[which_num] + q1[which_num])/ sum, holding_time);
        sum_prob += (q0[which_num] + q1[which_num])/ sum;
        
        //assert((q0[which_num] + q1[which_num])/ sum > 0 );
        
        //assert(q1[which_num] > 0);

        if( is_edge_free(e1))
        {
            elist.push_back(e1);

            v->edges_out.push_back(e1);
            v1->edges_in.push_back(e1);

        }
        else
            delete e1;

        which_num++;

        kd_res_next(res);
    }
    kd_res_free(res);
    
    Edge* e1 = new Edge(v, v, 1 - sum_prob, holding_time);
    elist.push_back(e1);
    v->edges_out.push_back(e1);
    v->edges_in.push_back(e1);

    //normalize_edges(v);

    //propagate_viterbi(v);

    //cout<<"getchar: "; getchar();

    return 0;
}

void Graph::get_best_path()
{

#if 0
    float pos[NUM_DIM] = {0};
    float to_query[NUM_DIM];
    cout<<"truth: "<< truth.back().x[0] <<" " << truth.back().x[1] << endl;
    system->get_key(truth.back(), to_query);

    kdres *res;
    res = kd_nearest_range(time_tree, &(to_query[0]), goal_width);

    float max_prob = -1;
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
#endif

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
    float curr_time = 0;
    float max_time = 1.0;
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

void Graph::put_init_samples()
{
#if 0
    for(int i=0; i < 50; i++)
    {
        add_sample();
    }

    // normalize the init_samples
    float totprob = 0;
    for(vector<Vertex*>::iterator i = vlist.begin(); i != vlist.end(); i++)
    {
        multivar_normal( system->init_state.x, system->init_var, v->s.x, NUM_DIM);
        v->prob_best_path = normal_val( system->init_state.x, system->init_var, v->s.x, NUM_DIM);      
        
        totprob += v->prob_best_path;
        Vertex* v = *i;
        float tmp = v->prob_best_path/totprob;
        v->prob_best_path = tmp;
    }
#endif

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

// returns 
int Graph::simulate_trajectory()
{
    list<State> seq;

    State stmp = system->init_state;

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
    float max_time = 1.0;

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
                cout<<"finished one sim" << endl;
                break;
            }
        }
        cout<<curr_time << " ";
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
    for(list<float>::iterator i = monte_carlo_probabilities.begin(); i!= monte_carlo_probabilities.end(); \
            i++)
    {
        totprob += (*i);
    }

    list<float>::iterator prob_iter = monte_carlo_probabilities.begin();
    for(list< list<State> >::iterator i= monte_carlo_trajectories.begin(); i != monte_carlo_trajectories.end(); \
            i++)
    {
        list<State> curr_traj = (*i);
        float curr_prob = (*prob_iter)/totprob;
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

