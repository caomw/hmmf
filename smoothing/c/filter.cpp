#include "common.h"

//#include "systems/vanderpol_parameter.h"
#include "systems/vanderpol.h"
//#include "systems/singleint.h"
//#include "systems/singleint1d.h"
//#include "systems/uhlmann.h"
//#include "systems/parameter.h"
//#include "systems/parameter_hard.h"
//#include "systems/ship.h"

class Node
{
    public:
        float x[ndim];
        float key[ndim];
        int index;
        float htime;
        float weight;
        int id;

        vector<int> eout;
        vector<int> ein;
        Node(int iin, float* xin, float win=0)
        {
            weight = win;
            id = 0;
            for(int i=0; i< ndim; i++)
                x[i] = xin[i];
            get_key(xin, key);
            index = iin;
        }
        int get_key(float* xin, float* keyout)
        {
            for(int i=0; i< ndim; i++)
            {
                keyout[i] = (xin[i]-zmin[i])/(zmax[i] - zmin[i]);
            }
            return 0;
        }
};
class Graph
{
    public:
        kdtree* node_tree;
        vector<Node*> nodes;
        int num_vert;
        float delta;
        float* P;
        float bowlr;

        Graph(int n, float delta_in)
        {
            num_vert = 2*n;
            P = new float[num_vert*num_vert];
            memset(P, 0, sizeof(float)*num_vert*num_vert);
            bowlr = 1; 
            delta = delta_in;
            node_tree = kd_create(ndim);
        }
        ~Graph()
        {
            for(int i=0; i< num_vert; i++)
                delete nodes[i];
            kd_free(node_tree);
            delete[] P;
        }
        int add_nodes(vector< vector<float> >& particles_prev, vector< vector<float> >& particles, vector<float>& win)
        {
            for(unsigned int i=0; i< particles_prev.size(); i++)
            {
                Node* n1 = new Node(i, &(particles_prev[i][0]), win[i]);
                n1->id = 0;
                nodes.push_back(n1);
            }
            for(unsigned int i=0; i< particles.size(); i++)
            {
                Node* n1 = new Node(particles_prev.size()+i, &(particles[i][0]), 0);
                n1->id = 1;
                nodes.push_back(n1);
                kd_insertf(node_tree, n1->key, n1);
            }
            return 0;
        }
        int connect_nodes()
        {
            float scaling = 1;
            //cout << "scaling: "<<scaling << endl;
            bowlr = 2.1*pow(1+1/(float)ndim, 1/(float)ndim)*scaling*pow(log(num_vert/2.0)/(float)(num_vert/2.0), 1/(float)ndim);
            for(int i=0; i< num_vert; i++)
            {
                if(nodes[i]->id == 0)
                    connect_nodes_iter(i);
            }
            return 0;
        }
        int connect_nodes_iter(int nodei)
        {
            Node* n = nodes[nodei];
            //n->htime = holding_time(n->x, bowlr);
            n->htime = delta;

            float next_state[ndim] ={0}, next_state_key[ndim]={0};
            copy(n->x, next_state);
            integrate_system(next_state, n->htime, true);
            n->get_key(next_state, next_state_key);
            float var[ndim];
            for(int j=0; j<ndim; j++)
                var[j] = pvar[j]*(n->htime);

            vector<int> neighbors;  
            kdres *res;
            res = kd_nearest_rangef(node_tree, next_state_key, bowlr);
            //cout<<"found: "<< kd_res_size(res) << endl;
            float pos[ndim] = {0};
            if(kd_res_size(res) > 0)
            {
                while( !kd_res_end(res) )
                {
                    Node* n1 = (Node*) kd_res_itemf(res, pos);
                    neighbors.push_back(n1->index);
                    kd_res_next(res);
                }
            }
            else
            {
                res = kd_nearestf(node_tree, next_state_key);
                while( !kd_res_end(res) )
                {
                    Node* n1 = (Node*) kd_res_itemf(res, pos);
                    neighbors.push_back(n1->index);
                    kd_res_next(res);
                }
            }
            kd_res_free(res);

            vector<float> probs(neighbors.size(), 0);
            float tprob = 0;
            for(unsigned int j=0; j<neighbors.size(); j++)
            {
                Node* n1 = nodes[neighbors[j]];
                probs[j] = normal_val(next_state, var, n1->x, ndim);
                tprob = tprob + probs[j];
            }
            for(unsigned int j=0; j<neighbors.size(); j++)
            {
                Node* n1 = nodes[neighbors[j]];
                probs[j] = probs[j]/tprob;
                P[n->index*num_vert + n1->index] = probs[j];
                n1->ein.push_back(n->index);
                //cout<<"pushed back "<<n->index<<" into "<< n1->index << endl;
            }
            n->eout = neighbors;
            return 0;
        }

        int calci_moment(float* retm, float* retv)
        {
            for(int i=0; i<ndim; i++)
            {
                retm[i] = 0;
                retv[i] = 0;
            }
            for(int i=0; i< num_vert; i++)
            {
                if (nodes[i]->id ==0)
                {
                    for(int j=0; j<ndim; j++)
                    {
                        retm[j] = retm[j] + nodes[i]->weight*nodes[i]->x[j];
                        retv[j] = retv[j] + nodes[i]->weight*sq(nodes[i]->x[j]);
                    }
                }
            }
            for(int i=0; i<ndim; i++)
            {
                retv[i] = abs(retv[i] - sq(retm[i]));
            }
            return 0;
        }
        int normalize_density()
        {
            float tprob = 0;
            for(int i=0; i< num_vert; i++)
            {
                if(nodes[i]->id == 0)
                    nodes[i]->weight = 0;
                else
                    tprob = tprob + nodes[i]->weight;
            }
            for(int i=0; i< num_vert; i++)
            {
                nodes[i]->weight = nodes[i]->weight/tprob;
            }
            return 0;
        }
        int update_density(vector<float>& obs_in)
        {
            float obs[ndim_obs] = {0};
            for(int i=0; i<ndim_obs; i++)
                obs[i] = obs_in[i];

            for(int i=0; i< num_vert; i++)
            {
                float to_add = 0;
                Node *n1 = nodes[i];
                if(n1->id == 1)
                {
                    for(unsigned int j=0; j< n1->ein.size(); j++)
                        to_add = to_add + nodes[n1->ein[j]]->weight*P[j*num_vert + i];
                    
                    float curr_obs[ndim_obs] = {0};
                    get_obs(n1->x, curr_obs, true);
                    n1->weight = to_add*normal_val(obs, ovar, curr_obs, ndim_obs);
                }
            }
            normalize_density();
            return 0;
        }
        int copy_weights(vector<float>& weights)
        {
            int np = weights.size();
            for(int i=0; i<np; i++)
            {
                weights[i] = nodes[np+i]->weight;
            }
            return 0;
        }
        int print()
        {
            cout<<"-----------------"<<endl;
            for(int i=0; i<num_vert; i++)
            {
                Node *n1 = nodes[i];
                cout<<i<<" : "<<n1->weight<<" - "; print_state(n1->x);
                /*
                cout<<"index: "<< n1->index<<" id: "<<n1->id<<" eins: "<<n1->ein.size()<<endl;
                for(unsigned int j=0; j< n1->ein.size(); j++)
                    cout<< nodes[n1->ein[j]]->index<<"-"<<nodes[n1->ein[j]]->weight<<" ";
                cout<<endl;
                */
            }
            /*
            cout<<"transition matrix"<<endl;
            for(int i=0;i<num_vert;i++)
            {
                for(int j=0;j<num_vert;j++)
                    cout<<P[i*num_vert + j]<<" ";
                cout<<endl;
            }
            */
            cout<<"-----------------"<<endl;
            getchar();
            return 0;
        }
};
class HMMF
{
    public:
        float delta;

        HMMF()
        {
            delta = 0.01;
        }

        vector< vector<float> > truth;
        vector< vector<float> > observations;
        int propagate_system(float max_time)
        {
            truth.clear();
            observations.clear();
            float curr_time = 0;
            float curr_state[ndim];
            float curr_obs[ndim_obs];
            copy(init_state_real, curr_state);
            while(curr_time <= max_time)
            {
                integrate_system(curr_state, delta);
                curr_time = curr_time + delta;

                vector<float> state_tmp; state_tmp.assign(curr_state, curr_state+ndim);
                truth.push_back(state_tmp);
                get_obs(curr_state, curr_obs);
                vector<float> obs_tmp; obs_tmp.assign(curr_obs, curr_obs+ndim);
                observations.push_back(obs_tmp);

                //cout<< fabs(curr_time/delta - (int)(curr_time/delta)) << endl;
                /*
                   cout<<curr_time<<" ";
                   for(int j=0; j<ndim; j++)
                   cout<<curr_state[j]<<" ";
                   cout<<endl;
                   */
            }

            return 0;
        }
        vector< vector<float> > festimates;
        int run_filter(int num_particles)
        {
            festimates.clear();
            int steps = observations.size();

            vector< vector<float> > particles(num_particles, vector<float>(ndim,0));
            vector< vector<float> > particles_prev(num_particles, vector<float>(ndim,0));
            vector<float> weights(num_particles, 0);
            for(int i=0; i< num_particles; i++)
            {
                multivar_normal(init_state, init_var, &(particles[i][0]), ndim);
                weights[i] = 1.0/(float)num_particles;
            }
            for(int i=0; i<steps; i++)
            {
                particles_prev = particles;

                for(int j=0; j<num_particles; j++)
                    integrate_system(&(particles[j][0]), delta, false);
#if 1
                Graph g(num_particles, delta);
                g.add_nodes(particles_prev, particles, weights);


                g.connect_nodes();
                //g.print();
                g.update_density(observations[i]);
                g.copy_weights(weights);
                
                vector<float> fmean(ndim,0);
                for(int j=0; j<num_particles; j++)
                {
                    for(int k=0; k<ndim; k++)
                        fmean[k] = fmean[k] + weights[j]*particles[j][k];
                }
                festimates.push_back(fmean); 
                
                if(i%10 == 0)
                    pfilter_resample(particles, weights);
#else
                float tot_prob = 0;
                vector<float> fmean(ndim,0);
                for(int j=0; j<num_particles; j++)
                {
                    float particle_obs[ndim_obs] ={0};
                    get_obs(&(particles[j][0]), particle_obs, true);
                    weights[j] = weights[j]*normal_val(&(observations[i][0]), ovar, particle_obs, ndim_obs);
                    tot_prob = tot_prob + weights[j];
                }
                for(int j=0; j<num_particles; j++)
                {
                    weights[j] = weights[j]/tot_prob;
                    for(int k=0; k<ndim; k++)
                        fmean[k] = fmean[k] + weights[j]*particles[j][k];
                }
                festimates.push_back(fmean); 
                
                // resample
                pfilter_resample(particles, weights);
#endif
            }
            return 0;
        }

        vector< vector<float> > pfestimates;
        int pfilter_resample(vector< vector<float> >& particles, vector<float>& weights)
        {
            int np = particles.size();
            float totweight = 0;
            for(int i=0; i<np; i++)
                totweight += weights[i];

            float* cum = new float[np];
            float curr_tot = 0;
            for(int i=0; i<np; i++)
            {
                weights[i] = weights[i]/totweight;
                curr_tot += weights[i];
                cum[i] = curr_tot;

                //reset to equal weights
                weights[i] = 1.0/np;
            }

            vector< vector<float> > pcopy = particles;
#if 1
            float u = RANDF/(float)np;
            int i = 0;
            for(int j=0; j<np; j++)
            {
                float tocheck = u + j/(float)np;
                while(tocheck > cum[i])
                    i++;

                pcopy[j] = particles[i];
            }
            particles = pcopy;
#endif
            delete[] cum;
            // weights are already set above
            return 0;
        }
        int run_pfilter(int num_particles)
        {
            pfestimates.clear();
            int steps = observations.size();

            vector< vector<float> > particles(num_particles, vector<float>(ndim,0));
            vector<float> weights(num_particles, 0);
            for(int i=0; i< num_particles; i++)
            {
                multivar_normal(init_state, init_var, &(particles[i][0]), ndim);
                weights[i] = 1.0/(float)num_particles;
            }

            for(int i=0; i<steps; i++)
            {
                for(int j=0; j<num_particles; j++)
                    integrate_system(&(particles[j][0]), delta);

                float tot_prob = 0;
                vector<float> pfmean(ndim,0), pfvar(ndim, 0);
                for(int j=0; j<num_particles; j++)
                {
                    float particle_obs[ndim_obs] ={0};
                    get_obs(&(particles[j][0]), particle_obs, true);
                    weights[j] = weights[j]*normal_val(&(observations[i][0]), ovar, particle_obs, ndim_obs);
                    tot_prob = tot_prob + weights[j];
                }
                for(int j=0; j<num_particles; j++)
                {
                    weights[j] = weights[j]/tot_prob;
                    for(int i=0; i<ndim; i++)
                    {
                        pfmean[i] = pfmean[i] + weights[j]*particles[j][i];
                        pfvar[i] = pfvar[i] + weights[j]*sq(particles[j][i]);
                    }
                }
                pfestimates.push_back(pfmean); 
                
                for(int i=0; i<ndim; i++)
                {
                    pfvar[i] = pfvar[i] - sq(pfmean[i]);
                }
                //print_state(&(pfvar[0]));
                // resample
                pfilter_resample(particles, weights);
            }
            return 0;
        }
        vector< vector<float> > kfestimates;
        int run_kfilter()
        {
            get_kalman_path(kfestimates, observations, delta);
            return 0;
        }
        int calculate_err(float& ferrt, float& pferrt, float& kferrt, float max_time)
        {
            if((festimates.size() == 0) && (pfestimates.size() == 0) )
                return 1;
            ferrt = 0;
            pferrt = 0;
            kferrt = 0;
            int steps = observations.size();
            for(int i=0; i<steps; i++)
            {
                float ferrc = 0;
                float pferrc = 0;
                float kferrc = 0;
                for(int j=0; j<ndim; j++)
                {
                    ferrc = ferrc + sq(truth[i][j] - festimates[i][j]);
                    pferrc = pferrc + sq(truth[i][j] - pfestimates[i][j]);
                    kferrc = kferrc + sq(truth[i][j] - kfestimates[i][j]);
                }
                ferrt = ferrt + ferrc;
                pferrt = pferrt + pferrc;
                kferrt = kferrt + kferrc;
            }
            ferrt = ferrt*max_time/(float)steps;
            pferrt = pferrt*max_time/(float)steps;
            kferrt = kferrt*max_time/(float)steps;
            return 0;
        }
        int output_trajectories()
        {
            ofstream ot("data/truth.dat");
            for(unsigned int i=0; i< truth.size(); i++)
            {
                for(int j=0; j<ndim; j++)
                    ot<<truth[i][j]<<" ";
                ot<<endl;
            }
            ot.close();
            ofstream ob("data/observations.dat");
            for(unsigned int i=0; i< observations.size(); i++)
            {
                for(int j=0; j<ndim; j++)
                    ob<<observations[i][j]<<" ";
                ob<<endl;
            }
            ob.close();
            ofstream of("data/filter.dat");
            for(unsigned int i=0; i< festimates.size(); i++)
            {
                for(int j=0; j<ndim; j++)
                    of<<festimates[i][j]<<" ";
                of<<endl;
            }
            of.close();
            ofstream pf("data/pfilter.dat");
            for(unsigned int i=0; i< pfestimates.size(); i++)
            {
                for(int j=0; j<ndim; j++)
                    pf<<pfestimates[i][j]<<" ";
                pf<<endl;
            }
            pf.close();
            ofstream kf("data/kfilter.dat");
            for(unsigned int i=0; i< kfestimates.size(); i++)
            {
                for(int j=0; j<ndim; j++)
                    kf<<kfestimates[i][j]<<" ";
                kf<<endl;
            }
            kf.close();

            return 0;
        }
};


int main(int argc, char** argv)
{
    // n per step
    int n = 100, e = 100;
    float max_time = 1;
    if(argc > 1)
        n = atoi(argv[1]);
    if(argc > 2)
        max_time = atof(argv[2]);
    if(argc > 3)
        e = max(n+1, atoi(argv[3]));

#if 1
    srand(time(NULL));
    //srand(0);
    HMMF h;
    h.propagate_system(max_time);
    
    tic();
    h.run_kfilter();
    cout<<"dt kf: "<< toc()<<endl;
    tic();
    h.run_filter(n);
    cout<<"dt hmmf: "<< toc()<<endl;
    tic();
    h.run_pfilter(n);
    cout<<"dt pf: "<< toc()<<endl;
    float ferr=0, pferr=0, kferr=0;
    h.calculate_err(ferr, pferr, kferr, max_time);
    h.output_trajectories();
    cout<<n<<" "<<ferr<<" "<<pferr<<" "<<kferr<<endl;
#endif

#if 0
    int s = n;
    int max_tries = 50;
    for(int i=s; i<e; i=i+10)
    {
        float favg=0, pfavg=0, ftavg=0, pftavg=0;
        for(int tries=1; tries<max_tries+1; tries++)
        {
            HMMF h;
            h.propagate_system(max_time);
            
            tic();
            h.run_pfilter(i);
            pftavg += toc();
            
            tic();
            h.run_filter(i);
            ftavg += toc();

            h.run_kfilter();

            float ferr=0, pferr=0, kferr=0;
            h.calculate_err(ferr, pferr, kferr, max_time);
            favg = favg + ferr;
            pfavg = pfavg + pferr;
            cout<<tries<<" "<<favg/(float)tries<<" "<<pfavg/(float)tries<<endl;
        }
        favg = favg/(float)max_tries;
        pfavg = pfavg/(float)max_tries;
        ftavg = ftavg/(float)max_tries;
        pftavg = pftavg/(float)max_tries;
        cout<<i<<" "<<ftavg<<" "<<pftavg<<" "<<favg<<" "<<pfavg<<endl;
    }
#endif

    return 0;
}

