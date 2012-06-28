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
        Node(int iin, float bowlr, int max_vert)
        {
#if 0
            for(int i=0; i< ndim; i++)
            {
                float r = RANDF;
                x[i] = zmin[i] + (zmax[i]-zmin[i])*r;
                key[i] = r;
            }
#else
            int n = pow(max_vert, 1/(float)ndim);
            int id[ndim]={0};
            int left = iin;
            for(int i=ndim-1; i>=0; i--)
            {
                id[i] = (int)(left/pow(n,i));
                left = left - (int)(id[i]*pow(n,i));
            }
            for(int i=0; i<ndim; i++)
            {
                x[i] = zmin[i] + (zmax[i] - zmin[i])*id[i]/(float)n;
                key[i] = id[i]/(float)n;
            }
#endif
            index = iin;
            htime = holding_time(x, bowlr);
        }
};

class Graph
{
    public:
        kdtree* node_tree;
        vector<Node*> nodes;
        int num_vert;
        float bowlr;
        float* P;
        vector< vector<int> > neighbor_list;
        float delta;
        Graph(int num_nodes)
        {
            num_vert = num_nodes;
            P = new float[num_nodes*num_nodes];
            memset(P, 0, sizeof(float)*num_nodes*num_nodes);
            float gamma = 2.2; //2.1*pow(1+1/(float)ndim, 1/(float)ndim);
            //bowlr = gamma*pow(log(num_vert+1)/(num_vert+1.0), 1/(float)ndim);
            bowlr = 1.4/pow(num_vert, 1/(float)ndim);
            node_tree = kd_create(ndim);
            for(int i=0; i< num_nodes; i++)
            {
                Node* n = new Node(i, bowlr, num_vert);
                nodes.push_back(n);
                kd_insertf(node_tree, n->key, n);
            }
            delta = -1;
        }
        ~Graph()
        {
            for(int i=0; i< num_vert; i++)
                delete nodes[i];
            kd_free(node_tree);
            delete[] P;
        }
        float min_htime()
        {
            float min_ht = 1000;
            for(int i=0; i< num_vert; i++)
            {
                if(nodes[i]->htime < min_ht)
                    min_ht = nodes[i]->htime;
            }
            return min_ht;
        }
        int connect_nodes()
        {
            float pos[ndim] = {0};
            delta = 0.99*min_htime();
            for(int i=0; i< num_vert; i++)
            {
                Node* n = nodes[i];
                vector<int> neighbors;        
                kdres *res;
                res = kd_nearest_rangef(node_tree, n->key, bowlr);
                while( !kd_res_end(res) )
                {
                    Node* n1 = (Node*) kd_res_itemf(res, pos);
                    if(n1 != n)
                        neighbors.push_back(n1->index);
                    kd_res_next(res);
                }
                kd_res_free(res);
                
                neighbor_list.push_back(neighbors);
                vector<float> probs(neighbors.size());
                float tprob = 0;

                float next_state[ndim] ={0};
                float var[ndim];
                copy(n->x, next_state);
                integrate_system(next_state, n->htime, true);
                for(int j=0; j<ndim; j++)
                    var[j] = pvar[j]*(n->htime);
                
                for(unsigned int j=0; j<neighbors.size(); j++)
                {
                    Node* n1 = nodes[neighbors[j]];
                    probs[j] = normal_val(next_state, var, n1->x, ndim);
                    tprob = tprob + probs[j];
                }
                float ps = 1 - delta/n->htime;
                P[n->index*num_vert + n->index] = ps;
                for(unsigned int j=0; j<neighbors.size(); j++)
                {
                    probs[j] = probs[j]/tprob;
                    P[n->index*num_vert + neighbors[j]] = (1-ps)*probs[j];
                }
            }
            
            return 0;
        }
        int integrate_system(float* curr_state, float dt, bool is_clean=false)
        {
            float integration_delta = min(1e-3, dt/2.0);
            float runner_time = 0;
            while(runner_time < dt)
            {
                float next_state_delta[ndim] ={0};
                drift(curr_state, next_state_delta, integration_delta, true);
                float noise[ndim] = {0};
                if(!is_clean)
                {
                    diffusion(curr_state, noise, integration_delta, true);
                    add(noise, next_state_delta);
                }
                add(next_state_delta, curr_state);
                runner_time = runner_time + integration_delta;
            }
            return 0;
        }
        vector< vector<float> > truth;
        vector< vector<float> > observations;
        int propagate_system(float max_time)
        {
            truth.clear();
            observations.clear();
            float curr_time = 0;
            float curr_state[ndim];
            float curr_obs[ndim];
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

        int calci_mean(vector< vector<float> >& density, int index, float* ret)
        {
            for(int i=0; i<ndim; i++)
                ret[i] = 0;

            float tprob = 0;
            for(int i=0; i<num_vert; i++)
            {
                for(int j=0; j<ndim; j++)
                {
                    ret[j] = ret[j] + density[index][i]*nodes[i]->x[j];
                }
                tprob = tprob + density[index][i];
            }
            for(int i=0; i<ndim; i++)
                ret[i] = ret[i]/tprob;

            return 0;
        }
        // use index alphas to get index+1 alphas
        int update_alphas(vector< vector<float> >& alphas, vector<float>& obs_in, int index)
        {
            float obs[ndim] = {0};
            for(int i=0; i<ndim; i++)
                obs[i] = obs_in[i];

            float tprob = 0;
            for(int i=0; i<num_vert; i++)
            {
                float toadd = 0;
                for(unsigned int j=0; j< neighbor_list[i].size(); j++)
                {
                    int neighbor_j = neighbor_list[i][j];
                    if(alphas[index][neighbor_j] > -1e-10)
                        toadd = toadd + P[neighbor_j*num_vert + i]*alphas[index][neighbor_j];
                }
                float curr_obs[ndim] = {0};
                get_obs(nodes[i]->x, curr_obs, true);
                alphas[index+1][i] = pow(10.0,9)*toadd*normal_val(obs, ovar, curr_obs, ndim_obs);
                tprob = tprob + alphas[index+1][i];
            }
            for(int i=0; i<num_vert; i++)
            {
                alphas[index+1][i] = alphas[index+1][i]/tprob;
            }
            return 0;
        }
        // use index betas to get index-1 betas
        int update_betas(vector< vector<float> >& betas, vector<float>& obs_in, int index)
        {
            float obs[ndim] = {0};
            for(int i=0; i<ndim; i++)
                obs[i] = obs_in[i];

            float tprob = 0;
            for(int i=0; i<num_vert; i++)
            {
                float toadd = 0;
                for(unsigned int j=0; j< neighbor_list[i].size(); j++)
                {
                    int neighbor_j = neighbor_list[i][j];
                    if(betas[index][neighbor_j] > -1e-40)
                    {
                        float curr_obs[ndim] = {0};
                        get_obs(nodes[neighbor_j]->x, curr_obs, true);
                        toadd = toadd + pow(10.0,9)*P[i*num_vert + neighbor_j]*betas[index][neighbor_j]*normal_val(obs, ovar, nodes[neighbor_j]->x, ndim_obs);
                    }
                }
                betas[index-1][i] = toadd;
                tprob = tprob + betas[index-1][i];
            }
            for(int i=0; i<num_vert; i++)
            {
                betas[index-1][i] = betas[index-1][i]/tprob;
            }
            return 0;
        }

        vector< vector<float> > festimates;
        vector< vector<float> > sestimates;
        int run_smoother()
        {
            int steps = observations.size();
            vector< vector<float> > alphas(steps+1, vector<float>(num_vert,0));
            vector< vector<float> > betas(steps, vector<float>(num_vert,0));
            vector< vector<float> > sdensity(steps, vector<float>(num_vert,0));

            float tprob = 0;
            for(int i=0; i<num_vert; i++)
            {
                alphas[0][i] = normal_val(init_state, init_var, nodes[i]->x, ndim);
                tprob = tprob + alphas[0][i];
                betas[steps-1][i] = 1.0;
            }
            for(int i=0; i<num_vert; i++)
                alphas[0][i] = alphas[0][i]/tprob;

            for(int i=0; i<steps; i++)
            {
                update_alphas(alphas, observations[i], i);
                
                //if(i%(int)(steps/10.0) == 0)
                //    cout<<"a: "<< i << endl;
            }
            for(int i=steps-1; i >= 1; i--)
            {
                //if(i%(int)(steps/10.0) == 0)
                //    cout<<"b: "<< i << endl;
                update_betas(betas, observations[i], i);
            }
            for(int i=0; i<steps; i++)
            {
                for(int j=0; j<num_vert; j++)
                    sdensity[i][j] = alphas[i+1][j]*betas[i][j];
            }
            for(int i=0; i< steps; i++)
            {
                float cmean[ndim];

                calci_mean(alphas, i+1, cmean);
                vector<float> tmp1; tmp1.assign(cmean, cmean+ndim);
                festimates.push_back(tmp1);

                calci_mean(sdensity, i, cmean);
                vector<float> tmp2; tmp2.assign(cmean, cmean+ndim);
                sestimates.push_back(tmp2);
            }

            return 0;
        }

        vector< vector<float> > pfestimates;
        int pfilter_resample()
        {
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
                vector<float> pfmean(ndim,0);
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
                        pfmean[i] = pfmean[i] + weights[j]*particles[j][i];
                }
                pfestimates.push_back(pfmean); 
            }
            return 0;
        }
        int calculate_err(float& ferrt, float& serrt, float& pferrt, float max_time)
        {
            ferrt = 0;
            serrt = 0;
            pferrt = 0;
            int steps = observations.size();
            for(int i=0; i<steps; i++)
            {
                float ferrc = 0;
                float serrc = 0;
                float pferrc = 0;
                for(int j=0; j<ndim; j++)
                {
                    ferrc = ferrc + sq(truth[i][j] - festimates[i][j]);
                    serrc = serrc + sq(truth[i][j] - sestimates[i][j]);
                    pferrc = pferrc + sq(truth[i][j] - pfestimates[i][j]);
                }
                ferrt = ferrt + ferrc;
                serrt = serrt + serrc;
                pferrt = pferrt + pferrc;
            }
            ferrt = ferrt*max_time/(float)steps;
            serrt = serrt*max_time/(float)steps;
            pferrt = pferrt*max_time/(float)steps;
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
            ofstream os("data/smoothing.dat");
            for(unsigned int i=0; i< sestimates.size(); i++)
            {
                for(int j=0; j<ndim; j++)
                    os<<sestimates[i][j]<<" ";
                os<<endl;
            }
            os.close();
            ofstream pf("data/pfilter.dat");
            for(unsigned int i=0; i< pfestimates.size(); i++)
            {
                for(int j=0; j<ndim; j++)
                    pf<<pfestimates[i][j]<<" ";
                pf<<endl;
            }
            pf.close();

            return 0;
        }

};

int main(int argc, char** argv)
{
    int n = 100;
    float max_time = 0.5;
    if(argc > 1)
        n = atoi(argv[1]);
    if(argc > 2)
        max_time = atof(argv[2]);

#if 0
    srand(0);

    float hmmf_time = 0;
    Graph g = Graph(n);
    tic();
    g.connect_nodes();
    hmmf_time += toc();
    cout<<"delta: "<< g.delta<<endl;
    
    g.festimates.clear();
    g.sestimates.clear();
    g.propagate_system(max_time);
    
    tic();
    g.run_pfilter(100);
    cout<<"pf dt: "<< toc()<<endl;

    tic();
    g.run_smoother();
    cout<<"hmmf dt: "<< hmmf_time + toc()<<endl;
    
    g.output_trajectories();
    float ferr, serr, pferr;
    g.calculate_err(ferr, serr, pferr, max_time);
    cout<<n<<" "<<ferr<<" "<<serr<<" "<<pferr<<endl;
#endif

#if 1
    int max_tries = 100;
    for(int n=5000; n<5500; n=n+1000)
    {
        float hmmf_time = 0, pf_time=0;
        srand(0);
        tic();
        Graph g = Graph(n);
        g.connect_nodes();
        hmmf_time += toc();
        float favg=0, savg=0, pfavg=0;
        for(int tries=0; tries<max_tries; tries++)
        {
            g.festimates.clear();
            g.sestimates.clear();
            g.propagate_system(1);
            
            tic();
            g.run_pfilter(100);
            pf_time += toc();
            tic();
            g.run_smoother();
            hmmf_time += toc();
            float ferr, serr, pferr;
            g.calculate_err(ferr, serr, pferr, 1);
            cout<<ferr<<" "<<serr<<endl;
            favg = favg + ferr;
            savg = savg + serr;
            pfavg = pfavg + pferr;
        }
        favg = favg/(float)max_tries;
        savg = savg/(float)max_tries;
        pfavg = pfavg/(float)max_tries;
        cout<<n<<" "<<favg<<" "<<savg<<" "<<pfavg<<" "<< hmmf_time<<" "<<pf_time<<endl;
    }
#endif

    return 0;
}

