#include "common.h"

#include "systems/rotation.h"

int add(double* src, double* dest)
{
    // puts result in dest
    for(int i=0; i< ndim; i++)
        dest[i] = dest[i] + src[i];
    return 0;
}
int copy(double* src, double* dest)
{
    for(int i=0; i< ndim; i++)
        dest[i] = src[i];
    return 0;
}
int print(double* src, int dim=ndim)
{
    for(int i=0; i<dim; i++)
        cout<<src[i]<<" ";
    cout<<endl;
    return 0;
}

class Node
{
    public:
        double x[ndim];
        double key[ndim];
        int index;
        double htime;
        Node(int iin, double *state)
        {
            for(int i=0; i< ndim; i++)
            {
                x[i] = state[i];
                key[i] = (state[i] +1)/2.0;
            }
            key[0] = state[0];
            //print(key);
            index = iin;
        }
};

class Graph
{
    public:
        kdtree* node_tree;
        vector<Node*> nodes;
        int num_vert;
        double bowlr;
        double* P;
        vector< vector<int> > neighbor_list;
        double delta;
        Graph()
        {
            num_vert = 0;
            node_tree = kd_create(ndim);
            ifstream qin("systems/quat.dat");
            while(1)
            {
                double quat[ndim] = {0};
                qin >> quat[1] >> quat[2] >> quat[3] >> quat[0];
                if(qin.eof())
                    break;
                //print(quat);

                Node* n = new Node(num_vert, quat);
                nodes.push_back(n);
                kd_insert(node_tree, n->key, n);
                num_vert++;
            }
            qin.close();
             
            bowlr = pow(2.1*(1+1/(double)ndim), 1/(double)ndim)*pow(log(num_vert+1)/(num_vert+1.0), 1/(double)ndim);
            cout<<"num_vert: "<< num_vert<<" bowlr: "<< bowlr << endl;
            P = new double[num_vert*num_vert];
            memset(P, 0, sizeof(double)*num_vert*num_vert);
            delta = -1;
        }
        ~Graph()
        {
            for(int i=0; i< num_vert; i++)
                delete nodes[i];
            kd_free(node_tree);
            delete[] P;
        }
        int calci_htime()
        {
            for(int i=0; i<num_vert; i++)
                nodes[i]->htime = holding_time(nodes[i]->x, bowlr);
            return 0;
        }
        double min_htime()
        {
            double min_ht = 1000;
            for(int i=0; i< num_vert; i++)
            {
                if(nodes[i]->htime < min_ht)
                    min_ht = nodes[i]->htime;
            }
            return min_ht;
        }
        int connect_nodes()
        {
            calci_htime();
            double pos[ndim] = {0};
            delta = min(0.01, 0.99*min_htime());
            for(int i=0; i< num_vert; i++)
            {
                Node* n = nodes[i];
                vector<int> neighbors;        
                kdres *res;
                res = kd_nearest_range(node_tree, n->key, bowlr);
                //cout<<"nodes in bowl: "<< kd_res_size(res) << endl;
                while( !kd_res_end(res) )
                {
                    Node* n1 = (Node*) kd_res_item(res, pos);
                    if(n1 != n)
                        neighbors.push_back(n1->index);
                    kd_res_next(res);
                }
                kd_res_free(res);
                
                neighbor_list.push_back(neighbors);
                vector<double> probs(neighbors.size());
                double tprob = 0;
                double fdt[ndim] ={0}; 
                drift(n->x, fdt, n->htime);
                double var[ndim];
                for(int j=0; j<ndim; j++)
                {
                    fdt[j] = fdt[j] + n->x[j];
                    var[j] = pvar[j]*(n->htime);
                }
                for(unsigned int j=0; j<neighbors.size(); j++)
                {
                    Node* n1 = nodes[neighbors[j]];
                    probs[j] = normal_val(fdt, var, n1->x, ndim);
                    tprob = tprob + probs[j];
                }
                double ps = 1 - delta/n->htime;
                P[n->index*num_vert + n->index] = ps;
                for(unsigned int j=0; j<neighbors.size(); j++)
                {
                    probs[j] = probs[j]/tprob;
                    P[n->index*num_vert + neighbors[j]] = (1-ps)*probs[j];
                }
            }
#if 0
            for(int i=0; i<num_vert; i++)
            {
                for(int j=0; j<num_vert; j++)
                    cout<<P[i*num_vert + j]<<" ";
                cout<<endl;
            }
            getchar();
#endif
            
            return 0;
        }

        vector< vector<double> > truth;
        vector< vector<double> > observations;
        int propagate_system(double max_time)
        {
            truth.clear();
            observations.clear();
            double curr_time = 0;
            double integration_delta = min(1e-3, delta/2.0);
            double curr_state[ndim];
            double curr_obs[ndim_obs];
            copy(init_state_real, curr_state);
            while(curr_time <= max_time)
            {
                double runner_time = 0;
                while(runner_time < delta)
                {
                    double next_state[ndim] ={0};
                    drift(curr_state, next_state, integration_delta, true);
                    double noise[ndim] = {0};
                    diffusion(curr_state, noise, integration_delta, true); 
                    add(noise, next_state);
                    add(next_state, curr_state);
                    if(curr_state[0] > 1.0)
                        curr_state[0] = 1.0 - (curr_state[0] - 1.0);
                    runner_time = runner_time + integration_delta;
                }
                curr_time = curr_time + delta;

                vector<double> state_tmp; state_tmp.assign(curr_state, curr_state+ndim);
                truth.push_back(state_tmp);
                get_obs(curr_state, curr_obs);
                vector<double> obs_tmp; obs_tmp.assign(curr_obs, curr_obs+ndim_obs);
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

        int calci_mean(vector< vector<double> >& density, int index, double* ret)
        {
            for(int i=0; i<ndim; i++)
                ret[i] = 0;

            double tprob = 0;
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
        int update_alphas(vector< vector<double> >& alphas, vector<double> obs_in, int index)
        {
            double obs[ndim] = {0};
            for(int i=0; i<ndim; i++)
                obs[i] = obs_in[i];

            double tprob = 0;
            for(int i=0; i<num_vert; i++)
            {
                double toadd = 0;
                for(unsigned int j=0; j< neighbor_list[i].size(); j++)
                {
                    int neighbor_j = neighbor_list[i][j];
                    if(alphas[index][neighbor_j] > -1e-40)
                        toadd = toadd + P[neighbor_j*num_vert + i]*alphas[index][neighbor_j];
                }
                double curr_obs[ndim_obs] = {0};
                get_obs(nodes[i]->x, curr_obs, true);
                alphas[index+1][i] = toadd*normal_val(obs, ovar, curr_obs, ndim_obs);
                tprob = tprob + alphas[index+1][i];
            }
            for(int i=0; i<num_vert; i++)
            {
                alphas[index+1][i] = alphas[index+1][i]/tprob;
            }
            return 0;
        }
        // use index betas to get index-1 betas
        int update_betas(vector< vector<double> >& betas, vector<double> obs_in, int index)
        {
            double obs[ndim] = {0};
            for(int i=0; i<ndim; i++)
                obs[i] = obs_in[i];

            double tprob = 0;
            for(int i=0; i<num_vert; i++)
            {
                double toadd = 0;
                double curr_obs[ndim_obs] = {0};
                for(unsigned int j=0; j< neighbor_list[i].size(); j++)
                {
                    int neighbor_j = neighbor_list[i][j];
                    get_obs(nodes[neighbor_j]->x, curr_obs, true);
                    if(betas[index][neighbor_j] > 1e-40)
                        toadd = toadd + P[i*num_vert + neighbor_j]*betas[index][neighbor_j]*normal_val(obs, ovar, curr_obs, ndim_obs);
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

        vector< vector<double> > festimates;
        vector< vector<double> > sestimates;
        int run_smoother()
        {
            int steps = observations.size();
            vector< vector<double> > alphas(steps+1, vector<double>(num_vert,0));
            vector< vector<double> > betas(steps, vector<double>(num_vert,0));
            vector< vector<double> > sdensity(steps, vector<double>(num_vert,0));

            double tprob = 0;
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
            /*
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
            */
            for(int i=0; i< steps; i++)
            {
                double cmean[ndim];

                calci_mean(alphas, i+1, cmean);
                vector<double> tmp1; tmp1.assign(cmean, cmean+ndim);
                festimates.push_back(tmp1);

                calci_mean(sdensity, i, cmean);
                vector<double> tmp2; tmp2.assign(cmean, cmean+ndim);
                sestimates.push_back(tmp2);
            }

            return 0;
        }

        vector< vector<double> > pfestimates;
        int pfilter_resample()
        {
            return 0;
        }
        int run_pfilter(int num_particles)
        {
            pfestimates.clear();
            int steps = observations.size();
            double integration_delta = min(1e-3, delta/2.0);

            vector< vector<double> > particles(num_particles, vector<double>(ndim,0));
            vector<double> weights(num_particles, 0);
            for(int i=0; i< num_particles; i++)
            {
                multivar_normal(init_state, init_var, &(particles[i][0]), ndim);
                weights[i] = 1.0/(double)num_particles;
            }
            for(int i=0; i<steps; i++)
            {
                for(int j=0; j<num_particles; j++)
                {
                    double runner_delta = 0;
                    while(runner_delta < delta)
                    {
                        double next_state_delta[ndim] ={0};
                        drift(&(particles[j][0]), next_state_delta, integration_delta, false);
                        double noise[ndim] = {0};
                        diffusion(&(particles[j][0]), noise, integration_delta, false);
                        add(noise, next_state_delta);
                        add(next_state_delta, &(particles[j][0]));
                        runner_delta += integration_delta;
                    }
                }
                double tot_prob = 0;
                vector<double> pfmean(ndim,0);
                for(int j=0; j<num_particles; j++)
                {
                    double particle_obs[ndim_obs] ={0};
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
        int calculate_err(double& ferrt, double& serrt, double& pferrt, double max_time)
        {
            ferrt = 0;
            serrt = 0;
            pferrt = 0;
            int steps = observations.size();
            for(int i=0; i<steps; i++)
            {
                double ferrc = 0;
                double serrc = 0;
                double pferrc = 0;
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
            ferrt = ferrt*max_time/(double)steps;
            serrt = serrt*max_time/(double)steps;
            pferrt = pferrt*max_time/(double)steps;
            return 0;
        }
        int output_trajectories()
        {
            ofstream ot("truth.dat");
            for(unsigned int i=0; i< truth.size(); i++)
            {
                for(int j=0; j<ndim; j++)
                    ot<<truth[i][j]<<" ";
                ot<<endl;
            }
            ot.close();
            ofstream ob("observations.dat");
            for(unsigned int i=0; i< observations.size(); i++)
            {
                for(int j=0; j<ndim_obs; j++)
                    ob<<observations[i][j]<<" ";
                ob<<endl;
            }
            ob.close();
            ofstream of("filter.dat");
            for(unsigned int i=0; i< festimates.size(); i++)
            {
                for(int j=0; j<ndim; j++)
                    of<<festimates[i][j]<<" ";
                of<<endl;
            }
            of.close();
            ofstream os("smoothing.dat");
            for(unsigned int i=0; i< sestimates.size(); i++)
            {
                for(int j=0; j<ndim; j++)
                    os<<sestimates[i][j]<<" ";
                os<<endl;
            }
            os.close();
            ofstream pf("pfilter.dat");
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
    double max_time = 0.5;
    if(argc > 1)
        max_time = atof(argv[1]);

#if 1
    srand(0);

    double hmmf_time = 0;
    Graph g = Graph();
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
    double ferr, serr, pferr;
    g.calculate_err(ferr, serr, pferr, max_time);
    cout<<ferr<<" "<<serr<<" "<<pferr<<endl;
#endif

#if 0
    int max_tries = 1;
    for(int n=5000; n<10000; n=n+100)
    {
        double hmmf_time = 0, pf_time=0;
        srand(0);
        tic();
        Graph g = Graph(n);
        g.connect_nodes();
        hmmf_time += toc();
        double favg=0, savg=0, pfavg=0;
        for(int tries=0; tries<max_tries; tries++)
        {
            g.festimates.clear();
            g.sestimates.clear();
            g.propagate_system(1);
            
            tic();
            g.run_pfilter(n);
            pf_time += toc();
            tic();
            g.run_smoother();
            hmmf_time += toc();
            double ferr, serr, pferr;
            g.calculate_err(ferr, serr, pferr, 1);
            //cout<<ferr<<" "<<serr<<endl;
            favg = favg + ferr;
            savg = savg + serr;
            pfavg = pfavg + pferr;
        }
        favg = favg/(double)max_tries;
        savg = savg/(double)max_tries;
        pfavg = pfavg/(double)max_tries;
        cout<<n<<" "<<favg<<" "<<savg<<" "<<pfavg<<" "<< hmmf_time<<" "<<pf_time<<endl;
    }
#endif

    return 0;
}

