#include "common.h"

//#include "systems/vanderpol_parameter.h"
//#include "systems/vanderpol.h"
#include "systems/singleint.h"
//#include "systems/uhlmann.h"
//#include "systems/parameter.h"
//#include "systems/parameter_hard.h"

#define large_num       (10000)

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
double max_norm(double* src, int dim=ndim)
{
    double m = -1;
    for(int i=0; i<dim; i++)
    {
        if(m < src[i])
            m = src[i];
    }
    return m;
}

class Node
{
    public:
        double x[ndim];
        double key[ndim];
        int index;
        double htime;
        
        set<int> ein;
        vector<int> eout;

        Node(int iin, double* xin)
        {
            for(int i=0; i< ndim; i++)
            {
                x[i] = xin[i];
                key[i] = (xin[i]-zmin[i])/(zmax[i] - zmin[i]);
            }
            index = iin;
        }
        ~Node()
        {
        }
};

class Graph
{
    public:
        kdtree* node_tree;
        vector<Node*> nodes;
        int num_vert;
        double delta;
        float* P;
        vector<double> density;

        Graph()
        {
            P = new float[large_num*large_num];
            num_vert = 0;
            delta = 0.001;
            node_tree = kd_create(ndim);
        }
        ~Graph()
        {
            for(int i=0; i< num_vert; i++)
                delete nodes[i];
            kd_free(node_tree);
            delete[] P;
        }
        int clear_graph()
        {
            for(int i=0; i< num_vert; i++)
                delete nodes[i];
            kd_free(node_tree);
            delete[] P;

            P = new float[large_num*large_num];
            num_vert = 0;
            delta = 0.001;
            node_tree = kd_create(ndim);
            return 0;
        }
        int add_node(Node* n1)
        {
            nodes.push_back(n1);
            num_vert++;
            kd_insert(node_tree, n1->key, n1);
            return 0;
        }
        int node_clear_eout(int ii)
        {
            Node* n = nodes[ii];
            for(unsigned int i=0; i< n->eout.size(); i++)
            {
                P[ii*large_num + n->eout[i]] = 0;
                nodes[n->eout[i]]->ein.erase(ii);
            }
            n->eout.clear();
            return 0;
        }
        Node* sample_node(double* mean, double* var)
        {
            double tmp[ndim];
            //multivar_normal(mean, var, tmp, ndim);
            for(int i=0; i<ndim; i++)
            {
                tmp[i] = mean[i] + (RANDF - 0.5)*6*sqrt(var[i]);
            }
            Node* n = new Node(num_vert, tmp);
            return n;
        }
        int connect_node(int nodei, double bowlr)
        {
            node_clear_eout(nodei);
            
            Node* n = nodes[nodei];
            n->htime = holding_time(n->x, bowlr);
            if(n->htime < delta)
            {
                cout<<"htime "<<n->htime<<" < delta "<<delta<<endl;
                exit(0);
            }
            vector<int> neighbors;  
            kdres *res;
            res = kd_nearest_range(node_tree, n->key, bowlr);
            double pos[ndim] = {0};
            while( !kd_res_end(res) )
            {
                Node* n1 = (Node*) kd_res_item(res, pos);
                if(n1 != n)
                    neighbors.push_back(n1->index);
                kd_res_next(res);
            }
            kd_res_free(res);

            vector<float> probs(neighbors.size(), 0);
            double tprob = 0;
            double fdt[ndim] ={0}; drift(n->x, fdt, n->htime);
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
            P[n->index*large_num + n->index] = ps;
            for(unsigned int j=0; j<neighbors.size(); j++)
            {
                Node* n1 = nodes[neighbors[j]];
                probs[j] = probs[j]/tprob;
                P[n->index*large_num + n1->index] = (1-ps)*probs[j];
                n1->ein.insert(n->index);
            }

            neighbors.push_back(n->index);
            n->eout = neighbors;
            n->ein.insert(n->index);
            return 0;
        }
        int reconnect_neighbors(int nodei, double bowlr)
        {
            Node* n = nodes[nodei];
            double pos[ndim] = {0};
            kdres *res;
            res = kd_nearest_range(node_tree, n->key, bowlr);
            while( !kd_res_end(res) )
            {
                Node* n1 = (Node*) kd_res_item(res, pos);
                if(n1 != n)
                    connect_node(n1->index, bowlr);
                kd_res_next(res);
            }
            kd_res_free(res);
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
            double curr_obs[ndim];
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
                    runner_time = runner_time + integration_delta;
                }
                curr_time = curr_time + delta;

                vector<double> state_tmp; state_tmp.assign(curr_state, curr_state+ndim);
                truth.push_back(state_tmp);
                get_obs(curr_state, curr_obs);
                vector<double> obs_tmp; obs_tmp.assign(curr_obs, curr_obs+ndim);
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

        int calci_moment(double* retm, double* retv)
        {
            for(int i=0; i<ndim; i++)
            {
                retm[i] = 0;
                retv[i] = 0;
            }
            double tprob = 0;
            for(unsigned int i=0; i<density.size(); i++)
            {
                for(int j=0; j<ndim; j++)
                {
                    retm[j] = retm[j] + density[i]*nodes[i]->x[j];
                    retv[j] = retv[j] + density[i]*sq(nodes[i]->x[j]);
                }
                tprob = tprob + density[i];
            }
            for(int i=0; i<ndim; i++)
            {
                retm[i] = retm[i]/tprob;
                retv[i] = retv[i]/tprob;
                retv[i] = retv[i] - sq(retm[i]);
            }
            return 0;
        }
        int normalize_density(vector<double>& d)
        {
            double tprob = 0;
            for(unsigned int i=0; i< d.size(); i++)
                tprob = tprob + d[i];
            for(unsigned int i=0; i< d.size(); i++)
                d[i] = d[i]/tprob;
            return 0;
        }
        int approximate_density(double bowlr, int ii)
        {
            double small_bowl = bowlr/pow(3, 1/(float)ndim);
            double mp = 0;
            int mp_num = 0;
            double pos[ndim] = {0};
            kdres *res;
            res = kd_nearest_range(node_tree, nodes[ii]->key, small_bowl);
            while( !kd_res_end(res) )
            {
                Node* n1 = (Node*) kd_res_item(res, pos);
                if( (n1->index != ii) && (density[n1->index] > 1e-30) )
                {
                    mp = mp + density[n1->index];
                    mp_num++;
                }
                kd_res_next(res);
            }
            kd_res_free(res);
            density[ii] = mp/(double)mp_num;
            return 0;
        }
        int update_density(double* next_mean, double bowlr, vector<double>& obs_in)
        {
            double next_mean_key[ndim] ={0};
            for(int i=0; i<ndim; i++)
                next_mean_key[i] = (next_mean[i] - zmin[i])/(zmax[i] - zmin[i]);
            
            vector<int> which_ones;
            double pos[ndim] = {0};
            kdres *res;
            res = kd_nearest_range(node_tree, next_mean_key, bowlr);
            while( !kd_res_end(res) )
            {
                Node* n1 = (Node*) kd_res_item(res, pos);
                which_ones.push_back(n1->index);
                kd_res_next(res);
            }
            kd_res_free(res);

            double obs[ndim] = {0};
            for(int i=0; i<ndim; i++)
                obs[i] = obs_in[i];

            vector<double> dcopy = density;
            memset(&(density[0]), 0, density.size()*sizeof(double));
            for(unsigned int i=0; i< which_ones.size(); i++)
            {
                Node* curr = nodes[which_ones[i]];
                double to_add = 0;
                for(set<int>::iterator j=curr->ein.begin(); j != curr->ein.end(); j++)
                {
                    int fromi = *j;
                    //cout<<curr->index<<" "<<fromi <<" "<<dcopy[fromi]<<" "<<P[fromi*large_num + curr->index]<<" "<< to_add<<endl;
                    to_add = to_add + dcopy[fromi]*(P[fromi*large_num + curr->index]);
                }
                double curr_obs[ndim] = {0};
                get_obs(curr->x, curr_obs, true);
                density[i] = to_add*normal_val(obs, ovar, curr_obs, ndim_obs);
            }
            normalize_density(density);
            return 0;
        }
        
        // returns a new estimate and changes density
        // samples near next estimated density and updates current density
        vector<double> filter_iterate(int per_step, vector<double>& curr_obs)
        {
            double curr_mean[ndim], curr_var[ndim];
            calci_moment(curr_mean, curr_var);
            double next_var[ndim], next_mean[ndim];
            for(int j=0; j<ndim; j++)
                next_var[j] = curr_var[j] + pvar[j]*delta;
            drift(curr_mean, next_mean, delta);
            add(curr_mean, next_mean);
            
            double scaling = (3*sqrt(max_norm(next_var))/max_norm(zdiff));
            double bowlr = 2.1*pow((1+1/(float)ndim), 1/(float)ndim)*scaling *pow(log(per_step)/(double)per_step, 1/(double)ndim);
            
            vector<int> which_ones;
            for(int i=0; i<per_step; i++)
            {
                Node* n2 = sample_node(next_mean, next_var);
                which_ones.push_back(num_vert);
                add_node(n2); density.push_back(0);
                //approximate_density(bowlr, n2->index);
            }
            for(int i=0; i<per_step; i++)
            {
                connect_node(which_ones[i], bowlr);
                //reconnect_neighbors(which_ones[i], bowlr/5);
            }


            update_density(next_mean, scaling, curr_obs);
            normalize_density(density);
            
            calci_moment(next_mean, next_var);
            cout<<next_mean[0]<<" "<<next_mean[1]<<endl;
            vector<double> tmp1; tmp1.assign(next_mean, next_mean+ndim);
            return tmp1;
        }

        vector< vector<double> > festimates;
        int run_filter(int per_step)
        {
            density.clear();
            festimates.clear();
            int steps = observations.size();

            vector<int> which_ones;
            for(int i=0; i<per_step; i++)
            {
                Node* n2 = sample_node(init_state, init_var);
                which_ones.push_back(num_vert);
                add_node(n2); density.push_back(0);
            }
            double scaling = (3*sqrt(max_norm(init_var))/max_norm(zdiff));
            double bowlr = 2.1*pow((1+1/(float)ndim), 1/(float)ndim)*scaling*pow(log(num_vert)/(float)num_vert, 1/(float)ndim);
            for(int i =0; i<per_step; i++)
            {
                connect_node(which_ones[i], bowlr);
                density[i] = normal_val(init_state, init_var, nodes[which_ones[i]]->x, ndim);
            }
            normalize_density(density);
            double next_var[ndim], next_mean[ndim];
            calci_moment(next_mean, next_var);
            cout<<"init mean: "<<next_mean[0]<<" "<<next_mean[1]<<endl;
            
            for(int i=0; i<steps; i++)
                festimates.push_back(filter_iterate(per_step, observations[i]));

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
        int calculate_err(double& ferrt, double& pferrt, double max_time)
        {
            ferrt = 0;
            pferrt = 0;
            int steps = observations.size();
            for(int i=0; i<steps; i++)
            {
                double ferrc = 0;
                double pferrc = 0;
                for(int j=0; j<ndim; j++)
                {
                    ferrc = ferrc + sq(truth[i][j] - festimates[i][j]);
                    pferrc = pferrc + sq(truth[i][j] - pfestimates[i][j]);
                }
                ferrt = ferrt + ferrc;
                pferrt = pferrt + pferrc;
            }
            ferrt = ferrt*max_time/(double)steps;
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
                for(int j=0; j<ndim; j++)
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
    // n per step
    int n = 100;
    double max_time = 0.1;
    if(argc > 1)
        n = atoi(argv[1]);
    if(argc > 2)
        max_time = atof(argv[2]);

#if 1
    srand(time(NULL));
    Graph g = Graph();
    
    g.festimates.clear();
    g.pfestimates.clear();
    g.propagate_system(max_time);
    
    g.run_pfilter(n);
    tic();
    g.run_filter(n);
    cout<<"dt: "<< toc()<<endl;
    g.output_trajectories();
    double ferr, pferr;
    g.calculate_err(ferr, pferr, max_time);
    cout<<n<<" "<<ferr<<" "<<pferr<<endl;
#endif

#if 0
    int max_tries = 50;
    for(int n=1000; n<3000; n=n+500)
    {
        srand(0);
        Graph g = Graph();
        double favg=0, pfavg=0;
        for(int tries=0; tries<max_tries; tries++)
        {
            g.festimates.clear();
            g.pfestimates.clear();
            g.propagate_system(1);
            
            g.run_pfilter(n);
            g.run_filter(n);
            double ferr, pferr;
            g.calculate_err(ferr, pferr, 1);
            //cout<<ferr<<" "<<pferr<<endl;
            favg = favg + ferr;
            pfavg = pfavg + pferr;
        }
        favg = favg/(double)max_tries;
        pfavg = pfavg/(double)max_tries;
        cout<<n<<" "<<favg<<" "<<pfavg<<endl;
    }
#endif

    return 0;
}

