#include "common.h"

#define ndim (2)

double zmin[ndim] = {0, 0.5};
double zmax[ndim] = {1, 2.5};
double init_var[ndim] = {1e-2, 1e-2};
double init_state[ndim] = {0, 1.0};
double pvar[ndim] = {1e-2, 1e-2};
double ovar[ndim] = {1e-3, 1e-3};
double zero[ndim] = {0, 0};

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
double norm(double* s)
{
    double prod=1;
    for(int i=0; i<ndim;i++)
        prod = prod*sq(s[i]);
    return sqrt(prod);
}
int drift(double* s, double *ret, double dt=1.0)
{
    ret[0] = s[1]*dt;
    ret[1] = (-s[0] + 2.0*s[1]*(1-sq(s[0])))*dt;
    return 0;
}
int diffusion(double* s, double* ret, double dt=1.0)
{
    double var[ndim] ={0};
    var[0] = pvar[0]*dt;
    var[1] = pvar[1]*dt;
    multivar_normal(zero, var, ret, ndim);
    return 0;
}
int get_obs(double* s, double* obs)
{
    for(int i=0; i< ndim; i++)
        obs[i] = 0;
    double noise[ndim] = {0};
    multivar_normal(zero, ovar, noise, 2); 
    obs[0] = s[0] + noise[0];
    obs[1] = s[1] + noise[1];
    return 0;
}
double holding_time(double* s, double r)
{
    double h = r*(zmax[1] - zmin[0]);
    double ret[ndim];
    drift(s, ret);
    return h*h/(pvar[0] + h*norm(ret));
}

class Node
{
    public:
        double x[ndim];
        double key[ndim];
        int index;
        double htime;
        Node(int iin, double bowlr)
        {
            for(int i=0; i< ndim; i++)
            {
                double r = RANDF;
                x[i] = zmin[i] + (zmax[i]-zmin[i])*r;
                key[i] = r;
            }
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
        double bowlr;
        double* P;
        double delta;
        Graph(int num_nodes)
        {
            num_vert = num_nodes;
            P = new double[num_nodes*num_nodes];
            memset(P, 0, sizeof(double)*num_nodes*num_nodes);
            bowlr = 2.2*pow(log(num_vert+1)/(num_vert+1.0), 1/(double)ndim);
            node_tree = kd_create(ndim);
            for(int i=0; i< num_nodes; i++)
            {
                Node* n = new Node(i, bowlr);
                nodes.push_back(n);
                kd_insert(node_tree, n->key, n);
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
            delta = 0.8*min_htime();
            for(int i=0; i< num_vert; i++)
            {
                Node* n = nodes[i];
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
                
                vector<double> probs(neighbors.size());
                double tprob = 0;
                double fdt[ndim]; drift(n->x, fdt, n->htime);
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
            copy(init_state, curr_state);
            while(curr_time <= max_time)
            {
                double runner_time = 0;
                while(runner_time < delta)
                {
                    double next_state[ndim];
                    drift(curr_state, next_state, integration_delta);
                    double noise[ndim] = {0};
                    diffusion(curr_state, noise, integration_delta); 
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
                for(int j=0; j<num_vert; j++)
                {
                    toadd = toadd + P[j*num_vert + i]*alphas[index][j];
                }
                alphas[index+1][i] = toadd*normal_val(obs, ovar, nodes[i]->x, 2);
                tprob = tprob + alphas[index+1][i];
            }
            for(int i=0; i<num_vert; i++)
            {
                alphas[index+1][i] = alphas[index+1][i]/tprob;
            }
            return 0;
        }
        // use index+1 betas to get index betas
        int update_betas(vector< vector<double> >& betas, vector<double> obs_in, int index)
        {
            double obs[ndim] = {0};
            for(int i=0; i<ndim; i++)
                obs[i] = obs_in[i];
            
            double tprob = 0;
            for(int i=0; i<num_vert; i++)
            {
                double toadd = 0;
                for(int j=0; j<num_vert; j++)
                {
                    toadd = toadd + P[i*num_vert + j]*betas[index+1][j]*normal_val(obs, ovar, nodes[j]->x, 2);
                }
                betas[index][i] = toadd;
                tprob = tprob + betas[index][i];
            }
            for(int i=0; i<num_vert; i++)
            {
                betas[index][i] = betas[index][i]/tprob;
            }
            return 0;
        }

        vector< vector<double> > festimates;
        vector< vector<double> > sestimates;
        int run_smoother(double max_time)
        {
            festimates.clear();
            sestimates.clear();
            propagate_system(max_time);
            
            int steps = observations.size();
            vector< vector<double> > alphas(steps+1, vector<double>(num_vert));
            vector< vector<double> > betas(steps+1, vector<double>(num_vert));
            vector< vector<double> > sdensity(steps, vector<double>(num_vert));
            int li = steps;

            double tprob = 0;
            for(int i=0; i<num_vert; i++)
            {
                alphas[0][i] = normal_val(init_state, init_var, nodes[i]->x, ndim);
                tprob = tprob + alphas[0][i];
                betas[li][i] = 1.0;
            }
            for(int i=0; i<num_vert; i++)
                alphas[0][i] = alphas[0][i]/tprob;
            
            for(int i=0; i<steps; i++)
                update_alphas(alphas, observations[i], i);
            for(int i=steps-1; i >= 0; i--)
                update_betas(betas, observations[i], i);
            for(int i=0; i<steps; i++)
            {
                for(int j=0; j<num_vert; j++)
                {
                    sdensity[i][j] = alphas[i+1][j]*betas[i][j];
                }
            }

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
                        
            // write output
            ofstream of("filter.dat");
            for(int i=0; i< steps; i++)
            {
                for(int j=0; j<ndim; j++)
                    of<<festimates[i][j]<<" ";
                of<<endl;
            }
            of.close();
            ofstream os("smoothing.dat");
            for(int i=0; i< steps; i++)
            {
                for(int j=0; j<ndim; j++)
                    os<<sestimates[i][j]<<" ";
                os<<endl;
            }
            os.close();
            return 0;
        }
        
};

int main(int argc, char** argv)
{
    int n = 100;
    double max_time = 0.5;
    if(argc > 1)
        n = atoi(argv[1]);
    if(argc > 2)
        max_time = atof(argv[2]);

    tic();
    Graph g = Graph(n);
    g.connect_nodes();
    cout<<"delta: "<< g.delta<<endl;
    g.run_smoother(max_time);
    
    cout<<"dt: "<< toc()<<endl;

    return 0;
}
