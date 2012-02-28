#include "common.h"

//#include "systems/vanderpol_parameter.h"
#include "systems/vanderpol.h"
//#include "systems/singleint.h"
//#include "systems/uhlmann.h"
//#include "systems/parameter.h"
//#include "systems/parameter_hard.h"

class Node
{
    public:
        double x[ndim];
        double key[ndim];
        int index;
        double htime;
        double weight;
        int id;

        vector<int> eout;
        vector<int> ein;
        Node(int iin, double* xin, double win=0)
        {
            weight = win;
            id = 0;
            for(int i=0; i< ndim; i++)
                x[i] = xin[i];
            get_key(xin, key);
            index = iin;
        }
        int get_key(double* xin, double* keyout)
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
        double delta;
        float* P;
        double bowlr;

        Graph(int n, double delta_in)
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
        int add_node_iter(Node* n1)
        {
            nodes.push_back(n1);
            kd_insert(node_tree, n1->key, n1);
            return 0;
        }
        int add_nodes(vector< vector<double> >& particles_prev, vector< vector<double> >& particles, vector<double>& win)
        {
            for(unsigned int i=0; i< particles_prev.size(); i++)
            {
                Node* n1 = new Node(i, &(particles_prev[i][0]), win[i]);
                n1->id = 0;
                add_node_iter(n1);
            }
            for(unsigned int i=0; i< particles.size(); i++)
            {
                Node* n1 = new Node(particles_prev.size()+i, &(particles[i][0]), 0);
                n1->id = 1;
                add_node_iter(n1);
            }
            return 0;
        }
        int connect_nodes()
        {
            double fmean[ndim]={0}, fvar[ndim] ={0};
            calci_moment(fmean, fvar);
            //double scaling = pow(max_norm(fvar), ndim);
            double scaling = 1;
            bowlr = 2.1*pow((1+1/(float)ndim), 1/(float)ndim)*scaling*pow(log(num_vert/2)/(double)num_vert/2, 1/(double)ndim);
            //cout<<"bowlr: "<<bowlr<<endl;
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
            
            double next_state[ndim] ={0}, next_state_key[ndim]={0};
            copy(n->x, next_state);
            integrate_system(next_state, n->htime, true);
            n->get_key(next_state, next_state_key);
            double var[ndim];
            for(int j=0; j<ndim; j++)
                var[j] = pvar[j]*(n->htime);

            vector<int> neighbors;  
            kdres *res;
            res = kd_nearest_range(node_tree, next_state_key, bowlr);
            double pos[ndim] = {0};
            while( !kd_res_end(res) )
            {
                Node* n1 = (Node*) kd_res_item(res, pos);
                if(n1->id == 1)
                    neighbors.push_back(n1->index);
                kd_res_next(res);
            }
            kd_res_free(res);

            vector<float> probs(neighbors.size(), 0);
            double tprob = 0;
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

        int calci_moment(double* retm, double* retv, int id=0)
        {
            int s = 0;
            int e = num_vert/2;
            if(id)
            {
                s = num_vert/2;
                e = num_vert;
            }
            for(int i=0; i<ndim; i++)
            {
                retm[i] = 0;
                retv[i] = 0;
            }
            double tprob = 0;
            for(int i=s; i< e; i++)
            {
                for(int j=0; j<ndim; j++)
                {
                    retm[j] = retm[j] + nodes[i]->weight*nodes[i]->x[j];
                    retv[j] = retv[j] + nodes[i]->weight*sq(nodes[i]->x[j]);
                }
                tprob = tprob + nodes[i]->weight;
            }
            for(int i=0; i<ndim; i++)
            {
                retm[i] = retm[i]/tprob;
                retv[i] = retv[i]/tprob;
                retv[i] = retv[i] - sq(retm[i]);
            }
            return 0;
        }
        int normalize_density()
        {
            for(int i=0; i< num_vert/2; i++)
            {
                if(nodes[i]->id == 0)
                    nodes[i]->weight = 0;
            } 

            double tprob = 0;
            for(int i=0; i< num_vert; i++)
                tprob = tprob + nodes[i]->weight;
            for(int i=0; i< num_vert; i++)
                nodes[i]->weight = nodes[i]->weight/tprob;
            return 0;
        }
        int update_density(vector<double>& obs_in)
        {
            double obs[ndim] = {0};
            for(int i=0; i<ndim; i++)
                obs[i] = obs_in[i];

            for(int i=0; i< num_vert; i++)
            {
                double to_add = 0;
                Node *n1 = nodes[i];
                if(n1->id == 1)
                {
                    for(unsigned int j=0; j< n1->ein.size(); j++)
                        to_add = to_add + nodes[n1->ein[j]]->weight*P[j*num_vert + i];
                    
                    double curr_obs[ndim_obs] = {0};
                    get_obs(n1->x, curr_obs, true);
                    nodes[i]->weight = to_add*normal_val(obs, ovar, curr_obs, ndim_obs);
                }
            }
            normalize_density();
            return 0;
        }
        int print()
        {
            cout<<"-----------------"<<endl;
            for(int i=num_vert/2; i<num_vert; i++)
            {
                Node *n1 = nodes[i];
                cout<<"index: "<< n1->index<<" id: "<<n1->id<<" eins: "<<n1->ein.size()<<endl;
                for(unsigned int j=0; j< n1->ein.size(); j++)
                    cout<< nodes[n1->ein[j]]->index<<"-"<<nodes[n1->ein[j]]->weight<<" ";
                cout<<endl;

            }
            cout<<"transition matrix"<<endl;
            for(int i=0;i<num_vert;i++)
            {
                for(int j=0;j<num_vert;j++)
                    cout<<P[i*num_vert + j]<<" ";
                cout<<endl;
            }
            cout<<"-----------------"<<endl;
            return 0;
        }
};
class HMMF
{
    public:
        double delta;

        HMMF()
        {
            delta = 0.01;
        }

        vector< vector<double> > truth;
        vector< vector<double> > observations;
        int propagate_system(double max_time)
        {
            truth.clear();
            observations.clear();
            double curr_time = 0;
            double curr_state[ndim];
            double curr_obs[ndim_obs];
            copy(init_state_real, curr_state);
            while(curr_time <= max_time)
            {
                integrate_system(curr_state, delta);
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
        vector< vector<double> > festimates;
        int run_filter(int num_particles)
        {
            festimates.clear();
            int steps = observations.size();

            vector< vector<double> > particles(num_particles, vector<double>(ndim,0));
            vector< vector<double> > particles_prev(num_particles, vector<double>(ndim,0));
            vector<double> weights(num_particles, 0);
            for(int i=0; i< num_particles; i++)
            {
                multivar_normal(init_state, init_var, &(particles[i][0]), ndim);
                weights[i] = 1.0/(double)num_particles;
            }
            for(int i=0; i<steps; i++)
            {
                particles_prev = particles;

                for(int j=0; j<num_particles; j++)
                    integrate_system(&(particles[j][0]), delta);

                Graph g(num_particles, delta);
                g.add_nodes(particles_prev, particles, weights);
                g.connect_nodes();
                //g.print();
                g.update_density(observations[i]);

                vector<double> fmean(ndim,0);
                for(int j=0; j<num_particles; j++)
                {
                    weights[j] = g.nodes[num_particles+j]->weight;
                    for(int i=0; i<ndim; i++)
                        fmean[i] = fmean[i] + weights[j]*particles[j][i];
                }
                festimates.push_back(fmean); 
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
                    integrate_system(&(particles[j][0]), delta);

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
            if((festimates.size() == 0) && (pfestimates.size() == 0) )
                return 1;
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

            return 0;
        }
};



int main(int argc, char** argv)
{
    // n per step
    int n = 100;
    double max_time = 1;
    if(argc > 1)
        n = atoi(argv[1]);
    if(argc > 2)
        max_time = atof(argv[2]);

#if 1
    srand(time(NULL));
    HMMF h;
    h.propagate_system(max_time);

    tic();
    h.run_pfilter(n);
    cout<<"dt pf: "<< toc()<<endl;
    tic();
    h.run_filter(n);
    cout<<"dt hmmf: "<< toc()<<endl;
    h.output_trajectories();
    double ferr=0, pferr=0;
    h.calculate_err(ferr, pferr, max_time);
    cout<<n<<" "<<ferr<<" "<<pferr<<endl;
#endif

#if 0
    int s = n;
    int max_tries = 50;
    for(int i=s; i<s+5; i=i+10)
    {
        srand(0);
        double favg=0, pfavg=0, ftavg=0, pftavg=0;
        for(int tries=0; tries<max_tries; tries++)
        {
            HMMF h;
            h.propagate_system(max_time);
            
            tic();
            h.run_pfilter(i);
            pftavg += toc();
            
            tic();
            h.run_filter(i);
            ftavg += toc();

            double ferr=0, pferr=0;
            h.calculate_err(ferr, pferr, max_time);
            favg = favg + ferr;
            pfavg = pfavg + pferr;
            cout<<tries<<" "<<favg/(float)tries<<" "<<pfavg/(float)tries<<endl;
        }
        favg = favg/(double)max_tries;
        pfavg = pfavg/(double)max_tries;
        ftavg = ftavg/(double)max_tries;
        pftavg = pftavg/(double)max_tries;
        cout<<i<<" "<<ftavg<<" "<<pftavg<<" "<<favg<<" "<<pfavg<<endl;
    }
#endif

    return 0;
}

