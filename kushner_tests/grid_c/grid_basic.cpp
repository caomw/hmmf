#include "common.h"
using namespace std;

double init_state = 0.7;
double init_var = 0.001;
double zmax = 1.0;
double zmin = 0.2;
double pvar = 0.001;
double delta = 0.001;

double drift(double s)
{
    return -10*s;
}

int find_closest(vector<double> &array, double var)
{
    int index = 0;
    double dist = 100;
    for(unsigned int i=0; i< array.size(); i++)
    {
        if( fabs(array[i] - var) <dist)
        {
            dist = fabs(array[i] - var);
            index = i;
        }
    }
    return index;
}

int main(int argc, char** argv)
{
    int num_states = 1000;
    if(argc > 1)
        num_states = atoi(argv[1]);

    for(int num = 1000; num<2000; num=num+10)
    {
        num_states = num;
        double h = (zmax - zmin)/(double)num_states;
        vector<double>states(num_states);
        vector< vector<double> >edge_probs(num_states, vector<double>(2));
        vector<double>htimes(num_states);
        
        for(int i=0; i<num_states; i++)
        {
            states[i] = zmin + i*h;
            //states[i] = zmin + RANDF*(zmax-zmin);
        }
        sort(states.begin(), states.end());
        for(int i=0; i<num_states; i++)
        {
            double c = pvar + h*fabs(drift(states[i]));
            htimes[i] = h*h/c;
        }
        vector<double>::iterator pos = min_element(htimes.begin(), htimes.end());
        delta = *pos;
        
        for(int i=0; i<num_states; i++)
        {
            double ps = 1*(1-delta/htimes[i]);
            double c = pvar + h*fabs(drift(states[i]));
            if ( (i!= 0) &&(i!=(num_states-1)))
            {
                
                /*
                // linear program
                double dz1 = states[i+1] - states[i];
                double dz2 = states[i-1] - states[i];
                h = (fabs(dz1)+fabs(dz2))/2.0;
                double a1 = 0; double a2 = drift(states[i])*htimes[i]/dz2;
                double b2 = pvar*htimes[i]/(dz2*dz2 - dz2*dz1);
                double b1 = -b2*dz2/dz1;
                edge_probs[i][0] = (a2+b2)/(a1+b1+a2+b2);
                edge_probs[i][1] = (a1+b1)/(a1+b1+a2+b2);
                */
                
                /*
                // Gaussian
                double xpfdt = states[i] + drift(states[i])*htimes[i];
                double newvar = pvar*htimes[i];
                double totprob = 0;
                edge_probs[i][0] = normal_val(&xpfdt, &newvar, &(states[i-1]),1);
                edge_probs[i][1] = normal_val(&xpfdt, &newvar, &(states[i+1]),1);
                totprob = edge_probs[i][0] + edge_probs[i][1];
                edge_probs[i][0] = edge_probs[i][0]/totprob;
                edge_probs[i][1] = edge_probs[i][1]/totprob;
                */
                
                // fixed formulae
                edge_probs[i][0] = (1-ps)*(pvar/2.0 + h*fabs(drift(states[i])))/c;
                edge_probs[i][1] = (1-ps)*pvar/2.0/c;
            }
            else if(i == 0)
            {
                edge_probs[i][0] = 0;
                edge_probs[i][1] = 1-ps;
            }
            else if(i== num_states-1)
            {
                edge_probs[i][0] = 1-ps;
                edge_probs[i][1] = 0;
            }
        }
        //cout<<"i already laid the grid"<<endl;
#if 1
        int num_traj = 50000; //50000;
        double max_time = 0.05;
        double end_state = init_state*exp(drift(1.0)/1.0*max_time);
        double em1 = 0;
        for(int i=0; i< num_traj; i++)
        {
            double ctime = 0;
            double scurr = init_state;
            //multivar_normal(&init_state, &pvar, &scurr, 1);
            if(scurr > zmax)
                scurr = zmax-0.001;
            else if(scurr < zmin)
                scurr = zmin+0.001;

            //int si = find_closest(states, scurr);
            int si = (int)((scurr - zmin)/h);
            //cout<<"-----------"<<endl;
            //cout<<"si: "<< si << endl;
            double cprob = 1;
            //double cprob = normal_val(&init_state, &init_var, &(states[si]), 1);
            while (ctime < max_time)
            {
                double cointoss = RANDF - 1*(1.0 - delta/htimes[si]);
                //cout<<"si: "<< si << " cointoss: "<< cointoss<<endl;
                //ctime += htimes[si];
                ctime += delta;
                if(cointoss < 0)
                {   
                    cprob = cprob*(1.0 - delta/htimes[i]);
                    si = si;
                }
                else if(cointoss < edge_probs[si][0])
                {
                    cprob = cprob*edge_probs[si][0];
                    si = si-1;
                }
                else
                {
                    cprob = cprob*edge_probs[si][0];
                    si = si+1;
                }
            }
            //cout<<si<<" ";
            em1 = em1 + states[si];
        }
        em1 = em1/float(num_traj);
        em1 = fabs(end_state - em1);
        cout<<num_states<<" "<<em1<<endl;
#endif
    }
    return 0;
}

