#include "random.h"
#include "kdtree.h"
#include "rrt.h"

struct kdtree *obstree;
State goal;
State robot;
Box box;
double robot_radius;

static double dist_sq( double *a1, double *a2, int dims ) 
{
    double dist_sq = 0, diff;
    while( --dims >= 0 ) {
        diff = (a1[dims] - a2[dims]);
        dist_sq += diff*diff;
    }
    return dist_sq;
}

void read_input()
{
    FILE *fpgoal, *fpobs, *fpbot, *fpbox;

    fpgoal = fopen("input/goal.txt", "r");
    fpobs = fopen("input/obstacles.txt", "r");
    fpbot = fopen("input/bot.txt", "r");
    fpbox = fopen("input/box.txt", "r");

    float dat[2*NUM_STATES];
    
    // Get goal
    assert(0 != fscanf(fpgoal, "%f, %f, %f", &dat[0], &dat[1], &dat[2]));
    vector<double> t (dat, dat+NUM_STATES);
    goal = State(t);
    printf("Got goal: %f, %f, %f\n", dat[0], dat[1], dat[2]);
    cout<<endl;
    fclose(fpgoal);

    // Get obstacles
    obstree = kd_create(2);
    while(fscanf(fpobs, "%f, %f, %f", &dat[0], &dat[1], &dat[2]) == 3)
    {
        double to_put[2] = {dat[0], dat[1]};
        kd_insert(obstree, to_put, &dat[2]);
        printf("Inserted %f, %f, %f\n", dat[0], dat[1], dat[2]);
    }
    cout<<endl;
    fclose(fpobs);
    
    // Get bot
    assert(0 != fscanf(fpbot, "%f, %f, %f", &dat[0], &dat[1], &dat[2]));
    t = vector<double>(dat, dat+NUM_STATES);
    robot = State(t);
    robot_radius = dat[2];
    printf("Got robot: %f, %f, %f\n", dat[0], dat[1], dat[2]);
    cout<<endl;
    fclose(fpbot);
     
    // Get box
    assert(0 != fscanf(fpbox, "%f, %f, %f, %f", &dat[0], &dat[1], &dat[2], &dat[3]));
    t = vector<double>(dat, dat+NUM_STATES);
    box.center = State(t);
    box.width = dat[NUM_STATES];
    box.height = dat[NUM_STATES+1];
    printf("Got center: %f, %f, %f, %f\n", dat[0], dat[1], dat[2], dat[3]);
    fclose(fpbox);
}

int main()
{
    init_rand();
    read_input();
    
    

    kd_free (obstree);
    return 0;
}
