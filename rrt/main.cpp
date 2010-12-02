#include "random.h"
#include "kdtree.h"
#include "rrt.h"

FILE *fpoints;

void read_input(char obs_file[])
{
    FILE *fpgoal, *fpobs, *fpbot, *fpbox;

    fpgoal = fopen("input/goal.txt", "r");
    fpobs = fopen(obs_file, "r");
    fpbot = fopen("input/bot.txt", "r");
    fpbox = fopen("input/box.txt", "r");

    double dat[2*NUM_STATES];
 
    // Get bot
    assert(0 != fscanf(fpbot, "%lf, %lf, %lf", &dat[0], &dat[1], &dat[2]));
    robot = State(dat);
    robot_radius = dat[2];
    //printf("Got robot: %f, %f, %f\n\n", dat[0], dat[1], dat[2]);
    fclose(fpbot);
   
    // Get goal
    assert(0 != fscanf(fpgoal, "%lf, %lf, %lf", &dat[0], &dat[1], &dat[2]));
    goal.state = State(dat);
    goal.size = dat[2];
    //printf("Got goal: %f, %f, %f\n\n", dat[0], dat[1], dat[2]);
    fclose(fpgoal);

    // Get obstacles
    obstree = kd_create(NUM_STATES);
    int c = 0;
    while(fscanf(fpobs, "%lf, %lf, %lf", &dat[0], &dat[1], &dat[2]) == 3)
    {
        double to_put[2] = {dat[0], dat[1]};
        // create config-space by adding bot_rad to obs_rad
        obs_rad[c] = dat[2] + robot_radius;
        kd_insert(obstree, to_put, &obs_rad[c]);
        
        if( obs_rad[c] > MAX_OBS_SIZE)
            MAX_OBS_SIZE = obs_rad[c];

        NUM_OBSTACLES = c;

        //printf("Inserted %f, %f, %f, %f\n", dat[0], dat[1], dat[2], obs_rad[c]);
        c++;
    }
    //for(int i=0; i<11; i++)
    //    printf("%f\n", obs_rad[i]);

    //printf("MAX_OBS_SIZE: %f\n\n", MAX_OBS_SIZE);
    fclose(fpobs);
    
     
    // Get box
    assert(0 != fscanf(fpbox, "%lf, %lf, %lf, %lf", &dat[0], &dat[1], &dat[2], &dat[3]));
    box.center = State(dat);
    for(int i=0; i<NUM_STATES; i++)
        box.size[i] = dat[NUM_STATES+i];
    
    //printf("Got center: %f, %f, %f, %f\n\n", dat[0], dat[1], dat[2], dat[3]);
    fclose(fpbox);
}

// print node-parent pair for plotting
void print_path(list<Node> whichpath)
{
    list<Node>::iterator i;
    for(i=whichpath.begin(); i != whichpath.end(); i++)
    {
        Node curr = *i;
        if(curr.parent != NULL)
        {
            for(unsigned int i=0; i<NUM_STATES; i++)
            {
                fprintf(fpoints, "%f ", curr.state.x[i]);
            }
            fprintf(fpoints, "\n");
            for(unsigned int i=0; i<NUM_STATES; i++)
            {
                fprintf(fpoints, "%f ", (curr.parent)->state.x[i]);
            }
            fprintf(fpoints, "\n");
        }
    }
}

int main(int argc, char* argv[])
{
    init_rand();
        
    int output_path = 0;
    if(argc == 3)
        output_path = atoi(argv[2]);
    
    if(argc >= 2)
        read_input(argv[1]);
    else
        read_input("input/obs1.txt");
    
    fpoints = fopen("points.dat", "w");
    double start = get_msec(), delt;
    /*
    int i = 0;
    while(i < 1)
    {
        start = get_msec();
        double cost = rrt_plan(15000);
        delt = get_msec() - start;
        printf("%d\t%f\t%f\n", i, optpath.back().cgoal, delt);
        i++;
    }
    */
    double cost = rrtstar_plan(10000);
    fprintf(fpoints, "%f \n", cost);
    fprintf(fpoints, "Duration: %.3f \n", get_msec() - start);
    print_path(tree);
    fprintf(fpoints, "optpath\n");
    print_path(optpath);
    fprintf(fpoints, "optpath_cost: %f \n", optpath.back().cgoal);

    fclose(fpoints);
    kd_free (obstree);
    return 0;
}
