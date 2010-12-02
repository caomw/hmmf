#include "random.h"
#include "kdtree.h"
#include "rrt.h"

#define INF             (1000)
#define EXTEND_DIST     (.2)
#define GOAL_PROB       (0.01)

typedef struct kdtree kdtree;
typedef struct kdres kdres;

list<Node> tree;
kdtree *obstree;
Goal goal;
State robot;
Box box;
double robot_radius;
double MAX_OBS_SIZE = 0;
double NUM_OBSTACLES = 0;
double obs_rad[20];

double dist(State s1, State s2) 
{
    double dist_sq = 0, diff;
    int dims = NUM_STATES;
    while( --dims >= 0 ) 
    {
        diff = (s1.x[dims] - s2.x[dims]);
        dist_sq += diff*diff;
    }
    return sqrt(dist_sq);
}

void read_input()
{
    FILE *fpgoal, *fpobs, *fpbot, *fpbox;

    fpgoal = fopen("input/goal.txt", "r");
    fpobs = fopen("input/obstacles.txt", "r");
    fpbot = fopen("input/bot.txt", "r");
    fpbox = fopen("input/box.txt", "r");

    double dat[2*NUM_STATES];
    
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
        obs_rad[c] = dat[2];
        kd_insert(obstree, to_put, &obs_rad[c]);
        
        if( obs_rad[c] > MAX_OBS_SIZE)
            MAX_OBS_SIZE = obs_rad[c];

        NUM_OBSTACLES = c;

        //printf("Inserted %f, %f, %f\n", dat[0], dat[1], obs_rad[c]);
        c++;
    }
    //for(int i=0; i<11; i++)
    //    printf("%f\n", obs_rad[i]);

    //printf("MAX_OBS_SIZE: %f\n\n", MAX_OBS_SIZE);
    fclose(fpobs);
    
    // Get bot
    assert(0 != fscanf(fpbot, "%lf, %lf, %lf", &dat[0], &dat[1], &dat[2]));
    robot = State(dat);
    robot_radius = dat[2];

    // add robot-size to obs radius
    for(int i=0; i<NUM_OBSTACLES; i++)
        obs_rad[i] += robot_radius;    

    //printf("Got robot: %f, %f, %f\n\n", dat[0], dat[1], dat[2]);
    fclose(fpbot);
     
    // Get box
    assert(0 != fscanf(fpbox, "%lf, %lf, %lf, %lf", &dat[0], &dat[1], &dat[2], &dat[3]));
    box.center = State(dat);
    for(int i=0; i<NUM_STATES; i++)
        box.size[i] = dat[NUM_STATES+i];
    
    //printf("Got center: %f, %f, %f, %f\n\n", dat[0], dat[1], dat[2], dat[3]);
    fclose(fpbox);
}

bool is_inside_goal(State c)
{
    if( dist(c, goal.state) <= goal.size)
        return true;
    else
        return false;
}

bool is_obstructed(State s)
{
    kdres *presults;
    double t[NUM_STATES];

    for(int i=0; i<NUM_STATES; i++)
        t[i] = s.x[i];

    double pos[NUM_STATES];
    presults = kd_nearest_range(obstree, t, MAX_OBS_SIZE);

    while( !kd_res_end(presults))
    {
        double rad = *((double *)kd_res_item(presults, pos)); 
        State temp = State(pos);
        
        //printf("Found: %f, %f, %.3f\n", pos[0], pos[1], rad);
        //printf("Dist is %f\n", dist(temp, s));

        if( dist(temp, s) <= (rad))
            return true;
        
        kd_res_next( presults );

    }
    kd_res_free(presults);
    return false;
}

State sample_state()
{
    double t = randdouble();
    if(t < GOAL_PROB)
    {
        return goal.state;
    }
    double sample[NUM_STATES];
    for(int i=0; i<NUM_STATES; i++)
    {
        double min = box.center.x[i] - box.size[i]/2;
        double max = box.center.x[i] + box.size[i]/2;
        //printf("min: %f, max: %f\n", min, max);

        double t = randdouble(min, max);
        sample[i] = t;
    }
    State test = State(sample);
    //test.print();
    return test;
}

Node* nearest(State s)
{
    list<Node>::iterator i;
    Node *min_node = NULL;
    double min = INF;

    for(i = tree.begin(); i != tree.end(); i++)
    {
        double t = dist(s, (*i).state);
        if( t < min)
        {
            min = t;
            min_node = &(*i);
        }
    }

    return min_node;
}

/*
Node* nearest(kdtree *node_tree, State s)
{
    kdres *pres;
    pres = kd_nearest(node_tree, s.x);

    if(kd_res_end(pres))
    {
        printf("Couldn't find nearest node, exiting..\n");
        exit(1);
    }
    else
    {
        Node *n = (Node *)kd_res_item_data(pres);
        return n;
    }
}
*/

Node extend(Node *near, State s)
{
    double totdist = dist(near->state, s);
    if(totdist < EXTEND_DIST)
    {
        Node newnode(s, near);
        return newnode;
    }
    else
    {
        double newdist[NUM_STATES];
        for(int i=0; i<NUM_STATES; i++)
        {
            newdist[i] = (near->state).x[i] + (s.x[i] - (near->state).x[i])/totdist*EXTEND_DIST;
        }
        State newstate(newdist);
        Node newnode(newstate, near);
        return newnode;
    }
}


/*
 * Given that p1 and p2 both lie outside the obstacle, the line segment between p1
 * and p2 hits an obstacle iff:
 * 1) The distance between the line joining p1, p2 is less than the obstacle radius.
 * 2) The projection of the center of the obstacle is lies between p1 and p2.
 */
bool does_line_hit(State s1, State s2, double *pos, double rad)
{
    State obs = State(pos);
    if( (dist(s1, obs) <= rad) || (dist(s2, obs) <= rad))
        return true;
    
    double xp, yp, dotProduct, obsDist;
    double xo = pos[0], yo = pos[1];        // Obstacle x, y co-ordinates
    
    double endy = s2.x[1], starty = s1.x[1], startx = s1.x[0], endx = s1.x[0];
    double a = endy - starty;
    double b = startx - endx;
    double c = (endx - startx)*starty - (endy - starty)*startx;


    xp = xo - a*(a*xo + b*yo + c)/(a*a + b*b);
    yp = yo - b*(a*xo + b*yo + c)/(a*a + b*b);

    obsDist = sqrt( (xp - xo)*(xp - xo) + (yp -yo)*(yp -yo) );

    if( obsDist <  rad)
    {
        dotProduct = (startx - xp)*(endx - xp) + (starty - yp)*(endy - yp);
        if(dotProduct < 0)
            return true;
    }

    return false;
}


bool can_join_nodes(Node n1, Node n2)
{
    double length = dist(n1.state, n2.state);
    double center[NUM_STATES];
    kdres *pres; 
    double pos[NUM_STATES];
    
    for(int i=0; i<NUM_STATES; i++)
        center[i] = (n1.state.x[i] + n2.state.x[i])/2;

    pres = kd_nearest_range(obstree, center, (MAX_OBS_SIZE + length/2) );

    while( !kd_res_end(pres))
    {
        double rad = *((double *)kd_res_item(pres, pos)); 
        
        State temp = State(pos);
        
        //printf("Found: %f, %f, %.3f\n", pos[0], pos[1], rad);

        if(does_line_hit(n1.state, n2.state, pos, rad))
            return false;

        kd_res_next( pres);

    }
    kd_res_free(pres);

    return true;
}

void rrt_plan()
{
    //kdtree *node_tree;
    //node_tree = kd_create(NUM_STATES);

    Node start;
    start.state = robot;
    start.parent = NULL;
    tree.push_back(start);

    //kd_insert(node_tree, start.state.x, tree.front().parent);

    bool reached = false;
    while(!reached)
    {
        State t = sample_state();
        Node *near = nearest(t);
        if(near != NULL)
        {
            Node curr = extend(near, t);        // extend near in the direction of t

            if(can_join_nodes(curr, *near))
            {
                //curr.state.print();
                //(curr.parent)->state.print();
                //printf("\n");

                tree.push_back(curr);
                //assert(0 == kd_insert(node_tree, curr.state.x, curr.parent));
                reached = is_inside_goal(curr.state);
            }
        }
    }
    //kd_free(node_tree);
    return;
};

void print_path(list<Node> path)
{
    list<Node>::iterator i;
    for(i=path.begin(); i != path.end(); i++)
    {
        Node curr = *i;
        curr.state.print();
        if(curr.parent != NULL)
        {
            (curr.parent)->state.print();
            //Node *p = curr.parent;
            //p->state.print();
        }
    }
}

int main()
{
    init_rand();
    read_input();
    
    rrt_plan();
    
    print_path(tree);

    kd_free (obstree);

    return 0;
}
