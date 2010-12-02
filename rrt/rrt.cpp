#include "random.h"
#include "kdtree.h"
#include "rrt.h"

#define INF                 (1000)
#define EXTEND_DIST         (1.0)
#define MAX_RRT_BOWL_SIZE   (10.0)
#define GAMMA               (30)
#define GOAL_PROB           (0.05)
#define USE_KDTREE          (1)
#define BRANCH_N_BOUND      (1)

list<Node> tree;        // stores the tree
list<Node> path;        // stores path that reaches goal
list<Node> optpath;     // stores softened path

kdtree *obstree;
kdtree *node_tree;

Goal goal;
State robot;
Box box;
double robot_radius;
double MAX_OBS_SIZE = 0;
double NUM_OBSTACLES = 0;
double obs_rad[100];
double RRT_BOWL_SIZE = MAX_RRT_BOWL_SIZE;
int num_nodes = 0;

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

Node* nearest_kdtree(kdtree *node_tree, State s)
{
    kdres *pres;
    //printf("query nearest ...."); getchar();
    //cout<<endl;
    pres = kd_nearest(node_tree, s.x);

    if(kd_res_end(pres))
    {
        printf("Couldn't find nearest node, exiting..\n");
        exit(1);
    }
    else
    {
        //printf("get item from kdtree ...."); getchar();
        //cout<<endl;

        Node *n = (Node *)kd_res_item_data(pres);
        return n;
    }
}

Node extend_rrt(Node *near, State s)
{
    double totdist = dist(near->state, s);
    Node newnode;
    if(totdist < EXTEND_DIST)
    {
        newnode = Node(s, near);
    }
    else
    {
        double newdist[NUM_STATES];
        for(int i=0; i<NUM_STATES; i++)
        {
            newdist[i] = (near->state).x[i] + (s.x[i] - (near->state).x[i])/totdist*EXTEND_DIST;
        }
        State newstate(newdist);
        newnode = Node(newstate, near);
    }

    return newnode;
}

int extend_rrtstar(kdtree *node_tree, Node *near, State s, Node &returned_node, double curr_min_cost)
{
    Node newnode = extend_rrt(near, s);

    if( can_join_nodes(newnode, *near) )
    {
        kdres *pres;
        RRT_BOWL_SIZE = GAMMA*sqrt(log((double)num_nodes)/(double)num_nodes);
        if(RRT_BOWL_SIZE > MAX_RRT_BOWL_SIZE)
            RRT_BOWL_SIZE = MAX_RRT_BOWL_SIZE;
        pres = kd_nearest_range(node_tree, newnode.state.x, RRT_BOWL_SIZE);

        double min_cost_neighbor = 1000;

        kdres *pres_temp = kd_nearest(node_tree, s.x);
        Node *near_node = (Node *)kd_res_item_data(pres_temp);
        kd_res_free(pres_temp);

        Node *temp_node;
        while( !kd_res_end(pres) )
        {
            temp_node = (Node *)kd_res_item_data(pres); 

            if (can_join_nodes(newnode, *temp_node) )
            {
                double t = temp_node->csrc + dist(temp_node->state, newnode.state);
                if( t < min_cost_neighbor)
                {
                    near_node = temp_node;
                    min_cost_neighbor = t;
                }
            }
            kd_res_next( pres );
        }

        // near_node now has the nearest code node
        // make it the parent of newnode and write the costs
        // insert in tree
        newnode.parent = near_node;
        double t = dist(newnode.state, near_node->state);
        newnode.csrc = near_node->csrc + t;
        newnode.cparent = t;
        newnode.cgoal = dist(newnode.state, goal.state);

        tree.push_back(newnode);

        // try to rewire parents of all the nodes in the bowl 
        while( !kd_res_end(pres) )
        {
            temp_node = (Node *)kd_res_item_data(pres); 

            if (can_join_nodes(newnode, *temp_node) )
            {
                double pardist = dist(temp_node->state, newnode.state);
                double t = newnode.csrc + pardist;
                if( t < temp_node->csrc)
                {
                    // you can reach temp_node quicker through newnode than its parent, rewire it
                    temp_node->parent = &(tree.back());
                    temp_node->csrc = t;
                    temp_node->cparent = pardist;
                }
            }
            kd_res_next( pres );
        }

        kd_res_free(pres);

        returned_node = newnode;
        return 0;
    }
    return 1;
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
    {
        //printf("Point inside obstacle\n");
        return true;
    }
    double xp, yp, dotProduct, obsDist;
    double xo = pos[0], yo = pos[1];        // Obstacle x, y co-ordinates
    
    double endy = s2.x[1], starty = s1.x[1], startx = s1.x[0], endx = s2.x[0];
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

void process_tree_rrt(Node goal_node)
{
    path.clear();
    optpath.clear();

    // calculate the smoothened path
    // curr is the goal state now -> travel back along it's parents
    Node *n = new Node(goal_node.state, goal_node.parent);

    path.push_back(*n);
    while( n->parent != NULL)
    {
        Node *n1 = n->parent;
        path.push_back(*n1);
        n = n1;
    }
    
    list<Node> temp(path);
    n = &(temp.front());
    while( n->state != robot)
    {
        Node *runner = n->parent;
        Node *runner_child = n;
        
        while( can_join_nodes(*n, *runner) )
        {
            if(runner != NULL)
            {
                /*
                printf("n ");
                n->state.print();
                printf("r ");
                runner->state.print();
                printf("can_join_flag: %d\n", can_join_flag);
                getchar();
                */    
                runner_child = runner;
                runner = runner->parent;
            }
            if(runner == NULL)
                break;
        }
        // update cost and wire it to new parent
        n->parent = runner_child;
        double t = dist( n->state, runner_child->state);
        runner_child->cgoal = n->cgoal + t;
        n->cparent = t;

        optpath.push_back(*n);
        n = n->parent;
    }
    optpath.push_back(*n);
    
    return;
}

void remove_bad_nodes(double min_cost)
{
    list<Node>::iterator s;
    list<Node>::iterator r;

    int num_deleted = 0; 
    for(s = tree.begin(); s != tree.end(); s++)
    {
        if( ((*s).csrc + (*s).cgoal) > min_cost)
        {
            for(r = tree.begin(); r != tree.end(); r++)
            {
                if( r != s )
                {
                    if( r->parent == (&(*s)) )
                    {
                        //r->state.print();
                        //getchar();
                        r = tree.erase(r);
                        num_deleted++;
                    }
                }
            }
            //s = tree.erase(s);
            //s->state.print();
            //getchar();
            num_deleted++;
        }
    }
    
    //remove all elements of old tree
    kd_clear(node_tree);
    
    //build tree again
    num_nodes = 0;
    for(s = tree.begin(); s!= tree.end(); s++)
    {
        assert( 0 == kd_insert(node_tree, (*s).state.x, &(*s)) );
        num_nodes++;
    }
    printf("Removed: %d, Left: %d\n", num_deleted, num_nodes);
    
    return;
}

double rrt_plan(unsigned int num_iter)
{
    double min_cost = 1000;
    tree.clear();
    
    node_tree = kd_create(NUM_STATES);

    Node start(robot, NULL);
    tree.push_back(start);

    kd_insert(node_tree, start.state.x, &(tree.front()) );

    Node curr;
    bool reached = false;
    unsigned int iter = 0;
    bool run_flag = true;

    while(run_flag)
    {
        State t = sample_state();
        if(! is_obstructed(t))
        {
#if USE_KDTREE
            Node *near = nearest_kdtree(node_tree, t);
#else
            Node *near = nearest(t);
#endif
            if(near != NULL)
            {
                //printf("near is not NULL\n");
                curr = extend_rrt(near, t);        // extend near in the direction of t

                double t = dist(curr.state, near->state);
                curr.csrc = near->csrc + t;
                curr.cparent = t;
                curr.cgoal = dist(curr.state, goal.state);

                if(can_join_nodes(curr, *near))
                {
                    //curr.state.print();
                    //(curr.parent)->state.print();
                    //printf("\n");

                    tree.push_back(curr);
                    assert(0 == kd_insert(node_tree, curr.state.x, &(tree.back()) ));
                    
                    if( (reached == 0) && (is_inside_goal(curr.state)) )
                    {
                        reached = true;
                        printf("Reached at iter: %d\n", iter);
                        if(min_cost > (curr.csrc + curr.cgoal) )
                            min_cost = (curr.csrc + curr.cgoal);
                    }
                    iter++;
                }
            }
        }
/*
#if BRANCH_N_BOUND
        if( (iter % 1000 == 0) && (reached) )
        {
            double start = get_msec();
            printf("Going into remove_bad\n");
            remove_bad_nodes(min_cost);
            printf("remove_bad took: %.3f [ms]\n", get_msec() - start);
        }
        run_flag = ( (reached == false) || (iter < num_iter) );          // run till end of steps for BRANCH AND BOUND
#else
*/
        run_flag = (!reached);
//#endif
    }
    kd_free(node_tree);

    process_tree_rrt(curr);
    
    return min_cost;
};

void process_tree_rrtstar(Node *goal_node)
{
    path.clear();
    optpath.clear();

    Node *n = goal_node;
    while( n->parent != NULL)
    {
        optpath.push_back(*n);
        n = n->parent;
    }
    return;
}


double rrtstar_plan(unsigned int num_iter)
{
    tree.clear();

    Node start(robot, NULL);
    tree.push_back(start);
    num_nodes++;

    node_tree = kd_create(NUM_STATES);
    kd_insert(node_tree, start.state.x, &(tree.front()) );

    Node curr;
    bool reached = false;
    Node *node_that_reached = NULL;
    double min_cost = 500;
    
    unsigned int iter = 0;
    while( (iter< num_iter) || (!reached) )
    {
        State t = sample_state();
        
        if(! is_obstructed(t))
        {
            Node *near = nearest_kdtree(node_tree, t);

            if(near != NULL)
            {
                //printf("near is not NULL\n");
                if( extend_rrtstar(node_tree, near, t, curr, min_cost) == 0)          // extend near in the direction of t
                {
                    //curr.state.print();
                    //(curr.parent)->state.print();
                    //printf("\n");
                    //getchar();

                    assert(0 == kd_insert(node_tree, curr.state.x, &(tree.back()) ));
                    num_nodes++;
                    if( is_inside_goal(curr.state) && (curr.csrc < min_cost) )
                    {
                        printf("iter_reached: %d\n", iter);
                        reached = true;
                        min_cost = curr.csrc;
                        node_that_reached = &(tree.back());     // just inserted this curr in the tree, hence valid
                    }
                    iter++;
                }
            }
        }
#if BRANCH_N_BOUND
        if ( (iter % 200 == 0) && reached )
        {
            printf("iter: %d ", iter);
            remove_bad_nodes(min_cost);
        }
#endif
    }
    kd_free(node_tree);

    process_tree_rrt( *node_that_reached );
    
    return min_cost;
};


