#include "random.h"
#include "kdtree.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

static double dist_sq( double *a1, double *a2, int dims ) 
{
    double dist_sq = 0, diff;
    while( --dims >= 0 ) {
        diff = (a1[dims] - a2[dims]);
        dist_sq += diff*diff;
    }
    return dist_sq;
}

int main()
{
    init_rand();
    
    struct kdtree *ptree;
    struct kdres *presults;

    ptree = kd_create(2);

    double pos[2];
    double org[] = {0.0, 0.0};
    
    printf("%.3f\n", get_msec());

    for(int i=0; i<10; i++)
    {
        double radius = 0.3;
        double data[2] = {randdouble(), randdouble()};
        printf("%.3f %.3f\n", data[0], data[1]);
        kd_insert( ptree, data, &radius);
    }
    printf("%.3f\n", get_msec());

    presults = kd_nearest_range(ptree, org, 2.0);

    int count = 0;
    while( !kd_res_end( presults ) ) 
    {
        double* rad = (double*)kd_res_item( presults, pos );

        float dist = sqrt( dist_sq( org, pos, 2 ) );

        printf( "%dth node at (%.3f, %.3f) is %.3f away and has size=%.3f\n", count++, 
                pos[0], pos[1], dist, *rad);

        kd_res_next( presults );
    }

    kd_res_free( presults );
    kd_free( ptree );
    
    return 0;
}
