#include "common.h"

double curr_time;
void tic()
{
    struct timeval start;
    gettimeofday(&start, NULL);
    curr_time = start.tv_sec*1000 + start.tv_usec/1000.0;
}

double toc()
{
    struct timeval start;
    gettimeofday(&start, NULL);
    double delta_t = start.tv_sec*1000 + start.tv_usec/1000.0 - curr_time;
    
    //cout<< delta_t/1000.0 << " [s]\n";
    return delta_t/1000.0;
}

double get_msec()
{
    struct timeval start;
    gettimeofday(&start, NULL);
    return start.tv_sec*1000 + start.tv_usec/1000.0;
}

