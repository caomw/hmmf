#include "random.h"

void init_rand()
{
    srand( time(0));
}

// random double between min and max
double randdouble(double min, double max)
{
    if (min>max)
    {
        return ((double)rand()/RAND_MAX)*(min-max) + max;
    }
    else
    {
        return ((double)rand()/RAND_MAX)*(max-min) + min;
    }
}

// random float between min and max
float randfloat(float min, float max)
{
    if (min>max)
    {
        return ((float)rand()/RAND_MAX)*(min-max) + min;
    }
    else
    {
        return ((float)rand()/RAND_MAX)*(max - min) + min;
    }
}

double get_msec()
{
    struct timeval start;
    gettimeofday(&start, NULL);
    return start.tv_sec*1000 + start.tv_usec/1000.0;
}

