#include "random.h"

void init_rand()
{
    srand( time(0));
}

// random int between min and max
int randint(int min, int max)
{
    if (min>max)
    {
        return max+int((min-max+1)*((double)rand())/(RAND_MAX+1.0));
    }
    else
    {
        return min+int((max-min+1)*((double)rand())/(RAND_MAX+1.0));
    }
}

// random double between min and max
double randdouble(double min, double max)
{
    if (min>max)
    {
        return rand()/(double(RAND_MAX)+1)*(min-max)+max;
    }
    else
    {
        return rand()/(double(RAND_MAX)+1)*(max-min)+min;
    }
}

// random float between min and max
float randfloat(float min, float max)
{
    if (min>max)
    {
        return rand()/(float(RAND_MAX)+1)*(min-max)+max;
    }
    else
    {
        return rand()/(float(RAND_MAX)+1)*(max-min)+min;
    }
}

double get_msec()
{
    struct timeval start;
    gettimeofday(&start, NULL);
    return start.tv_sec*1000 + start.tv_usec/1000.0;
}

