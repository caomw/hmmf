#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN         (-50)
#define MAX         (50)

double rd()
{
    return (double)rand()/RAND_MAX;
}

double rdlim(double min, double max)
{
    return ( (max - min)*rd() + min);
}

int is_valid(double x, double y, double r)
{
    double dist = sqrt(x*x + y*y);
    if (dist < (r+1))
        return 0;
    else
        return 1;
}

int main(int argc, char *argv[])
{
    srand(time(0));
    if(argc != 3)
    {
        printf("Usage: ./genobs <num_of_obs>  <max_size>\n");
        exit(0);
    }

    int num_obs = atoi(argv[1]);
    double max_size = atoi(argv[2]);

    int c = 0;
    while(c < num_obs)
    {
        double x = rdlim(MIN, MAX);
        double y = rdlim(MIN, MAX);
        double r;
        if(c < 0.9*num_obs)
            r = rdlim(0.5, 0.5*max_size);
        else
            r = rdlim(0.5*max_size, max_size);

        if( is_valid(x, y, r))
        {
            c++;
            printf("%f, %f, %f\n", x, y, r);
        }
    }
    
    return 0;
}
