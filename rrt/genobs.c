#include <stdio.h>
#include <stdlib.h>

#define MIN         (-10)
#define MAX         (10)

double rd()
{
    return (double)rand()/RAND_MAX;
}

double rdlim(double min, double max)
{
    return (max - min)*rd() + min;
}

bool is_valid(double x, double y, double r)
{
    double dist = x*x + y*y;
    if (dist < (r+1))
        return false;
    else
        return true;
}

int main(int argc, char *argv[])
{
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
        double r = rdlim(1, max_size);

        if( is_valid(x, y, r))
            printf("%f, %f, %f\n", x, y, r);
    }
    
    return 0;
}
