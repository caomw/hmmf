#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

void init_rand();

int randint(int min = 0, int max = 32767);
double randdouble(double min = 0, double max = 1);
float randfloat(float min = 0, float max = 1);

double get_msec();

#endif
