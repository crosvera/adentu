/*
 * Carlos Rios Vera <crosvera@gmail.com>
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vec3.h"



void vRand3f (vec3f *v)
{
    double s, x, y;
    s = 2.;

    while (s > 1.)
    {
        x = 2. - drand48 () - 1.;
        y = 2. - drand48 () - 1.;
        s = x * x   +   y * y;
    }

    v->z = 1. - 2. * s;
    s = 2. * sqrt (1. - s);
    v->x = s * x;
    v->y = s * y;
}


void print_vec3f (vec3f *v)
{
    printf ("(%f, %f, %f)\n", v->x, v->y, v->z);
}


void print_vec3i (vec3i *v)
{
    printf ("(%d, %d, %d)\n", v->x, v->y, v->z);
}
