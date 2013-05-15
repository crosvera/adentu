/*
    Adentu: An hybrid molecular dynamic software.
    https://github.com/crosvera/adentu
    
    Copyright (C) 2013 Carlos Ríos Vera <crosvera@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
