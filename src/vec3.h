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

#ifndef __ADENTU_VEC3_H__
#define __ADENTU_VEC3_H__

#define real double 


/* Vector structures */
typedef struct _vec3f vec3f;
struct _vec3f {
    real x, y, z;
};


typedef struct _vec3i vec3i;
struct _vec3i {
    int x, y, z;
};


/* Vector functions */
#define vecAdd(v1, v2, v3)         \
        (v1).x = (v2).x + (v3).x,   \
        (v1).y = (v2).y + (v3).y,   \
        (v1).z = (v2).z + (v3).z

#define vecSub(v1, v2, v3)         \
        (v1).x = (v2).x - (v3).x,   \
        (v1).y = (v2).y - (v3).y,   \
        (v1).z = (v2).z - (v3).z

#define vecDot(v1, v2)             \
        ( (v1).x * (v2).x  +  (v1).y * (v2).y  +  (v1).z * (v2).z )

#define vecMod(v)                   \
        sqrt (vecDot(v, v))

#define vecSet(v1, sx, sy, sz)     \
        (v1).x = sx,                \
        (v1).y = sy,                \
        (v1).z = sz

#define vecScale(v1, v2, s)        \
        (v1).x = (v2).x * s,        \
        (v1).y = (v2).y * s,        \
        (v1).z = (v2).z * s

#define vecCross(v1, v2, v3)        \
        (v1).x = (v2).y * (v3).z  -  (v2).z * (v3).y,   \
        (v1).y = (v2).z * (v3).x  -  (v2).x * (v3).z,   \
        (v1).z = (v2).x * (v3).y  -  (v2).y * (v3).x


void vRand3f (vec3f *v);


void print_vec3f (vec3f *v);
void print_vec3i (vec3i *v);


#define print3i(v)  printf ("(%d, %d, %d)", (v).x, (v).y, (v).z);
#define print3f(v)  printf ("(%f, %f, %f)", (v).x, (v).y, (v).z);



#endif /* __ADENTU_VEC3_H__ */
