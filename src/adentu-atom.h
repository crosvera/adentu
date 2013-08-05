/*
    Adentu: An hybrid molecular dynamic software.
    https://github.com/crosvera/adentu
    
    Copyright (C) 2013 Carlos Ríos Vera <crosvera@gmail.com>
    Universidad del Bío-Bío.

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

#ifndef __ADENTU_ATOM_H__
#define __ADENTU_ATOM_H__


#include "adentu-types.h"

typedef enum {
    ADENTU_ATOM_GRAIN,
    ADENTU_ATOM_FLUID
} AdentuAtomType;


typedef struct _AdentuAtom {
    AdentuAtomType type;
    int n;
    adentu_real *h_pos;
    adentu_real *h_vel;
    adentu_real *h_velRel;
    //int *h_lastTime;
    int *h_nCol;
    float *h_mass;
    float *h_radius;

    adentu_real *d_pos;
    adentu_real *d_vel;
    adentu_real *d_velRel;
    //int *d_lastTime;
    int *d_nCol;
    float *d_mass;
    float *d_radius;
} AdentuAtom;


typedef enum {
    ADENTU_PROP_CONSTANT,
    ADENTU_PROP_NORMAL,
    ADENTU_PROP_DELTA} AdentuRangeType;

typedef struct _AdentuPropRange {
    double from, to;
    AdentuRangeType rangeType;
} AdentuPropRange;

typedef struct _AdentuAtomConfig {
    int nAtoms;
    AdentuAtomType type;
    AdentuPropRange mass;
    AdentuPropRange radii;
} AdentuAtomConfig;


void adentu_atom_create_from_config (AdentuAtom *atoms, AdentuAtomConfig *conf);
//void adentu_atom_set_init_vel (AdentuAtom *atoms, AdentuModel *model);
//void adentu_atom_set_random_pos (AdentuAtom *atoms, AdentuGrid *grid)



#endif  /* __ADENTU_ATOM_H__ */
