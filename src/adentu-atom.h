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

#ifndef __ADENTU_ATOM_H__
#define __ADENTU_ATOM_H__


#include "adentu.h"

typedef enum {
    ADENTU_ATOM_GRAIN,
    ADENTU_ATOM_FLUID
} AdentuAtomType;


typedef struct _AdentuAtom {
    AdentuAtomType type;
    int n;
    double *h_pos;
    double *h_vel;
    double *h_velRel;
    //int *h_lastTime;
    int *h_nCol;
    double *h_mass;
    double *h_radius;

    double *d_pos;
    double *d_vel;
    double *d_velRel;
    //int *d_lastTime;
    int *d_nCol;
    double *d_mass;
    double *d_radius;
} AdentuAtom;


typedef struct _AdentuPropRange {
    double from, to;
    enum {ADENTU_PROP_CONSTANT, ADENTU_PROP_NORMAL, ADENTU_PROP_DELTA} rangeType;
} AdentuPropRange;

typedef struct _AdentuAtomConfig {
    int nAtoms;
    AdentuAtomType type;
    AdentuPropRange mass;
    AdentuPropRange radii;
} AdentuAtomConfig;


void adentu_atom_create_from_config (AdentuAtom *atoms, AdentuAtomConfig *conf);
//void adentu_atom_set_init_vel (AdentuAtom *atoms, AdentuModel *model);
//void adentu_atom_set_init_pos (AdentuAtom *atoms, AdentuGrid *grid)



#endif  /* __ADENTU_ATOM_H__ */
