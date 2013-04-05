/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#ifndef __ADENTU_ATOM_H__
#define __ADENTU_ATOM_H__


#include "vec3.h"
//#include "adentu-model.h"
//#include "adentu-grid.h"

typedef enum {
    ADENTU_ATOM_GRAIN,
    ADENTU_ATOM_FLUID
} AdentuAtomType;

/*const char *AdentuAtomTypeStr[] = { 
                        [ADENTU_ATOM_GRAIN] = "GRAIN ATOM",
                        [ADENTU_ATOM_FLUID] = "FLUID ATOM"
}; */



typedef struct _AdentuAtom {
    AdentuAtomType type;
    int n;
    vec3f *pos;
    vec3f *vel;
    int *lastTime;
    int *nCol;

    double *mass;
    double *radius;

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
