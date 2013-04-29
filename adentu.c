/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "adentu-model.h"
#include "adentu-atom.h"
#include "adentu-grid.h"
#include "adentu-event.h"
#include "adentu-neighbourhood.h"
#include "vec3.h"


/* Graphics */
#ifdef ADENTU_GRAPHICS
    #include "adentu-graphic.h"
#endif /* ADENTU_GRAPHICS*/


/* events */
#include "adentu-event-mpc.h"
#include "adentu-event-bc.h"


const char *AdentuBoundaryTypeStr[] = {
    [ADENTU_BOUNDARY_PBC] = "PBC",  
    [ADENTU_BOUNDARY_BBC] = "BBC",  
    [ADENTU_BOUNDARY_RBC] = "RBC",  
    [ADENTU_BOUNDARY_FBC] = "FBC"
}; 



AdentuEventHandler *handler[] = {
    [ADENTU_EVENT_START] = NULL,
    [ADENTU_EVENT_MPC] = &AdentuMPCEventHandler,
    [ADENTU_EVENT_BC_GRAIN] = NULL, //&AdentuBCEventHandler,
    [ADENTU_EVENT_BC_FLUID] = &AdentuBCEventHandler,
    [ADENTU_EVENT_GGC] = NULL,
    [ADENTU_EVENT_GFC] = NULL, 
    [ADENTU_EVENT_END] = NULL
};


int main (int argc, char *argv[])
{

    //set seeds
    srand (time (NULL));
    srand48 (time (NULL));
    /* leer configuración */


    /* crear modelo */
    AdentuModel m;
    vecSet (m.accel, 3.0, 2.0, 1.0);
    m.totalTime = 100;
    m.dT = 10;
    m.alpha = 3.141593;
    m.gTemp = 4.0;
    m.fTemp = 0.0;
    vecSet (m.gVel, 5.0, 1.0, 3.0);
    vecSet (m.fVel, 1.0, 1.0, 1.0);
    vecSet (m.vcmGrain, 1., 1., 1.);
    vecSet (m.vcmFluid, 0., 0., 0.);
    vecSet (m.bCond, ADENTU_BOUNDARY_PBC, 
            ADENTU_BOUNDARY_PBC, ADENTU_BOUNDARY_PBC);
    m.grain = NULL;
    m.fluid = NULL;
    m.gGrid = NULL;
    m.fGrid = NULL;
    m.mpcGrid = NULL;

    
    
    /* creating grain grid */
    AdentuGridConfig gc;
    vecSet (gc.origin, 0.0, 0.0, 0.0);
    vecSet (gc.length, 10.0, 20.0, 30.0);
    vecSet (gc.cells, 5, 2, 3);
    gc.type = ADENTU_GRID_DEFAULT;

  //  AdentuGrid g;
    //adentu_grid_set_from_config (&g, &gc);

    /*set grit into the model*/
   // m.gGrid = &g;

    /* creating fliud and MPC grid */
    AdentuGrid fg;
    adentu_grid_set_from_config (&fg, &gc);

    gc.type = ADENTU_GRID_MPC;
    AdentuGrid mpcg;
    adentu_grid_set_from_config (&mpcg, &gc);

    m.fGrid = &fg;
    m.mpcGrid = &mpcg;








    /* create grains */
    AdentuAtomConfig ac;
    ac.nAtoms = 8;
    ac.type = ADENTU_ATOM_GRAIN;
    ac.mass.from = ac.mass.to = 5.0;
    ac.mass.rangeType = ADENTU_PROP_CONSTANT;
    ac.radii.from = ac.radii.to = 2.0;
    ac.radii.rangeType = ADENTU_PROP_CONSTANT;

  /*  AdentuAtom a;
    adentu_atom_create_from_config (&a, &ac);
    adentu_atom_set_init_vel (&a, &m);
    adentu_atom_set_init_pos (&a, &g); */

    /*set grains into the model*/
    //m.grain = &a;

    /* set atoms into grid */
    //adentu_grid_set_atoms (&g, &a, &m);


    /****************************************************/
    /* creating fluid*/
    ac.nAtoms = 32;
    ac.type = ADENTU_ATOM_FLUID;
    ac.mass.from = ac.mass.to = 0.5;
    ac.mass.rangeType = ADENTU_PROP_CONSTANT;
    ac.radii.from = ac.radii.to = 0.1;
    ac.radii.rangeType = ADENTU_PROP_CONSTANT;

    AdentuAtom f;
    adentu_atom_create_from_config (&f, &ac);
    adentu_atom_set_init_vel (&f, &m);
    adentu_atom_set_init_pos (&f, &fg);
    m.fluid = &f;
    adentu_grid_set_atoms (&fg, &f, &m);
    adentu_grid_set_atoms (&mpcg, &f, &m);




    /* General debug Info */
    vec3f half, center;
    vecScale (half, fg.length , 0.5);
    center.x = fg.origin.x + half.x;
    center.y = fg.origin.y + half.y;
    center.z = fg.origin.z + half.z;
               
    printf ("Origin: ");
    print_vec3f (&m.fGrid->origin);
    printf ("Length: ");
    print_vec3f (&fg.length);
    printf ("Half:   ");
    print_vec3f (&half);
    printf ("Center: ");
    print_vec3f (&center);
    printf ("Acceleration: ");
    print_vec3f (&m.accel);
    printf ("gVel: ");
    print_vec3f (&m.gVel);
    printf ("Boundary Conditions: ");
    printf ("(%s, %s, %s)\n", AdentuBoundaryTypeStr[m.bCond.x], 
                              AdentuBoundaryTypeStr[m.bCond.y], 
                              AdentuBoundaryTypeStr[m.bCond.z]);


   /* printf ("Created %d grains, mass: %f, radii: %f\n", m.grain->n, 
            m.grain->mass[0], m.grain->radius[0]);
    for (int i = 0, j = m.grain->n; i < j; ++i)
        printf ("%3d    %f, %f, %f    %f, %f, %f\n", i,
                m.grain->pos[i].x, m.grain->pos[i].y, m.grain->pos[i].z, 
                m.grain->vel[i].x, m.grain->vel[i].y, m.grain->vel[i].z);
*/


    /* setup event engine */
    GSList *eList = NULL;

    
#ifndef ADENTU_GRAPHICS
    eList = adentu_event_init (eList, handler, &m);
    puts (""); 
    /* start event engine loop */
    eList = adentu_event_loop (eList, handler, &m);
#elif
    /* graphics (: */

#endif /* ADENTU_GRAPHICS*/









    /* Destroy everything */


    return 0;
}

