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
#include <sys/time.h>

#include "adentu-model.h"
#include "adentu-atom.h"
#include "adentu-grid.h"
#include "adentu-event.h"
#include "adentu-neighbourhood.h"
#include "adentu-runnable.h"
#include "vec3.h"


/* Graphics */
#include "adentu-graphic.h"
#include <GLUT/glut.h>

/* events */
#include "adentu-event-mpc.h"
#include "adentu-event-bc.h"
#include "adentu-event-gfc.h"
#include "adentu-event-ggc.h"
#include "adentu-event-usr.h"


/* usr modules */
#include "usr/atoms-pos-cuda.h"
#include "usr/print-event-info.h"


AdentuEventHandler *handler[] = {
    [ADENTU_EVENT_START] = NULL,
    [ADENTU_EVENT_MPC] = &AdentuMPCEventHandler,
    [ADENTU_EVENT_BC_GRAIN] = &AdentuBCGrainEventHandler,
    [ADENTU_EVENT_BC_FLUID] = &AdentuBCFluidEventHandler,
    [ADENTU_EVENT_GGC] = &AdentuGGCEventHandler, /* NULL,*/
    [ADENTU_EVENT_GFC] = &AdentuGFCEventHandler,/* NULL,*/
    [ADENTU_EVENT_USR] = &AdentuUSREventHandler,
    [ADENTU_EVENT_END] = NULL
};




int main (int argc, char *argv[])
{
    g_message ("Reseting CUDA Device.");
    adentu_usr_cuda_reset_device ();

    g_message ("Initializing adentu.");
    //set seeds
    srand (time (NULL));
    srand48 (time (NULL));
    //srand (1234567);
    //srand48 (1234567);
    /* leer configuración */


    /* crear modelo */
    AdentuModel m;
    vecSet (m.accel, 0.0, -1.0, 0.0);
    m.totalTime = 50;
    m.dT = 1;
    m.alpha = 3.141593;
    m.gTemp = 0.0;
    m.fTemp = 0.0;
    vecSet (m.gVel, 0.0, 0.0, 0.0);
    vecSet (m.fVel, 0.0, 0.0, 0.0);
    vecSet (m.vcmGrain, 5., 0., 0.);
    vecSet (m.vcmFluid, 2., 0., 0.);
    vecSet (m.bCond, ADENTU_BOUNDARY_BBC, 
            ADENTU_BOUNDARY_BBC, ADENTU_BOUNDARY_BBC);
    m.grain = NULL;
    m.fluid = NULL;
    m.gGrid = NULL;
    m.fGrid = NULL;
    m.mpcGrid = NULL;
    m.pre_event_func = NULL;
    m.post_event_func = NULL;

    
    
    /* creating grain grid */
    AdentuGridConfig gc;
    vecSet (gc.origin, 0.0, 0.0, 0.0);
    vecSet (gc.length, 3.1, 3.1, 3.1);
    vecSet (gc.cells, 2, 1, 1);
    gc.type = ADENTU_GRID_DEFAULT;

    AdentuGrid g;
    adentu_grid_set_from_config (&g, &gc);

    /*set grid into the model*/
    m.gGrid = &g;

    /* creating fliud and MPC grid */
    AdentuGrid fg;
    adentu_grid_set_from_config (&fg, &gc);
    m.fGrid = &fg;

    gc.type = ADENTU_GRID_MPC;
    AdentuGrid mpcg;
    adentu_grid_set_from_config (&mpcg, &gc);

    m.mpcGrid = &mpcg;








    /* create grains */
    AdentuAtomConfig ac;
    ac.nAtoms = 1;
    ac.type = ADENTU_ATOM_GRAIN;
    ac.mass.from = ac.mass.to = 1.0;
    ac.mass.rangeType = ADENTU_PROP_CONSTANT;
    ac.radii.from = ac.radii.to = 0.500000;
    ac.radii.rangeType = ADENTU_PROP_CONSTANT;

    AdentuAtom a;
    adentu_atom_create_from_config (&a, &ac);
    adentu_atom_set_init_vel (&a, &m);
    //adentu_atom_set_init_pos (&a, &g); 

    /*set grains into the model*/
    m.grain = &a;

    /* set atoms into grid */
    //adentu_grid_set_atoms (&g, &a, &m.bCond);


    /****************************************************/
    /* creating fluid*/
    ac.nAtoms = 1;
    ac.type = ADENTU_ATOM_FLUID;
    ac.mass.from = ac.mass.to = 0.5;
    ac.mass.rangeType = ADENTU_PROP_CONSTANT;
    ac.radii.from = ac.radii.to = 0.00000000000001;
    ac.radii.rangeType = ADENTU_PROP_CONSTANT;

    AdentuAtom f;
    adentu_atom_create_from_config (&f, &ac);
    adentu_atom_set_init_vel (&f, &m);
    //adentu_atom_set_init_pos (&f, &fg);
    m.fluid = &f;
    //adentu_grid_set_atoms (&fg, &f, &m.bCond);
    //adentu_grid_set_atoms (&mpcg, &f, &m.bCond);


    adentu_usr_cuda_set_atoms_pos (&m);


    /* General debug Info */
    vec3f half, center;
    vecScale (half, m.gGrid->length , 0.5);
    center.x = m.gGrid->origin.x + half.x;
    center.y = m.gGrid->origin.y + half.y;
    center.z = m.gGrid->origin.z + half.z;
               
    printf ("gGrid Origin: ");
    print_vec3f (&m.gGrid->origin);
    printf ("gGrid Length: ");
    print_vec3f (&m.gGrid->length);
    printf ("gGrid Half:   ");
    print_vec3f (&half);
    printf ("gGrid Center: ");
    print_vec3f (&center);
/*
    vecScale (half, m.fGrid->length , 0.5);
    center.x = m.fGrid->origin.x + half.x;
    center.y = m.fGrid->origin.y + half.y;
    center.z = m.fGrid->origin.z + half.z;
               
    printf ("fGrid Origin: ");
    print_vec3f (&m.fGrid->origin);
    printf ("fGrid Length: ");
    print_vec3f (&m.fGrid->length);
    printf ("fGrid Half:   ");
    print_vec3f (&half);
    printf ("fGrid Center: ");
    print_vec3f (&center);

    vecScale (half, m.mpcGrid->length , 0.5);
    center.x = m.mpcGrid->origin.x + half.x;
    center.y = m.mpcGrid->origin.y + half.y;
    center.z = m.mpcGrid->origin.z + half.z;
               
    printf ("mpcGrid Origin: ");
    print_vec3f (&m.mpcGrid->origin);
    printf ("mpcGrid Length: ");
    print_vec3f (&m.mpcGrid->length);
    printf ("mpcGrid Half:   ");
    print_vec3f (&half);
    printf ("mpcGrid Center: ");
    print_vec3f (&center);

*/




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


    adentu_runnable_add_pre_func (&m, print_pre_event);
    //adentu_runnable_add_post_func (&m, print_event);
    adentu_runnable_add_post_func (&m, print_post_event);

    /* setup event engine */
    //GSList *eList = NULL;

    
    m.eList = adentu_event_init (handler, &m);
    puts ("");

    adentu_event_usr_set_dt (.5);
//    eList = adentu_event_loop (m.eList, handler, &m);
    /* graphics (: */
    adentu_graphic_init (argc, argv, &m, &handler);//&eList, &handler);
    adentu_graphic_set_time_sleep (50000);
    glutMainLoop ();









    /* Destroy everything */


    return 0;
}

