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
#include <string.h>
#include <sys/time.h>

#include "adentu-model.h"
#include "adentu-atom.h"
#include "adentu-grid.h"
#include "adentu-event.h"
#include "adentu-neighbourhood.h"
#include "adentu-runnable.h"
#include "adentu-types.h"


/* Graphics */
#include "adentu-graphic.h"

/* events */
#include "event/adentu-event-mpc.h"
#include "event/adentu-event-bc.h"
#include "event/adentu-event-gfc.h"
#include "event/adentu-event-ggc.h"
#include "event/adentu-event-usr.h"


/* usr modules */
#include "usr/atoms-pos-cuda.h"
#include "usr/print-event-info.h"


AdentuEventHandler *handler[] = {
    [ADENTU_EVENT_START] = NULL,
    [ADENTU_EVENT_MPC] = NULL, //&AdentuMPCEventHandler,
    [ADENTU_EVENT_BC_GRAIN] = &AdentuBCGrainEventHandler,
    [ADENTU_EVENT_BC_FLUID] = &AdentuBCFluidEventHandler,
    [ADENTU_EVENT_GGC] = &AdentuGGCEventHandler, /* NULL,*/
    [ADENTU_EVENT_GFC] = NULL, //&AdentuGFCEventHandler,/* NULL,*/
    [ADENTU_EVENT_USR] = &AdentuUSREventHandler,
    [ADENTU_EVENT_END] = NULL
};




int main (int argc, char *argv[])
{
    g_message ("Reseting CUDA Device.");
    adentu_usr_cuda_reset_device ();

    g_message ("Initializing adentu.");
    //set seeds
    //srand (time (NULL));
    //srand48 (time (NULL));
    srand (1234567);
    srand48 (1234567);
    /* leer configuración */


    /* crear modelo */
    AdentuModel m;
    vecSet (m.accel, 0.0, 0.0, 0.0);
    m.totalTime = 50;
    m.gTemp = 1.0;
    m.fTemp = 0.0;
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


    /* if using ADENTU_BOUNDARY_FBC as BC, 
     * set the grain and fluid velocity */
    vec3f gvel, fvel;
    vecSet (gvel, 0.0, 0.0, 0.0);
    vecSet (fvel, 0.0, 0.0, 0.0);
    adentu_event_bc_set_fbc_vel (gvel, fvel);
   

    /* if using MPC events, set the interval time
     * which it will generate events. Also set the
     * rotation angle.*/
    adentu_event_mpc_set_dt (0.001);
    adentu_event_mpc_set_alpha (3.141592653589793238462);

    
    /* creating grain grid */
    AdentuGridConfig gc;
    vecSet (gc.origin, 0.0, 0.0, 0.0);
    vecSet (gc.length, 3.1, 3.1, 3.1);
    vecSet (gc.cells, 3, 3, 3);
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
    ac.nAtoms = 2;
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
    ac.nAtoms =  2540;
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


    /* oink test */
    vecSet (m.grain->pos[0], 1.5500, 2.3500, 1.5500);
    vecSet (m.grain->vel[0], 0, -1, 0);
    vecSet (m.grain->pos[1], 1.5500, 1.2500, 1.5500);
    vecSet (m.grain->vel[1], 0, 1, 0);



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





    printf ("Acceleration: ");
    print_vec3f (&m.accel);
    printf ("Boundary Conditions: ");
    printf ("(%s, %s, %s)\n", AdentuBoundaryTypeStr[m.bCond.x], 
                              AdentuBoundaryTypeStr[m.bCond.y], 
                              AdentuBoundaryTypeStr[m.bCond.z]);




    //adentu_runnable_add_pre_func (&m, print_pre_event);
    adentu_runnable_add_post_func (&m, print_event);
    //adentu_runnable_add_post_func (&m, print_post_event);

    /* setup event engine */
    m.eList = adentu_event_init (handler, &m);
    puts ("");



    /* run in graphic mode or text mode */
    if (argc == 2 && !strncmp (argv[1], "-g", 2))
        {
            adentu_graphic_init (argc, argv, &m, &handler);
            adentu_graphic_set_time_sleep (0000);
            adentu_graphic_start ();
        }
    else
        m.eList = adentu_event_loop (handler, &m);


    //adentu_event_usr_set_dt (.5);
    //m.eList = adentu_event_loop (handler, &m);
    /* graphics (: */
    









    /* Destroy everything */


    return 0;
}

