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

const char *AdentuBoundaryTypeStr[] = {
    [ADENTU_BOUNDARY_CBP] = "CBP",  
    [ADENTU_BOUNDARY_CBB] = "CBB",  
    [ADENTU_BOUNDARY_CBR] = "CBR",  
    [ADENTU_BOUNDARY_CBF] = "CBF"
};



int main (int argc, char *argv[])
{

    //set seeds
    srand (time (NULL));
    srand48 (time (NULL));
    /* leer configuración */


    /* crear modelo */
    AdentuModel m;
    vecSet (m.accel, 0.0, 1.0, 0.0);
    m.totalTime = 500;
    m.deltaTime = 50;
    m.tempGrain = 4.0;
    m.tempFluid = 0.0;
    vecSet (m.vcmGrain, 1., 1., 1.);
    vecSet (m.vcmFluid, 0., 0., 0.);
    vecSet (m.bCond, ADENTU_BOUNDARY_CBP, 
            ADENTU_BOUNDARY_CBP, ADENTU_BOUNDARY_CBP);

    
    
    /* creating grid */
    AdentuGridConfig gc;
    vecSet (gc.origin, 0.0, 0.0, 0.0);
    vecSet (gc.length, 10.0, 20.0, 30.0);
    vecSet (gc.cells, 5, 2, 3);
    gc.type = ADENTU_GRID_DEFAULT;

    AdentuGrid g;
    adentu_grid_set_from_config (&g, &gc);

    /* create grains */
    AdentuAtomConfig ac;
    ac.nAtoms = 32;
    ac.type = ADENTU_ATOM_GRAIN;
    ac.mass.from = ac.mass.to = 5.0;
    ac.mass.rangeType = ADENTU_PROP_CONSTANT;
    ac.radii.from = ac.radii.to = 2.0;
    ac.radii.rangeType = ADENTU_PROP_CONSTANT;

    AdentuAtom a;
    adentu_atom_create_from_config (&a, &ac);
    adentu_atom_set_init_vel (&a, &m);
    adentu_atom_set_init_pos (&a, &g);

    /* set atoms into grid */
    adentu_grid_set_atoms (&g, &a, &m);


    /* create event structures */
    AdentuEvent *ge = NULL; 
    adentu_event_set_init_events (ge, &g);


    /* get neighbours */
    int nNeighbours;
    int *neighbours = adentu_neighbourhood_get_atom_neighbourhood (20, 
                                                                   &nNeighbours, 
                                                                   &a, &g, &m);


    /* check */
    vec3f half, center;
    vecScale (half, g.length , 0.5);
    center.x = g.origin.x + half.x;
    center.y = g.origin.y + half.y;
    center.z = g.origin.z + half.z;
    
    printf ("Origin: ");
    print_vec3f (&g.origin);
    printf ("Length: ");
    print_vec3f (&g.length);
    printf ("Half:   ");
    print_vec3f (&half);
    printf ("Center: ");
    print_vec3f (&center);
    printf ("Boundary Conditions: ");
    printf ("(%s, %s, %s)\n", AdentuBoundaryTypeStr[m.bCond.x],
                              AdentuBoundaryTypeStr[m.bCond.y],
                              AdentuBoundaryTypeStr[m.bCond.z]);

    for (int i = 0; i < a.n; ++i)
    {
        //printf ("Atom %u - Vel: ", i);
        //print_vec3f (a.vel + i);
        //printf ("\t\tPos: ");
        printf ("Atom %u - Pos: ", i);
        print_vec3f (a.pos + i);
    }

    
    puts ("nAtoms at Cells:");
    for (int i = 0; i < g.tCell; ++i)
        printf ("cell[%d] = %d atoms\n", i, g.cells.nAtoms[i]);
    
    
    puts ("HEAD:");
    for (int i = 0; i < g.tCell; ++i)
        printf ("head[%d] = %d\n", i, g.head[i]);

    puts ("LINKED:");
    for (int i =0; i < a.n; ++i)
        printf ("linked[%d] = %d\n", i, g.linked[i]);


/*    puts ("GRID:");
    for (int i = 0; i < g.tCell; ++i)
    {
        printf ("Celda %3d", i);
        if (g.head[i] == -1)
            puts ("\tCell is empty");
        else{
            int x = g.head[i];
            int y = g.linked[x];
            printf ("\tAtoms: %d - y: %d ", x, y);
            while (y != -1){
                x = g.linked[y];
                printf ("%d ", x);
                y = g.linked[x];
            }
            printf ("\n");

        }
    }
*/

    printf ("\nAtom:20's neighbours are %d in total\n", nNeighbours);

    for (int i=0; i < nNeighbours; ++i)
        printf ("%2d, ", neighbours[i]);
    puts ("");

    return 0;
}
