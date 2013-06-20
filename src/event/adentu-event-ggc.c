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
#include <math.h>
#include <glib.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "adentu-event.h"
#include "vec3.h"

#include "adentu-event-ggc.h"
#include "adentu-event-ggc-cuda.h"
#include "adentu-event-mpc-cuda.h"

#include "adentu-event-bc.h"

#include "adentu-graphic.h"


AdentuEventHandler AdentuGGCEventHandler = {adentu_event_ggc_init,
                                            adentu_event_ggc_is_valid,
                                            adentu_event_ggc_attend,
                                            adentu_event_ggc_get_next};



GSList *adentu_event_ggc_init (AdentuModel *model)//,
                               //GSList *eList)
{
    model->eList = adentu_event_schedule (model->eList, 
                            adentu_event_ggc_get_next (model));


    return model->eList;
}


AdentuEvent *adentu_event_ggc_get_next (AdentuModel *model)
{
    adentu_grid_set_atoms (model->gGrid, 
                           model->grain, 
                           &model->bCond);

    AdentuEvent *ev = adentu_event_ggc_cuda_get_next (model);
    ev->time += model->elapsedTime;

    /*g_message ("New GGC event, time: %f", ev->time);
    getchar ();
*/
    return ev;
}


int adentu_event_ggc_is_valid (AdentuModel *model,
                               AdentuEvent *event)
{
    if (event->type != ADENTU_EVENT_GGC)
        return 0;

    int owner = event->owner;
    int partner = event->partner;
    int nEvents = event->nEvents;
    int eventData = *(int *)event->eventData;
    AdentuAtom *grain = model->grain;

    if (grain->nCol[owner] == nEvents &&
        grain->nCol[partner] == eventData)
            return 1;

    return 0;
}


void adentu_event_ggc_attend (AdentuModel *model,
                               AdentuEvent *event)
{
    AdentuAtom *grain = model->grain;

    double dT = event->time - model->elapsedTime;

    adentu_event_mpc_cuda_integrate (grain,
                                     model->gGrid,
                                     model->accel, dT);
    adentu_grid_set_atoms (model->gGrid, 
                           grain, 
                           &model->bCond);


    /* if graphics are set, update it*/
    if (adentu_graphic_is_set)
        adentu_graphic_display ();


    int owner = event->owner;
    int partner = event->partner;
    
    vec3f *gvel = &grain->vel[owner];
    vec3f *fvel = &grain->vel[partner];
    vec3f *gpos = &grain->pos[owner];
    vec3f *fpos = &grain->pos[partner];
    double gmass = grain->mass[owner];
    double fmass = grain->mass[partner];
    double gradius = grain->radius[owner];
    double fradius = grain->radius[partner];

    double radius = gradius + fradius;
    vec3f pos, vel;
    vecSub (pos, *fpos, *gpos);
    vecSub (vel, *fvel, *gvel);

    /* testing */
    vec3f p1, p2, diffp;
    vec3f p1p, p2p;
    vecScale (p1, gvel[0], gmass);
    vecScale (p2, fvel[0], fmass);
    vec3f pAntes;
    vecAdd (pAntes, p1, p2);

    double e1, e2, diffe, eBefore;
    e1 = gmass * vecMod (gvel[0]);
    e2 = fmass * vecMod (fvel[0]);
    eBefore = e1 + e2;



    double pu = vecMod (pos);
    
    if (fabs (pu - radius) > 10e-6)
        {
            g_error ("Bad Prediction! - PU: %f != Radius: %f", pu, radius);
            //getchar ();
            //model->eList = adentu_event_schedule (model->eList,
            //                                adentu_event_ggc_get_next (model));
            //getchar ();
            //return ;
        }

    vec3f n;
    vecScale (n, pos, pu);
   
    g_message ("radius: %f, pu: %f, pu-radius: %.9f\npos: (%f, %f, %f), n: (%f, %f, %f), vecDot (n, n) = %f", 
                radius, pu, pu-radius, pos.x, pos.y, pos.z, n.x, n.y, n.z, vecDot (n, n));
    //getchar ();


    /* update states */
    double VN = vecDot (n, vel);

    double dP = (2 * gmass * fmass * VN) / (gmass + fmass);
    vec3f j;
    
    j.x = dP * n.x;
    j.y = dP * n.y;
    j.z = dP * n.z;
   
    gvel->x = gvel->x + j.x / gmass;
    gvel->y = gvel->y + j.y / gmass;
    gvel->z = gvel->z + j.z / gmass;

    fvel->x = fvel->x - j.x / fmass;
    fvel->y = fvel->y - j.y / fmass;
    fvel->z = fvel->z - j.z / fmass;
    grain->nCol[owner]++;
    grain->nCol[partner]++;



    /* testing */
    vecScale (p1p, gvel[0], gmass);
    vecScale (p2p, fvel[0], fmass);
    vec3f pDespues;
    vecAdd (pDespues, p1p, p2p);

    vecSub (diffp, pDespues, pAntes);

    double e1p, e2p, eAfter;
    e1p = gmass * vecMod (gvel[0]);
    e2p = fmass * vecMod (fvel[0]);
    eAfter = e1p + e2p;
    diffe = eAfter - eBefore;
    
    g_message ("Momentum vecMod(diffp): %.8f, Kinetic Energy diffe: %.8f", vecMod(diffp), diffe);

    //getchar ();


    /* check next BC event */
    model->elapsedTime = event->time;
    model->eList = adentu_event_schedule (model->eList,
                                    adentu_event_bc_grain_get_next (model));
    
    /* model->eList = adentu_event_schedule (model->eList,
                                    adentu_event_bc_fluid_get_next (model));
    */

}



