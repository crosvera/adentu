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

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <glib.h>
#include <execinfo.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "adentu-event.h"
#include "adentu-types.h"

#include "adentu-event-gfc.h"
#include "adentu-event-gfc-cuda.h"


#include "adentu-graphic.h"

const char *ADENTU_EVENT_GFC = "ADENTU_EVENT_GFC";

AdentuEventHandler AdentuGFCEventHandler = {adentu_event_gfc_init,
                                            adentu_event_gfc_is_valid,
                                            adentu_event_gfc_attend,
                                            adentu_event_gfc_get_next};



GSList *adentu_event_gfc_init (AdentuModel *model)
{
    model->eList = adentu_event_schedule (model->eList, 
                                adentu_event_gfc_get_next (model));


    return model->eList;
}


AdentuEvent *adentu_event_gfc_get_next (AdentuModel *model)
{
    adentu_grid_set_atoms (model->gGrid, 
                           model->grain, 
                           &model->bCond);
    adentu_grid_set_atoms (model->fGrid, 
                           model->fluid, 
                           &model->bCond);

    AdentuEvent *ev = adentu_event_gfc_cuda_get_next (model);
    ev->time += model->elapsedTime;

    /* testing */
    //g_message ("New GFC event, time: %f, partner: %d", ev->time, ev->partner);
    /* testing, get the backtrace of caller functions */
#ifdef Debug
    void *buffer[100];
    char **strings;
    int nptrs = backtrace (buffer, 100);
    strings = backtrace_symbols (buffer, nptrs);
    for (int i = 0; i < 4; ++i)
        g_message ("%s", strings[i]);
    free (strings);
#endif 

    return ev;
}


int adentu_event_gfc_is_valid (AdentuModel *model,
                               AdentuEvent *event)
{
    if (event->type != ADENTU_EVENT_GFC)
        return 0;

    int owner = event->owner;
    int nEvents = event->nEvents;
    AdentuAtom *grain = model->grain;

    if (grain->nCol[owner] == nEvents)
        return 1;

    return 0;
}



void adentu_event_gfc_attend (AdentuModel *model,
                              AdentuEvent *event)
{
    AdentuAtom *grain = model->grain;
    AdentuAtom *fluid = model->fluid;

    double dT = event->time - model->elapsedTime;

    adentu_cuda_integrate_atoms (grain,
                                 model->gGrid,
                                 model->accel, dT);

    adentu_cuda_integrate_atoms (fluid,
                                 model->fGrid,
                                 model->accel, dT);


    adentu_grid_set_atoms (model->gGrid, 
                           grain, 
                           &model->bCond);
    adentu_grid_set_atoms (model->fGrid, 
                           fluid, 
                           &model->bCond);

    /* if graphics are set, update it*/
    if (adentu_graphic_is_set)
        adentu_graphic_display ();

    int owner = event->owner;
    int partner = event->partner;
    
    vec3f gvel = get_vec3f_from_array4f(grain->h_vel, owner);
    vec3f fvel = get_vec3f_from_array4f(fluid->h_vel, partner);
    vec3f gpos = get_vec3f_from_array4f(grain->h_pos, owner);
    vec3f fpos = get_vec3f_from_array4f(fluid->h_pos, partner);
    double gmass = grain->h_mass[owner];
    double fmass = fluid->h_mass[partner];
    double gradius = grain->h_radius[owner];
    double fradius = fluid->h_radius[partner];

    double radius = gradius + fradius;
    vec3f pos, vel;
    vecSub (pos, fpos, gpos);
    vecSub (vel, fvel, gvel);

    /* testing */
    vec3f p1, p2, diffp;
    vec3f p1p, p2p;
    vecScale (p1, gvel, gmass);
    vecScale (p2, fvel, fmass);
    vec3f pAntes;
    vecAdd (pAntes, p1, p2);

    double e1, e2, diffe, eBefore;
    e1 = gmass * vecDot (gvel, gvel);
    e2 = fmass * vecDot (fvel, fvel);
    eBefore = e1 + e2;



    double pu = vecMod (pos);
    
    vec3f n;
    vecScale (n, pos, (1/pu));
   
    g_message ("radius: %f, pu: %f, pu-radius: %.9f\npos: (%f, %f, %f), n: (%f, %f, %f), vecDot (n, n) = %f", 
                radius, pu, pu-radius, pos.x, pos.y, pos.z, n.x, n.y, n.z, vecDot (n, n));
    if (fabs (pu - radius) > 10e-6)
        {
            g_error ("Bad Prediction! - PU: %f != Radius: %f", pu, radius);
            return ;
        }



    double VN = vecDot (n, vel);

    double dP = (2 * gmass * fmass * VN) / (gmass + fmass);
    vec3f j;
    
    j.x = dP * n.x;
    j.y = dP * n.y;
    j.z = dP * n.z;
   
    /* update states */
    gvel.x = gvel.x + j.x / gmass;
    gvel.y = gvel.y + j.y / gmass;
    gvel.z = gvel.z + j.z / gmass;

    fvel.x = fvel.x - j.x / fmass;
    fvel.y = fvel.y - j.y / fmass;
    fvel.z = fvel.z - j.z / fmass;
    
    /* in cpu */
    array4_set_vec3 (grain->h_vel, owner, gvel);
    array4_set_vec3 (fluid->h_vel, partner, fvel);
    grain->h_nCol[owner]++;
    fluid->h_nCol[partner]++;

    /* in gpu */
    ADENTU_CUDA_MEMCPY_H2D (array4_get_ptr_at (grain->d_vel, owner), 
                            array4_get_ptr_at (grain->h_vel, owner),
                            4 * sizeof (adentu_real));
    
    ADENTU_CUDA_MEMCPY_H2D (array4_get_ptr_at (fluid->d_vel, partner), 
                            array4_get_ptr_at (fluid->h_vel, partner),
                            4 * sizeof (adentu_real));

    ADENTU_CUDA_MEMCPY_H2D ((grain->d_nCol + owner),
                            (grain->h_nCol + owner),
                            sizeof (int));

    ADENTU_CUDA_MEMCPY_H2D ((fluid->d_nCol + partner),
                            (fluid->h_nCol + partner),
                            sizeof (int));
    /* if the above copy code doesn't work, try this: */
    /* 
    int gn = grain->n;
    int fn = fluid->n;
    ADENTU_CUDA_MEMCPY_H2D (grain->d_vel, grain->h_vel, 
                            4 * gn * sizeof (adentu_real));
    ADENTU_CUDA_MEMCPY_H2D (fluid->d_vel, fluid->h_vel, 
                            4 * fn * sizeof (adentu_real));
    ADENTU_CUDA_MEMCPY_H2D (grain->d_nCol, grain->h_nCol, 
                            gn * sizeof (int));
    ADENTU_CUDA_MEMCPY_H2D (fluid->d_nCol, fluid->h_nCol, 
                            fn * sizeof (int));
    */   



    /* testing */
    vecScale (p1p, gvel, gmass);
    vecScale (p2p, fvel, fmass);
    vec3f pDespues;
    vecAdd (pDespues, p1p, p2p);

    vecSub (diffp, pDespues, pAntes);

    double e1p, e2p, eAfter;
    e1p = gmass * vecDot (gvel, gvel);
    e2p = fmass * vecDot (fvel, fvel);
    eAfter = e1p + e2p;
    diffe = eAfter - eBefore;
    
    g_message ("Momentum vecMod(diffp): %.8f, Kinetic Energy diffe: %.8f", vecMod(diffp), diffe);
}


