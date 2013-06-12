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

#include "adentu-event-gfc.h"
#include "adentu-event-gfc-cuda.h"
#include "adentu-event-mpc-cuda.h"


AdentuEventHandler AdentuGFCEventHandler = {adentu_event_gfc_init,
                                            adentu_event_gfc_is_valid,
                                            adentu_event_gfc_attend2,
                                            adentu_event_gfc_get_next};



GSList *adentu_event_gfc_init (AdentuModel *model,
                               GSList *eList)
{
    eList = adentu_event_schedule (eList, adentu_event_gfc_get_next (model));
    //exit (1);


    return eList;
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


void adentu_event_gfc_attend2 (AdentuModel *model,
                               AdentuEvent *event)
{
    AdentuAtom *grain = model->grain;
    AdentuAtom *fluid = model->fluid;

    double dT = event->time - model->elapsedTime;

    adentu_event_mpc_cuda_integrate (grain,
                                     model->gGrid,
                                     model->accel, dT);
    adentu_grid_set_atoms (model->gGrid, 
                           grain, 
                           &model->bCond);


    adentu_event_mpc_cuda_integrate (fluid,
                                     model->fGrid,
                                     model->accel, dT);
    adentu_grid_set_atoms (model->fGrid, 
                           fluid, 
                           &model->bCond);


    int owner = event->owner;
    int partner = event->partner;
    
    vec3f *gvel = &grain->vel[owner];
    vec3f *fvel = &fluid->vel[partner];
    vec3f *gpos = &grain->pos[owner];
    vec3f *fpos = &fluid->pos[partner];
    double gmass = grain->mass[owner];
    double fmass = fluid->mass[partner];
    double gradius = grain->radius[owner];
    double fradius = fluid->radius[partner];

    double radius = gradius + fradius;
    vec3f pos, vel;
    vecSub (pos, *fpos, *gpos);
    vecSub (vel, *fvel, *gvel);

    double pu = vecMod (pos);
    
    if (fabs (pu - radius) > 10e-6)
        {
            printf ("Bad Prediction! - PU: %f != Radius: %f\n", pu, radius);
            //return ;
        }

    vec3f n;
    vecScale (n, pos, pu);
   
    printf ("radius: %f, pu: %f, pos: (%f, %f, %f), n: (%f, %f, %f), vecDot (n, n) = %f\n", 
                radius, pu, pos.x, pos.y, pos.z, n.x, n.y, n.z, vecDot (n, n));

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
    fluid->nCol[partner]++;
}


void adentu_event_gfc_attend (AdentuModel *model,
                              AdentuEvent *event)
{
    AdentuAtom *grain = model->grain;
    AdentuAtom *fluid = model->fluid;

    double dT = event->time - model->elapsedTime;

    adentu_event_mpc_cuda_integrate (grain,
                                     model->gGrid,
                                     model->accel, dT);

    adentu_event_mpc_cuda_integrate (fluid,
                                     model->fGrid,
                                     model->accel, dT);


    adentu_grid_set_atoms (model->gGrid, 
                           grain, 
                           &model->bCond);
    adentu_grid_set_atoms (model->fGrid, 
                           fluid, 
                           &model->bCond);


    int owner = event->owner;
    int partner = event->partner;
    
    vec3f *gvel = &grain->vel[owner];
    vec3f *fvel = &fluid->vel[partner];
    vec3f *gpos = &grain->pos[owner];
    vec3f *fpos = &fluid->pos[partner];
    double gmass = grain->mass[owner];
    double fmass = fluid->mass[partner];



    vec3f pBefore, eBefore;
    vec3f pAfter, eAfter;
    vec3f dpos, dvel, n, dp;
    double eBeforeTotal, eAfterTotal, mod, dpAux;

    pBefore.x = gmass * gvel->x + fmass * fvel->x;
    pBefore.y = gmass * gvel->y + fmass * fvel->y;
    pBefore.z = gmass * gvel->z + fmass * fvel->z;

    eBefore.x = 0.5 * gmass * (gvel->x * gvel->x) +
                0.5 * fmass * (fvel->x * fvel->x);
    eBefore.y = 0.5 * gmass * (gvel->y * gvel->y) +
                0.5 * fmass * (fvel->y * fvel->y);
    eBefore.x = 0.5 * gmass * (gvel->z * gvel->z) +
                0.5 * fmass * (fvel->z * fvel->z);

    eBeforeTotal = eBefore.x + eBefore.y + eBefore.z;

    vecSub (dpos, *fpos, *gpos);
    n = dpos;
    
    mod = vecMod (n);
    n.x /= mod;
    n.y /= mod;
    n.z /= mod;
    
    vecSub (dvel, *fvel, *gvel);

    dpAux = -2.0 * ((gmass * fmass) / (gmass + fmass)) * vecDot (dvel, n);
    dp.x = dpAux * n.x;
    dp.y = dpAux * n.y;
    dp.z = dpAux * n.z;

    gvel->x -= (dp.x / gmass);
    gvel->y -= (dp.y / gmass);
    gvel->z -= (dp.z / gmass);
    fvel->x -= (dp.x / fmass);
    fvel->y -= (dp.y / fmass);
    fvel->z -= (dp.z / fmass);

    pAfter.x = gmass * gvel->x + fmass * fvel->x;
    pAfter.y = gmass * gvel->y + fmass * fvel->y;
    pAfter.z = gmass * gvel->z + fmass * fvel->z;

    eAfter.x = 0.5 * gmass * (gvel->x * gvel->x) + 0.5 * fmass * (fvel->x * fvel->x);
    eAfter.y = 0.5 * gmass * (gvel->y * gvel->y) + 0.5 * fmass * (fvel->y * fvel->y);
    eAfter.z = 0.5 * gmass * (gvel->z * gvel->z) + 0.5 * fmass * (fvel->z * fvel->z);

    eAfterTotal = eAfter.x + eAfter.y + eAfter.z;

    vecSub (dvel, *fvel, *gvel);
    vecSub (dpos, *fpos, *gpos);

    grain->nCol[owner]++;
    fluid->nCol[partner]++;

}

