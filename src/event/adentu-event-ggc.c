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

#include "event/adentu-event-ggc.h"
#include "event/adentu-event-ggc-cuda.h"

#include "adentu-graphic.h"

char ADENTU_EVENT_GGC[] = "ADENTU_EVENT_GGC";

AdentuEventHandler AdentuGGCEventHandler = {ADENTU_EVENT_GGC,
                                            adentu_event_ggc_init,
                                            adentu_event_ggc_is_valid,
                                            adentu_event_ggc_cuda_attend,
                                            adentu_event_ggc_get_next};



GSList *adentu_event_ggc_init (AdentuModel *model)
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

    /* testing */
    //g_message ("New GGC event, time: %f", ev->time);
    /* testing, get the backtrace of caller functions */
#ifdef Debug
    void *buffer[100];
    char **strings;
    int nptrs = backtrace (buffer, 100);
    strings = backtrace_symbols (buffer, nptrs);
    for (int i = 0; i < 3; ++i)
        g_message ("%s", strings[i]);
    free (strings);
#endif 
    //getchar ();

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

    if (grain->h_nCol[owner] == nEvents &&
        grain->h_nCol[partner] == eventData)
            return 1;

    return 0;
}




