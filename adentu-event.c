/*
 * Carlos Ríos Vera <crosvera@gmail.com>
 */

#include <stdlib.h>
#include <glib.h>

#include "adentu-event.h"
#include "adentu-grid.h"

void adentu_event_set_init_events (AdentuEvent *events, AdentuGrid *grid)
{
    events = calloc (grid->tCell, sizeof (AdentuEvent));
}
