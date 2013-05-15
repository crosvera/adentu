/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#ifndef __ADENTU_EVENT_H__
#define __ADENTU_EVENT_H__

#include <glib.h>

#include "vec3.h"
#include "adentu-model.h"


typedef enum {
    ADENTU_EVENT_START,
    ADENTU_EVENT_MPC,   // MPC events
    ADENTU_EVENT_BC_GRAIN,    // Boundary Condition
    ADENTU_EVENT_BC_FLUID,    // Boundary Condition
    ADENTU_EVENT_GGC,   // Grain-Grain Collision
    ADENTU_EVENT_GFC,   // Grain-Fluid Collision
    ADENTU_EVENT_END
} AdentuEventType;


typedef struct _AdentuEvent {
    AdentuEventType type;
    double time;

    int owner;
    int partner;
    int nEvents;
    void *eventData;
} AdentuEvent;


typedef struct _AdentuEventHandler {
    GSList      *(*event_init)      (AdentuModel *, GSList *);
    int         (*event_is_valid)   (AdentuModel *, AdentuEvent *);
    void        (*event_attend)     (AdentuModel *, AdentuEvent *);
    AdentuEvent *(*event_get_next)  (AdentuModel *);
} AdentuEventHandler;


extern const char *AdentuEventTypeStr[];


AdentuEvent *adentu_event_get_next (GSList **eList);
GSList *adentu_event_schedule (GSList *eList, AdentuEvent *event);
int adentu_event_compare (gconstpointer a, gconstpointer b);


GSList * adentu_event_init (GSList *eList,
                            AdentuEventHandler *handler[],
                            AdentuModel *model);

GSList *adentu_event_loop (GSList *eList,
                        AdentuEventHandler *handler[],
                        AdentuModel *model);



#endif /* __ADENTU_EVENT_H__ */
