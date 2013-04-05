/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#ifndef __ADENTU_EVENT_H__
#define __ADENTU_EVENT_H__

#include <glib.h>

#include "vec3.h"

typedef enum {
    ADENTU_EVENT_START,
    ADENTU_EVENT_END,
    ADENTU_EVENT_MPC,   // MPC events
    ADENTU_EVENT_BC,    // Boundary Condition
    ADENTU_EVENT_CC,    // Cell Crossing
    ADENTU_EVENT_GGC,   // Grain-Grain Collision
    ADENTU_EVENT_GFC,   // Grain-Fluid Collision
    ADENTU_EVENT_USER,  // User-defined event
    ADENTU_EVENT_DUMMY  // Dummy event
} AdentuEventType;


typedef struct _AdentuEvent {
    AdentuEventType *type;
    double *time;

    int *object;
    int *partner;
    double *lastTimePartner;

    int nEvents;

} AdentuEvent;



#endif /* __ADENTU_EVENT_H__ */
