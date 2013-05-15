/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#ifndef __ADENTU_RUNNABLE_H__
#define __ADENTU_RUNNABLE_H__

#include <glib.h>

#include "vec3.h"
#include "adentu-model.h"
#include "adentu-event.h"

/*
typedef struct _AdentuRunnable
{
    GSList *pre_attend_func;
    GSList *post_attend_func;
} AdentuRunnable; */

typedef void (*AdentuRunnableFunc)(AdentuModel *, AdentuEvent *);

void adentu_runnable_add_pre_func (AdentuModel *model,
                                   AdentuRunnableFunc func);

void adentu_runnable_add_post_func (AdentuModel *model,
                                    AdentuRunnableFunc func);

void adentu_runnable_exec_pre_func (AdentuModel *model,
                                    AdentuEvent *event);

void adentu_runnable_exec_post_func (AdentuModel *model,
                                     AdentuEvent *event);

#endif /* __ADENTU_RUNNABLE_H__  */
