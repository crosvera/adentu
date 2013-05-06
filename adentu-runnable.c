/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#include <stdio.h>
#include <stdlib.h>
#include <glib.h>

#include "adentu-runnable.h"

#include "vec3.h"
#include "adentu-atom.h"
#include "adentu-grid.h"
#include "adentu-model.h"

void adentu_runnable_add_pre_func (AdentuModel *model,
                                   AdentuRunnableFunc func)
{
    model->pre_event_func = g_slist_append (model->pre_event_func,
                                            (void *)func);
}

void adentu_runnable_add_post_func (AdentuModel *model,
                                    AdentuRunnableFunc func)
{
    model->post_event_func = g_slist_append (model->post_event_func,
                                             (void *)func);
}

void adentu_runnable_exec_pre_func (AdentuModel *model,
                                    AdentuEvent *event)
{
    GSList *list = model->pre_event_func;
    AdentuRunnableFunc func = NULL;
    while (list)
        {
            GSList *next = list->next;
            func = list->data;
            (*func) (model, event);
            list = next;
        }
}

void adentu_runnable_exec_post_func (AdentuModel *model,
                                     AdentuEvent *event)
{
    GSList *list = model->post_event_func;
    AdentuRunnableFunc func = NULL;
    while (list)
        {
            GSList *next = list->next;
            func = list->data;
            (*func) (model, event);
            list = next;
        }
}
