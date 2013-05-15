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
