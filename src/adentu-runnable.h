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

#ifndef __ADENTU_RUNNABLE_H__
#define __ADENTU_RUNNABLE_H__

#include <glib.h>

#include "adentu-event.h"
#include "adentu-model.h"

typedef void (*AdentuRunnableFunc)(const AdentuModel *, const AdentuEvent *);

void adentu_runnable_add_pre_func (AdentuModel *model,
                                   AdentuRunnableFunc func);

void adentu_runnable_add_post_func (AdentuModel *model,
                                    AdentuRunnableFunc func);

void adentu_runnable_exec_pre_func (AdentuModel *model,
                                    AdentuEvent *event);

void adentu_runnable_exec_post_func (AdentuModel *model,
                                     AdentuEvent *event);

#endif /* __ADENTU_RUNNABLE_H__  */
