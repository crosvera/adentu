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

#ifndef __ADENTU_GRAPHIC_H__
#define __ADENTU_GRAPHIC_H__

#include <stdlib.h>
#include <glib.h>

#include "vec3.h"
#include "adentu-atom.h"
#include "adentu-grid.h"
#include "adentu-model.h"
#include "adentu-event.h"


extern double phi;
extern double theta;

extern AdentuModel *adentu_model;
extern GSList **adentu_eList;
extern AdentuEventHandler **adentu_handler;


void adentu_graphic_init (int argc,
                          char **argv,
                          AdentuModel *model, 
                          GSList **eList,
                          AdentuEventHandler **handler);


void adentu_graphic_key (unsigned char key, 
                         int x, int y);

void adentu_graphic_special_key (int key, 
                                 int x, int y);

void adentu_graphic_reshape (int width, int height);

void adentu_graphic_display (void);

void adentu_graphic_event_loop (void);







#endif /* __ADENTU_GRAPHIC_H__ */
