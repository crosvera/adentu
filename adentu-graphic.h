/*
    Carlos Ríos Vera <crosvera@gmail.com>
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
