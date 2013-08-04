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

#ifndef __ADENTU_H__
#define __ADENTU_H__

#include <stdlib.h>


#ifndef ADENTU_VERSION
#define ADENTU_VERSION  "0.2"
#endif


/* Seed for random numbers */
extern unsigned int adentu_srand;
#define ADENTU_SET_SRAND(SEED) \
    srand (SEED); srand48 (SEED); \
    adentu_srand = SEED

/* Adentu core libs */
#include "adentu-model.h"
#include "adentu-atom.h"
#include "adentu-grid.h"
#include "adentu-event.h"
//#include "adentu-neighbourhood.h"
#include "adentu-runnable.h"
#include "adentu-types.h"

/* Graphics */
#include "adentu-graphic.h"

#endif /* __ADENTU_h__ */
