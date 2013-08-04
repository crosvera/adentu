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
#include <unistd.h>
#include <math.h>
#include <glib.h>
#include <string.h>
#include <GLUT/glut.h>

#include "adentu-types.h"
#include "adentu-model.h"
#include "adentu-event.h"
#include "adentu-graphic.h"
#include "adentu-runnable.h"


double phi = 0.0;
double theta = 90.0;
useconds_t _graphic_pause = 10000;
AdentuModel *adentu_model = NULL;
AdentuEventHandler **adentu_handler = NULL;
int adentu_graphic_is_set = 0;



void adentu_graphic_init (int argc,
                          char **argv, 
                          AdentuModel *model, 
                          AdentuEventHandler **handler)
{
    phi = 0.0;
    theta = 90.0;

    adentu_model = model;
    adentu_handler = handler;

    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_SINGLE|GLUT_DOUBLE|GLUT_DEPTH);
    glutInitWindowPosition (10,10);
    glutInitWindowSize (1024, 768);
    glutCreateWindow("Adentu Simulation v0.1");

    glClearColor(0.0, 0.0, 0.0, 1.0); // Black Background

    glShadeModel(GL_SMOOTH);

    GLfloat lightdiffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat lightposition[] = {1.0, 1.0, 1.0, 0.0};
    GLfloat lightambient[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat lightspecular[] = { 1.0, 1.0, 1.0, 1.0 };


    // Light 1
    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightdiffuse);
    glLightfv(GL_LIGHT1, GL_POSITION, lightposition);
    glLightfv(GL_LIGHT1, GL_AMBIENT, lightambient);

    // Material
    glMaterialfv(GL_FRONT, GL_DIFFUSE, lightdiffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, lightspecular);

    glEnable(GL_LIGHT1);    // Enable Light 1
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    

    glutDisplayFunc (adentu_graphic_display);
    glutReshapeFunc (adentu_graphic_reshape);
    glutIdleFunc (adentu_graphic_event_loop);
    glutKeyboardFunc (adentu_graphic_key);
    glutSpecialFunc (adentu_graphic_special_key);
    //glutMainLoop();

    adentu_graphic_is_set = 1;
}


void adentu_graphic_key (unsigned char key,
                         int x, int y)
{
    switch (key)
        {
        case 'q':
            /* exit, terminate all*/
            g_message ("Forcing end of simulation...");
            exit (0);
            break ;

        default:
            break ;
        }
}


void adentu_graphic_special_key (int key,
                                 int x, int y)
{
    switch (key)
        {
        case GLUT_KEY_LEFT:
            phi -= 1.0;
            break ;

        case GLUT_KEY_UP:
            if (theta > 2.0)
                theta -= 1.0;
            break ;

        case GLUT_KEY_RIGHT:
            phi += 1.0;
            break ;

        case GLUT_KEY_DOWN:
            if (theta < 179.0)
                theta += 1.0;
            break ;

        default:
            break;
        }
}


void adentu_graphic_reshape (int width, int height)
{
    glViewport (0, 0, width, height);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluPerspective (45.0, 1.0, 1.0, 1000.0);    
    glMatrixMode (GL_MODELVIEW);   
}


void adentu_graphic_display (void)
{
    AdentuAtom *grain = adentu_model->grain;
    AdentuAtom *fluid = adentu_model->fluid;

    vec3f origin, length, half, center;
    origin = adentu_model->gGrid->origin;
    length = adentu_model->gGrid->length;

    vecScale (half, length , 0.5);
    center.x = origin.x + half.x;
    center.y = origin.y + half.y;
    center.z = origin.z + half.z;
    
    glClear (GL_COLOR_BUFFER_BIT);


    //camera
    glLoadIdentity ();
    gluLookAt (half.x * 5 * sin (M_PI*theta / 180.0) * sin (M_PI*phi / 180.0),
               half.y * 5 * cos (M_PI*theta / 180.0),
               half.z * 5 * sin (M_PI*theta / 180.0) * cos (M_PI*phi / 180.0),
               0, 0, 0, 0, 1, 0);

    glTranslatef (origin.x, origin.y, origin.z);
    glColor3ub (255, 255, 255);
    glutWireCube (length.x);

    vec3f pos;
    for (int i = 0; i < grain->n; ++i)
    {
        pos = get_vec3f_from_array4f (grain->h_pos, i);
        glPushMatrix ();
        if (i%2)
            glColor3ub (140, 140, 140);
        else
            glColor3ub (254, 140, 100);
        glTranslatef (pos.x - half.x,
                      pos.y - half.y,
                      pos.z - half.z);
        glutSolidSphere (grain->radius[i], 50, 50);
        glPopMatrix ();
    }


    for (int i = 0; i < fluid->n; ++i)
    {
        pos = get_vec3f_from_array4f (fluid->h_pos, i);
        glPushMatrix ();
        if (i%2)
            glColor3ub (104, 227, 215);
        else
            glColor3ub (0, 227, 215);

        glTranslatef (pos.x - half.x,
                      pos.y - half.y,
                      pos.z - half.z);
        glutSolidSphere (0.01, 20, 20);
        glPopMatrix ();
    }






    glutSwapBuffers ();
}


void adentu_graphic_event_loop (void)
{
    AdentuEvent *event = NULL;
    char *t;
    AdentuEventHandler *_handler;
  
    /* if (adentu_model->elapsedTime == 0)
        getchar ();
    */

    event = adentu_event_get_next (&(adentu_model->eList));
    t = event->type;
    /* testing */
    /* g_message ("%s: Next event to attend: %s, time: %f", __FILE__, 
                AdentuEventTypeStr[t], event->time);
    */

    
    if (t != ADENTU_EVENT_END)
        {
            handler = adentu_event_get_handler (adentu_handler, t);
            if (_handler == NULL)
                {
                    free (event->eventData);
                    free (event);
                }
            else
            if ((*_handler).event_is_valid (adentu_model, event))
                {
                    adentu_runnable_exec_pre_func (adentu_model, event);
                    (*_handler).event_attend (adentu_model, event);

                    adentu_model->elapsedTime = event->time;
                    
                    adentu_runnable_exec_post_func (adentu_model, event);
                }

            /* predict new events */
            for (int i = 0; adentu_handler[i] != NULL; ++i) 
                {
                    adentu_model->eList = adentu_event_schedule (adentu_model->eList,
                               (*adentu_handler[i]).event_get_next (adentu_model));
                }


            free (event->eventData);
            free (event);
        }
    else
        {
            g_message ("END_EVENT reached, terminating Adentu Simulation Loop.");
            free (event);
            exit (0);
        }


    glutPostRedisplay();

    usleep (_graphic_pause);
}

void adentu_graphic_start (void)
{
    glutMainLoop ();
}

void adentu_graphic_set_time_sleep (useconds_t useconds)
{
    _graphic_pause = useconds;
}

