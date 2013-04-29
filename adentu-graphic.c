/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#include <stdio.h>
#include <math.h>
#include <glib.h>
#include <GL/glut.h>

#include "vec3.h"
#include "adentu-model.h"
#include "adentu-event.h"
#include "adentu-graphic.h"

void adentu_graphic_init (char **argc,
                          int argv, 
                          AdentuModel *model, 
                          GSList *eList,
                          AdentuEventHandler **handler)
{
    phi = 0.0;
    theta = 90.0;

    adentu_model = model;
    *adentu_eList = eList;
    adentu_handler = handler;

    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_SINGLE|GLUT_DOUBLE|GLUT_DEPTH);
    glutInitWindowPosition (10,10);
    glutInitWindowSize (1024, 768);
    glutCreateWindow("Adentu Simulation v0.1");
    //init();

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
}


void adentu_graphic_key (unsigned char key,
                         int x, int y)
{
    switch (key)
        {
        case 'q':
            /* exit, terminate all*/
            g_message ("Frocing end of simulation...");
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

    for (int i = 0; i < grain->n; ++i)
    {
        glPushMatrix ();
        glColor3ub (140, 140, 140);
        glTranslatef (grain->pos[i].x - half.x,
                      grain->pos[i].y - half.y,
                      grain->pos[i].z - half.z);
        glutSolidSphere (grain->radius[i], 50, 50);
        glPopMatrix ();
    }


    for (int i = 0; i < fluid->n; ++i)
    {
        glPushMatrix ();
        if (i%2)
            glColor3ub (104, 227, 215);
        else
            glColor3ub (0, 227, 215);

        glTranslatef (fluid->pos[i].x - half.x,
                      fluid->pos[i].y - half.y,
                      fluid->pos[i].z - half.z);
        glutSolidSphere (length.x * 2e-3, 20, 20);
        glPopMatrix ();
    }

    glutSwapBuffers ();
}


void adentu_graphic_event_loop (void)
{
    AdentuEvent *event = NULL;
    AdentuEventType t;
    
    event = adentu_event_get_next (adentu_eList);
    t = event->type;
    
    if (t != ADENTU_EVENT_END)
        {
            if ((*adentu_handler[t]).event_is_valid (adentu_model, event))
                {
                    (*adentu_handler[t]).event_attend (adentu_model, event);
                    adentu_model->elapsedTime += (event->time - 
                                                  adentu_model->elapsedTime);
                    *eList = adentu_event_schedule (*eList,
                            (*adentu_handler[t]).event_get_next (adentu_model));
                }
            else
                *eList = adentu_event_schedule (*eList, 
                            (*adentu_handler[t]).event_get_next (adentu_model));

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
}
