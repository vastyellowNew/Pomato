/*************************************************************************
POMATO: POincare MAp TOpology: Extraction and visualization of the
        topological structure of the circular restricted three-body 
        problem for orbital mechanics applications. 

Authors: Wayne Schlei and Xavier Tricoche 

Copyright (c) 2013-2018, Purdue University 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**************************************************************************/


#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <iomanip>
#include <functional>

#include <boost/format.hpp>
#include <boost/limits.hpp>
#include <math/rational.hpp>

#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>

// display
#include <graphics/colors.hpp>
#include <graphics/GLUT_helper.hpp>
#include <graphics/GUI/GLUI_Wrapper.h>

using namespace xavier;

// -------------------------
//
//              UI
//
// -------------------------
int     width, height;
int     res[2];
double  scale = 1;

std::string filename;

// -------------------------
//
//          Display
//
// -------------------------

int main_window;
nvis::vec2 wc;
void keyboard(unsigned char key, int x, int y);
void display();

void idle(void)
{
    // switch context back to main window after GLUI activity
    glutSetWindow(main_window);
}

double current_size()
{
    nvis::vec3 min = GLUT_helper::world_coordinates(0,0);
    nvis::vec3 d = GLUT_helper::world_coordinates(width, 0) - min;
    d[0] = fabs(d[0]);
    d[1] = fabs(d[1]);
    d[2] = fabs(d[2]);
    return *std::max_element(&d[0], &d[3]);
}

void display_matrix(const GLfloat m[16])
{
    std::cerr << "(" << m[0] << ", " << m[1] << ", " << m[2] << ", " << m[3] << ")\n"
              << "(" << m[4] << ", " << m[5] << ", " << m[6] << ", " << m[7] << ")\n"
              << "(" << m[8] << ", " << m[9] << ", " << m[10] << ", " << m[11] << ")\n"
              << "(" << m[12] << ", " << m[13] << ", " << m[14] << ", " << m[15] << ")\n";
}

std::vector< std::pair<nvis::vec2, nvis::vec2> > vectors;

nvis::vec2 _last;

void draw(void)
{
    static int count = 0;
    
    glDisable(GL_DEPTH_TEST);
    
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    
    // draw vector associated with each vertex
    GLUT_helper::draw_vectors(vectors, nvis::fvec3(1,0,0));
}

void display(void)
{
    GLUT_helper::setup_display(draw);
    glutSwapBuffers();
}

void guiCallback(int)
{
    display();
}

bool save_to_file;
void mouse(int button, int state, int x, int y)
{
    std::cerr << "mouse, button = " << button << ", state = " << state
              << ", x = " << x << ", y = " << y << std::endl;
              
    // standard_map pmap(_k);
    
    nvis::vec3 _wc = GLUT_helper::world_coordinates(x, y);
    wc = nvis::vec2(_wc[0], _wc[1]);
    
    if (nvis::all(wc == _last)) {
        return;
    } else {
        _last = wc;
    }
    
    std::ostringstream os;
    os << "Cursor at " << wc;
    glutSetWindowTitle(os.str().c_str());
    // if (button != GLUT_MIDDLE_BUTTON)
    GLUT_helper::glut_helper_mouse(button, state, x, y);
    display();
}

bool arrow_pressed = false;
int __x0, __y0;
void keySpecial(int key, int x, int y)
{

    if (!arrow_pressed) {
        arrow_pressed = true;
        __x0 = x;
        __y0 = y;
        mouse(GLUT_MIDDLE_BUTTON, GLUT_DOWN, __x0, __y0);
    } else {
        switch (key) {
            case GLUT_KEY_LEFT: {
                __x0 -= 1;
                break;
            }
            case GLUT_KEY_RIGHT: {
                __x0 += 1;
                break;
            }
            case GLUT_KEY_UP: {
                __y0 += 1;
                break;
            }
            case GLUT_KEY_DOWN: {
                __y0 -= 1;
            }
            default: {
            }
        }
        nvis::vec3 _wc = GLUT_helper::world_coordinates(__x0, __y0);
        wc = nvis::vec2(_wc[0], _wc[1]);
        if (nvis::all(wc == _last)) {
            return;
        } else {
            _last = wc;
        }
        std::ostringstream os;
        os << "Cursor at " << wc;
        glutSetWindowTitle(os.str().c_str());
        
        GLUT_helper::glut_helper_motion(__x0, __y0);
    }
}

void keySpecialUp(int key, int x, int y)
{
    arrow_pressed = false;
    mouse(GLUT_MIDDLE_BUTTON, GLUT_UP, __x0, __y0);
}

void keyboard(unsigned char key, int x, int y)
{
    static double sensitivity = 1000;
    if(key == 'r') {
        GLUT_helper::resetCamera();
    }
    glutPostRedisplay();
}
// --------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    _last = nvis::vec2(-1,-1);
    
    std::string arg;
    for (int i=1 ; i<argc ; ++i) {
        arg = argv[i];
        if (arg == "-i") {
            filename = argv[++i];
        } else if (arg == "-s") {
            scale = atof(argv[++i]);
        }
    }
    
    std::fstream in(filename.c_str(), std::ios::in);
    nvis::vec2 x, v(0);
    nvis::bbox2 bounds;
    while (!in.eof()) {
        char c;
        in >> x[0] >> c >> x[1] >> c >> v[0] >> c >> v[1];
        if (!nvis::norm(v)) {
            break;
        }
        vectors.push_back(std::make_pair(x, scale*v));
        bounds.add(x);
        bounds.add(x+scale*v);
        v[0] = v[1] = 0;
    }
    in.close();
    
    width = 512;
    height = 512;
    
    // initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayString("samples rgba double alpha");
    glutInitWindowSize(width, height);
    glutInitWindowPosition(20, 20);
    main_window = glutCreateWindow(argv[0]);
    
    // configure OpenGL for aesthetic results
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH, GL_NICEST);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_POINT_SMOOTH, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    
    GLUT_helper::box = bounds;
    
    // set GLUT callback functions
    glutDisplayFunc(display);
    glutReshapeFunc(GLUT_helper::glut_helper_reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(GLUT_helper::glut_helper_motion);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    glutSpecialFunc(keySpecial);
    glutSpecialUpFunc(keySpecialUp);
    
    GLUT_helper::resetCamera();
    display();
    glutPostRedisplay();
    
    // adjust camera to mesh
    GLUT_helper::update_panning_sensitivity(100);
    
    // Enter GLUT event loop
    glutMainLoop();
    
    return 0;
}
