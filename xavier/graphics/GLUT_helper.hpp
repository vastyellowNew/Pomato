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


#ifndef __GLUT_HELPER_HPP__
#define __GLUT_HELPER_HPP__

#include <iostream>
#include <map>
#include <list>
#include <vector>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <glut.h>
#include "Camera/CameraWrapper_Misc.h"

namespace GLUT_helper {

typedef nvis::fvec3             color_type;
typedef nvis::vec2              point2d;
typedef nvis::vec3              point3d;

static float                screen_ratio;
static int                  width, height;
static nvis::bbox2          box;
static CameraWrapper        camera;

// API to set user defined callback functions
extern void set_my_display(void (*f)(void));      /* drawing instructions*/
extern void set_my_idle(void (*f)(void));         /* idle instructions */

inline double max_dim(const nvis::bbox2& abox)
{
    nvis::vec2 span = abox.size();
    return std::max(span[0], span[1]);
}

inline nvis::bbox2 current_bounds()
{
    nvis::bbox2 cur_box;
    cur_box.min()[0] = camera.getLeftClipPlane();
    cur_box.max()[0] = camera.getRightClipPlane();
    cur_box.min()[1] = camera.getBottomClipPlane();
    cur_box.max()[1] = camera.getTopClipPlane();
}

inline double current_size()
{
    return max_dim(current_bounds());
}

inline void update_panning_sensitivity(float s = 2.)
{
    nvis::bbox2 cur_box = current_bounds();
    
    float sensitivityScale = max_dim(cur_box) * s;
    
    if (sensitivityScale == 0) {
        return;
    }
    
    camera.setPanSensitivityX(sensitivityScale);
    camera.setPanSensitivityY(sensitivityScale);
}

inline nvis::vec3 world_coordinates(int x, int y)
{
    camera.getActualPixelCoordinates(x,y);
    Vector3f v = camera.unProject(x,y);
    return nvis::vec3(v(0), v(1), v(2));
}

inline void resetCamera()
{
    camera.setOrthographic(box.min()[0], box.max()[0], box.min()[1], box.max()[1],
                           camera.getNearZ(), camera.getFarZ());
                           
                           
    // I should have a function that turns off the rotate complete, but this little hack does the trick
    camera.setRotateAboutLookAtCenterSensitivityX(0);
    camera.setRotateAboutLookAtCenterSensitivityY(0);
    camera.setViewport(0, 0, width, height);
    
    glutPostRedisplay();
}

// setup_display takes as input pointer to actual drawing instructions function
inline void setup_display(void (*f)(void))
{
    camera.setProjectionMatrix();
    camera.setModelviewMatrix();
    
    f();
    
    camera.markUnchanged();
}

// following functions must be static since their address will be passed to GLUT
static void glut_helper_reshape(int w, int h)
{
    width = w;
    height = h;
    
    if (h == 0) {
        h = 1;
    }
    
    if (screen_ratio > ((float)w) / h) {
        glViewport(0, 0, w, w / screen_ratio);
    } else if (screen_ratio < ((float)w) / h) {
        glViewport(0, 0, h*screen_ratio, h);
    } else {
        glViewport(0, 0, w, h);
    }
    
    update_panning_sensitivity();
    glutPostRedisplay();
}

static void glut_helper_mouse(int button, int state, int x, int y)
{
    CameraWrapper_GUI::glutMouseCallback(&camera, button, state, x, y);
    
    // example of getting actual coordinates of pixel
    //   when third parameter is false, the function uses gluProject to determine the z "pixel" value
    //   when it is true, it reads the depth buffer to get the z value
    Vector3f pos = camera.unProject(x, y, false);
    
    glutPostRedisplay();
}

static void glut_helper_motion(int x, int y)
{
    CameraWrapper_GUI::glutMouseMotionCallback(&camera, x, y);
    glutPostRedisplay();
}

static void glut_helper_keyboard(unsigned char key, int x, int y)
{
    if (key == 'r') {
        resetCamera();
    }
    
    glutPostRedisplay();
}

inline void __draw_vector(const nvis::vec2& x, const nvis::vec2& v)
{
    const double cos_alpha = cos(M_PI/12.);
    const double sin_alpha = sin(M_PI/12.);
    
    if (nvis::norm(v) == 0) {
        return;
    }
    nvis::vec2 y = x + v;
    glVertex2f(x[0], x[1]);
    glVertex2f(y[0], y[1]);
    nvis::vec2 e0 = -0.2 * v;
    nvis::vec2 e1(-e0[1], e0[0]);
    glVertex2f(y[0], y[1]);
    nvis::vec2 z = y + cos_alpha * e0 + sin_alpha * e1;
    glVertex2f(z[0], z[1]);
    glVertex2f(y[0], y[1]);
    z = y + cos_alpha * e0 - sin_alpha * e1;
    glVertex2f(z[0], z[1]);
}

inline void draw_vector(const nvis::vec2& x, const nvis::vec2& v, const nvis::fvec3& col, float width=1)
{
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(width);
    glColor3f(col[0], col[1], col[2]);
    glBegin(GL_LINES);
    __draw_vector(x, v);
    glEnd();
}

inline void draw_vectors(const std::vector<nvis::vec2>& xs,
                         const std::vector<nvis::vec2>& vs, const nvis::fvec3& col, float width=1)
{

    assert(xs.size() == vs.size());
    
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(width);
    glColor3f(col[0], col[1], col[2]);
    glBegin(GL_LINES);
    for (int i=0 ; i<xs.size() ; ++i) {
        const nvis::vec2& x = xs[i];
        const nvis::vec2& v = vs[i];
        __draw_vector(x, v);
    }
    glEnd();
}

inline void draw_vectors(const std::vector<std::pair<nvis::vec2, nvis::vec2> >& ps, const nvis::fvec3& col, float width=1)
{
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(width);
    glColor3f(col[0], col[1], col[2]);
    glBegin(GL_LINES);
    for (int i=0 ; i<ps.size() ; ++i) {
        const nvis::vec2& x = ps[i].first;
        const nvis::vec2& v = ps[i].second;
        __draw_vector(x, v);
    }
    glEnd();
}

inline void draw_curve(const std::vector<nvis::vec2>& xs, const nvis::fvec3& col, float width=1)
{
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glColor3f(col[0], col[1], col[2]);
    glLineWidth(width);
    glBegin(GL_LINE_STRIP);
    for (int i=0 ; i<xs.size() ; ++i) {
        glVertex2f(xs[i][0], xs[i][1]);
    }
    glEnd();
}

inline void draw_quad(const nvis::bbox2& box, const nvis::fvec3& col, float width=1)
{
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glColor3f(col[0], col[1], col[2]);
    glLineWidth(width);
    const nvis::vec2& minp = box.min();
    const nvis::vec2& maxp = box.max();
    glBegin(GL_LINE_STRIP);
    glVertex2f(minp[0], minp[1]);
    glVertex2f(maxp[0], minp[1]);
    glVertex2f(maxp[0], maxp[1]);
    glVertex2f(minp[0], maxp[1]);
    glVertex2f(minp[0], minp[1]);
    glEnd();
}

inline void draw_dots(const std::vector<nvis::vec2>& xs, const nvis::fvec3& col, float sz=1)
{
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
    glPointSize(sz);
    glColor3f(col[0], col[1], col[2]);
    glBegin(GL_POINTS);
    for (int i=0 ; i<xs.size() ; ++i) {
        nvis::vec2 x = xs[i];
        glVertex2f(x[0], x[1]);
    }
    glEnd();
}

} // GLUT_helper

#endif



























