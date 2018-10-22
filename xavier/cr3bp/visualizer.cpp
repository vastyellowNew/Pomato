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
#include <list>
#include <vector>
#include <set>
#include <sstream>
#include <assert.h>
#include <glut.h>
#include <math.h>

struct point {
    point(double x, double y)
    {
        data[0] = x;
        data[1] = y;
    }
    point(const point& p)
    {
        data[0] = p[0];
        data[1] = p[1];
    }
    double& operator[](unsigned int i)
    {
        assert(i<2);
        return data[i];
    }
    const double& operator[](unsigned int i) const
    {
        assert(i<2);
        return data[i];
    }
    point& operator-=(const point& p)
    {
        data[0] -= p[0];
        data[1] -= p[1];
        return *this;
    }
    point& operator+=(const point& p)
    {
        data[0] += p[0];
        data[1] += p[1];
        return *this;
    }
    point& operator*=(double a)
    {
        data[0] *= a;
        data[1] *= a;
        return *this;
    }
    
    double data[2];
};

inline point operator+(const point& p0, const point& p1)
{
    point p(p0);
    p += p1;
    return p;
}

inline point operator-(const point& p0, const point& p1)
{
    point p(p0);
    p -= p1;
    return p;
}

inline point operator*(double a, const point& p)
{
    point q(p);
    q *= a;
    return q;
}

typedef point           point_type;
std::vector<point_type> points;
std::vector<int>        interior;
std::vector<int>        hull;

inline double whichSide(const point_type& p0, const point_type& p1, const point_type& p2)
{
    point_type v0 = p1-p0;
    point_type v1 = p2-p1;
    return v0[0]*v1[1] - v0[1]*v1[0];
}

struct lexico_less {
    bool operator()(int i0, int i1)
    {
        const point_type& p0 = points[i0];
        const point_type& p1 = points[i1];
        return (p0[0] < p1[0] || (p0[0] == p1[0] && p0[1] < p1[1]));
    }
};

struct lexico_more {
    bool operator()(int i0, int i1)
    {
        const point_type& p0 = points[i0];
        const point_type& p1 = points[i1];
        return (p0[0] > p1[0] || (p0[0] == p1[0] && p0[1] > p1[1]));
    }
};

void grow(std::list<int>& L, const std::vector<int>& ids)
{
    L.push_back(ids[0]);
    L.push_back(ids[1]);
    for (int i=2 ; i<ids.size() ; ++i) {
        L.push_back(ids[i]);
        while(L.size() >= 3) {
            std::list<int>::iterator i0, i1, i2;
            i0 = L.end();
            i0--;
            i2 = i0;
            i0--;
            i1 = i0;
            i0--;
            if (whichSide(points[*i0], points[*i1], points[*i2]) <= 0) {
                L.erase(i1);
            } else {
                break;
            }
        }
    }
}

int main_window;
void idle(void)
{
    glutSetWindow(main_window);
}

inline point_type ctrans(const point_type& x)
{
    return point_type(-0.9,-0.9) + 1.8*x;
}

inline void draw_circle(const point_type& c, double r, double dz=0)
{
    float dt = M_PI/20;
    glBegin( GL_TRIANGLE_FAN );
    glVertex3f(c[0], c[1], dz);
    for( float t = 0; t <=2*M_PI+dt; t += dt ) {
        glVertex3f( c[0] + r*cos(t), c[1] + r*sin(t), dz );
    }
    glEnd();
}

inline void print_char(const point_type& x, const std::string& s)
{
    glRasterPos2f(x[0], x[1]);
    for (int i = 0; i < s.size(); i++) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, s.c_str()[i]);
    }
}

const point_type shift(0.02, 0.02);
void display()
{
    glDisable(GL_DEPTH_TEST);
    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glColor3f(0,0,0);
    for (int i=0 ; i<interior.size() ; ++i) {
        draw_circle(ctrans(points[interior[i]]), 0.02, 0);
        std::ostringstream os;
        os << interior[i];
        print_char(ctrans(points[interior[i]])+shift, os.str());
    }
    
    if (hull.size()) {
        glColor3f(0.5,0,0);
        glBegin(GL_LINE_STRIP);
        for (int i=0 ; i<hull.size() ; ++i) {
            point_type x = ctrans(points[hull[i]]);
            glVertex3f(x[0], x[1], 0);
        }
        point_type x = ctrans(points[hull[0]]);
        glVertex3f(x[0], x[1], 0);
        glEnd();
        
        glColor3f(1,0,0);
        for (int i=0 ; i<hull.size() ; ++i) {
            draw_circle(ctrans(points[hull[i]]), 0.02, 0);
            std::ostringstream os;
            os << hull[i];
            print_char(ctrans(points[hull[i]])+shift, os.str());
        }
    }
    
    glutSwapBuffers();
}

int main(int argc, char* argv[])
{
    int n;
    std::cin >> n;
    points.reserve(n);
    for (int i=0 ; i<n ; ++i) {
        double x, y;
        std::cin >> x >> y;
        points.push_back(point_type(x, y));
    }
    
    std::vector<int> ids(n);
    for (int i=0 ; i<n ; ++i) {
        ids[i]=i;
    }
    
    std::sort(ids.begin(), ids.end(), lexico_more());
    std::list<int> Lupper;
    grow(Lupper, ids);
    
    std::sort(ids.begin(), ids.end(), lexico_less());
    std::list<int> Llower;
    grow(Llower, ids);
    
    Llower.pop_back();
    Llower.pop_front();
    Lupper.insert(Lupper.end(), Llower.begin(), Llower.end());
    
    std::list<int>::iterator lowest = std::min_element(Lupper.begin(), Lupper.end());
    std::list<int> _back(Lupper.begin(), lowest);
    Lupper.erase(Lupper.begin(), lowest);
    Lupper.insert(Lupper.end(), _back.begin(), _back.end());
    
    for (std::list<int>::const_iterator it=Lupper.begin(); it!=Lupper.end(); ++it) {
        std::cout << *it << " (" << points[*it][0] << ", " << points[*it][1] << ")";
        std::list<int>::const_iterator tmp = it;
        if (++tmp != Lupper.end()) {
            std::cout << ", ";
        }
    }
    std::cout << '\n';
    std::set<int> __hull;
    for (std::list<int>::const_iterator it=Lupper.begin(); it!=Lupper.end(); ++it) {
        __hull.insert(*it);
    }
    std::copy(Lupper.begin(), Lupper.end(), std::back_inserter(hull));
    for (int i=0 ; i<n ; ++i) {
        if (__hull.find(i) == __hull.end()) {
            interior.push_back(i);
        }
    }
    
    glutInit(&argc, argv);
    glutInitDisplayString("samples rgba double alpha");
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(20, 20);
    main_window = glutCreateWindow(argv[0]);
    
    // configure OpenGL for aesthetic results
    glHint(GL_LINE_SMOOTH, GL_NICEST);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_POINT_SMOOTH, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    display();
    glutPostRedisplay();
    
    glutMainLoop();
    
    return 0;
}
