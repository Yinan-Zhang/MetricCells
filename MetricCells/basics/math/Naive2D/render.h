//
//  render.h
//  Naive2D
//
//  Created by Yinan Zhang on 3/11/15.
//  Copyright (c) 2015 Yinan Zhang. All rights reserved.
//

#ifndef Naive2D_render_h
#define Naive2D_render_h


#if __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#elif __linux
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#include "geometry.h"
#include "polygon.h"
#include <vector>


namespace N2D {
    
    namespace render {
        
        struct Color
        {
            int R, G, B, A;
            Color( int R, int G, int B, int A = 255 ):R(std::max(std::min(R, 255), 0)), G(std::max(std::min(G, 255), 0)), B(std::max(std::min(B, 255), 0)), A(std::max(std::min(A, 255), 0)){}
        };
        
        static bool window_initialized   = false;
        static bool display_set          = false;
        
        static int WIDTH  = 0;
        static int HEIGHT = 0;
        
        /* create a window
         * @param width: the width of the window
         * @param height: the height of the window
         * @param title: title of the window
         * @param position: where you want the window to be shown in the screen.
         */
        void create_window( int width, int height, const char* title, v2 position, Color bgcolor = Color(255,255,255,255), bool animation = false )
        {
            if( window_initialized )
                return;
            int argc = 1;
            char *argv[1] = {(char*)"whatever"};
            glutInit(&argc, argv);
            animation? glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA) : glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
            glutInitWindowPosition(position.x, position.y); // Position the window's initial top-left corner
            glutInitWindowSize(width, height);   // Set the window's initial width & height
            glutCreateWindow(title);  // Create window with the given title
            WIDTH = width; HEIGHT = height;
            window_initialized = true;
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glEnable( GL_BLEND );
            glClearColor( bgcolor.R/255.0, bgcolor.G/255.0, bgcolor.B/255.0, bgcolor.A/255.0); // Black and opaque
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
        }
        
        void set_display_func( void (*func)(void) )
        {
            if(!window_initialized)
                throw "Please initialize rendering then initalize window.";
            glutDisplayFunc(func);
            display_set = true;
        }
        
        void set_keyboard_func( void (*func)(unsigned char, int, int) )
        {
            if(!window_initialized)
                throw "Please initialize rendering then initalize window.";
            glutKeyboardFunc(func);
        }
        
        void main_loop()
        {
            if(!window_initialized)
                throw "Please initialize rendering then initalize window.";
            if(!display_set)
                throw "Please set display function first";
            glutMainLoop();
        }
        
        void clean_screen()
        {
            glClear(GL_COLOR_BUFFER_BIT);   // Clear the color buffer with current clearing color
        }
        
        void flush()
        {
            glFlush();
        }
        
        // Scale the world
        void scale_world( double x, double y )
        {
            glScaled(x, y, 1);
            //glTranslated( (1.0/x-1)/2.0, (1.0/y-1)/2.0, 1);
        }
        
        // Scale the world
        void translate_world( double x, double y )
        {
            glTranslated( x, y, 1);
        }
        
        
        /* Set line with of rendering
         * @param width: the width of lines
         */
        void set_line_width(double width)
        {
            glLineWidth(width);
        }
        
        /* Render a Line segment
         * @param line: the line segment we are rendering
         * @param color: color of the polygon
         */
        void line_seg( const Line_segment& line, Color color)
        {
            glBegin(GL_LINE_STRIP);
            glColor4f(color.R/255.0f, color.G/255.0f, color.B/255.0f, color.A/255.0f);
            glVertex2f( line.start.x/WIDTH, line.start.y/HEIGHT);
            glVertex2f( line.end.x/WIDTH, line.end.y/HEIGHT );
            glEnd();
        }
        
        // Render lines that connect points[0] to points[-1]
        void lines( const std::vector<v2>& points, Color color )
        {
            glBegin(GL_LINES);
            glColor4f(color.R/255.0f, color.G/255.0f, color.B/255.0f, color.A/255.0f);
            for (int i = 0; i < points.size(); i++) {
                glVertex2f( points[i].x/WIDTH, points[i].y/HEIGHT );
            }
            glEnd();
        }
        
        // Render lines.
        void lines( const std::vector<Line_segment>& lines, Color color )
        {
            for (Line_segment line : lines) {
                line_seg(line, color);
            }
        }
        
        /* Render a polygon
         * @param polygon: the polygon we are rendering
         * @param color: color of the polygon
         * @param fill: if fill in the polygon.
         */
        void polygon( const Polygon& polygon, Color color, bool fill = true )
        {
            if( fill )  glBegin ( GL_POLYGON );
            else glBegin(GL_LINE_LOOP);
            glColor4f(color.R/255.0f, color.G/255.0f, color.B/255.0f, color.A/255.0f);
            for( int i = 0; i < polygon.vertices.size(); i++ )
            {
                v2 vert = polygon.vertices[i];
                glVertex2f( vert.x/WIDTH, vert.y/HEIGHT);
            }
            glEnd();
        }
        
        
        /* Render all polygons */
        void polygons( std::vector<Polygon>& polygons, Color color = Color(100, 100, 100, 150), bool fill = true )
        {
            for( int i = 0; i < polygons.size(); i++  )
            {
                Polygon poly = polygons[i];
                polygon(poly, color, fill);
            }
        }
        
        /* Render a sphere
         * @param sphere: the sphere we are rendering
         * @param color: color of the polygon
         * @param fill: if fill in the polygon.
         */
        void sphere(const sphere& sphere, Color color, bool fill = true)
        {
            double PI = 3.14159265f;
            v2 center = sphere.center();
            if( fill )
            {
                glBegin(GL_TRIANGLE_FAN);
                glColor4f(color.R/255.0f, color.G/255.0f, color.B/255.0f, color.A/255.0f);
            }
            else
            {
                glBegin(GL_LINE_LOOP);
                glColor4f(color.R/255.0f, color.G/255.0f, color.B/255.0f, color.A/255.0f);
            }
            if (sphere.metric == SPHEREMETRIC::L1) {
                double hr = sphere.radius();
                v2 p1( (center.x-hr), (center.y) );
                v2 p2( (center.x), (center.y+hr) );
                v2 p3( (center.x+hr), (center.y) );
                v2 p4( (center.x), (center.y-hr) );
                glVertex2f((p1.x/WIDTH), (p1.y/WIDTH));
                glVertex2f((p2.x/WIDTH), (p2.y/WIDTH));
                glVertex2f((p3.x/WIDTH), (p3.y/WIDTH));
                glVertex2f((p4.x/WIDTH), (p4.y/WIDTH));
            }
            else if (sphere.metric == SPHEREMETRIC::L2) {
                int numSegments = 50;
                GLfloat angle;
                for (int i = 0; i <= numSegments; i++) { // Last vertex same as first vertex
                    angle = i * 2.0f * PI / numSegments;  // 360 deg for all segments
                    v2 temp( center.x + cos(angle)*sphere.radius(), center.y + sin(angle)*sphere.radius() );
                    float x = (temp.x)/(WIDTH);
                    float y = (temp.y)/(HEIGHT);
                    glVertex2f( x, y );
                }
            }
            else // l-infinity
            {
                double hr = sphere.radius();
                v2 p1( (center.x-hr), (center.y+hr) );
                v2 p2( (center.x+hr), (center.y+hr) );
                v2 p3( (center.x+hr), (center.y-hr) );
                v2 p4( (center.x-hr), (center.y-hr) );
                glVertex2f((p1.x/WIDTH), (p1.y/WIDTH));
                glVertex2f((p2.x/WIDTH), (p2.y/WIDTH));
                glVertex2f((p3.x/WIDTH), (p3.y/WIDTH));
                glVertex2f((p4.x/WIDTH), (p4.y/WIDTH));
            }
            glEnd();
        }
        
        
        /* Render all spheres */
        void spheres( std::vector<struct sphere>& spheres, Color color = Color(100, 100, 100, 150), bool fill = true )
        {
            for( int i = 0; i < spheres.size(); i++  )
            {
                struct sphere s = spheres[i];
                sphere(s, color, fill);
            }
        }
    }
}

#endif
