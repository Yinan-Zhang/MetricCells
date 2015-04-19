//
//  GJK_utility.h
//  GJK
//
//  Reference:
//      http://in2gpu.com/2014/05/12/gjk-algorithm-collision-detection-2d-in-c/
//
//      http://www.codezealot.org/archives/153
//
//  Created by Yinan Zhang on 1/22/15.
//  Copyright (c) 2015 Yinan Zhang. All rights reserved.
//
//  Modified by Yu-han Lyu

#ifndef GJK_GJK_utility_h
#define GJK_GJK_utility_h
#include "geometry.h"
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
using namespace std;

namespace N2D {
    namespace GJK
    {
        constexpr double EPSILON = 1e-7;
        static inline void print( const std::vector<v2>& array )
        {
            std::cout << "vector looks like this:\n";
            for (v2 point : array) {
                std::cout << point << std::endl;
            }
            std::cout << "----------------------------------\n";
        }
        
        static inline const std::vector<v2> mink_diff( const std::vector<v2>& poly_points1, const std::vector<v2>& poly_points2 )
        {
            std::vector<v2> results;
            results.reserve(poly_points1.size()*poly_points2.size());
            for (v2 p1 : poly_points1) {
                for (v2 p2 : poly_points2) {
                    results.push_back(p1-p2);
                }
            }
            
            return results;
        }
        
        static inline v2 farest_point_in_dir( const std::vector<v2>& points, const v2& dir )
        {
            int size = (unsigned)points.size();
            int index = 0;
            double max_dot = points[0].dot(dir);
            double dot = points[1].dot(dir);
            for(unsigned int i = 2; i < size; ++i)
            {
                if(dot > max_dot)
                {
                    max_dot = dot;
                    index = i - 1;
                }
                dot = points[i].dot(dir);
            }
            return points[dot > max_dot ? size - 1 : index];
        }
        
        static inline v2 support_func( const std::vector<v2>& poly_points1, const std::vector<v2>& poly_points2, const v2& dir )
        {
            return farest_point_in_dir(poly_points1, dir) - farest_point_in_dir(poly_points2, -dir);
        }
        
        static inline v2 triple_product(const v2& a, const v2& b, const v2& c)
        {
            double dotac = a.dot(c), dotbc = b.dot(c);
            return v2{b.x * dotac - a.x * dotbc, b.y * dotac - a.y * dotbc};
        }
        
        static inline v2 closest_to_origin( const v2& a, const v2& b )
        {
            v2 ab(b - a);
            return ab * std::max(0.0, std::min(-a.dot(ab) / ab.dot(ab), 1.0)) + a;
        }
        
        static inline bool contains_origin( v2 (&simplex)[3], int& n, v2& dir )
        {
            // Triangle case
            if (n == 3)
            {
                // get b and c
                v2 ab{simplex[1] - simplex[2]};
                //direction perpendicular to AB
                dir = v2(-ab.y,ab.x);
                
                //away from C
                if(simplex[0].dot(dir) > 0.0)// if same direction, make d opposite
                    dir = -dir;
                
                //If the new vector (dir) perpenicular on AB is in the same direction with the origin (A0)
                //it means that C is the furthest from origin and remove to create a new simplex
                if(simplex[2].dot(dir) < 0.0)//same direction
                {
                    simplex[0] = simplex[2]; --n;
                    return false;
                }
                
                v2 ac{simplex[0] - simplex[2]};
                //direction to be perpendicular to AC
                dir = v2(-ac.y, ac.x);
                
                //away form B
                if(simplex[1].dot(dir) > 0.0)
                {
                    dir = -dir;
                }
                
                //If the new vector (d) perpenicular on AC edge, is in the same direction with the origin (A0)
                //it means that B is the furthest from origin and remove to create a new simplex
                if(simplex[2].dot(dir) < 0.0)
                {
                    simplex[1] = simplex[2]; --n;
                    return false;
                }
                //origin must be inside the triangle, so this is the simplex
                return true;
            }
            // line segment case
            else
            {
                // then its the line segment case
                // compute AB
                v2 ab{simplex[0] - simplex[1]};
                
                //direction perpendicular to ab, to orgin: ABXAOXAB
                dir = v2(-ab.y, ab.x);
                if(simplex[1].dot(dir) > 0.0)
                {
                    dir = -dir;
                }
            }
            return false;
        }

        
        static bool intersects( const std::vector<v2>& poly1, const std::vector<v2>& poly2 )
        {
            v2 simplex[3];
            v2 dir{1, -1};
            int count = 0;
            simplex[count++] = support_func(poly1, poly2, -dir);
            
            while (true) {
                simplex[count++] = support_func(poly1, poly2, dir);
                // make sure that the last point we added actually passed the origin
                if (simplex[count-1].dot(dir) <= 0.0)
                {
                    // if the point added last was not past the origin in the direction of d
                    // then the Minkowski Sum cannot possibly contain the origin since
                    // the last point added is on the edge of the Minkowski Difference
                    return false;
                }
                if (contains_origin( simplex, count, dir ) )//also change direction
                {
                    // if it does then we know there is a collision
                    return true;
                }
            }
        }
        
        static inline double distance( const std::vector<v2>& poly1, const std::vector<v2>& poly2 )
        {
            v2 dir{1, -1};
            v2 a{support_func(poly1, poly2, dir)};
            v2 b{support_func(poly1, poly2, -dir)};
            dir = -closest_to_origin(a, b);
            if ( dir.rsq() <= EPSILON )
                return 0.0;
            while (true) {
                v2 c{support_func(poly1, poly2, dir)};
                double sa = a.cross(b);
                double da = a.dot(dir);
                double db = b.dot(dir);
                double sb = b.cross(c);
                double sc = c.cross(a);
                double dc = c.dot(dir);
                // Test whether origin is in the triangle of abc
                if (std::min(sa * sb, sa * sc) > 0.0)
                    return 0.0;

                // Didn't make progress, c is the closest point to the origin
                if (std::min(dc - da, dc - db) <= EPSILON)
                    return std::sqrt(-dc);
                
                v2 p1{closest_to_origin(a, c)};
                v2 p2{closest_to_origin(b, c)};
                double p1_mag = p1.rsq();
                double p2_mag = p2.rsq();
                if (std::min(p1_mag, p2_mag) <= EPSILON)
                    return 0.0;
                
                if (p1_mag <= p2_mag) {
                    b = c;
                    dir = -p1;
                } else {
                    a = c;
                    dir = -p2;
                }
            }
        }       
    }
}

#endif