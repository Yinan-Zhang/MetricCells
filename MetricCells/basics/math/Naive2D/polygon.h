//
//  polygon.h
//  RSS2015
//
//  Created by Yinan Zhang on 12/25/14.
//  Copyright (c) 2014 Yinan Zhang. All rights reserved.
//

#ifndef RSS2015_polygon_h
#define RSS2015_polygon_h

#include <cassert>
#include <vector>
#include <cmath>
#include "geometry.h"

#include "GJK_utility.h"

namespace N2D {
    // Or maybe I should changed the name to convex.
    struct Polygon
    {
        std::vector<v2> vertices;
        
        // Please make sure the order of these points are clockwise.
        Polygon(std::vector<v2> points) : vertices(std::move(points)) { vertices.shrink_to_fit(); }
       
        // Please make sure the order of these points are clockwise.
        Polygon(const v2* points, int n) : vertices(points, points+n) { vertices.shrink_to_fit(); }
        
        // translate.
        void self_translate( const v2& vect )
        {
            int size = (int)vertices.size();
            for(int i = 0; i < size; i++)
            {
                vertices[i] += vect;
            }
        }
        
        // rotate self.
        void self_rotate( double dtheta, const v2& center )
        {
            double cos_dtheta = cos(dtheta);
            double sin_dtheta = sin(dtheta);
            int size = (int)vertices.size();
            for(int i = 0; i < size; i++)
            {
                v2 temp = vertices[i]-center;
                vertices[i].x = ( temp.x * cos_dtheta - temp.y * sin_dtheta ) + center.x;
                vertices[i].y = ( temp.x * sin_dtheta + temp.y * cos_dtheta ) + center.y;
            }
        }
        
        // returns true if the polygon contains the point
        bool contains(const v2& point) const
        {
            Line_segment ray( point, v2(300000, 300000) );
            int intersectTimes = 0;
            unsigned size = (unsigned)vertices.size();
            for(unsigned int i = 0; i < size; i++)
            {
                Line_segment line(this->vertices[i], this->vertices[(i + 1)% size]);
                if( line.intersects(ray) )
                    intersectTimes++;
            }
            return intersectTimes & 1;
        }
        
        // returns true if the polygon intersects with the line
        bool intersects( const Line_segment& line ) const
        {
            if( this->contains(line.start) || this->contains(line.end) )
                return true;
            unsigned size = (unsigned)vertices.size();
            for(unsigned int i = 0; i < size; i++)
            {
                Line_segment li(this->vertices[i], this->vertices[(i + 1)% size]);
                if( li.intersects(line) )
                    return true;
            }
            return false;
        }
        
        // Returns true if it intersects with another polygon.
        // (Given that self and other are both convices, if not, please use naive_intersects)
        bool intersects( const Polygon& other ) const
        {
            return GJK::intersects( this->vertices, other.vertices );
        }
        
        bool naive_intersects( const Polygon& other ) const
        {
            unsigned size = (unsigned)vertices.size();
            for(unsigned int i = 0; i < size; i++)
            {
                Line_segment line(this->vertices[i], this->vertices[(i + 1) % size]);
                if(this->intersects(line) ) return true;
            }
            return false;
        }
        
        // How much deep is a point inside the polygon?
        double penetration( const v2 pt ) const
        {
            double min = std::numeric_limits<double>::infinity();
            unsigned size = (unsigned)vertices.size();
            for(unsigned int i = 0; i < size; i++)
            {
                Line_segment li(this->vertices[i], this->vertices[(i + 1) % size]);
                //double dist = li.dist_to(pt);
                min = std::min(min, li.dist_to(pt));
            }
            return min;
        }
        
        // Returns the distance to a point
        double distance_to(const v2& pt) const
        {
            if(this->contains(pt))
                return 0.0;
            
            double min = MAX_DOUBLE;
            unsigned size = (unsigned)vertices.size();
            for(unsigned int i = 0; i < size; i++)
            {
                Line_segment li(this->vertices[i], this->vertices[(i + 1) % size]);
                min = std::min(min, li.dist_to(pt));
            }
            return min;
        }
        
        // The distance to a line segment.
        double distance_to(const Line_segment& line) const
        {
            double min = std::numeric_limits<double>::infinity();
            unsigned size = (unsigned)vertices.size();
            for(unsigned int i = 0; i < size; i++)
            {
                Line_segment li(this->vertices[i], this->vertices[(i + 1) % size]);
                min = std::min(min, li.dist_to_line_seg(line));
            }
            return min;
        }
        
        // If your polygon is a convex, so is other, this function applies GJK algorithm which can be very fast.
        // Otherwise, please use "naive_distance_to" method
        double distance_to( const Polygon& other ) const
        {
            if( this->intersects(other) ) return 0.0;
            return GJK::distance(this->vertices, other.vertices);
        }
        
        // This method loop over all line segments of the polygon and other to test min distance
        // That's why it is naive.
        double naive_distance_to(const Polygon& other ) const
        {
            if( this->naive_intersects(other) )
                return 0.0;
            double min = MAX_DOUBLE;
            unsigned size = (unsigned)vertices.size();
            for(unsigned int i = 1; i < size; i++)
            {
                Line_segment line(this->vertices[i], this->vertices[(i + 1) % size]);
                min = std::min( min, this->distance_to(line) );
            }
            return min;
        }
        
        /* return this closest point to the given point
         */
        v2 closest_pt_to( const v2& point ) const
        {
            v2 nearest;
            double min_dist = std::numeric_limits<double>::infinity();
            
            unsigned size = (unsigned)vertices.size();
            for(unsigned int i = 0; i < size; i++)
            {
                Line_segment line(this->vertices[i], this->vertices[(i + 1) % size]);
                v2 temp = line.project_in(point);
                double temp_dist = (point - temp).r();
                if( min_dist > temp_dist )
                {
                    min_dist = temp_dist;
                    nearest = temp;
                }
            }
            return nearest;
        }
    };
    
}


#endif
