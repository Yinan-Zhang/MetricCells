//
//  geometry.h
//  RSS2015
//
//  Created by Yinan Zhang on 12/24/14.
//  Copyright (c) 2014 Yinan Zhang. All rights reserved.
//

#ifndef NAIVE2D_GEOMETRY_h
#define NAIVE2D_GEOMETRY_h

#include <cassert>
#include <vector>
#include <cmath>       /* sqrt */
#include <ostream>
#include <iostream>
#include <limits>
#include <array>

namespace N2D
{
    double N2D_MAX_DOUBLE = std::numeric_limits<double>::infinity();
    
    /*****************************
     * 2D Vector
     *****************************/
    struct v2
    {
        double x, y;
        
        // Constructor
        explicit v2(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
        
        // Plus
        const v2 operator+(const v2 &b) const { return v2(this->x+b.x, this->y+b.y); }
        
        // +=
        v2& operator+=(const v2 &b) { this->x+= b.x; this->y+=b.y; return *this; }
        
        // Minus
        const v2 operator-(const v2 &b) const {return v2(this->x-b.x, this->y-b.y);}
        
        // Unary Minus
        const v2 operator-() const {return v2(-(this->x), -(this->y));}
        
        // -=
        v2& operator-=(const v2 &b) { this->x-= b.x; this->y-=b.y; return *this;}
        
        // times a T "b".  a*b = v2( a.x*b, a.y*b, a.z*b )
        const v2 operator*(double b) const {return v2(this->x*b,  this->y*b);}
        
        //v2 operator*(v2 b) const {return v2(x*b.x, y*b.y );}
        
        // *=
        v2& operator*=(double b) { this->x*= b; this->y*=b; return *this;}
        
        // devided by a T "b".  a/b = v2( a.x/b, a.y/b, a.z/b )
        const v2 operator/(double b) const {return v2(*this) *= 1/b; }
        
        // /=
        v2& operator/=(double b) { return *this *= (1/b); }
        
        bool operator==( const v2 &other ) const { return this->x == other.x && this->y == other.y; }
        
        bool operator!=( const v2 &other ) const { return this->x != other.x || this->y != other.y; }
        
        // multiply a.mult(b) = ( a.x*b.x, a.y*b.y, a.z*b.z )
        //v2 mult( const v2 &b ) const {return v2(x*b.x, y*b.y, z*b.z );}
        
        // Normalize
        void normalize(){ *this *= 1 / r(); }
        v2 norm() const { return v2(*this) * (1/r()); }
        
        // dot product
        double dot(const v2 &b) const { return x*b.x + y*b.y; }
        
        double cross(const v2&b) const { return x * b.y - y * b.x; }
        
        // corss product
        //v2 cross( v2 &b ){return v2(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);} // Cross;
        
        // get the length of the vector
        double r() const { return sqrt( this->dot(*this) ); }
        
        double rsq() const { return this->dot(*this); }
        
        double l1() const { return fabs(x)+fabs(y); }
        
        double l2() const { return r(); }
        
        double linfty() const { return std::max(fabs(x),fabs(y)); }
    };
    
#define INFINITE_POINT N2D::v2(N2D_MAX_DOUBLE, N2D_MAX_DOUBLE)
    
    std::ostream& operator<< (std::ostream& str, const v2& v);
    
    /*****************************
     * 2D Line Segment
     *****************************/
    struct Line_segment
    {
        v2 start, end;
        
        explicit Line_segment() = default;
        explicit Line_segment(const v2& arg_start, const v2& arg_end) : start(arg_start), end(arg_end) {}
        
        // return the vector from start to end
        v2 vec() const { return this->end - this->start; }
        
        // length of the line segment
        double length() const {
            return (end - start).r();
        }
        
        // determine if two line segments are equal
        bool operator==( const Line_segment &other ){ return this->start == other.start && this->start == other.start; }
        
        // Unary Minus. return a line segment with opposite direction
        const Line_segment operator-() const {return Line_segment(this->end, this->start);}
        
        /*Project `pt` onto the infinite line through `self` and return the
        value `t` where the projection = ``start + (end-start)*t``.
        This has the convenient property that t \in [0,1] if the projection is on
        the line segment.*/
        double project_t( const v2& pt ) const
        {
            //v2 v = this->vec();
            double r = this->length();
            if( r == 0.0 )
                return 0.0;
            else
            {
                return (pt-this->start).dot(this->vec()) / r;
            }
        }
        
        /* Project `pt` onto the infinite line through `self` and return the
         * projected point.
         */
        v2 project(const v2& pt) const
        {
            double t = this->project_t(pt);
            return this->start + this->vec()*t;
        }
        
        /* Project `pt` onto the line segment through `self` and return the
         * projected point. This method stricts that
         */
        v2 project_in(const v2& pt ) const
        {
            v2 ab = end - start;
            return start + ab * std::max( 0.0, std::min((pt - start).dot(ab) / ab.rsq(), 1.0));
        }
        
        // determine if three points are listed in a counterclockwise order
        bool __ccw__( const v2& A, const v2& B, const v2& C ) const {return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x);}
        
        // Determine if two line segments intersects with each other
        bool intersects( const Line_segment& other ) const
        {
            v2 A = this->start; v2 B = this->end;
            v2 C = other.start; v2 D = other.end;
            return __ccw__(A,C,D) != __ccw__(B,C,D) && __ccw__(A,B,C) != __ccw__(A,B,D);
        }
        
        //Return the closest point on the line segment to `pt`.
        v2 closest_point(const v2& pt) const
        {
            return this->project_in(pt);
        }
        
        // Returns the closest distance from a point to the line segment
        double dist_to( const v2& pt ) const
        {
            return (pt-project_in(pt)).r();
        }
        
        // Returns the distance from another line to self.
        double dist_to( const Line_segment& other ) const
        {
            return this->dist_to_line_seg(other);
        }
        
        // return this closest distance between two line segments
        double dist_to_line_seg(const Line_segment& other) const
        {
            if(intersects(other)) return 0.0;
            
            v2 other_proj_start    = this->project_in( other.start );
            v2 other_proj_end      = this->project_in( other.end );
            v2 proj_to_otehr_start = other.project_in( this->start );
            v2 proj_to_otehr_end   = other.project_in( this->end );
            
            double d1 = (other_proj_start - other.start).r();
            double d2 = (other_proj_end - other.end).r();
            double d3 = (proj_to_otehr_start - this->start).r();
            double d4 = (proj_to_otehr_end - this->end).r();
            
            return  std::min( std::min(d1, d2), std::min(d3, d4) );
        }
        
        /* whether two segments in the plane intersect:
         * one segment is (x11, y11) to (x12, y12)
         * the other is   (x21, y21) to (x22, y22)
         * Source: http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
         */
        v2 intersection_point( const Line_segment& other ) const
        {
            double x11 = this->start.x;
            double y11 = this->start.y;
            double x12 = this->end.x;
            double y12 = this->end.y;
            double x21 = other.start.x;
            double y21 = other.start.y;
            double x22 = other.end.x;
            double y22 = other.end.y;
            double dx1 = x12 - x11;
            double dy1 = y12 - y11;
            double dx2 = x22 - x21;
            double dy2 = y22 - y21;
            
            double delta = dx2 * dy1 - dy2 * dx1;
            if (fabs(delta) < 0.0000001) return INFINITE_POINT; // parallel segments, no intersection
            double s = (dx1 * (y21 - y11) + dy1 * (x11 - x21)) / delta;
            double t = (dx2 * (y11 - y21) + dy2 * (x21 - x11)) / (-1 * delta);
            if ( !((0 <= s && s <= 1) && (0 <= t && t <= 1)))
                return INFINITE_POINT; // No intersection
            double i_x = x11 + t * dx1;
            double i_y = y11 + t * dy1;
            return v2(i_x, i_y);
        }
    };
    
    std::ostream& operator << (std::ostream& str, const Line_segment& line);
    
    /*****************************
     * 2D Sphere
     *****************************/
    enum class SPHEREMETRIC{ L1, L2, LINFTY };
    
    struct sphere
    {
        v2 c_;
        double r_;
        SPHEREMETRIC metric;
        
        explicit sphere(const v2& center=v2(0.0, 0.0), double radius=0.0, SPHEREMETRIC metric_ = SPHEREMETRIC::L2) : c_(center), r_(radius), metric(metric_) {}
        sphere( const sphere& copy ):c_(copy.c_), r_(copy.r_), metric(copy.metric){}
        
        // Returns the center of the sphere
        const v2& center() const {return this->c_;}
        // Returns this radius of the sphere
        double radius() const { return this->r_; }
        
        // Determines if a point is on the boundary of the sphere.
        // That is if | |point - center| - radius| <= tolerance
        bool on_boundary( const v2& point, double tolerance ) const
        {
            double dist;
            switch (metric) {
                case SPHEREMETRIC::L1:
                    dist = (this->c_ - point).l1();
                    break;
                case SPHEREMETRIC::L2:
                    dist = (this->c_ - point).r();
                    break;
                case SPHEREMETRIC::LINFTY:
                    dist = (this->c_ - point).linfty();
                    break;
                default:
                    throw "metric has to be l1, l2 or li.";
            }
            return (dist - this->r_) < tolerance;
        }
        
        // Returns true if the sphere contains the point
        bool contains( const v2& point ) const
        {
            double dist;
            switch (metric) {
                case SPHEREMETRIC::L1:
                    dist = (this->c_ - point).l1();
                    break;
                case SPHEREMETRIC::L2:
                    dist = (this->c_ - point).r();
                    break;
                case SPHEREMETRIC::LINFTY:
                    dist = (this->c_ - point).linfty();
                    break;
                default:
                    throw "metric has to be l1, l2 or li.";
            }
            return dist < this->r_;
        }
        
        // returns true if two spheres intersects
        bool intersects( sphere& other ) const
        {
            double dist;
            switch (metric) {
                case SPHEREMETRIC::L1:
                    dist = (this->c_ - other.c_).l1();
                    break;
                case SPHEREMETRIC::L2:
                    dist = (this->c_ - other.c_).r();
                    break;
                case SPHEREMETRIC::LINFTY:
                    dist = (this->c_ - other.c_).linfty();
                    break;
                default:
                    throw "metric has to be l1, l2 or li.";
            }
            return dist < this->r_ + other.r_;
        }
        
        /*Returns true if the sphere intersects with a given line.*/
        bool intersects(Line_segment& line) const
        {
            return line.dist_to(center()) - radius() < 0;
        }
        
        /* Returns the distance between the sphere and point
         */
        double dist_to(v2& point) const
        {
            double dist;
            switch (metric) {
                case SPHEREMETRIC::L1:
                    dist = (this->c_ - point).l1();
                    break;
                case SPHEREMETRIC::L2:
                    dist = (this->c_ - point).r();
                    break;
                case SPHEREMETRIC::LINFTY:
                    dist = (this->c_ - point).linfty();
                    break;
                default:
                    throw "metric has to be l1, l2 or li.";
            }
            return dist-radius();
        }
        
        /* determines if this sphere and other are neighbors by checking their distance( < tolerance ). */
        bool neighbor( sphere& other, double tolerance ) const
        {
            double center_dist;
            switch (metric) {
                case SPHEREMETRIC::L1:
                    center_dist = (this->c_ - other.c_).l1();
                    break;
                case SPHEREMETRIC::L2:
                    center_dist = (this->c_ - other.c_).r();
                    break;
                case SPHEREMETRIC::LINFTY:
                    center_dist = (this->c_ - other.c_).linfty();
                    break;
                default:
                    throw "metric has to be l1, l2 or li.";
            }
            return (center_dist - this->r_ - other.r_) <= tolerance;
        }
        
        sphere& operator=(const sphere& copy)
        {
            this->c_ = copy.c_;
            this->r_ = copy.r_;
            this->metric = copy.metric;
            return *this;
        }
    };
    
    inline std::ostream& operator<< (std::ostream& str, const v2& v)
    { return str<< "v2(" << v.x<< ',' <<v.y << ')'; }
    
    
    inline std::ostream& operator << (std::ostream& str, const Line_segment& line)
    {
        return str << "line segment | " << (line.start) << "------" << line.end;
    }
}

#endif
