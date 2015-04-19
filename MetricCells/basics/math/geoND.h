//
//  geometry.h
//  RSS2015
//
//  Created by Yinan Zhang on 12/24/14.
//  Copyright (c) 2014 Yinan Zhang. All rights reserved.
//

#ifndef RSS2015_geometry_h
#define RSS2015_geometry_h

#include <cassert>
#include <vector>
#include <cmath>       /* sqrt */
#include <ostream>
#include <iostream>
#include <limits>
#include <array>

static constexpr double MAX_DOUBLE = std::numeric_limits<double>::infinity();

namespace ND
{
    /*****************************
     * N-D Vector
     *****************************/
    template <size_t N> struct vec
    {
        std::array<double, N> data;
        
        explicit vec() = default;
        
        explicit vec( const std::array<double, N>& data_ ) : data(data_) {}
        
        //vec(const double* data_){ std::copy(data, data+N, this->data); }
        
        const vec<N>& vector() const { return *this; }
        
        int size() const { return N; }
        
        // Plus
        const vec<N> operator+(const vec<N> &other) const
        {
            return vec<N>(*this) += other;
        }
        
        // +=
        vec<N>& operator+=(const vec<N> &other)
        {
            for(int i = 0; i < N; i++)
            {
                this->data[i] += other.data[i];
            }
            return *this;
        }
        
        const double& operator[] (std::size_t idx) const { return (this->data[idx]); }
        
        double& operator[] (std::size_t idx) { return (this->data[idx]); }
        
        // Minus
        const vec<N> operator-(const vec<N> &other) const
        {
            return vec<N>(*this) -= other;
        }
        
        // Unary Minus
        const vec<N> operator-() const
        {
            double data_[N];
            for(int i = 0; i < N; i++)
            {
                data_[i] = - this->data[i];
            }
            return vec<N>(data_);
        }
        
        // -=
        vec<N>& operator-=(const vec<N> &other)
        {
            for(int i = 0; i < N; i++)
            {
                this->data[i] -= other.data[i];
            }
            return *this;
        }
        
        // times a T "b".  a*b = v2( a.x*b, a.y*b, a.z*b )
        const vec<N> operator*(double b) const
        {
            return vec<N>(*this) *= b;
        }
        
        // *=
        vec<N>& operator*=(double b)
        {
            for(int i = 0; i < N; i++)
            {
                this->data[i] *= b;
            }
            return *this;
        }
        
        // devided by a T "b".  a/b = v2( a.x/b, a.y/b, a.z/b )
        const vec<N> operator/(double b) const
        {
            return vec<N>(*this) *= 1/b;
        }
        
        // /=
        vec<N>& operator/=(double b)
        {
            return *this *= 1 / b;
        }
        
        bool operator==( const vec<N>& other ) const
        {

            for(int i = 0; i < N; i++)
            {
                if( this->data[i] != other.data[i] )
                    return false;
            }
            return true;
        }
        
        // get the length of the vector
        double r() const
        {
            double sum = 0;
            for(int i = 0; i < N; i++)
            {
                sum += this->data[i] * this->data[i];
            }
            return std::sqrt(sum);
        }
        
        double rsq() const
        {
            double sum = 0;
            for(int i = 0; i < N; i++)
            {
                sum += this->data[i] * this->data[i];
            }
            return sum;
        }
        
        double l2() const { return r(); }
        
        double l1() const
        {
            double sum = 0;
            for(int i = 0; i < N; i++)
            {
                sum += fabs( this->data[i] );
            }
            return sum;
        }
        
        double linfty() const
        {
            double max = 0;
            for(int i = 0; i < N; i++)
            {
                if( fabs(data[i]) > max )
                    max = fabs(data[i]);
            }
            return max;
        }
        
        // Normalize
        void normalize() { *this=*this*( 1/r()); }
        vec<N> norm() const
        {
            double r = this->r();
            vec<N> result;
            for( int i = 0; i < N; i++  )
            {
                result[i] = this->data[i]/r;
            }
            return result;
        }
        
        // dot product
        double dot(const vec<N> &other) const
        {
            double sum = 0;
            for(int i = 0; i < N; i++)
            {
                sum += this->data[i] * other.data[i];
            }
            return sum;
        }
        
        // corss product
        //v2 cross( v2 &b ){return v2(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);} // Cross;
    };
    
    template <size_t N> std::ostream& operator<< (std::ostream& str, const vec<N>& vector);
    
    /*****************************
     * N-D Sphere
     *****************************/
    enum class SPHEREMETRIC{ L1, L2, LINFTY };
    
    template <size_t N> struct sphere
    {
        vec<N> c_;
        double r_;
        SPHEREMETRIC metric;
        
        explicit sphere(const vec<N>& center, double radius, SPHEREMETRIC metric_ = SPHEREMETRIC::LINFTY) : c_(center), r_(radius), metric(metric_) {}
        
        const vec<N>& center() const {return this->c_;}
        double radius() const { return this->r_; }
        
        bool on_boundary( const vec<N>& point, double tolerance ) const
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
        
        bool contains( const vec<N>& point ) const
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
        
        bool intersects( sphere<N>& other ) const
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
        
        /* determine if this sphere and other are neighbors by checking their distance( < tolerance ). */
        bool neighbor( sphere<N>& other, double tolerance ) const
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
    };

    
    template <size_t N> std::ostream& operator<< (std::ostream& str, const vec<N>& vector)
    {
        str << "vec(";
        for(int i = 0; i < N; i++)
            if(i < N-1)
                str << vector[i] << ',';
            else
                str << vector[i];
        str << ')';
        return str;
    }
}

#endif
