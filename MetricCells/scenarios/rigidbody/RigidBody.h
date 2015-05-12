//
//  RigidBody.h
//  MetricCells
//
//  Created by Yinan Zhang on 5/11/15.
//  Copyright (c) 2015 Yinan Zhang. All rights reserved.
//

#ifndef MetricCells_RigidBody_h
#define MetricCells_RigidBody_h

#include <memory>
#include <vector>
#include <limits>
#include <array>
#include "../../basics/math/geoND.h"
#include "../../basics/math/Naive2D/geometry.h"
#include "../../basics/math/Naive2D/polygon.h"
#include "../../basics/robotics/obstacles.h"

#include <cmath>       /* sin */
#include <vector>
#include <array>

#include <memory>
#include <vector>
#include <limits>
#include <array>
#include "../../basics/math/geoND.h"
#include "../../basics/math/Naive2D/geometry.h"
#include "../../basics/math/Naive2D/polygon.h"
#include "../../basics/robotics/obstacles.h"

/*
 A Reed-Shepp car looks like this:
 A
 |-o-|
 |   |
 |   |
 |-o-|
 B
 (simplified to a rectangle).
 It's definition can be found here: http://planning.cs.uiuc.edu/node822.html
 */
using namespace std;
class RigidBody
{
public:
    static constexpr int DIM = 3;
    static constexpr int PENETRATION_CONSTANT = 2;
    static constexpr bool MERGE_CELLS = true;
    using CONFIG = ND::vec<DIM>;
    double eps = 0.00001;
    RigidBody( double width_, double length_, double turning_radius )
    :width(width_),length(length_), radius(turning_radius), config(ND::vec<3>({0,0,0}))
    {
        max_speed = radius * (1 + std::max(width, length)/radius);
        parameter_speeds = std::array<double, DIM>{1,1,1};
        range = std::array<robotics::Range, DIM>{robotics::Range( 0+eps, M_PI-eps ), robotics::Range( 0+eps, M_PI-eps ), robotics::Range( 0+eps, M_PI-eps ) };
        this->shapes.emplace_back(this->get_points());
        this->shapes.shrink_to_fit();
    }
    
    void set_config( const CONFIG& config_ )
    {
        N2D::v2 tran( config_[0]-config[0], config_[1]-config[1] );
        double dtheta = config_[2] - config[2];
        
        this->config[0] = config_[0];
        this->config[1] = config_[1];
        this->config[2] = config_[2];
        
        this->shapes[0].self_translate(tran);
        this->shapes[0].self_rotate(dtheta, N2D::v2( this->config[0], config[1] ));
    }
    
    const CONFIG& get_config() const { return config; }
    
    const std::array<robotics::Range, DIM>& get_config_ranges() const
    {
        return range;
    }
    
    double dist_to( const N2D::Polygon& obst ) const
    {
        return this->shapes[0].distance_to(obst);
    }
    
    bool intersects( const N2D::Polygon& obst ) const
    {
        return this->shapes[0].intersects(obst);
    }
    
    const std::vector<N2D::Polygon>& get_shape() const
    {
        return this->shapes;
    }
    
    const std::array<double, DIM>& get_parameter_speeds() const
    {
        return parameter_speeds;
    }
    
    double free_space_oracle(const robotics::ObstManager& obstacle_manager, double lower = 0.0, double upper = M_PI) const {
        return lower / std::sqrt(2.0);
    }
    
    double h(double d) const {
        return d / 6;
    }
    
    double get_s_small(double r) const {
        return h(r / std::sqrt(2.0));
    }
    
    bool subconvex_radius(const robotics::ObstManager& obstacle_manager, double r, double s) const
    {
        return r >= s;
    }
    
    double get_radius() const {
        return radius;
    }
    
    double get_max_speed() const {
        return max_speed;
    }
    
    double get_length() const {
        return length;
    }
    
    double get_width() const {
        return width;
    }
    
    static inline bool cmp(const CONFIG* lhs, const CONFIG* rhs) {
        double l = (*lhs)[0], r = (*rhs)[0];
        if (l != r)
            return l < r;
        l = (*lhs)[1], r = (*rhs)[1];
        if (l != r)
            return l < r;
        l = (*lhs)[2], r = (*rhs)[2];
        if (l != r)
            return l < r;
        return lhs < rhs;
    }
    
    static inline bool approximate(const CONFIG* lhs, const CONFIG* rhs)
    {
        return (fabs((*lhs)[0] - (*rhs)[0]) <= 1e-15)
        && (fabs((*lhs)[1] - (*rhs)[1]) <= 1e-15)
        && (fabs((*lhs)[2] - (*rhs)[2]) <= 1e-15);
    }
    
    static inline double heuristic(const RigidBody& body, const CONFIG& start, const CONFIG& goal)
    {
        return std::max((N2D::v2( start[0], start[1] )-N2D::v2(goal[0], goal[1])).r(), std::fabs(goal[2]-start[2]));
    }
    
    static inline double metric(const RigidBody& body, const CONFIG& start, const CONFIG& goal)
    {
        return heuristic(body, start, goal);
    }
    
    
private:
    
    /* Get position of point A. (see the graph above.) */
    N2D::v2 A() const
    {
        double hl = this->length / 2.0;
        return N2D::v2(this->config[0] + hl*cos(this->config[2]), this->config[1] + hl*sin(this->config[2]));
    }
    
    /* Get position of point B. (see the graph above.) */
    N2D::v2 B() const
    {
        double hl = this->length / 2.0;
        return N2D::v2(this->config[1] - hl*cos(this->config[2]), this->config[1] - hl*sin(this->config[2]));
    }
    
    /* The center of A and B. */
    N2D::v2 center() const { return N2D::v2(config[0], config[1]); }
    
    /* Get points of four corners.
     * p1 |--| p2
     *    |  |
     * p4 |--| p3
     */
    std::vector<N2D::v2> get_points() const
    {
        double theta = this->config[2];
        double hw = this->width/2.0;
        //N2D::v2 A = this->A();
        N2D::v2 B = this->B();
        std::vector<N2D::v2> result;
        result.emplace_back(B.x-hw*sin(theta), B.y+hw*cos(theta));
        result.emplace_back(A());
        result.emplace_back(B.x+hw*sin(theta), B.y-hw*cos(theta));
        return result;
    }
    
private:
    std::array<double, DIM> parameter_speeds;
    std::array<robotics::Range, DIM> range;
    CONFIG config;
    double max_speed;
    double width;
    double length;
    double radius;
    mutable std::vector<N2D::Polygon> shapes;
};


#endif
