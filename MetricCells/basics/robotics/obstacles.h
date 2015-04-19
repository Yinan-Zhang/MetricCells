//
//  obstacles.h
//  RSS2015
//
//  Created by Yinan Zhang on 12/26/14.
//  Copyright (c) 2014 Yinan Zhang. All rights reserved.
//

#ifndef RSS2015_obstacles_h
#define RSS2015_obstacles_h

#include <vector>
#include "../math/Naive2D/polygon.h"
#include "robot.h"

namespace robotics {
    /****************************************
     * Obstacle class
     ****************************************/
    class Obstacle: public N2D::Polygon
    {
    public:
        Obstacle( std::vector<N2D::v2> points ):N2D::Polygon(std::move(points)){}
        Obstacle( const N2D::Polygon& poly ):N2D::Polygon( poly.vertices ){}
    };
    
    /****************************************
     * Obstacle Manager class
     ****************************************/
    class ObstManager
    {
    public:
        ObstManager(){};
        
        ObstManager( const std::vector<Obstacle>& obstacles )
        {
            this->set_obstacles(obstacles);
        }
        
        void set_obstacles(const std::vector<Obstacle>& obstacles)
        {
            this->obstacles = obstacles;
            this->obstacles.shrink_to_fit();
        }
        
        /* determine if the robot intersects with any obstacles in the space
         */
        template<typename IROBOT>
        bool intersects( const IROBOT& robot ) const
        {
            const std::vector<N2D::Polygon>& robot_shapes = robot.get_shape();
            for( N2D::Polygon poly : robot_shapes )
                if(this->intersects(poly))
                        return true;
            return false;
        }
        
        /* determine if the polygon intersects with any obstacles in the space
         */
        bool intersects( const N2D::Polygon& poly ) const
        {
            for(unsigned int i = 0; i < this->obstacles.size(); i++)
            {
                if(this->obstacles[i].intersects(poly))
                    return true;
            }
            return false;
        }
        
        
        /* return the min dist from a point to any obstacles. */
        double dist_to_obsts(const N2D::v2& point) const
        {
            double min = std::numeric_limits<double>::infinity();
            for(unsigned int i = 0; i < this->obstacles.size(); i++)
            {
                min = std::min(min, this->obstacles[i].distance_to(point));
            }
            return min;
        }
        
        /* return the min dist from the robot to any obstacles. */
        template<typename IROBOT>
        double dist_to_obsts(const IROBOT& robot) const
        {
            double min = std::numeric_limits<double>::infinity();
            const std::vector<N2D::Polygon>& robot_shapes = robot.get_shape();
            //std::cout << robot->get_config()[0] << '\t' << robot->get_config()[1] << '\t' << robot->get_config()[2] << '\n';
            
            for(unsigned int i = 0; i < this->obstacles.size(); i++)
            {
                //std::cout << i << '\n';
                for( const N2D::Polygon& robot_shape : robot_shapes )
                    min = std::min(min, robot_shape.distance_to(this->obstacles[i]));
			}
            return min;
        }
        
        /* return the min dist from the robot to any obstacles. */
        template<typename IROBOT>
        double penetration(const IROBOT& robot) const
        {
            bool flag = false;
            const std::vector<N2D::Polygon>& robot_shapes = robot.get_shape();
            
            double max = 0.0;
            for( const N2D::Polygon& robot_shape : robot_shapes )
            {
                for(unsigned int i = 0; i < this->obstacles.size(); i++)
                {
                    if( obstacles[i].intersects(robot_shape) )
                    {
                        double min = std::numeric_limits<double>::infinity();
                        for (const N2D::v2& point : robot_shape.vertices)
                        {
                            if (this->obstacles[i].contains(point))
                            {
                                min = obstacles[i].penetration(point);
                            }
                            if (std::isfinite(min))
                                max = std::max(max, min);
                        }
                    }
                }
            }
            
            return max;
            //throw "This should not happen";
        }
        
        /* return the min time for a robot to collide with any obstacles */
        template<typename IROBOT>
        double time_leave_obsts(IROBOT& robot, double velocity) const
        {
            return this->penetration(robot) / velocity;
        }
        
        /* return the min time for a robot to collide with any obstacles */
        template<typename IROBOT>
        double time_to_obsts(const IROBOT& robot, double velocity) const
        {
            return this->dist_to_obsts(robot) / velocity;
        }
        
        /* get the nearest point in obstacles to a point */
        N2D::v2 closest_pt_to_obst( const N2D::v2& point ) const
        {
            double min = std::numeric_limits<double>::infinity();
            N2D::v2 nearest;
            
            for(unsigned int i = 0; i < this->obstacles.size(); i++)
            {
                N2D::v2 temp = this->obstacles[i].closest_pt_to(point);
                double temp_dist = (temp - point).r();
                if( temp_dist < min )
                {
                    min = temp_dist;
                    nearest = temp;
                }
            }
            
            return nearest;
        }
        
        /* test if a point is inside any polygon */
        bool inside_obsts( const N2D::v2& point ) const
        {
            for(unsigned int i = 0; i < this->obstacles.size(); i++)
            {
                if(this->obstacles[i].contains(point))
                    return true;
            }
            return false;
        }

		const std::vector<Obstacle>& get_obsts() const {
		    return obstacles;
		}
        
    private:
        std::vector<Obstacle> obstacles;
    };
}



#endif
