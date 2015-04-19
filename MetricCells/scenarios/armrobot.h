#ifndef ARM_ROBOT_H
#define ARM_ROBOT_H

#include <memory>
#include <vector>
#include <limits>
#include <array>
#include "../basics/math/geoND.h"
#include "../basics/math/Naive2D/geometry.h"
#include "../basics/math/Naive2D/polygon.h"
#include "../basics/robotics/obstacles.h"

template<size_t DIMENSION>
class ArmRobot {
public:
    static constexpr int DIM = DIMENSION;
    
    static constexpr int PENETRATION_CONSTANT = 1;
    
    static constexpr bool MERGE_CELLS = true;
    
    using CONFIG = ND::vec<DIM>;
    
    // Constructor
    ArmRobot(double anchor_x,                               // base x position
             double anchor_y,                               // base y position
             const std::array<double, DIM>& lengths,        // lengthe of each link
             const std::array<double, DIM>& max_speeds);    // max speed of each joint
    
    
    // Get the range of each config
    const std::array<robotics::Range, DIMENSION>& get_config_ranges() const {
        return config_ranges_;
    }
    
    
    // Get each link as an array of line segments
    std::array<N2D::Line_segment, DIM> get_lines() const {
        std::array<N2D::Line_segment, DIM> segments;
        double prev_joint_x = anchor_x_;
        double prev_joint_y = anchor_y_;
        double sum_of_angles = 0.0;
        for (unsigned int i = 0; i < DIM; i++) {
            sum_of_angles += config[i];
            double curr_joint_x = chain_lengths_[i] * std::sin(sum_of_angles) + prev_joint_x;
            double curr_joint_y = -chain_lengths_[i] * std::cos(sum_of_angles) + prev_joint_y;
            segments[i] = N2D::Line_segment(N2D::v2(prev_joint_x, prev_joint_y),
                                                 N2D::v2(curr_joint_x, curr_joint_y));
            prev_joint_x = curr_joint_x;
            prev_joint_y = curr_joint_y;
        }
        return segments;
    }
    
    // determine if the robot intersects with an obstacle
    bool intersects(const N2D::Polygon& obstacle) const {
        std::array<N2D::Line_segment, DIM> curr_lines = get_lines();
        for (int i = 0; i < DIM; i++)
            if (obstacle.intersects(curr_lines[i]))
                return true;
        return false;

    }
    
    // @return the distance to an obstacle
    double dist_to_obstacle(const N2D::Polygon& obstacle) const {
        dist_to_obst(obstacle, std::numeric_limits<int>::max());
        return dist_to_obst(obstacle, std::numeric_limits<int>::max());
    }
    
    
    // Free space oracle
    double free_space_oracle(robotics::ObstManager& obstacle_manager, double lower = 0.0, double upper = M_PI) const {
        double threshold = 0.01;
        while (upper - lower > threshold) {
            double midpt = (lower + upper) / 2;
            std::array<double, DIM> max_reaches;
            for (int i = 0; i < DIM; ++i)
                max_reaches[i] = std::numeric_limits<double>::max();
            CalcExtremalReaches(midpt, max_reaches);
            // Get "effective lengths" by taking differences in extremal lengths
            for (unsigned int i = 1; i < DIM; i++)
                max_reaches[i] = std::max(max_reaches[i] - max_reaches[i-1], 0.0);
            
            double time_to_collision = MinTimeToCollision(obstacle_manager, max_reaches);
            // TODO: Check max angular speeds
            if (time_to_collision > midpt / parameter_speeds[0]) {
                lower = midpt;
            } else {
                upper = midpt;
            }
        }
        return lower;
    }
    
    
    // @return: true if the cell is sub-convex of s.
    bool subconvex_radius(const robotics::ObstManager& obstacle_manager, double r, double s) const {
        return r >= s;
    }
    
    double get_s_small(double epsilon) {
        if (s_small_ < 0.00000001) {
            double quad_l = 0.0;
            for (int i = 0; i < DIM; i++) {
                for (int j = i; j < DIM; j++) {
                    quad_l += chain_lengths_[j];
                }
            }
            double numerator = pow(epsilon, 2);
            double denominator = 8 * pow(quad_l, 2);
            s_small_ = acos(1-(numerator / denominator));
        }
        return s_small_;
    }
    
    
    void set_config(const CONFIG& rhs) {
        config = rhs;
        
        std::array<N2D::v2, DIM+1> points = this->get_points();
        
        for( int i = 0; i < DIM; i++ )
        {
            this->shapes[i].vertices[0] = points[i];
            this->shapes[i].vertices[1] = points[i+1];
        }

    }
    
    const std::array<double, DIM>& get_parameter_speeds() const { return parameter_speeds; }
    double get_max_speed() const { return max_speed; }
    const std::vector<N2D::Polygon>& get_shape() const {
        return this->shapes;
    }
    
    static bool cmp(const CONFIG* lhs, const CONFIG* rhs) {
        for (int i = 0; i < DIM; ++i) {
            if ((*lhs)[i] < (*rhs)[i])
                return true;
            else if ((*lhs)[i] > (*rhs)[i])
                return false;
        }
        return lhs < rhs;
    }
    
    static bool approximate(const CONFIG* lhs, const CONFIG* rhs)
    {
        for (int i = 0; i < DIM; ++i) {
            if (fabs((*lhs)[i] - (*rhs)[i]) > 1e-15)
                return false;
        }
        return true;
    }
    
    static double heuristic(const ArmRobot<DIMENSION>& robot, const CONFIG& start, const CONFIG& goal)
    {
        return (start-goal).linfty();
    }
    
    static double metric(const ArmRobot<DIMENSION>& car, const CONFIG& start, const CONFIG& goal)
    {
        return (start-goal).linfty();
    }
    
private:
    double dist_to_obst(const N2D::Polygon& obstacle, int chain_limit) const {
        double dist = std::numeric_limits<int>::max();
        std::array<N2D::Line_segment, DIM> curr_lines = get_lines();
        
        int limit;
        if (static_cast<unsigned int>(chain_limit) > curr_lines.size()) {
            limit = curr_lines.size();
        }
        else {
            limit = chain_limit;
        }
        for (int i = 0; i < limit; i++) {
            double curr_dist = obstacle.distance_to(curr_lines[i]);
            if (curr_dist < dist)
                dist = curr_dist;
        }
        return dist;

    }
    void CalcExtremalReaches(double radius, std::array<double, DIM>& max_reach_upper_bounds /* out param */) const {
        // Upper bounds the extent of each joint on the arm.
        // To do this we add the max length of each successive pair of segments
        // as if they were laid out straight (sum of pair lengths)
        for (int attempt = 0; attempt < 2; attempt++) {
            // Two attempts: one set of results is from pairs starting at the first segment,
            // the second set of results is from pairs starting from starting at the second segment
            double curr_reach = 0.0;
            int num_pairs = DIM / 2;
            if (attempt == 1) {
                curr_reach = chain_lengths_[0];
                max_reach_upper_bounds[0] = curr_reach;
                num_pairs = (DIM - 1) / 2;
            }
            
            for (int i = 0; i < num_pairs; i++) {
                int s1_i = i * 2 + attempt;
                int s2_i = s1_i + 1;
                double ideal_angle_adj = config[s2_i]; // clockwise
                if (ideal_angle_adj <= radius) {
                    // We can make the ideal angle adjustment, segments are straight
                    curr_reach += chain_lengths_[s1_i] + chain_lengths_[s2_i];
                } else {
                    // Angle adjustment is not in the ball, find the closest to ideal angle using law of cosines
                    double opp_angle = M_PI - (ideal_angle_adj - radius);
                    double max_length = (pow(chain_lengths_[s1_i], 2) + pow(chain_lengths_[s2_i], 2)
                                         - 2 * chain_lengths_[s1_i] * chain_lengths_[s2_i] * cos(opp_angle));
                    max_length = pow(max_length, 0.5);
                    curr_reach += max_length;
                }
                
                if (s1_i == 0 || max_reach_upper_bounds[s1_i - 1] + chain_lengths_[i] < max_reach_upper_bounds[s1_i]) {
                    // Update max reach for first segment
                    if (s1_i == 0)
                        max_reach_upper_bounds[s1_i] = chain_lengths_[0];
                    else
                        max_reach_upper_bounds[s1_i] = max_reach_upper_bounds[s1_i - 1] + chain_lengths_[i];
                }
                if (curr_reach < max_reach_upper_bounds[s2_i]) {
                    // Update max reach for second segment
                    max_reach_upper_bounds[s2_i] = curr_reach;
                }
            }
            
            // Add last segment's length if not in pair already
            if (num_pairs * 2 + attempt != DIM)
                if (max_reach_upper_bounds[DIM - 1] > curr_reach + chain_lengths_[DIM - 1])
                    max_reach_upper_bounds[DIM - 1] = curr_reach + chain_lengths_[DIM - 1];
        }
    }
    double MinTimeToCollision(robotics::ObstManager& obstacle_manager, const std::array<double, DIM>& max_reaches) const {
        const std::vector<robotics::Obstacle>& obsts = obstacle_manager.get_obsts();
        std::array<double, DIM> min_dists{0.0};
        for (int chain_i = 0; chain_i < DIM; chain_i++) {
            double curr_min_dist = std::numeric_limits<double>::infinity();
            for (unsigned int obs_i = 0; obs_i < obsts.size(); obs_i++) {
                double curr_dist = dist_to_obst(obsts[obs_i], chain_i + 1 /* chain limit */);
                if (curr_dist < curr_min_dist) {
                    curr_min_dist = curr_dist;
                }
            }
            min_dists[chain_i] = curr_min_dist;
        }
        
        std::array<double, DIM> max_speeds{0.0};
        for (unsigned int i = 0; i < chain_lengths_.size(); i++) {
            for (unsigned int j = 0; j < i + 1; j++) {
                double sum_len = 0.0;
                for (unsigned int k = j; k < i + 1; k++) {
                    sum_len += max_reaches[k];
                }
                max_speeds[i] += parameter_speeds[j] * sum_len;
            }
        }
        
        double min_time_to_collision = std::numeric_limits<double>::max();
        for (int i = 0; i < DIM; i++) {
            double curr_time = min_dists[i] / max_speeds[i];
            if (curr_time < min_time_to_collision)
                min_time_to_collision = curr_time;
        }
        return min_time_to_collision;
    }
    
    inline std::array<N2D::v2, DIM+1> get_points() const
    {
        std::array<N2D::v2, DIM+1> points;
        double prev_joint_x = anchor_x_;
        double prev_joint_y = anchor_y_;
        points[0] = N2D::v2(prev_joint_x, prev_joint_y);
        double sum_of_angles = 0.0;
        for (unsigned int i = 0; i < DIM; i++) {
            sum_of_angles += config[i];
            double curr_joint_x = chain_lengths_[i] * std::sin(sum_of_angles) + prev_joint_x;
            double curr_joint_y = prev_joint_y - chain_lengths_[i] * std::cos(sum_of_angles);
            prev_joint_x = curr_joint_x;
            prev_joint_y = curr_joint_y;
            points[i+1] = N2D::v2(prev_joint_x, prev_joint_y);
        }
        return points;
    }

    std::array<robotics::Range, DIM> config_ranges_;
    const std::array<double, DIM> parameter_speeds;
    const std::array<double, DIM> chain_lengths_;
    CONFIG config;
    const double anchor_x_;
    const double anchor_y_;
    double s_small_;
    double max_speed;
    mutable std::vector<N2D::Polygon> shapes;
};
#include "armrobot.cpp"
#endif /* ARM_ROBOT_H*/
