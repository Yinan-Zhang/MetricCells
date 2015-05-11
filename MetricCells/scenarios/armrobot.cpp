#ifndef ARM_ROBOT_CPP
#define ARM_ROBOT_CPP
#include <limits>
#include <cmath>
#include <unistd.h>
#include <algorithm>
#include <fstream>
#include <array>
#include "armrobot.h"


template <size_t DIM>
ArmRobot<DIM>::ArmRobot(double anchor_x,
                        double anchor_y,
                        const std::array<double, DIM>& lengths,
                        const std::array<double, DIM>& max_speeds)
                        : anchor_x_ (anchor_x),
                          anchor_y_ (anchor_y),
                          chain_lengths_ (lengths),
                          parameter_speeds(max_speeds)
{
    //parameter_speeds = max_speeds;
    // Calculate max linear speed
    std::array<double, DIM> max_linear_speeds{0};
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < i + 1; j++) {
            double chain_length_sum = 0.0;
            for (int k = j; k < i + 1; k++) {
                chain_length_sum += chain_lengths_[k];
            }
            max_linear_speeds[i] += parameter_speeds[j] * chain_length_sum;
        }
    }
    max_speed = max_linear_speeds[DIM - 1];
    
    // Set config ranges
    double e = 0.000000001;
    for (int i = 0; i < DIM; i++)
        config_ranges_[i] = robotics::Range(-M_PI + e, M_PI - e);
    
    // build shapes of the robot
    std::array<N2D::v2, DIM+1> points = this->get_points();
    for (int i = 1; i < DIM+1; i++) {
        N2D::Polygon poly({points[i-1],points[i]});
        this->shapes.emplace_back( poly );
    }
}
#endif
