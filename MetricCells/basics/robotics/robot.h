//
//  robot.h
//  MetricCells
//
//  Created by Yinan Zhang on 4/19/15.
//  Copyright (c) 2015 Yinan Zhang. All rights reserved.
//

#ifndef MetricCells_robot_h
#define MetricCells_robot_h

namespace robotics
{
    
    /* use tuple to represent range
     * - to create a range instance, use std::make_tuple( lower, upper )
     * - to get the lower bound, use std::get<0>(range)
     * - to get the upper bound, use std::get<1>(ranges)
     */
    typedef std::tuple <double, double> Range;
    
}
#endif
