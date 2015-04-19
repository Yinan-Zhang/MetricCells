//
//  cspace_scanner.h
//  RSS2015
//
//  Created by Yinan Zhang on 2/18/15.
//  Copyright (c) 2015 Yinan Zhang. All rights reserved.
//

#ifndef RSS2015_cspace_scanner_h
#define RSS2015_cspace_scanner_h

#include "../math/Naive2D/geometry.h"

#include "obstacles.h"
#include "robot.h"

template <typename IROBOT>
class CSpaceScanner
{
public:
    CSpaceScanner( IROBOT& robot_, robotics::ObstManager& obstMgr_ ):robot (robot_), obstMgr (obstMgr_)
    {
    }
    
    // Scan the whole space
    // @return: invalid configs (as 2d spheres)
    std::vector<N2D::sphere> scan( robotics::Range xrange, robotics::Range yrange,  double interval )
    {
        for( double x = std::get<0>(xrange); x < std::get<1>(xrange); x += interval )
        {
            for (double y = std::get<0>(yrange); y < std::get<1>(yrange); y += interval) {
                ND::vec<2> cfg({x,y});
                robot.set_config(cfg);
                if(this->obstMgr.intersects(robot))
                {
                    obst_cfgs.emplace_back(N2D::v2(x,y), interval/2.0, N2D::SPHEREMETRIC::LINFTY);
                }
            }
        }
        
        return obst_cfgs;
    }
    
private:
    IROBOT& robot;
    robotics::ObstManager& obstMgr;
    
    std::vector<N2D::sphere> obst_cfgs;
};

#endif
