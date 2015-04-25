//
//  PRM.h
//  MetricCells
//
//  Created by Yinan Zhang on 4/20/15.
//  Copyright (c) 2015 Yinan Zhang. All rights reserved.
//

#ifndef MetricCells_PRM_h
#define MetricCells_PRM_h

#include <vector>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include "../math/geoND.h"
#include "../robotics/robot.h"
#include "../robotics/obstacles.h"

namespace algorithms
{
    namespace RoadMap
    {
        typedef std::unordered_map<int, unordered_set<int>> MAP;
        
        template <class IROBOT>
        class PRM {
            //using IROBOT::CONFIG;
            using MapNode = typename IROBOT::CONFIG;
        public:
            PRM( std::vector<robotics::Range>& ranges_, IROBOT robot_, robotics::ObstManager& obstMgr_ )
            {
                this->ranges = ranges_;
                this->obstMgr = obstMgr_;
                this->robot = robot_;
            }
            
            MAP construct(int max_num)
            {
                MAP map;
                this->nodes.reserve( max_num );
                for( int i = 0; i < max_num; i++ )
                {
                    MapNode node;
                    while( true )
                    {
                        for (int d = 0; d< IROBOT::DIM; d++) {
                            double low = std::get<0>( ranges[d] );
                            double up  = std::get<1>( ranges[d] );
                            node[d] = low +  double(std::rand())/RAND_MAX * (up-low);
                        }
                        if(valid_node(node))
                            break;
                    }
                    
                    this->nodes.push_back( node );
                }
                this->nodes.shrink_to_fit();
                
                for( int i = 0; i < max_num; i++ )
                {
                    for( int j = 0; j < max_num; j++ )
                    {
                        if( local_planner(nodes[i], nodes[j]) )
                        {
                            map[i].insert(j);
                            map[j].insert(i);
                        }
                    }
                }
            }
            
            bool connects( PRM<IROBOT> other_prm )
            {
                int other_size = other_prm.get_nodes().size();
                for( int i=0; i < nodes.size(); i++ )
                {
                    for( int j = 0; j < other_size; j++ )
                    {
                        if( local_planner(get_node(i), other_prm.get_node(j)) )
                            return true;
                    }
                }
                return false;
            }
            
            std::vector<MapNode>& get_nodes() { return this->nodes; }
            
            MapNode& get_node( int index ){ return this->nodes[index]; }
            
        private:
            bool valid_node( MapNode& node )
            {
                robot.set_config(node);
                if( obstMgr.intersects(robot) )
                    return false;
                return true;
            }
            
            bool local_planner( MapNode& node1, MapNode& node2 )
            {
                double t = 0.0;
                double interval = 0.05;
                while (t<1.0) {
                    
                    MapNode intermedial = node1 + t*(node2-node1);
                    if( valid_node(intermedial) )
                    t += interval;
                }
                return false;
            }
            
        private:
            std::vector<robotics::Range>& ranges;
            std::vector<MapNode> nodes;
            IROBOT robot;
            robotics::ObstManager& obstMgr;
        };
    }
}


#endif
