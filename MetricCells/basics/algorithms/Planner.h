//
//  Planner.h
//  RSS2015
//
//  Created by Yinan Zhang on 1/2/15.
//  Copyright (c) 2015 Yinan Zhang. All rights reserved.
//

#ifndef __METRICCELLS__Planner__
#define __METRICCELLS__Planner__

#include <vector>
#include <random>
#include <iostream>
#include <unordered_map>

#include "../robotics/robot.h"
#include "../math/geoND.h"

#include "cspace_decomposer.h"

namespace algorithms {
    template <typename IROBOT> class Cell;
    template <typename IROBOT>
    class Planner
    {
    public:
        static constexpr int DIM = IROBOT::DIM;
        using CONFIG = typename IROBOT::CONFIG;
        Planner( algorithms::KDDecomposer<IROBOT>& decomposer )  : cells(decomposer.move_free_cells()), all_boundaries(decomposer.move_boundaries()) {}
        
        //void initialize( void );
        
        void initialize(algorithms::KDDecomposer<IROBOT>& decomposer, const CONFIG& goal, const IROBOT& robot);
        
        /* @param start: the start state
         * @param goal : goal state
         * @param heuristic_cost: this parameter is a function that takes two parameters: (current_state, goal_state)
         * @param metric: returns metric distance between two states. It takes three parameters ( robot, state1, stat2 )
         * @param robot: will be used as a parameter in metric
         */
        std::vector<CONFIG> AStar(const IROBOT& robot, const CONFIG& start, const CONFIG& goal);
        
        const std::vector<Cell<IROBOT>> get_cells() const { return cells; }
        
    private:
        void get_unique_crn_cfgs(std::vector<const CONFIG*>& result, std::vector<CONFIG>& all_corners) const;
        void build_edges(algorithms::KDDecomposer<IROBOT>& decomposer, const IROBOT& robot);
        void compute_potential(const IROBOT& robot, const CONFIG& goal);
        //bool check_cell_conflict( CONFIG current_cfg, Cell* current_cell, std::unordered_map<robotics::Config, robotics::Config>& from_map, std::unordered_map<robotics::Config, algorithms::Cell*>& cell_map );
        int compute_fix_dim(int lhs, int rhs) const ;
        // Given a boundary config, it is also in a neighbor cell
        int find_next_owner(int cell, const CONFIG& bnd_cfg ) const;
        
        // Given a boundary config, it is also in a neighbor cell
        void find_next_owners(int cell, const CONFIG& bnd_cfg, std::vector<int>& result ) const;
        
        bool is_adjacent(int from, int to) const;
        
        bool compute_bound(int boundary);
        
        double get_smallest_radius() const;
        
        class Info final {
        public:
            static constexpr int CLOSED_MASK = 1;
            static constexpr int OPEN_MASK = 2;
            explicit Info(const CONFIG& arg_cfg) : cfg(arg_cfg) {}
            bool is_open() const { return info & OPEN_MASK; }
            bool is_closed() const { return info & CLOSED_MASK; }
            void set_open() { info |= OPEN_MASK; }
            void set_closed() { info |= CLOSED_MASK; }
            
            CONFIG cfg;
            double g_score = std::numeric_limits<double>::infinity();
            double f_score = std::numeric_limits<double>::infinity();
            double heuristic = std::numeric_limits<double>::infinity(); // g_score + heuristic = f_score
            int    cell = NOT_FOUND;
            int    came_from = NOT_FOUND;
            unsigned char info = 0;
        };
        std::vector<Info> nodes;
        std::vector<Cell<IROBOT>> cells;
        std::vector<Boundary<IROBOT>> all_boundaries;
                                                                                                                       
        std::vector<CONFIG> trace_back( const Info& info, const CONFIG& start, const CONFIG& goal );
        
    public:
        int get_new_info(const CONFIG& cfg) {
            nodes.emplace_back(cfg);
            return (int)nodes.size() - 1;
        }
        Info& get_info(int index) {
            return nodes[index];
        }
        const Info& get_info(int index) const {
            return nodes[index];
        }
        Boundary<IROBOT>& get_boundary(int index) {
            return all_boundaries[index];
        }
        const Boundary<IROBOT>& get_boundary(int index) const {
            return all_boundaries[index];
        }
        Cell<IROBOT>& get_cell(int index) {
            return cells[index];
        }
        const Cell<IROBOT>& get_cell(int index) const {
            return cells[index];
        }
    };
}
//#include "Planner.cpp"
#endif /* defined(CSPACE_DECOMPOSER_H) */
