#ifndef CSPACE_DECOMPOSER_CPP
#define CSPACE_DECOMPOSER_CPP

#include <limits>
#include <ctime>
#include <utility>
#include <fstream>
#include <queue>
#include <stack>
#include <array>
#include "cspace_decomposer.h"
#include "../math/geoND.h"
#include "cell_analysis.h"

#define DEBUG_NODE 1281

template <typename IROBOT>
algorithms::KDDecomposer<IROBOT>::KDDecomposer(IROBOT& robot,
                                               robotics::ObstManager& obstacle_manager,
                                               double epsilon)
: robot_ (robot), obstacle_manager_ (obstacle_manager), method_ (algorithms::SamplingMethod::KD), epsilon_ (epsilon)
{
    const std::array<double, DIM>& param_speeds = robot_.get_parameter_speeds();
    min_param_speed_ = std::numeric_limits<int>::max();
    for (unsigned int i = 0; i < param_speeds.size(); i++)
        if (param_speeds[i] < min_param_speed_)
            min_param_speed_ = param_speeds[i];
}


/**
 * Decompose the space
 */
template <typename IROBOT>
void algorithms::KDDecomposer<IROBOT>::ShallowDecompose( double min_radius )
{
    const std::array<robotics::Range, DIM>& ranges = robot_.get_config_ranges();
    ND::vec<DIM> center;
    for (int i = 0; i < DIM; i++) {
        center[i] = (std::get<0>(ranges[i]) + std::get<1>(ranges[i])) * 0.5;
    }
    radius_array[0] = (std::get<1>(ranges[0]) - std::get<0>(ranges[0])) * 0.5;
    for (int i = 1; i < sizeof(radius_array) / sizeof(radius_array[0]); ++i)
        radius_array[i] = radius_array[i-1] * 0.5;
    // NOTE: Assumes hypercube config space!! (Commented code more general)
    nodes.reserve(360000000);
    cells.reserve(3000000);
    root_index = get_new_node(center, 0);
    std::stack<int> stack;
    stack.emplace(root_index);
    while (!stack.empty())
    {
        int node_index = stack.top();
        stack.pop();
        
        double cell_radius = radius_array[get_node(node_index).depth()];
        
        robot_.set_config(get_node(node_index).center());
        double dist_to_obsts = obstacle_manager_.dist_to_obsts(robot_);
        if (dist_to_obsts > 0)
        {
            // Create free space ball
            double initial_guess = CalcFreeCellRadius(dist_to_obsts);
            double radius = robot_.free_space_oracle(obstacle_manager_, initial_guess /* lower bound */, 2 * initial_guess /* upper bound */);
            if (robot_.subconvex_radius(obstacle_manager_, radius, cell_radius) )
            {
                get_node(node_index).set_covered();
                get_node(node_index).set_free();
                if (!MERGE_CELLS)
                    get_node(node_index).set_cell(get_new_cell(get_node(node_index).center(), cell_radius, node_index));
            }
            else if ( cell_radius > min_radius )
            {
                SplitCell(node_index);
                for (int i = 0; i < NUM_SPLITS; ++i)
                    stack.emplace(get_node(node_index).get_children() + i);
            }
            else
            {
                get_node(node_index).set_covered();
                //get_node(node_index).set_mix();
                get_node(node_index).set_free();
            }
        }
        else
        {
            const double penetration = obstacle_manager_.penetration(robot_);
            if ( penetration >= cell_radius  )
            {
                continue;
            }
            /*
            else if (cell_radius > penetration )
            {
                get_node(node_index).set_covered();
                get_node(node_index).set_free();
            }*/
            else if (cell_radius > min_radius )
            {
                SplitCell(node_index);
                for (int i = 0; i < NUM_SPLITS; ++i)
                    stack.emplace(get_node(node_index).get_children() + i);
            }
        }
        
    }
    
    nodes.shrink_to_fit();
    cells.shrink_to_fit();
    clean_tree();
    build_edges(false);
}

/**
 * Adaptively decompose the space
 */
template <typename IROBOT>
void algorithms::KDDecomposer<IROBOT>::AdaptiveDecompose(double larg_radius, double min_radius)
{
    this->ShallowDecompose(larg_radius);
    algorithms::Analyzer<IROBOT> analyzer( *this );
    //importance = analyzer.build_simple_weighted_centrality_matrix(M_PI/8, M_PI/64);
    std::vector<double> importance = analyzer.build_path_importance_matrix(larg_radius * 7, larg_radius);
    // node_index --> importance
    std::unordered_map<int, double> impt_map;
    std::stack<int> stack;
    for( int i = 0; i < this->cells.size(); i++ )
    {
        double cell_radius = get_cell(i).radius();
        if(cell_radius > larg_radius*2)
            continue;
        stack.push(get_cell(i).node_id);
        impt_map[get_cell(i).node_id] = importance[i];
    }
    this->cells.clear();
    
    this->DecomposeSubspaces( stack, larg_radius, min_radius, impt_map );
    this->cells.reserve(nodes.size());
    nodes.shrink_to_fit();
    cells.shrink_to_fit();
    clean_tree();
    build_edges();
}
template <typename IROBOT>
void algorithms::KDDecomposer<IROBOT>::DecomposeSubspaces(std::stack<int>& stack, double large_radius, double min_radius, std::unordered_map<int, double>& impt_map)
{
    while (!stack.empty())
    {
        int node_index = stack.top();
        stack.pop();
        
        double cell_radius = radius_array[get_node(node_index).depth()];
        // 3R Arm:
        //double local_min_radius = std::min( large_radius/2.0, std::max(min_radius, large_radius/(impt_map[node_index]*10)));
        // 4R Arm:
        double local_min_radius = std::min( large_radius/2.0, std::max(min_radius, large_radius/(impt_map[node_index])));
        
        robot_.set_config(get_node(node_index).center());
        double dist_to_obsts = obstacle_manager_.dist_to_obsts(robot_);
        if (dist_to_obsts > 0)
        {
            get_node(node_index).reset_info();
            // Create free space ball
            double initial_guess = CalcFreeCellRadius(dist_to_obsts);
            double radius = robot_.free_space_oracle(obstacle_manager_, initial_guess /* lower bound */, 2 * initial_guess /* upper bound */);
            if (robot_.subconvex_radius(obstacle_manager_, radius, cell_radius))
            {
                get_node(node_index).set_covered();
                get_node(node_index).set_free();
                if (!MERGE_CELLS)
                    get_node(node_index).set_cell(get_new_cell(get_node(node_index).center(), cell_radius, node_index));
            }
            else if( cell_radius > local_min_radius )
            {
                get_node(node_index).set_children( NOT_FOUND );
                SplitCell(node_index);
                for (int i = 0; i < NUM_SPLITS; ++i)
                {
                    impt_map[get_node(node_index).get_children()+i] = impt_map[node_index];
                    stack.emplace(get_node(node_index).get_children()+i);
                }
            }
        }
        else
        {
            get_node(node_index).reset_info();
            const double penetration = obstacle_manager_.penetration(robot_);
            if ( penetration / robot_.get_max_speed() >= (PENETRATION_CONSTANT / min_param_speed_) * cell_radius  )
            {
                continue;
            }
            else if (cell_radius > local_min_radius )
            {
                get_node(node_index).set_children( NOT_FOUND );
                SplitCell(node_index);
                for (int i = 0; i < NUM_SPLITS; ++i)
                {
                    impt_map[get_node(node_index).get_children() + i] = impt_map[node_index];
                    stack.emplace(get_node(node_index).get_children() + i);
                }

            }
        }
        
    }
}

/**
 * Decompose the space
 */
template <typename IROBOT>
void algorithms::KDDecomposer<IROBOT>::DecomposeSpace() {
    std::clock_t    t0 = std::clock();
    const std::array<robotics::Range, DIM>& ranges = robot_.get_config_ranges();
    ND::vec<DIM> center;
    for (int i = 0; i < DIM; i++) {
        center[i] = (std::get<0>(ranges[i]) + std::get<1>(ranges[i])) * 0.5;
    }
    radius_array[0] = (std::get<1>(ranges[0]) - std::get<0>(ranges[0])) * 0.5;
    for (int i = 1; i < sizeof(radius_array) / sizeof(radius_array[0]); ++i)
        radius_array[i] = radius_array[i-1] * 0.5;
    // NOTE: Assumes hypercube config space!! (Commented code more general)
    nodes.reserve(360000000);
    cells.reserve(3000000);
    root_index = get_new_node(center, 0);
    std::stack<int> stack;
    stack.emplace(root_index);
    while (!stack.empty())
    {
        int node_index = stack.top();
        stack.pop();

        robot_.set_config(get_node(node_index).center());
        double dist_to_obsts = obstacle_manager_.dist_to_obsts(robot_);
        if (dist_to_obsts > 0)
        {
            // Create free space ball
            double initial_guess = CalcFreeCellRadius(dist_to_obsts);
            double radius = robot_.free_space_oracle(obstacle_manager_, initial_guess /* lower bound */, 2 * initial_guess /* upper bound */);
            if (robot_.subconvex_radius(obstacle_manager_, radius, radius_array[get_node(node_index).depth()]) )
            {
                get_node(node_index).set_covered();
                get_node(node_index).set_free();
                if (!MERGE_CELLS)
                    get_node(node_index).set_cell(get_new_cell(get_node(node_index).center(), radius_array[get_node(node_index).depth()], node_index));
            }
            else if (dist_to_obsts >= epsilon_ * 0.5 || radius_array[get_node(node_index).depth()] >= robot_.get_s_small(epsilon_ * 0.5) )
            {
                SplitCell(node_index);
                for (int i = 0; i < NUM_SPLITS; ++i)
                    stack.emplace(get_node(node_index).get_children() + i);
            }
        }
        else if (obstacle_manager_.penetration(robot_) / robot_.get_max_speed() >= (PENETRATION_CONSTANT / min_param_speed_) * radius_array[get_node(node_index).depth()] )
        {
            continue;
        }
        else if (radius_array[get_node(node_index).depth()] >= robot_.get_s_small(epsilon_ * 0.5))
        {
            SplitCell(node_index);
            for (int i = 0; i < NUM_SPLITS; ++i)
                stack.emplace(get_node(node_index).get_children() + i);
        }
    }
    
    nodes.shrink_to_fit();
    cells.shrink_to_fit();
    clean_tree();
    build_edges();
}

template<typename IROBOT>
void algorithms::KDDecomposer<IROBOT>::clean_tree()
{
    std::clock_t    t0 = std::clock();
    delete_unsafe_cells(root_index);
    if (MERGE_CELLS) {
        merge_free_cells(root_index);
        create_free_cells(root_index);
    }
    std::clock_t    t1 = std::clock();
}

template<typename IROBOT>
bool algorithms::KDDecomposer<IROBOT>::delete_unsafe_cells(int node)
{
    if (get_node(node).get_children() == NOT_FOUND)
        return !get_node(node).is_free();
    bool contain_free_cell = false;
    for (int i = get_node(node).get_children(); i < get_node(node).get_children() + NUM_SPLITS; ++i)
        if (!delete_unsafe_cells(i))
            contain_free_cell = true;
    if (!contain_free_cell) {
        get_node(node).set_children(NOT_FOUND);
        return true;
    }
    return false;
}

template<typename IROBOT>
bool algorithms::KDDecomposer<IROBOT>::merge_free_cells(int node)
{
    if (get_node(node).is_free())
        return true;
    if (get_node(node).get_children() == NOT_FOUND)
        return false;
    bool is_free = true;
    for (int i = get_node(node).get_children(); i < get_node(node).get_children() + NUM_SPLITS; ++i)
        if (!merge_free_cells(i))
        {
            is_free = false;
        }
    if (is_free) {
        get_node(node).set_free();
        get_node(node).set_children(NOT_FOUND);
    }
    return is_free;
}

template<typename IROBOT>
void algorithms::KDDecomposer<IROBOT>::create_free_cells(int node)
{
    if (get_node(node).is_free()) {
        get_node(node).set_cell(get_new_cell(get_node(node).center(), radius_array[get_node(node).depth()], node));
        return;
    }
    if (get_node(node).get_children() != NOT_FOUND)
        for (int i = get_node(node).get_children(); i < get_node(node).get_children() + NUM_SPLITS; ++i)
            create_free_cells(i);
}

/**
 * Build the connection between cells
 * Time cost for creating corners of each cell:9327.53ms
 * Time cost for finding neighbors of each cell with 8014392 boundaries:24505ms
 */
template<typename IROBOT>
void algorithms::KDDecomposer<IROBOT>::build_edges(bool clear_nodes)
{
    all_boundaries.clear();
    std::vector<CONFIG> all_corners;
    std::vector<const CONFIG*> corner_configs;
    all_corners.reserve(1000000);
    corner_configs.reserve(200000);
    get_unique_crn_cfgs(corner_configs, all_corners);
    
    std::clock_t    t0 = std::clock();
    all_boundaries.reserve(300000);
    // Build all boundaries
    std::vector<int> containing_cells;
    for (const CONFIG* cfg : corner_configs) {
        containing_cells.clear();
        ContainedCells(cfg->vector(), containing_cells);
        if(containing_cells.size() == 0)
            continue;
        for (int i = 0; i < containing_cells.size() - 1; ++i) {
            for (int j = i + 1; j < containing_cells.size(); ++j) {
                if (is_adjacent(containing_cells[i], containing_cells[j]))
                    continue;
                int fix_dim = compute_fix_dim(containing_cells[i], containing_cells[j]);
                // If the fixed dimension is not one, then it is degenerate case
                if (fix_dim == -1)
                    continue;
                int boundary_index = get_new_boundary(containing_cells[i], containing_cells[j]);
                cells[containing_cells[i]].add_boundary(boundary_index);
                cells[containing_cells[j]].add_boundary(boundary_index);
                double weight = IROBOT::heuristic(robot_, CONFIG(cells[containing_cells[i]].center()), CONFIG(cells[containing_cells[j]].center()));
                all_boundaries[boundary_index].set_weight(weight);
                all_boundaries[boundary_index].set_fix_dim(fix_dim);
            }
        }
    }
    for (int i = 0; i < all_boundaries.size(); ++i) {
        compute_bound(i);
    }

    
    std::clock_t    t2 = std::clock();
    all_boundaries.shrink_to_fit();
    if(clear_nodes) nodes.clear();
    nodes.shrink_to_fit();
    for (Cell<IROBOT>& cell : cells)
        cell.shrink_to_fit();
    std::clock_t    t3 = std::clock();
}

template<typename IROBOT>
void algorithms::KDDecomposer<IROBOT>::get_unique_crn_cfgs(std::vector<const CONFIG*>& result, std::vector<CONFIG>& all_corners) const
{
    std::clock_t    t0 = std::clock();
    for (const Cell<IROBOT>& cell : cells) {
        cell.get_corner_configs(all_corners);
    }
    all_corners.shrink_to_fit();
    const CONFIG** temp = new const CONFIG*[all_corners.size()];
    
    for (int i = 0; i < all_corners.size(); ++i)
        temp[i] = &all_corners[i];
    
    std::sort(temp, temp + all_corners.size(),  IROBOT::cmp);
    const CONFIG** end = unique( temp, temp + all_corners.size(), IROBOT::approximate);
    result.insert(result.begin(), temp, end);
    delete[] temp;
    std::clock_t    t1 = std::clock();
}
#endif