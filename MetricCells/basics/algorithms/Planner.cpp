

#ifndef __RSS2015__Planner__CPP
#define __RSS2015__Planner__CPP
#include <ctime>
#include <vector>
#include <iostream>
#include <queue>
#include <limits>
#include <array>
#include <algorithm>

#include "Planner.h"
#include "../robotics/robot.h"
#include "../math/geoND.h"
#include "priority_queue.h"

struct Potential_node {
    Potential_node(int arg_cell, double arg_priority) : cell(arg_cell), priority(arg_priority) {}
    bool operator<(const Potential_node& rhs) const
    {
        return priority > rhs.priority;
    }
    int cell;
    double priority;
    
};

//
//Compute the potential function
//
template<typename IROBOT>
void algorithms::Planner<IROBOT>::compute_potential(const IROBOT& robot, const CONFIG& goal)
{
    std::clock_t    t0 = std::clock();
    int goal_cell  = NOT_FOUND;
    for (int i = 0; i < cells.size(); ++i) {
        if( goal_cell == NOT_FOUND && cells[i].contains(goal.vector())) {
            goal_cell = i;
            break;
        }
    }
    if (goal_cell == NOT_FOUND)
        std::cout << "Can't find goal cell\n";
    // Build potential field
    cells[goal_cell].set_potential(0.0);
    std::priority_queue<Potential_node> queue;
    
    for (int index : cells[goal_cell].get_boundaries()) {
        int neighbor = all_boundaries[index].otherside(goal_cell);
        cells[neighbor].set_potential(0.0);
        queue.emplace(neighbor, 0.0);
    }
    // Dijkstra's algorithm
    while( !queue.empty() )
    {
        Potential_node node = queue.top();
        queue.pop();
        int current = node.cell;
        if (node.priority > cells[current].get_potential())
            continue;
        Cell<IROBOT>& current_cell = cells[current];
        for (int index : current_cell.get_boundaries()) {
            Boundary<IROBOT>& boundary = all_boundaries[index];
            int neighbor = boundary.otherside(current);
            Cell<IROBOT>& neighbor_cell = cells[neighbor];
            double tentative_potential = current_cell.get_potential() + boundary.get_weight();
            if (tentative_potential < neighbor_cell.get_potential()) {
                neighbor_cell.set_potential(tentative_potential);
                queue.emplace(neighbor, tentative_potential);
            }
        }
    }
    std::clock_t    t1 = std::clock();
    std::cout << "Time cost for building potential field: \n\t" << (t1 - t0) / (double)(CLOCKS_PER_SEC / 1000) << " ms\n";
}

// Initialize the planner
template<typename IROBOT>
void algorithms::Planner<IROBOT>::initialize(algorithms::KDDecomposer<IROBOT>& decomposer, const CONFIG& goal,
                                             const IROBOT& robot)
{
    compute_potential(robot, goal);
}

struct AStar_node {
    AStar_node(int arg_info_index, double arg_priority) : info_index(arg_info_index), priority(arg_priority) {}
    bool operator<(const AStar_node& rhs) const
    {
        return priority > rhs.priority;
    }
    int info_index;
    double priority;
    
};

//
//A star search for finding path
//
template<typename IROBOT>
std::vector<typename IROBOT::CONFIG> algorithms::Planner<IROBOT>::AStar(const IROBOT& robot, const CONFIG& start, const CONFIG& goal)
{
    
    double smallest_radius = get_smallest_radius();
    // Find the cells containing the start and the goal
    int start_cell = NOT_FOUND, goal_cell = NOT_FOUND;
    for (int i = 0; i < cells.size(); ++i) {
        if( cells[i].contains(goal.vector()))
            goal_cell = i;
        if (cells[i].contains(start.vector())) {
            start_cell = i;
            std::cout << "Start cell potential : " << cells[i].get_potential() << '\n';
        }
    }
    for (Boundary<IROBOT>& boundary : all_boundaries)
        boundary.reset();
    
    if( start_cell == NOT_FOUND)
        throw "start config is not in any cell!";
    else if (goal_cell == NOT_FOUND )
        throw "goal config is not in any cell!";
    
    // Initialize for the nodes
    nodes.clear();
    nodes.reserve(1000000);
    int start_node = get_new_info(start);
    get_info(start_node).came_from = -1;
    int goal_node = get_new_info(goal);
    get_info(goal_node).cell = goal_cell;
    int metric_query = 0;
    std::cout << "Initialize the queue of AStar\n";
    // Initialize for the queue
    std::priority_queue<AStar_node> openset;
    {
        double g_score_start = 0.0;
        int current_cell_index = start_cell;
        double interval = (std::log(cells[current_cell_index].radius() / smallest_radius) + 1) * smallest_radius;
        for (int boundary_index : cells[current_cell_index].get_boundaries()) {
            int begin, end;
            get_boundary(boundary_index).get_boundary_configs(*this, interval, begin, end);
            for (int vertex = begin; vertex < end; ++vertex) {
                double tentative_g_score = g_score_start + IROBOT::metric(robot, start, nodes[vertex].cfg);
                nodes[vertex].came_from = start_node;
                nodes[vertex].g_score = tentative_g_score;
                nodes[vertex].f_score = tentative_g_score + cells[current_cell_index].get_potential();
                nodes[vertex].heuristic = cells[current_cell_index].get_potential();
                nodes[vertex].cell = start_cell;
                nodes[vertex].set_open();
                openset.emplace(vertex, nodes[vertex].f_score);
                ++metric_query;
            }
            
        }
    }
    std::cout << "Start AStar\n";
    // Starting A*
    while (!openset.empty()) {
        int current_node = openset.top().info_index;
        openset.pop();
        if (nodes[current_node].is_closed())
            continue;
        get_info(current_node).set_closed();
        if (current_node == goal_node) {
            std::cout << "Find Goal by expolring " << nodes.size() << " nodes with " << metric_query << " queries\n";
            return this->trace_back(nodes[current_node], start, goal);
        }
        
        int prev_cell;
        int current_cell;
        // Special case for start's neighbor
        if (nodes[current_node].came_from == start_node) {
            prev_cell = nodes[current_node].cell;
            current_cell = prev_cell;
        } else {
            prev_cell = nodes[nodes[current_node].came_from].cell;
            current_cell = find_next_owner(prev_cell, nodes[nodes[current_node].came_from].cfg); // search the neighbors of prev cell to see which cell contians current_cfg, that is the current_cell. ( there is only one such cell )
            
        }
        get_info(current_node).cell = current_cell;
        
        // time to check if the path has common cells.
        int neighbor_cell_index = find_next_owner(current_cell, nodes[current_node].cfg); // get the cell that also contains current config
        if( neighbor_cell_index == NOT_FOUND ) {
            continue;
        }
        cells[current_cell].set_visited();
        cells[neighbor_cell_index].set_visited();
        // Special case for goal's neighbor
        if (cells[neighbor_cell_index].contains(goal.vector())) {
            if( !nodes[goal_node].is_open() ) {
                nodes[goal_node].set_open();
                nodes[goal_node].g_score = get_info(current_node).g_score + IROBOT::metric(robot, nodes[current_node].cfg, goal);
                nodes[goal_node].came_from = current_node;
                openset.emplace(goal_node, nodes[goal_node].g_score);
                ++metric_query;
            } else {
                double distance_to_goal = IROBOT::metric(robot, nodes[current_node].cfg, goal);
                ++metric_query;
                if (get_info(current_node).g_score + distance_to_goal < nodes[goal_node].g_score) {
                    // 1. remove the neighbor_cfg from openset
                    // 2. re-add neighbor_cfg to openset with new priority.
                    nodes[goal_node].g_score = get_info(current_node).g_score + distance_to_goal;
                    nodes[goal_node].came_from = current_node;
                    openset.emplace(goal_node, nodes[goal_node].g_score);
                }
            }
            continue;
        }
        
        int neighbor_count = 0;
        // Compute successors
        for (int boundary_index : cells[neighbor_cell_index].get_boundaries()) {
            double interval = (std::log(cells[neighbor_cell_index].radius() / smallest_radius) + 1) * smallest_radius;
            int begin, end;
            get_boundary(boundary_index).get_boundary_configs(*this, interval, begin, end);
            for (int neighbor = begin; neighbor < end; ++neighbor) {
                if (nodes[neighbor].is_closed())    // already in closed set
                    continue;
                if (nodes[neighbor].is_open() && nodes[neighbor].g_score <= nodes[current_node].g_score)
                    continue;
                ++neighbor_count;
                double tentative_g_score = get_info(current_node).g_score + IROBOT::metric(robot, nodes[current_node].cfg, nodes[neighbor].cfg);
                ++metric_query;
                if (!nodes[neighbor].is_open() || tentative_g_score < nodes[neighbor].g_score)
                {
                    nodes[neighbor].came_from = current_node;
                    nodes[neighbor].g_score = tentative_g_score;
                    
                    nodes[neighbor].f_score = nodes[neighbor].g_score + cells[neighbor_cell_index].get_potential();
                    nodes[neighbor].heuristic = cells[neighbor_cell_index].get_potential();
                    if( !nodes[neighbor].is_open() ) {
                        //get_info(neighbor).f_score = get_info(neighbor).g_score + metric(robot, get_info(neighbor).cfg, goal);
                        nodes[neighbor].set_open();
                        openset.emplace(neighbor, nodes[neighbor].f_score);
                    } else{
                        //get_info(neighbor).f_score = get_info(neighbor).g_score + get_info(neighbor).heuristic;
                        // 1. remove the neighbor_cfg from openset
                        // 2. re-add neighbor_cfg to openset with new priority.
                        openset.emplace(neighbor, nodes[neighbor].f_score); // (two steps are done inside)
                    }
                }
            }
        }
        //std::cout << neighbor_count << '\n';
    }
    throw "The path doesn't exist";
}

template<typename IROBOT>
std::vector<typename IROBOT::CONFIG> algorithms::Planner<IROBOT>::trace_back( const Info& info, const CONFIG& start, const CONFIG& goal )
{
    Info curr_info = info;
    std::vector<CONFIG> reverse;
    reverse.emplace_back(goal);
    while ( curr_info.came_from != NOT_FOUND) {
        reverse.emplace_back(curr_info.cfg);
        curr_info = get_info(curr_info.came_from);
    };
    reverse.emplace_back(start);
    return reverse;
}

template<typename IROBOT>
int algorithms::Planner<IROBOT>::find_next_owner(int cell, const CONFIG& bnd_cfg ) const
{
    for (int index : get_cell(cell).get_boundaries()) {
        int neighbor = get_boundary(index).otherside(cell);
        if(get_cell(neighbor).on_boundary(bnd_cfg.vector(), 1e-15))
            return neighbor;
    }
    return NOT_FOUND;
}

template<typename IROBOT>
void algorithms::Planner<IROBOT>::find_next_owners(int cell, const CONFIG& bnd_cfg, std::vector<int>& result) const
{
    result.clear();
    for (int index : get_cell(cell).get_boundaries()) {
        int neighbor = get_boundary(index).otherside(cell);
        if(get_cell(neighbor).on_boundary(bnd_cfg, 1e-15))
            result.emplace_back( neighbor );
        
    }
}

template<typename IROBOT>
double algorithms::Planner<IROBOT>::get_smallest_radius() const
{
    double smallest_radius = std::numeric_limits<double>::infinity();
    for (const Cell<IROBOT>& cell : cells) {
        if(cell.radius() < smallest_radius)
            smallest_radius = cell.radius();
    }
    return smallest_radius;
}
#endif
