//
//  4armrobot_planning.cpp
//  MetricCells
//
//  Created by Yinan Zhang on 5/11/15.
//  Copyright (c) 2015 Yinan Zhang. All rights reserved.
//

#include <unistd.h>
#include <algorithm>
#include <fstream>

#include "armrobot.h"
#include "../basics/robotics/robot.h"
#include "../basics/math/Naive2D/geometry.h"
#include "../basics/math/Naive2D/render.h"
#include "../basics/algorithms/file_reader.h"
#include "../basics/algorithms/Planner.h"
#include "../basics/algorithms/cspace_decomposer.h"
#include "../basics/algorithms/cell_analysis.h"
#include "../basics/algorithms/PRM.h"

// Basic configuration set up
int     num_samples     = 1000;
bool    outside         = false;
double  epsilon         = 0.4;
bool    write_to_file   = false;
bool    use_oracle      = true;
int     num_joints      = 2;
bool    calc_collision_cells = false;
bool    render_cells    = true;
bool    scan_space      = false;
bool    render_workspace= true;
bool    search_path     = true;
#define DIM 4

// Set up the robot
std::array<double, DIM> arm_lengths{0.6, 0.6, 0.6, 0.6 };
std::array<double, DIM> max_angular_speeds{1.0, 1.0, 1.0, 1.0};
double eps = 0.000001;
std::array<robotics::Range, DIM> cfg_ranges = std::array<robotics::Range, DIM>{robotics::Range( -0+eps, M_PI-eps ), robotics::Range( -0+eps, M_PI-eps ), robotics::Range( -0+eps, M_PI-eps ), robotics::Range( -0+eps, M_PI-eps ) };
ArmRobot<DIM> robot ( 0.1, M_PI/2, arm_lengths, max_angular_speeds, cfg_ranges);
std::vector<ArmRobot<DIM>::CONFIG> path;
std::vector<N2D::sphere> obst_cfgs;
std::vector<algorithms::Cell<ArmRobot<DIM>>> planning_cells;
std::vector<double> importance;


void render_robot(ArmRobot<DIM>::CONFIG& config, int r, int g, int b, int a = 250)
{
    robot.set_config(config);
    const std::vector<N2D::Polygon> shapes2 = robot.get_shape();
    for (N2D::Polygon shape : shapes2)
    {
        N2D::render::polygon( shape,  N2D::render::Color(r, g, b, a), false);
    }
}

void keyPressed (unsigned char key, int x, int y) {
    if (key == 'r') { // If they ‘a’ key was pressed
        // Perform action associated with the ‘a’ key
        render_cells = !render_cells;
        glutPostRedisplay();
    }
    if (key == 's') { // If they ‘a’ key was pressed
        // Perform action associated with the ‘a’ key
        scan_space = !scan_space;
        glutPostRedisplay();
    }
    if (key == 'w') { // If they ‘a’ key was pressed
        // Perform action associated with the ‘a’ key
        render_workspace = !render_workspace;
        glutPostRedisplay();
    }
    if (key == 'p') { // If they ‘a’ key was pressed
        // Perform action associated with the ‘a’ key
        search_path = !search_path;
        glutPostRedisplay();
    }
}


void display()
{
    N2D::render::clean_screen();
    
    ArmRobot<DIM>::CONFIG start( {0.7, 0.5, 0.2, 0.30} );
    ArmRobot<DIM>::CONFIG goal( {M_PI/2 - .01, 0.11, 0.31, 0.31} );
    //ArmRobot<DIM>::CONFIG goal( {M_PI - .1, 0.15, 0.1} );
    
    // Set up the world
    //std::vector<N2D::Polygon> polies = parse_file("/Users/Yinan/workspace/RSS2014/C++/RSS2015/ji.txt", N2D::v2(M_PI, M_PI));
    std::vector<N2D::Polygon> polies = parse_file("/Users/IanZhang/Documents/workspace/RSS2015/C++/RSS2015/ji.txt", N2D::v2(M_PI, M_PI));
    //std::vector<N2D::Polygon> polies = parse_file("/Users/yuhanlyu/Documents/RSS2015/C++/RSS2015/ji.txt", geometry::v2(M_PI, M_PI));
    std::vector<robotics::Obstacle> obsts;
    /*
     for (N2D::Polygon poly : polies) {
     obsts.emplace_back(poly);
     }*/
    
    // Simple World
    //std::vector<N2D::v2> p1_points( { N2D::v2( 0.5, 1.2 ), N2D::v2(0.8, 1.2), N2D::v2(0.8, 0.9), N2D::v2( 0.5, 0.9) } );
    std::vector<N2D::v2> p2_points( { N2D::v2( 1.4, 1.4 ), N2D::v2(1.7, 1.4), N2D::v2(1.7, 1.1), N2D::v2( 1.4, 1.1) } );
    std::vector<N2D::v2> p3_points( { N2D::v2( 2.2, 0.6 ), N2D::v2(2.5, 0.6), N2D::v2(2.5, 0.3), N2D::v2( 2.2, 0.3) } );
    std::vector<N2D::v2> p4_points( { N2D::v2( 1.2, 2.45 ), N2D::v2(1.5, 2.45), N2D::v2(1.5, 2.15), N2D::v2( 1.2, 2.15) } );
    //N2D::Polygon p1(p1_points);
    N2D::Polygon p2(p2_points);
    N2D::Polygon p3(p3_points);
    N2D::Polygon p4(p4_points);
    //obsts.push_back(p1);
    obsts.push_back(p2);
    obsts.push_back(p3);
    obsts.push_back(p4);
    
    robotics::ObstManager obstacle_manager(obsts);
    
    if( planning_cells.size() == 0 && !render_workspace )
    {
        // Decompose C-Space
        algorithms::KDDecomposer<ArmRobot<DIM>> sampler(robot, obstacle_manager, epsilon);
        std::clock_t    t0 = std::clock();
        //sampler.DecomposeSpace();
        sampler.AdaptiveDecompose(M_PI/64, M_PI/4096);
        std::clock_t    t1 = std::clock();
        std::cout << "Time cost for decomposing C-space:\n\t" << (t1-t0) / (double)(CLOCKS_PER_SEC / 1000) << "ms\n";
        printf("Free cells: %d\n", (int)sampler.get_free_cells().size());
        
        planning_cells = sampler.get_free_cells();
        
        // Start Planning
        // Set up planner
        
        algorithms::Planner<ArmRobot<DIM>> planner( sampler );
        planner.initialize(sampler, goal, robot);
        
        // Search for path
        std::clock_t    t4 = std::clock();
        path = planner.AStar(robot, start, goal);
        std::clock_t    t5 = std::clock();
        std::cout << "Time cost for finding the shortest path:\n\t" << (t5-t4) / (double)(CLOCKS_PER_SEC / 1000) << "ms\n";
        
        planning_cells = planner.get_cells();
    }
    
    // Render Cells
    if( !render_workspace && render_cells )
    {
        const double max_weight = 0;//*std::max_element(std::begin(importance), std::end(importance));
        const double min_weight = 0;//*std::min_element(std::begin(importance), std::end(importance));
        
        int count = 0;
        for (int i = 0; i < planning_cells.size(); i++)
        {
            ND::sphere<DIM> temp_sphere = planning_cells[i].sphere();
            N2D::sphere temp2d_sphere( N2D::v2( temp_sphere.center()[0], temp_sphere.center()[1] ), temp_sphere.radius(), N2D::SPHEREMETRIC::LINFTY );
            
            if( planning_cells[i].visited() )
            {
                N2D::render::sphere(temp2d_sphere, N2D::render::Color( 0,220,0, 20 ), true);
                N2D::render::sphere(temp2d_sphere, N2D::render::Color( 0,250,0, 250 ),  false);
            }
            else
            {
                N2D::render::sphere(temp2d_sphere, N2D::render::Color( 50,50,50, 40 ), false);
            }
        }
        std::cout << count << '\n';
    }
    
    // Render obstacles
    if(render_workspace)
    {
        const std::vector<robotics::Obstacle> obsts = obstacle_manager.get_obsts();
        for (robotics::Obstacle obst : obsts)
        {
            N2D::Polygon poly( obst.vertices );
            N2D::render::polygon( poly, N2D::render::Color(100, 100, 100, 100));
        }
        robot.set_config(start);
        const std::vector<N2D::Polygon> shapes = robot.get_shape();
        for (N2D::Polygon shape : shapes)
        {
            N2D::render::polygon( shape,  N2D::render::Color(250, 0, 0), false);
        }
        robot.set_config(goal);
        const std::vector<N2D::Polygon> shapes2 = robot.get_shape();
        for (N2D::Polygon shape : shapes2)
        {
            N2D::render::polygon( shape,  N2D::render::Color(250, 0, 0), false);
        }
    }
    
    if(path.size() != 0 && render_workspace)
    {
        render_robot(path[0], 250, 0, 0);
        for( int i = 1; i < path.size(); i++ )
        {
            render_robot(path[i],135, 206, 250, std::max( (int)(float(i)/float(path.size()) * 60), 10));
        }
    }
    
    render_robot(start,250, 0, 0);
    render_robot(goal,0, 0, 250);
    
    N2D::render::flush();
}


int main()
{
    N2D::render::create_window(314, 314, "ArmRobot", N2D::v2(100,100));
    N2D::render::translate_world(0.08,0.08);
    N2D::render::scale_world(90, 90);
    
    N2D::render::set_display_func(display);
    N2D::render::set_keyboard_func(keyPressed);
    
    N2D::render::main_loop();
    return 0;
    
}
