//
//  armrobot_planning.cpp
//  RSS2015
//
//  Created by Yu-Han Lyu on 1/26/15.
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
#include "../basics/robotics/cspace_scanner.h"
#include "../basics/algorithms/cspace_decomposer.h"
#include "../basics/algorithms/cell_analysis.h"
#include "../basics/algorithms/PRM.h"

// Basic configuration set up
int     num_samples     = 1000;
bool    outside         = false;
double  epsilon         = 0.2;
bool    write_to_file   = false;
bool    use_oracle      = true;
int     num_joints      = 2;
bool    calc_collision_cells = false;
bool    render_cells    = true;
bool    scan_space      = false;
bool    render_workspace= false;
bool    search_path     = true;

// Set up the robot
std::array<double, 2> arm_lengths{1.0, 1.0};
std::array<double, 2> max_angular_speeds{1.0, 1.0};
ArmRobot<2> robot ( 0.1, M_PI/2, arm_lengths, max_angular_speeds);
std::vector<ArmRobot<2>::CONFIG> path;
std::vector<N2D::sphere> obst_cfgs;
std::vector<algorithms::Cell<ArmRobot<2>>> planning_cells;
std::vector<double> importance;


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
    
    //ArmRobot<2>::CONFIG start( {0.05, 0.05} );
    //ArmRobot<2>::CONFIG goal( {M_PI/2 + .25, 0.11} );
    
    ArmRobot<2>::CONFIG start( {1.1, 0.4} );
    ArmRobot<2>::CONFIG goal( {M_PI/2 + .2, 0.1} );
    
    
    // Set up the world
    std::vector<N2D::Polygon> polies = parse_file("/Users/Yinan/workspace/RSS2014/C++/RSS2015/ji.txt", N2D::v2(M_PI, M_PI));
    //std::vector<N2D::Polygon> polies = parse_file("/Users/IanZhang/Documents/workspace/RSS2015/C++/RSS2015/ji.txt", N2D::v2(M_PI, M_PI));
    //std::vector<N2D::Polygon> polies = parse_file("/Users/yuhanlyu/Documents/RSS2015/C++/RSS2015/ji.txt", geometry::v2(M_PI, M_PI));
    std::vector<robotics::Obstacle> obsts;
    for (N2D::Polygon poly : polies) {
        obsts.emplace_back(poly);
    }
    
    robotics::ObstManager obstacle_manager(obsts);
    
    // CSpace Scanner
    CSpaceScanner<ArmRobot<2>> scanner( robot, obstacle_manager );
    
    if(scan_space)
    {
        if(obst_cfgs.size() == 0)
            obst_cfgs = scanner.scan(robotics::Range(0.0, M_PI - 0.0000001), robotics::Range(0.0, M_PI - 0.0000001), 0.005);
        if(!render_workspace)
        {
            for ( N2D::sphere cfg: obst_cfgs) {
                N2D::render::sphere(cfg, N2D::render::Color( 0,0,0, 160 ), true);
            }
        }
    }
    
    if( planning_cells.size() == 0 && !render_workspace )
    {
        // Decompose C-Space
        algorithms::KDDecomposer<ArmRobot<2>> sampler(robot, obstacle_manager, epsilon);
        //sampler.DecomposeSpace();
        sampler.ShallowDecompose(M_PI/128);
        printf("Free cells: %d\n", (int)sampler.get_free_cells().size());
        
        algorithms::Analyzer<ArmRobot<2>> analyzer( sampler );
        //importance = analyzer.build_simple_weighted_centrality_matrix(M_PI/8, M_PI/64);
        importance = analyzer.build_path_importance_matrix(M_PI/8, M_PI/64);
        printf("Finished building importance matrix!\n");
        
        planning_cells = sampler.get_free_cells();
    }
    
    // Render Cells
    if( render_cells && !render_workspace )
    {
        const double max_weight = *std::max_element(std::begin(importance), std::end(importance));
        const double min_weight = *std::min_element(std::begin(importance), std::end(importance));
        
        for (int i = 0; i < planning_cells.size(); i++)
        {
            ND::sphere<2> temp_sphere = planning_cells[i].sphere();
            N2D::sphere temp2d_sphere( N2D::v2( temp_sphere.center()[0], temp_sphere.center()[1] ), temp_sphere.radius(), N2D::SPHEREMETRIC::LINFTY );
            
            if( planning_cells[i].visited() )
            {
                N2D::render::sphere(temp2d_sphere, N2D::render::Color( 0,220,0, 20 ), true);
                N2D::render::sphere(temp2d_sphere, N2D::render::Color( 0,250,0, 250 ),  false);
            }
            else
            {
                N2D::render::sphere(temp2d_sphere, N2D::render::Color( 50,50,50, 40 ), false);
                if( temp2d_sphere.radius() <= M_PI/64)
                {
                    double weight = importance[i] / ( max_weight/8.0 - min_weight );
                    N2D::render::sphere(temp2d_sphere, N2D::render::Color( 255*weight,0,0, 40 ), true);
                }
            }
        }
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
