//
//  cell_analysis.h
//  MetricCells
//
//  Created by Yinan Zhang on 4/20/15.
//  Copyright (c) 2015 Yinan Zhang. All rights reserved.
//

#ifndef MetricCells_cell_analysis_h
#define MetricCells_cell_analysis_h

#include <vector>

#include "priority_queue.h"
using namespace data_structure;

#define INFTY 100000000

namespace algorithms {
    template <typename IROBOT> class Cell;
    
    template <typename IROBOT>
    class Analyzer
    {
    public:
        Analyzer(algorithms::KDDecomposer<IROBOT>& decomposer) :
                            cells(decomposer.get_free_cells()),
                            all_boundaries(decomposer.get_boundaries())
        {
            min_radius = INFTY;
            for (Cell<IROBOT>& cell : cells) {
                if( cell.radius() < min_radius )
                    min_radius = cell.radius();
            }
        }
        
        
        std::vector<double> build_simple_weighted_centrality_matrix( double large_size_threshold, double bnd_size_threshold )
        {
            vector<double> importance = vector<double>( cells.size(), 0 );
            vector<int> row( cells.size(), INFTY );
            vector<vector<int>> origin_matrix( cells.size(), row );
            vector<vector<int>> matrix( cells.size(), row );
            build_simple_dist_matrix( large_size_threshold, origin_matrix );
            
            for( int i = 0; i < cells.size(); i++ )
            {
                if(get_cell(i).radius() > bnd_size_threshold) continue;
                
                build_simple_dist_matrix( large_size_threshold, matrix, i );
                importance[i] = matrix_diff_sum( origin_matrix, matrix, i );
            }
            
            return importance;
        }
        

    private:
        std::vector<Cell<IROBOT>> cells;
        std::vector<Boundary<IROBOT>> all_boundaries;
        double min_radius;
        
        
        
    private:
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
        
        
        void build_simple_dist_matrix( double cell_size_threshold, vector<vector<int>>& matrix, int exception = -1 )
        {
            for( int i = 0; i < matrix.size(); i++ )
            {
                for( int j=0; j<matrix[0].size(); j++ )
                {
                    matrix[i][j] = INFTY;
                }
                
                if( i != exception ) matrix[i][i] = 0;
            }
            
            for( int i= 0; i < cells.size(); i++ )
            {
                if(get_cell(i).radius() >= cell_size_threshold)
                {
                    Dijkstra( i, matrix[i], exception );
                }
            }
        }
        
        
        void Dijkstra( int source_idx, std::vector<int>& dists, int exception )
        {
            for (int i=0; i<dists.size(); i++) {
                dists[i] = INFTY;
            }
            dists[source_idx] = 0;
            
            data_structure::priority_queue<int> Q;
            
            for( int i = 0; i<cells.size(); i++ )
            {
                Q.push( i, dists[i] );
            }
            
            while (!Q.is_empty()) {
                int current = Q.pop();
                if( current == exception ) continue;
                
                
                for( int boundary_index : cells[current].get_boundaries() )
                {
                    int neighbor = get_boundary(boundary_index).otherside( current );
                    if( neighbor == exception ) continue;
                    int alt = dists[current] + int(length(current, neighbor)/min_radius);
                    if( alt < dists[neighbor] )
                    {
                        dists[neighbor] = alt;
                        Q.push( neighbor, alt );
                    }
                }
            }
        }
        
        double length( int u, int v )
        {
            return (get_cell(u).radius()+get_cell(v).radius());
        }
        
        int matrix_diff_sum(vector< vector<int> >& mtx1,
                            vector<vector<int> >& mtx2,
                            int exception)
        {
            int sum = 0;
            for (int i = 0; i < mtx1.size(); i++)
            {
                if(i==exception)
                    continue;
                for (int j = 0; j < mtx1[0].size(); j++)
                {
                    if( j == exception || i == j )
                        continue;
                    if( mtx1[i][j] < INFTY )
                    {
                        double weight =  std::pow(get_cell(i).radius() * get_cell(j).radius(),2);
                        double diff   = std::fabs( mtx1[i][j]-mtx2[i][j] );
                        sum += diff*weight;
                        if(sum < 0)
                        {
                            sum = INFTY+1;
                        }
                    }
                    if(sum > INFTY) return sum;
                }
            }
            return sum;
        }
        
    };
    
}

#endif
