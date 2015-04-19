//
//  file_reader.h
//  RSS2015
//
//  Created by Yinan Zhang on 1/8/15.
//  Copyright (c) 2015 Yinan Zhang. All rights reserved.
//

#ifndef RSS2015_file_reader_h
#define RSS2015_file_reader_h

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "../math/geoND.h"
#include "../math/Naive2D/geometry.h"

using namespace std;

N2D::Polygon parse_line( string line, N2D::v2 scale = N2D::v2(1,1) );
N2D::Polygon parse_polygon(const std::vector<std::string> tokens, N2D::v2 scale = N2D::v2(1,1));
std::vector<N2D::Polygon> parse_file( const string filename, N2D::v2 scale = N2D::v2(1,1) );

N2D::Polygon parse_line( string line, N2D::v2 scale )
{
    std::istringstream buf(line);
    std::istream_iterator<std::string> beg(buf), end;
    std::vector<std::string> tokens(beg, end); // done!
    return parse_polygon(tokens, scale);
}


// input example: "rect -50 50 100 100"
N2D::Polygon parse_polygon(const std::vector<std::string> tokens, N2D::v2 scale)
{
    std::vector<N2D::v2> ptlist;
    for (int i = 0; i < tokens.size(); i+= 2) {
        double x = std::stod(tokens[i], NULL) * scale.x;   // scale from the normal size to world world size
        double y = std::stod(tokens[i+1], NULL) * scale.y;
        int x_int = x*100;
        int y_int = y*100;
        x = x_int / 100.0;
        y = y_int / 100.0;
        ptlist.push_back( N2D::v2(x,y) );
    }
    return N2D::Polygon( ptlist );
}

std::vector<N2D::Polygon> parse_file( const string filename, N2D::v2 scale )
{
    std::vector<N2D::Polygon> polies;
    ifstream infile(filename);
    string line;
    while (std::getline(infile, line))
    {
        if(line[0] != '#')
            polies.push_back(parse_line(line, scale));
    }
    return polies;
}

template <size_t N>
static void write2file( const string filename, const std::vector<ND::vec<N>> data )
{
    std::ofstream output_file(filename);
    if (!output_file.is_open())
        std::cout << "Can't open the file\n";
    for(ND::vec<N> cfg : data) {
        for( int i = 0; i < N; i++ )
            output_file << cfg[i] << '\t';
        output_file << '\n';
    }
    
    output_file.flush();
    output_file.close();
}

template <size_t N>
static std::vector<ND::vec<N>> read_vec( const string filename )
{
    std::vector<ND::vec<N>> cfgs;
    ifstream infile(filename);
    string line;
    while (std::getline(infile, line))
    {
        if(line[0] == '#')
            continue;
        std::istringstream buf(line);
        std::istream_iterator<std::string> beg(buf), end;
        std::vector<std::string> tokens(beg, end); // done!
        ND::vec<N> cfg;
        for (int i = 0; i < tokens.size(); i+= 1) {
            double x = std::stod(tokens[i], NULL);
            cfg[i] = x;
        }
        cfgs.push_back(cfg);
        
    }
    return cfgs;
}

/*
int main()
{
    std::vector<N2D::Polygon> polies = parse_file("Users/IanZhang/Documents/workspace/RSS2015/C++/RSS2015/scene.txt");
    std::cout << polies.size() << std::endl;
    return 0;
}*/


#endif
