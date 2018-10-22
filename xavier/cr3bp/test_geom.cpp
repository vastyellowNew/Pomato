/*************************************************************************
POMATO: POincare MAp TOpology: Extraction and visualization of the
        topological structure of the circular restricted three-body 
        problem for orbital mechanics applications. 

Authors: Wayne Schlei and Xavier Tricoche 

Copyright (c) 2013-2018, Purdue University 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**************************************************************************/


#include "geom.hpp"
#include <math/fixed_vector.hpp>
#include <string>
#include <fstream>

const double invalid_double = std::numeric_limits<double>::max();

int main(int argc, char* argv[])
{
    std::fstream in(argv[1], std::ios::in);
    std::vector<nvis::vec2> points;
    std::string tmp;
    std::getline(in, tmp);
    while (!in.eof()) {
        double x, y=invalid_double;
        char c;
        in >> x >> c >> y;
        std::cerr << "read: " << x << c << y << '\n';
        if (y == invalid_double) {
            break;
        }
        points.push_back(nvis::vec2(x,y));
    }
    in.close();
    
    nvis::vec2 p = orbital::innermost<nvis::vec2>(points);
    std::cout << "innermost point is " << p << '\n';
    
}
