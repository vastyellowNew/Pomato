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


#include <iostream>
#include <vector>
#include <list>
#include <util/wall_timer.hpp>
#include <complex>
#include <sstream>
#include <limits>

// math
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/intersection.hpp>


using namespace xavier;

//Type definitions
typedef nvis::fixed_vector<nvis::vec2,2>  SegType;
typedef LineSegmentIntersection2D<SegType>  Intersector;

//Function to gather a string for an output option from an Intersector
std::string getSolutionType(const Intersector& lsx)
{
    std::string output;
    switch(lsx.getSolutionType()) {
        case Intersector::NONE :
            output = "NONE";
            break;
        case Intersector::PARALLEL :
            output = "PARALLEL";
            break;
        case Intersector::COLLINEAR_DISJOINT :
            output = "COLLINEAR_DISJOINT";
            break;
        case Intersector::COLLINEAR :
            output = "COLLINEAR";
            break;
        case Intersector::FOUND :
            output = "FOUND";
            break;
    }
    return output;
}
//Printing the test output
void printTest(SegType& s0, SegType& s1, const Intersector& lsx)
{
    std::string solType;
    solType = getSolutionType(lsx);
    std::cout << "-----------------------------------------------------------------\n";
    std::cout << "  Seg A: " << s0[0] << " , " << s0[1] << "\n";
    std::cout << "  Seg B: " << s1[0] << " , " << s1[1] << "\n";
    std::cout << "  SOLUTION :: " << solType << "\n";
    std::pair<bool,double> out(false,1000.0);
    switch(lsx.getSolutionType()) {
        case Intersector::NONE :
            out = lsx.getMinDistance();
            std::cout << "    Overlap = " << out.first << "  Min Distance = " << out.second << "\n";
            break;
        case Intersector::PARALLEL :
            //Write the parallel distance
            std::cout << "     Parallel Distance = " << lsx.getParallelDistance() << "\n";
            out = lsx.getMinDistance();
            std::cout << "    Overlap = " << out.first << "  Min Distance = " << out.second << "\n";
            break;
        case Intersector::COLLINEAR_DISJOINT :
            //Output some values
            out = lsx.getMinDistance();
            std::cout << "    Overlap = " << out.first << "  Min Distance = " << out.second << "\n";
            break;
        case Intersector::COLLINEAR :
            //Output some values
            out = lsx.getMinDistance();
            std::cout << "    Overlap = " << out.first << "  Min Distance = " << out.second << "\n";
            break;
        case Intersector::FOUND :
            //Write solution
            std::cout << "     Intersection = " << lsx.getIntersection() << "\n";
            std::cout << "       tau(SegA) = " << lsx.getIntersectionParameter() << "\n";
            std::cout << "       tau(SegB) = " << lsx.getIntersectionParameter(1) << "\n";
            break;
    }
    std::cout << "-----------------------------------------------------------------\n";
    
    
}
int main()
{

    //Test Points
    std::vector<nvis::vec2> points;
    points.push_back(nvis::vec2(-1.0,-1.0));
    points.push_back(nvis::vec2(1.0,1.0));
    points.push_back(nvis::vec2(0.0,1.0));
    points.push_back(nvis::vec2(1.0,2.0));
    points.push_back(nvis::vec2(0.5,-1.0));
    
    points.push_back(nvis::vec2(0.1,0.1));
    points.push_back(nvis::vec2(0.5,0.5));
    points.push_back(nvis::vec2(2.0,2.0));
    points.push_back(nvis::vec2(5.0,5.0));
    
    
    
    //Test Case 1: Parallel
    SegType seg0(points[0],points[1]);
    SegType seg1(points[2],points[3]);
    Intersector lsx1(seg0,seg1);
    std::cout << " Test 1:  Parallel\n";
    printTest(seg0,seg1,lsx1);
    std::cout << "\n\n";
    
    //Test Case 2: Collinear w/Overlap
    SegType seg2(points[5],points[6]);
    Intersector lsx2(seg0,seg2);
    std::cout << " Test 2:  Collinear w/Overlap\n";
    printTest(seg0,seg2,lsx2);
    std::cout << "\n\n";
    
    //Test Case 3: None
    SegType seg3(points[2],points[8]);
    Intersector lsx3(seg0,seg3);
    std::cout << " Test 3:  NONE\n";
    printTest(seg0,seg3,lsx3);
    std::cout << "\n\n";
    
    //Test Case 4: Collinear (disjoint)
    SegType seg4(points[7],points[8]);
    Intersector lsx4(seg0,seg4);
    std::cout << " Test 4:  Collinear Disjoint\n";
    printTest(seg0,seg4,lsx4);
    std::cout << "      Bounds overlap = " << lsx4.isIntersectionPossible() << "\n";
    std::cout << "      Apply compute():\n";
    lsx4.compute();
    printTest(seg0,seg4,lsx4);
    std::cout << "\n\n";
    
    //Test Case 5: Intersection
    SegType seg5(points[2],points[4]);
    Intersector lsx5(seg0,seg5);
    std::cout << " Test 5:  Intersection\n";
    printTest(seg0,seg5,lsx5);
    std::cout << "\n\n";
    
    //Test Case 6: Purely vertical segment that should have an intersection
    SegType vertSeg(nvis::vec2(0.0,4.0),nvis::vec2(0.0,0.0) );
    SegType seg6(points[0],points[3]);
    Intersector lsx6(seg6,vertSeg);
    std::cout << " Test 6:  Intersection with Vertical Segment\n";
    printTest(seg6,vertSeg,lsx6);
    std::cout << "      Bounds overlap = " << lsx6.isIntersectionPossible() << "\n";
    std::cout << "\n\n";
    
    //Test Case 7: Purely horizontal segment that should have an intersection
    SegType horzSeg(nvis::vec2(-1.0,0.0),nvis::vec2(2.0,0.0) );
    Intersector lsx7(seg6,horzSeg);
    std::cout << " Test 7:  Intersection with Horizontal Segment\n";
    printTest(seg6,horzSeg,lsx7);
    std::cout << "      Bounds overlap = " << lsx7.isIntersectionPossible() << "\n";
    std::cout << "\n\n";
    
    
    return 0;
}
