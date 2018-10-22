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


#ifndef __VTK_UTILS_HPP__
#define __VTK_UTILS_HPP__

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

namespace nvis {
typedef bounding_box<vec1> bbox1;
}

#include <vis/vtk_macros.hpp>
#include <vis/vtk_data_helper.hpp>
#include <vis/vtk_camera_helper.hpp>
#include <vis/vtk_image_helper.hpp>

#endif