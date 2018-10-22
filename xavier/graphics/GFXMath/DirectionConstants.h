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


#ifndef DIRECTION_CONSTANTS
#define DIRECTION_CONSTANTS

enum AxisDirections { DIRECTION_X = 0, DIRECTION_Y, DIRECTION_Z, NUM_AXIS_DIRECTIONS };

enum FaceDirections { FACE_RIGHT = 1, FACE_LEFT = 2, FACE_UP = 4, FACE_DOWN = 8, FACE_FRONT = 16, FACE_BACK = 32 };

enum CornerDirections {
	CORNER_LBB = 0, CORNER_LBF, CORNER_LTB, CORNER_LTF,
	CORNER_RBB, CORNER_RBF, CORNER_RTB, CORNER_RTF, NUM_CORNER_DIRECTIONS
};

enum PosNegDirection { DIRECTION_NEGATIVE = 0, DIRECTION_POSITIVE, NUM_POS_NEG_DIRECTION };

#endif