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


#ifndef CUDA_CONSTANTS_CUH
#define CUDA_CONSTANTS_CUH

#define COMPILE_WITH_CUDA

#define CUDA_MSS_COMPILE_AS_OCTREE

enum CudaMatrixIDs { CUDA_MATRIX_0 = 0, CUDA_MATRIX_1, CUDA_MATRIX_2 };

#define NODE_MATRIX_ID CUDA_MATRIX_0
#define CORNER_MATRIX_ID CUDA_MATRIX_1
#define NODE_EXISTS_MATRIX_ID CUDA_MATRIX_2


#define MAX_TEXTURE_2D_WIDTH_BYTES					65536
#define MAX_TEXTURE_2D_WIDTH_4_BYTE_VALUE			MAX_TEXTURE_2D_WIDTH_BYTES/4
#define MAX_TEXTURE_2D_WIDTH_8_BYTE_VALUE			MAX_TEXTURE_2D_WIDTH_BYTES/8
#define MAX_TEXTURE_2D_WIDTH_16_BYTE_VALUE			MAX_TEXTURE_2D_WIDTH_BYTES/16

#define MAX_TEXTURE_2D_WIDTH_BYTES_LG				16
#define MAX_TEXTURE_2D_WIDTH_4_BYTE_VALUE_LG		14
#define MAX_TEXTURE_2D_WIDTH_8_BYTE_VALUE_LG		13
#define MAX_TEXTURE_2D_WIDTH_16_BYTE_VALUE_LG		12

#define MAX_TEXTURE_2D_WIDTH_BYTES_MINUS_1			MAX_TEXTURE_2D_WIDTH_BYTES-1
#define MAX_TEXTURE_2D_WIDTH_4_BYTE_VALUE_MINUS_1	MAX_TEXTURE_2D_WIDTH_4_BYTE_VALUE-1
#define MAX_TEXTURE_2D_WIDTH_8_BYTE_VALUE_MINUS_1	MAX_TEXTURE_2D_WIDTH_8_BYTE_VALUE-1
#define MAX_TEXTURE_2D_WIDTH_16_BYTE_VALUE_MINUS_1	MAX_TEXTURE_2D_WIDTH_16_BYTE_VALUE-1





#ifndef DIRECTION_CONSTANTS
#define DIRECTION_CONSTANTS

	enum AxisDirections { DIRECTION_X = 0, DIRECTION_Y, DIRECTION_Z };

	enum FaceDirections { FACE_RIGHT = 0, FACE_LEFT, FACE_UP, FACE_DOWN, FACE_FRONT, FACE_BACK, NUM_FACE_DIRECTIONS };

	enum CornerDirections {
		CORNER_LBB = 0, CORNER_LBF, CORNER_LTB, CORNER_LTF,
		CORNER_RBB, CORNER_RBF, CORNER_RTB, CORNER_RTF, NUM_CORNER_DIRECTIONS
	};

	enum PosNegDirection { DIRECTION_NEGATIVE = 0, DIRECTION_POSITIVE, NUM_POS_NEG_DIRECTION };

#endif


#ifndef M_PI
	#define M_PI 3.14159265358979
#endif

#ifndef M_2_PI
	#define M_2_PI 6.28318530718
#endif

#ifndef M_PI_2
	#define M_PI_2 1.57079632679
#endif

#ifndef M_PI_4
	#define M_PI_4 0.785398163397
#endif

#ifndef M_PI_8
	#define M_PI_8 0.392699081699
#endif

#ifndef M_e
	#define M_e 2.71828182846
#endif


#endif