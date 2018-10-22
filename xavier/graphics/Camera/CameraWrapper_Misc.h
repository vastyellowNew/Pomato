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


#ifndef CAMERA_WRAPPER_MISC_H
#define CAMERA_WRAPPER_MISC_H

#include "CameraWrapper.h"

namespace CameraWrapper_GUI {
	void glutMouseCallback(CameraWrapper *camera, int button, int state, int x, int y) {
		if (state == GLUT_DOWN) {
			if (button == GLUT_LEFT_BUTTON) camera->recordLeftMouseDown(x, y);
			if (button == GLUT_RIGHT_BUTTON) camera->recordRightMouseDown(x, y);
			if (button == GLUT_MIDDLE_BUTTON) camera->recordMiddleMouseDown(x, y);
		}

		else if (state == GLUT_UP) {
			if (button == GLUT_LEFT_BUTTON) camera->recordLeftMouseUp(x, y);
			if (button == GLUT_RIGHT_BUTTON) camera->recordRightMouseUp(x, y);
			if (button == GLUT_MIDDLE_BUTTON) camera->recordMiddleMouseUp(x, y);
		}
	}

	void glutMouseMotionCallback(CameraWrapper *camera, int x, int y) {
		camera->recordMouseMotion(x, y);
	}

	void glutResizeCallback(CameraWrapper *camera, int width, int height) {
		camera->resizeViewport(width, height);
	}
}

#endif