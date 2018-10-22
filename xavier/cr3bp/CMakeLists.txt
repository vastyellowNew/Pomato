#*************************************************************************
# POMATO: POincare MAp TOpology: Extraction and visualization of the
#         topological structure of the circular restricted three-body
#         problem for orbital mechanics applications.
#
# Authors: Wayne Schlei and Xavier Tricoche
#
# Copyright (c) 2013-2018, Purdue University
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#**************************************************************************
set(EXTRA_SRC
	${XAVIER_DIR}/maps/map_analysis
)
add_library(topo ${EXTRA_SRC})

set(TEST_SRC
	cr3bp_plot
	computeMapOC
	computeMapOC_Apse
	computeMapOC_yL4
	computeMapOC_thetaL4
	pmateMergeFPData
	pmateTabulateFPX
	pmateAdaptiveGrid
	pmateTopologyExtraction
	pmateManifoldCreator
	test_MultipleShooting
	test_fpICgrid
	plot_winding2
	plot_winding_Avizo
	test_Monodromy
	test_SaddleMan
	test_AstroManifold
	test_TaskScheduling
	testIntersection
	writeOCLatexFigs
)
#Others - need to be updated!
#cr3bp_analysis
#plot_winding
#test_wnAdaptiveGrid

foreach(target ${TEST_SRC})
	add_executable(${target} ${target}.cpp)
	target_link_libraries(${target} topo  
        ${GLUT_LIBRARIES}
        ${OPENGL_LIBRARIES} 
        ${GLUI_LIBRARIES} 
        ${Teem_LIBRARIES} )
endforeach(target ${TEST_SRC})
