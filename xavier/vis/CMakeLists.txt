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
set( VIS_SRCS
    display_topology
)

if( ${VTK_FOUND} )
    include( ${VTK_USE_FILE} )
    foreach( name ${VIS_SRCS} )
        add_executable( ${name} ${EXECUTABLE_FLAG} ${name}.cpp )
        target_link_libraries( ${name} ${VTK_LIBRARIES} )
        install( TARGETS ${name} RUNTIME DESTINATION bin )
    endforeach()
else()
    message( WARNING "Vis executables: " ${VIS_SRCS} " cannot be built."
        "\nVTK was not found."
    )
endif()
