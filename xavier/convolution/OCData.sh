#!/bin/sh

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

#Bash script for computing Orbit convolution data
# Wayne Schlei & Xavier Tricoche (Purdue University)

#Intial options
fn="./bin/computeMapOC"
eps=1.e-8
niter=100
#Hamiltonian Values (Jacobi)
#Earth-Moon Set
#hLow=1.8
#hHi=3.2
#Saturn-Titan Set
#hLow=1.55
#hHi=3.15
#Saturn-Enceladus Set
hLow=1.55
hHi=3.05
hSteps=50
#Res
rx=1024
ry=512
#Bounds - Main Box
#xmin=-0.4
#ymin=-2.5
#xmax=1.1
#ymax=2.5
#Bounds - L3Side Box
#xmin=-1.9
#xmax=-0.4
#ymin=-2.5
#ymax=2.5
#Bounds - L3Side Box2
#xmin=-3.4
#xmax=-1.9
#ymin=-2.5
#ymax=2.5
#Bounds for Saturn-Enceladus
#MainBox
#xmin=0.0
#xmax=1.3
#ymin=-2.5
#ymax=2.5
#L3SideBox
xmin=-2.3
xmax=-0.8
ymin=-2.5
ymax=2.5

#EM System - Default
# ->Be sure to remove 'mu' option in function call for EM
#ST System
#mup=0.00023658080508871
#SEn System
mup=1.8984152807945e-07

#Global variables for name mangling
#baseName=emMapData.1024x512.x.xx.nrrd
#baseName=stMapData.1024x512.x.xx.nrrd
#baseName=senMapData.1024x512.x.xx.nrrd
#baseName=emMapData.L3Box.1024x512.x.xx.nrrd
#baseName=emMapData.L3Box2.1024x512.x.xx.nrrd
#baseName=stMapData.L3Box.1024x512.x.xx.nrrd
baseName=senMapData.L3Box.1024x512.x.xx.nrrd
cStrBase=x.xx


#Steps
dh=$(echo "${hHi} - ${hLow}" | bc)
dh=$(echo "scale=6; $dh / $hSteps " | bc)

#Loop through each value of Jacobi
for i in $(seq 0 $hSteps )
do
	cVal=$(echo "scale=6;${hLow} + $i*$dh" | bc)
	#Construct file names
	filename=$(echo ${baseName} | sed -e "s/$cStrBase/$cVal/" )
	#Call mapping function
	pc=$(echo "scale=6; 100 * $i / ${hSteps} " | bc)
	echo "Call $i / ${hSteps} = ${pc} %:"
	echo "$fn -e ${eps} -b ${xmin} ${ymin} ${xmax} ${ymax} -r ${rx} ${ry} "

	#With mup input
	echo "    -m ${mup} -C ${cVal} -p ${niter} -o ${filename}"
	$fn -e ${eps} -b ${xmin} ${ymin} ${xmax} ${ymax} -r ${rx} ${ry} \
	    -m ${mup} \
	    -C ${cVal} -p ${niter} -o ${filename}

	#Without mup input (default to Earth-Moon system)
	#echo "    -C ${cVal} -p ${niter} -o ${filename}"
	#$fn -e ${eps} -b ${xmin} ${ymin} ${xmax} ${ymax} -r ${rx} ${ry} \
	#    -C ${cVal} -p ${niter} -o ${filename}



	#Storage

done
