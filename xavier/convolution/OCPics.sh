
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

#Bash script for executing convolution process
# on connected data sets
# Wayne Schlei & Xavier Tricoche (Purdue University)

#Intial options
niter=100
#Hamiltonian Values (Jacobi)
#-------------------------------
#Earth-Moon Set
#hLow=1.5
hLow=1.8
hHi=3.2
#Earth-Moon Set 2
hLow=2.91
hHi=3.1932
hSteps=59
#Saturn-Titan Set
#hLow=1.55
#hHi=3.15
#Saturn-Enceladus Set
#hLow=1.55
#hHi=3.15
#hSteps=50

#Res
rx=1024
ry=512

#mup=1 #Implement?
numPasses=4

#Global variables for name mangling
#baseName=emMapData.1024x512.x.xx.nrrd
#ocName=emMap.x.xx.nrrd
baseName=emMapData2.1024x512.x.xx.nrrd
ocName=emMap2.x.xx.nrrd

#baseName=stMapData.1024x512.x.xx.nrrd
#ocName=stMap.x.xx.nrrd
#baseName=senMapData.1024x512.x.xx.nrrd
#ocName=senMap.x.xx.nrrd
#L3 Side
#baseName=emMapData.L3Box.1024x512.x.xx.nrrd
#ocName=emMap.L3Box.x.xx.nrrd


cStrBase=x.xx
nStrBase=nrrd
pStrBase=png
wpStrBase=wo.png
#Noise image
noiseName=noiseEM.nrrd

#Steps
dh=$(echo "${hHi} - ${hLow}" | bc)
dh=$(echo "scale=6; $dh / $hSteps " | bc)

#Loop through each value of Jacobi
for i in $(seq 0 $hSteps )
do
	cVal=$(echo "scale=6;${hLow} + $i*$dh" | bc)
	#Construct file names
	dataFile=$(echo ${baseName} | sed -e "s/$cStrBase/$cVal/" )
	imageFile=$(echo ${ocName} | sed -e "s/$cStrBase/$cVal/" )

	#Call convolution function
	pc=$(echo "scale=6; 100 * $i / ${hSteps} " | bc)
	echo "Call $i / ${hSteps} = ${pc} %:"
	sh OCProcess.sh $dataFile $imageFile $noiseName $numPasses

	#Create a White-out image with ImageMagick
	imagePNG=$(echo ${imageFile} | sed -e "s/$nStrBase/$pStrBase/")
	wImagePNG=$(echo ${imageFile} | sed -e "s/$nStrBase/$wpStrBase/")
	echo "Generating light version..."
	convert $imagePNG -fill white -colorize 60% $wImagePNG

done

#Storage
mkdir -p ocImages
mkdir -p ocData
mv emMapData.*.*.*.nrrd ./ocData/
mv *.png ./ocImages/
