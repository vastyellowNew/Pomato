
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
# Wayne Schlei & Xavier Tricoche (Purdue University)

#Example call:
#  sh OCProcess.sh mapData.nrrd output.nrrd noiseImage.nrrd numPasses


#Convlution function
# Note: 's' option is method to indicate shortfall
#  0 : Drop completely (if not enough hits)
#  1 : Use a gray Penalty
#  2 : Use a gray Penalty but skew based on completed returns
#  3 : Accept whatever is available for full color
convolution () {
	#Call the convolution function
	echo "singleConvolutionPass -d $1 -i $2 -k $3 -np $4 -s 0 -o $5"
	./bin/singleConvolutionPass -d $1 -i $2 -k $3 -np $4 -s 0 -o $5
	return 0
}

#Inputs
mapData=${1}
image=${2}
noise=${3}
N=${4}
#maxp=${5}

#Copy noise into image
cp ${noise} ${image}

#Find image name and create png handle
imagePNG=$(echo ${image} | sed -e "s/.nrrd/.png/")

#for ((a=1; a <= N ; a++)) #bash
for a in $(seq 0 $(($N-1)))
do
	convolution ${mapData} ${image} $a ${N} "tmp.nrrd"
	if [ $a -lt $(($N-1)) ]; then
	 #High-pass Filter with TEEM
	 #unu resample -i tmp.nrrd -k gauss:2,2 -s = x1 x1 | \
	 #unu 2op - tmp.nrrd - | unu 2op x 0.5 - -o ${image}

	 #Others convert images (png)
	 unu quantize -b 8 -i tmp.nrrd -o tmp.png

	 #Imagemagick: nrrd->png->hpfilter->png->nrrd
	 #Best thus far
	 #convert tmp.png -bias 50% -morphology Convolve LoG:0x2 tmp_hpf.png
	 #Others:
	 #convert tmp.png -define convolve:scale='100,100%' \
	 #-morphology Convolve 'LoG:0x2' tmp_hpf.png
	 #convert tmp.png -bias 50% -morphology Convolve DoG:0,0,2 tmp_hpf.png
	 #Laplacian filter
	 #convert tmp.png -define convolve:scale='!' -bias 50% \
	 #-morphology Convolve Laplacian:0 tmp_hpf.png
	 #convert tmp.png -define convolve:scale='!' -bias 50% \
	 #-morphology Convolve Laplacian:15 tmp_hpf.png

	 #Gimp
	 cp tmp.png tmp_hpf.png
	 gimp -i -b "(simple-unsharp-mask \"tmp_hpf.png\" 1 5 0)" -b "(gimp-quit 0)"

	 #Done
	 unu unquantize -i tmp_hpf.png -o ${image}
	 rm -f tmp.png
	 rm -f tmp_hpf.png
	 #cp tmp.nrrd ${image}
	else
	 cp tmp.nrrd ${image}
	fi
done
rm -f tmp.nrrd

#Create picture
unu quantize -b 8 -i ${image} -o ${imagePNG}
