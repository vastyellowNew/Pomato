
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

#Bash script for computing orbit convolutino data for apse sections
#Wayne Schlei

#Hard-coded dual call
fn="./bin/computeMapOC_Apse"
options=" -e 1.e-8 -b 0.8 -0.2 1.2 0.2 -r 512 512 -C 3.172 -P 2 -p 50 "
p0=" -s 0 -d 0 -o apseMapData.PeriPro.r512.p50.3.172.nrrd"
p1=" -s 0 -d 1 -o apseMapData.PeriRetro.r512.p50.3.172.nrrd"
a0=" -s 1 -d 0 -o apseMapData.ApoPro.r512.p50.3.172.nrrd"
a1=" -s 1 -d 1 -o apseMapData.ApoRetro.r512.p50.3.172.nrrd"

echo "$fn $options $p0"
$fn $options $p0

echo "$fn $options $p1"
$fn $options $p1

echo "$fn $options $a0"
$fn $options $a0

echo "$fn $options $a1"
$fn $options $a1
