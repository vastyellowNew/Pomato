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
#!/bin/sh

#Bash script for PMATE Results in CR3BP
# Wayne Schlei & Xavier Tricoche (Purdue University)
#
# Note: This is intended for large-scale runs on a large domain.
#       While using large domains, PMATE may experience some memory
#       issues that are temporarily circumvented by calling
#       PMATE for only half the domain at a time.  Remove relevant
#       split parts if you wish to run some smaller cases.

#####################################
#                                   #
#  Basic path & directory settings  #
#                                   #
#####################################
home=${HOME}
code="${home}/code/github"
pmate="${code}/pmate"
binDir="${pmate}/build/bin"
dataDir="${home}/visdata/Orbital_Mechanics/CR3BP/Jupiter-Europa"
paramDir="${pmate}/script/params"
fn="${binDir}/pmateTopologyExtraction"
fpxMerge="${binDir}/pmateMergeFPData"
fpxTable="${binDir}/pmateTabulateFPX"
manCompute="${binDir}/pmateManifoldCreator"
outputPath=${dataDir}
mapParam="${paramDir}/map_params_JE.txt"
manParam="${paramDir}/jeTopoMan.C3.0034.mset"

#############################
#                           #
#  Intial Sampling options  #
#                           #
#############################
eps=1.e-8
niter=50
depth=3
#PMATE Parameters
#lmin=2e-5 => This value will take a long time: ~4.5hrs for a half
#lmVal="2m5"
lmin=1e-3
lmVal="1m3"
#pmax=12 => This is a high value
pmax=10
#Res
rx=24
ry=80
maxa=2.356194490192345

#################################
#                               #
#  Hamiltonian Values (Jacobi)  #
#                               #
#################################


#Bounds - Main Box
bName="MainBox"
mup=2.528017682687079e-05
sysVal="Jupiter-Europa"
xmin=1.002300453692939
xmax=1.02
ymin=-0.135109616283352
ymax=0.135109616283352

######################
#                    #
#  Manifold options  #
#                    #
######################

# NB: These options are defined in a file to be imported by the manifold constuction algorithm with -is <params>

# min sampling distance between consecutive manifold points
# tau_min=2e-5 # "any lower and you will be waiting forever for the manifold to compute"
#tau_min=5e-5

# Maximum stability index
# Many periodic orbits have huge stability index magnitudes >1e5 which can jam up the computation process.

# The second parameter that I often enabled was the maximum segment count (“-xg 1 –ns 1200”)

#Global variables for name mangling
baseName=${outputPath}/sysPMATE.bbox.Cx.xxPART
sysStrBase=sys
cStrBase=x.xx
bboxStrBase=bbox
partStrBase=PART
outputName=${outputPath}/sysPMATE_bbox_Cx.xx.fpx
latexName=${outputPath}/sysPMATE_bbox_Cx.xx.tex

cVal=${1}
fpxFiles=""
#Break runs into 2 halves that split x-bounds

for j in $(seq 1 2 ); do
    #Evaluate the half-size x component
    xRange=$(echo "${xmax} - ${xmin}" | bc)
    dx=$(echo "scale=6;${xRange}/2.0" | bc)
    xMinPart=$(echo "scale=6;${xmin} + ($j-1)*${dx}" | bc)
    xMaxPart=$(echo "scale=6;${xMinPart} + ${dx}" | bc)
    rxHalf=$(echo "${rx}/2" | bc)
    partName=$(echo ".part$j")

    #Construct file names
    filename=$(echo ${baseName} | sed -e "s/$sysStrBase/$sysVal/" )
    filename=$(echo ${filename} | sed -e "s/$cStrBase/$cVal/" )
    filename=$(echo ${filename} | sed -e "s/$bboxStrBase/$bName/" )
    filename=$(echo ${filename} | sed -e "s/$partStrBase/$partName/" )
    fpxFiles=$(echo "${fpxFiles} ${filename}.fpx")

    #Call mapping function - Note that this splits
    # pc=$(echo "scale=6; 100 * $i / ${hSteps} " | bc)
    echo "******************************************************************************************"
    echo "PMATE: Call Part $j:"
    echo "******************************************************************************************"
    echo "$fn -e ${eps} -b ${xMinPart} ${ymin} ${xMaxPart} ${ymax} -r ${rxHalf} ${ry} \ "
    echo "    -m ${mup} -C ${cVal} -md ${depth} -xp ${pmax} -p ${niter} -o ${filename} \ "
    echo "    -wp 1 -wg 1 -of 1 -og 1 -lm ${lmin} -ma ${maxa}"

    if [ -n $2 ]; then
        echo "$fn -e ${eps} -b ${xMinPart} ${ymin} ${xMaxPart} ${ymax} -r ${rxHalf} ${ry} \ " >> ${2}
        echo "    -m ${mup} -C ${cVal} -md ${depth} -xp ${pmax} -p ${niter} -o ${filename} \ " >> ${2}
        echo "    -wp 1 -wg 1 -of 1 -og 1 -lm ${lmin} -ma ${maxa}" >> ${2}
    fi

    #With mup input
    time -va -o timestats.txt $fn -e ${eps} -b ${xmin} ${ymin} ${xmax} ${ymax} -r ${rxHalf} ${ry} \
        -m ${mup} -C ${cVal} -md ${depth} -xp ${pmax} \
        -p ${niter} -wp 1 -wg 1 -of 1 -og 1 -lm ${lmin} -ma ${maxa} -o ${filename}
done

#Add any other localized runs on smaller regions

outputFile=$(echo ${outputName} | sed -e "s/$sysStrBase/$sysVal/" )
outputFile=$(echo ${outputFile} | sed -e "s/$bboxStrBase/$bName/" )
cModVal=$(echo ${cVal} | sed -e "s/\./_/" ) #CHECK THIS!
outputFile=$(echo ${outputFile} | sed -e "s/$cStrBase/$cModVal/" )

latexFile=$(echo ${latexName} | sed -e "s/$sysStrBase/$sysVal/" )
latexFile=$(echo ${latexFile} | sed -e "s/$bboxStrBase/$bName/" )
latexFile=$(echo ${latexFile} | sed -e "s/$cStrBase/$cModVal/" )

#Reconnect the two resulting fixed point sets into one file

echo "******************************************************************************************"
echo "Merge fixed points files"
echo "******************************************************************************************"

echo "$fpxMerge -i ${fpxFiles} -o ${outputFile}"
$fpxMerge -i ${fpxFiles} -o ${outputFile}

#Tabulate Results
echo "******************************************************************************************"
echo "Tabulate results"
echo "******************************************************************************************"
echo "$fpxTable -i ${outputFile} -o ${latexFile}"
$fpxTable -i ${outputFile} -o ${latexFile}

if [ -n $2 ]; then
    echo "$fpxMerge -i ${fpxFiles} -o ${outputFile}" >> ${2}
    echo "$fpxTable -i ${outputFile} -o ${latexFile}" >> ${2}
fi

#Run the manifold generator
manoutputFile=$(echo ${outputFile} | sed -e "s/\.fpx/_manifolds/" )
echo "******************************************************************************************"
echo " MANIFOLD Creator:"
echo "******************************************************************************************"
echo "$manCompute -C ${cVal} -m ${mup} -b ${xmin} ${ymin} ${xmax} ${ymax} -ip ${mapParam} -ix ${outputFile} -is ${manParam} -o ${manoutputFile}"

if [ -n $2 ]; then
    echo "$manCompute -C ${cVal} -m ${mup} -b ${xmin} ${ymin} ${xmax} ${ymax} -ip ${mapParam} -ix ${outputFile} -is ${manParam} -o ${manoutputFile}" >> $2
    echo "Note:" >> $2
    echo "cat ${mapParam}" >> $2
    cat ${mapParam} >> $2
    echo "cat ${manParam}" >> $2
    cat ${manParam} >> $2
fi

time -va -o timestats.txt $manCompute -C ${cVal} -m ${mup} -b ${xmin} ${ymin} ${xmax} ${ymax} -ip ${mapParam} -ix ${outputFile} -is ${manParam} -o ${manoutputFile}
