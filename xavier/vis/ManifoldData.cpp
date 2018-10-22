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


/*
 *  ManifoldData.cpp - a Data class
 *
 *  Author: Wayne Schlei
 *          Purdue University
 *
 *  Date:  03/26/2015
 *
 */

//Avizo Headers
#include <QApplication>
#include <hxcore/HxMessage.h>
#include <poincare/FixedPointData.h>
#include <poincare/ManifoldData.h>

//Additional Headers
#include <ctime>
#include <QString>
#include <QByteArray>

using namespace std;

//Required Macro for Hx Classes
HX_INIT_CLASS(ManifoldData,HxData)


///Constructor
ManifoldData::ManifoldData() :
    fpxDataConnection(this,"FPData",QApplication::translate("ManifoldData","FPData"),FixedPointData::getClassTypeId()),
    portInfo(this,"Info",QApplication::translate("ManifoldData","Info")),
    portSeps(this,"SepXs",QApplication::translate("ManifoldData","SepXs")),
    portFPDataFile(this,"FpxFile",QApplication::translate("ManifoldData","FpxFile"))
{
    //Port Initializing
    portInfo.setValue("Invariant Manifold Storage");
    portSeps.setValue("---");
    portFPDataFile.registerFileType("FixedPointData","fpx");

}

ManifoldData::~ManifoldData()
{ }

///Avizo compute()
void ManifoldData::compute()
{
  int numSep = getNumManifolds();
  QString numStr;
  numStr.setNum(numSep);
  portSeps.setValue(numStr.toLocal8Bit().data());
  
  //Continuously update the icon position to below the FixedPointData icon?
}

