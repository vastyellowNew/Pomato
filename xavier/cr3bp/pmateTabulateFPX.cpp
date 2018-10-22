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


/// Generate a Latex table for FixedPointData
/// Author:  Wayne Schlei
/// Date:    7/23/2016


#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <regex>
#include <iterator>
// math
#include <util/wall_timer.hpp>
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>

//API - maps
#include <maps/DP45wrapper.hpp>
#include <maps/definitions.hpp>
#include <maps/map_analysis.hpp>
#include <maps/poincare_map.hpp>
//API - cr3bp
#include <cr3bp/cr3bp.hpp>
#include <cr3bp/planar_section.hpp>
#include <cr3bp/tracker.hpp>
#include <cr3bp/multipleAngleTracker.hpp>
#include <maps/fixpoints.hpp>
//API - pmate
#include <pmate/FixedPointData.hpp>
#include <pmate/FixedPointDataFilters.hpp>


using namespace nvis;
typedef pmate::FixedPointData<xavier::fixpoint>  FPData;
#if defined(_WIN32) || defined(C_0X)
typedef FPData::FixedPointChain         FixedPointChain;
typedef FPData::FPChainsContainer       FPChainsContainer;
typedef FixedPointChain::iterator       FPChainIter;
typedef FPChainsContainer::iterator     FPChainsContainerIter;
#else
typedef typename FPData::FixedPointChain         FixedPointChain;
typedef typename FPData::FPChainsContainer       FPChainsContainer;
typedef typename FixedPointChain::iterator       FPChainIter;
typedef typename FPChainsContainer::iterator     FPChainsContainerIter;
#endif
typedef xavier::dp5wrapper<double, 42>                                      ode_solver;
typedef orbital::cr3bp                                                      rhs_type;
typedef orbital::planar_section<rhs_type, 6, 42>                            section_type;
typedef orbital::MultipleAngleTracker<42, 3>                                tracker_type;
typedef xavier::poincare_map<rhs_type, ode_solver, section_type >           MapType;
typedef xavier::map_analysis_param                                          MapParams;

std::string me("");
std::string fpxFile("none");
std::string texFile("none");

//Define a custom numeric facet
class WithU: public std::numpunct<char> { // class for decimal numbers using _ instead of point
protected:
    char do_decimal_point() const
    {
        return '_';    // change the decimal separator
    }
};

//Extract the Jacobi Constant to be C='2_96000' (Could use regex)
std::string getJCString(const double& C)
{
    char buffer[10];
    int n = sprintf(buffer,"%.5f",C);
    std::string temp(buffer);
    /* Regex is being mean!!!
    std::string temp("2.96000");
    std::cout << "Jacobi string printf: " << temp << "\n";
    std::regex reg1("\\.");
    std::string replacement("_");
    //Check for match
    if( std::regex_search(temp,reg1) ) {
      std::cout << " Regex pattern found\n";
    } else {
      std::cout << " Regex NO MATCH\n";
    }
    std::string result =  std::regex_replace(temp,reg1,replacement);
    std::cout << "After regex_replace: " << result << " and temp = " << temp << "\n";
    return result;*/
    
    std::ostringstream Convert;
    std::locale MyLocale(  std::locale(), new WithU);// Crate customized locale
    Convert.imbue(MyLocale);       // Imbue the custom locale to the stringstream
    Convert << std::fixed << std::setprecision(5) << C; // Use some manipulators
    std::string result =  Convert.str(); // Give the result to the string
    //std::cout << " We want " << temp << "\n";
    //std::cout << " Actual : " << result << "\n";
    return result;
    
}

std::string getSystemString(const int sys)
{
    std::string systemStr;
    switch(sys) {
        case 0 :
            systemStr = "EM";
            break;
        case 1 :
            systemStr = "JE";
            break;
        case 2 :
            systemStr = "ST";
            break;
        case 3 :
            systemStr = "SEnc";
            break;
    }
    return systemStr;
}

//Write the start of FIRST table
void writeFirstTableStart(FILE* f,const int tableCount)
{
    fprintf(f,"%% Table number %d\n\n",tableCount);
    fprintf(f,"\\begin{table}[ht]\n");
    fprintf(f,"  \\scalebox{0.88} {\n");
    fprintf(f,"  \\centering\n");
    fprintf(f,"  \\begin{tabular}{|r|c| d{2} d{8} d{8}| c|c|}\n");
    fprintf(f,"    \\hline\n");
    fprintf(f,"    ID & Parameters & \n");
    fprintf(f,"    \\multicolumn{1}{c}{$t_{k}$} & \\multicolumn{1}{c}{$x_k$} & \\multicolumn{1}{c}{$\\dot{x}_k$}  & Orbit (Map) & Orbit $(xy)$  \\\\ \n");
    fprintf(f,"    \\hline \\hline\n");
}

//Write start of a new table
void writeTableStart(FILE* f,const int tableCount)
{
    fprintf(f,"%% Table number %d\n\n",tableCount);
    fprintf(f,"\\begin{table}[ht]\n");
    fprintf(f,"  \\scalebox{0.8} {\n");
    fprintf(f,"  \\centering\n");
    fprintf(f,"  \\begin{tabular}{|r|c| d{2} d{8} d{8}| c|c|}\n");
    fprintf(f,"    \\hline\n");
    fprintf(f,"    ID & Parameters & \n");
    fprintf(f,"    \\multicolumn{1}{c}{$t_{k}$} & \\multicolumn{1}{c}{$x_k$} & \\multicolumn{1}{c}{$\\dot{x}_k$}  & Orbit (Map) & Orbit $(xy)$  \\\\ \n");
    fprintf(f,"    \\hline \\hline\n");
}

//Write the end of a table (C is from fpx file)
void writeTableEnd(FILE* f, const double& mup, const double& C, const int sys, const int tableNumber)
{
    fprintf(f,"   \\end{tabular}\n");
    fprintf(f,"   }\n");
    fprintf(f,"   \\caption{Fixed points extracted from $\\varSigma:y=0$ in the Earth-Moon system\n");
    fprintf(f,"             ($\\mu = %.10e $) at $C = %.5f $ ($T$ and $t_k$ in days). } \n",mup,C);
    
    
    //Table reference
    std::string jcStr = getJCString(C);
    std::string sysStr = getSystemString(sys);
    std::stringstream  os("");
    os << "   \\label{tab:" << sysStr << "C" << jcStr << "_" << tableNumber << "}";
    std::string labelStr = os.str();
    fprintf(f,"%s\n",labelStr.c_str());
    fprintf(f,"\\end{table}\n\n");
    fprintf(f,"\\pagebreak\n\n");
}

//Write a blank return entry
void writeBlankEntry(FILE* f)
{
    fprintf(f,"     & &  &  &  & & \\\\ \n");
}

//Write an iterate
void writeEntry(FILE* f, const int i, const FixedPointChain& fpChain, const double& tstar)
{
    double tValue = fpChain[i].t * tstar; //In days
    double xdot = fpChain[i].pos[1];
    if (fabs(xdot) <= 1e-8) {
        xdot = fabs(xdot);
    }
    fprintf(f,"     & & %.2f & %.8f & %.8f & & \\\\ \n", tValue, fpChain[i].pos[0], xdot );
}
//Evaluate tstar in days
double eval_tstar(const double& gm1, const double& gm2, const double& lstar)
{
    return sqrt(lstar*lstar*lstar/(gm1+gm2))/(3600.0*24.0);// days
}

//Write fixed point entry given a chain
void writeFixedPointChainToTable(
    FILE* f,
    const int& sys,
    const double& C,
    const int& orbitID,
    const FixedPointChain& fpChain,
    const double& tstar
)
{
    int p = fpChain[0].K;
    const xavier::fixpoint& fp0 = fpChain[0];
    double periodDays = fp0.timePeriod*tstar;
    fprintf(f,"    %% A p = %d orbit\n",p);
    
    char buffer[5];
    int j=0;
    if (p>=10) {
        j=sprintf(buffer,"%2d%s",p,(fp0.saddle)? "S" : "C");
    } else {
        j=sprintf(buffer,"%1d%s",p,(fp0.saddle)? "S" : "C");
    }
    std::string pTypeString(buffer);
    
    int rowVal = std::max(p,4);
    
    std::string systemStr = getSystemString(sys);
    
    std::string jcStr = getJCString(C);
    std::stringstream orbitPicSS("");
    //orbitPicSS << "OrbitID" << orbitID << "_" << systemStr << "C" << jcStr << ".png";
    orbitPicSS << "OrbitID" << orbitID << "_" << systemStr << "C" << jcStr << ".jpg";
    std::stringstream mapPicSS("");
    //mapPicSS << "MapID" << orbitID << "_" << systemStr << "C" << jcStr << ".png";
    mapPicSS << "MapID" << orbitID << "_" << systemStr << "C" << jcStr << ".jpg";
    
    
    //Always writing the initial Table entries : OrbitID, pType, etc... as an aligned array of equations
    fprintf(f,"    \\multirow{%2d}{*}{%4d} & \\multirow{%2d}{*}{ \n",
            rowVal,orbitID,rowVal);
    fprintf(f,"        {$\\begin{aligned}\n");
    fprintf(f,"           p &= %s \\\\ \n",pTypeString.c_str());
    fprintf(f,"           T &= %.2f \\\\ \n",periodDays);
    fprintf(f,"           \\nu_{SI} &= %.2e \\\\ \n",fp0.si);
    fprintf(f,"           \\nu_{z,SI} &= %.2f \n",fp0.otherSI);
    fprintf(f,"         \\end{aligned}$ } } & \n");
    
    //TimePeriod, Stability Index, Out-of-plane Stability index
    //if (fp0.si < 99.99) { //Switch between standard and scientific
    //fprintf(f,"       \\multirow{%2d}{*}{%.2f} & \\multirow{%2d}{*}{%.2f} & \\multirow{%2d}{*}{%.2f} & \n",
    //rowVal,periodDays,rowVal,fp0.si,rowVal,fp0.otherSI);
    //} else {
    //  fprintf(f,"       \\multirow{%2d}{*}{%.2f} & \\multirow{%2d}{*}{%.2e} & \\multirow{%2d}{*}{%.2f} & \n",
    //    rowVal,periodDays,rowVal,fp0.si,rowVal,fp0.otherSI);
    //}
    
    //Write the first point (t[days],x,xdot) and Images
    double xdot = fp0.pos[1];
    if (fabs(xdot) <= 1e-8) {
        xdot = fabs(xdot);
    }
    fprintf(f,"         %.2f &  %.8f &  %.8f & \n", fp0.t*tstar, fp0.pos[0],xdot);
    //Map Pic
    fprintf(f,"         \\multirow{%2d}{*}{ \\includegraphics[scale=0.09]{fpAppendix/orbits/%s} } &  \n",
            rowVal,mapPicSS.str().c_str() );
    //rowVal,"MakeaFigure.pdf"); //mapPicSS.str().c_str());
    //Orbit Pic
    fprintf(f,"         \\multirow{%2d}{*}{ \\includegraphics[scale=0.075]{fpAppendix/orbits/%s} } \\\\  \n",
            rowVal,orbitPicSS.str().c_str() );
    //rowVal,"MakeaFigure.pdf"); //orbitPicSS.str().c_str());
    
    //Special print for less than 4 returns to section
    if (p==1) {
        //Write 3 blank entries to create 4 rows
        for(int i=0; i<3; i++) {
            writeBlankEntry(f);
        }
    } else if(p==2) {
        //Write the second return
        writeEntry(f,1,fpChain,tstar);
        //Then write two blanks to complete 4 rows
        writeBlankEntry(f);
        writeBlankEntry(f);
    } else if(p==3) {
        //Write the second and third returns
        writeEntry(f,1,fpChain,tstar);
        writeEntry(f,2,fpChain,tstar);
        //Write a blank to complete 4 rows
        writeBlankEntry(f);
    } else {
        for(int k=1; k<p; k++) {
            writeEntry(f,k,fpChain,tstar);
        }
    }
    
    //Wrap up whole orbit with a horizontal line
    fprintf(f,"    \\hline\n");
}


void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Tabulate fixed point data\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -i  | --input <string>           Input .fpx File\n"
            << " -n  | --returnsPerTable <int>    Returns allowed in one page [default = 42]\n"
            << " -s  | --system <int>             CR3BP system:\n"
            << "                                  0:Earth-Moon   | 1:Jupiter-Europa  |\n"
            << "                                  2:Saturn-Titan | 3:Saturn-Enceladus \n"
            << " -d  | --dir  <int>               Direction (>=0:Positive [d],<0:Negative)\n"
            << " -e  | --eps <double>             Tolerance on 'same point' test [d=5e-6]\n"
            << " -o  | --output <string>          Output LaTeX File\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    me = argv[0];
    int tableCounter = 0;
    int maxReturnsInTable = 51;
    int maxReturnsInFIRSTTable = 20;
    int currentPageReturnCount = 0;
    bool positiveDir = true;
    double eps = 5.e-6;
    int sysID = 0;
    
    //System information for conversion to real units:
    //Standard GMs
    double gmE=3.986004418e5, gmS = 1.327122e11,
           gmJ = 1.26686535e8,     gmSat = 37940585.; //km^3/s^2
    //Standard Radii
    double rE = 6378.14, rM = 1738.20, rS = 695990.00,
           rJ = 71492.0, rSat = 60268.00; //km - (mean)
    //Moons
    double gmM=4902.7949, gmEu = 3202.739, gmTitan = 8978.1382, gmEnceladus = 7.2027;  //km^3/s^2
    // gmRhea = 153.9426, gmHyperion = 0.3727;
    
    //Earth-Moon System
    double gm1 = gmE, gm2 = gmM;
    double lstar = 384388.174; //km
    double radius1 = rE; //km
    double radius2 = rM; //km
    
    //Generic values
    double tstar = eval_tstar(gm1,gm2,lstar);
    double mup = gm2/(gm1+gm2);
    
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("");
        } else if (arg == "-i" || arg == "--input") {
            if (i == argc-1) {
                printUsageAndExit("missing input name");
            }
            fpxFile = argv[++i];
        } else if (arg == "-n" || arg == "--returnsPerTable") {
            if (i == argc-1) {
                printUsageAndExit("missing returns per table");
            }
            maxReturnsInTable = atoi(argv[++i]);
        } else if (arg == "-s" || arg == "--system") {
            if (i == argc-1) {
                printUsageAndExit("missing system integer");
            }
            sysID = atoi(argv[++i]);
            switch(sysID) {
                case 0 :
                    //Earth-Moon
                    gm1 = gmE;
                    gm2 = gmM;
                    lstar = 384388.174; //km
                    radius1 = rE;
                    radius2 = rM;
                    break;
                case 1 :
                    //Jupiter-Europa
                    gm1 = gmJ;
                    gm2 = gmEu;
                    lstar = 670900.0; // km
                    radius1 = rJ; //km - equatorial
                    radius2 = 1560.8;//km
                    break;
                case 2 :
                    //Saturn-Titan
                    gm1 = gmSat;
                    gm2 = gmTitan;
                    lstar = 1221865.; //km
                    radius1 = rSat; //km - equatorial
                    radius2 = 2574.73;
                    break;
                case 3 :
                    //Saturn-Enceladus
                    gm1 = gmSat;
                    gm2 = gmEnceladus;
                    lstar = 238042.; //km
                    radius1 = rSat; //km - equatorial
                    radius2 = 252.1;
                    break;
                    
                default :
                    printUsageAndExit("Invalid System Value");
                    break;
            }
            //Evaluate common parameters
            mup = gm2 / (gm1+gm2);
            tstar = eval_tstar(gm1,gm2,lstar);
        } else if (arg == "-d" || arg == "--dir") {
            if (i == argc-1) {
                printUsageAndExit("missing positive section direction");
            }
            int d = atoi(argv[++i]);
            positiveDir = (d>=0);
        } else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) {
                printUsageAndExit("missing tolerance");
            }
            eps = (double) atof(argv[++i]);
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit("missing output name");
            }
            texFile = argv[++i];
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if ( texFile=="none" || fpxFile=="none" ) {
        printUsageAndExit("");
    }
    
    //Open a FixedPointData object from a file
    FPData* theData = new FPData(fpxFile.c_str());
    
    //Setup the Map object
    double C = theData->getJacobiConstant();
    //Build Poincare Map object
    const double LARGE = std::numeric_limits<double>::max();
    //Right-hand Side (CR3BP EOMs)
    rhs_type rhs(C, mup);
    //Section
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    section.isPositive = positiveDir;
    
    //The Poincare Map Engine
    MapType theMap(rhs, section);
    theMap.setPrecision(1.e-12);
    
    //Setup the Map Parameters
    MapParams   theMapParams; //Assigns defaults
    theMapParams.max_depth = 1;
    std::vector<double> wTols, wDist;
    for(int i=0; i<3; i++) {
        wTols.push_back( 0.5 );
        wDist.push_back( 1.0 );
    }
    wTols[1] = 1.0;
    wDist[1] = 2000.0; //w_xydot
    theMapParams.winding_convexity_tols = wTols;
    theMapParams.winding_cell_maxDist = wDist;
    theMapParams.the_metric.bounds().min()[0] = -LARGE;
    theMapParams.the_metric.bounds().min()[1] = -2.55;
    theMapParams.the_metric.bounds().max()[0] = LARGE;
    theMapParams.the_metric.bounds().max()[1] = 2.55;
    //Probably the most important part for map params:
    theMapParams.the_metric.periodic()[0] = false; //Not periodic bounds
    theMapParams.the_metric.periodic()[1] = false; //Not periodic bounds
    
    //Run a filter on the fixed points to make sure there is no copies of orbits in data
    pmate::filterFixedPointData( (*theData), theMap, theMapParams, eps );
    
    
    //Open a new file for writing data to LaTex tables
    FILE* f = fopen(texFile.c_str(), "w");
    
    
    
    //Gather each chain within this object:
    FPChainsContainer data;
    theData->getData( data );
    FPChainsContainerIter fpcIT;
    int tableCount = 1, dataCounter = 0;
    
    //Start the first table
    writeFirstTableStart(f,tableCount);
    
    //Loop through all fixed point chains to write to file
    for (int i=0; i<(int)data.size(); i++) {
    
        int p = data[i][0].K;
        dataCounter += std::max(p,4); //Always adding 4 lines
        
        //Note:  There is a maximum number of returns that can be collected
        //       on one page (i.e., one table)
        
        //The first table has a smaller number (say 24 ~= 6 p=1,2,3,4 orbits)
        if ( (tableCount == 1 && dataCounter >= maxReturnsInFIRSTTable) ||
                (dataCounter >= maxReturnsInTable) ) {
            //Close out the table
            writeTableEnd(f,mup,C,sysID,tableCount);
            //Start a new table
            tableCount++;
            writeTableStart(f,tableCount);
            dataCounter = std::max(p,4);
        }
        
        //Write orbit to table
        writeFixedPointChainToTable(f,sysID,C,i,data[i],tstar);
        
    } //End for each chain
    
    //Close out the current table when finished with data
    writeTableEnd(f,mup,C,sysID,tableCount);
    
    //May need Close-Approach distance in the future (for P1 and P2)
    
    
    
    //Close file
    fclose(f);
    //Free info
    delete theData;
    
    return 0;
}
