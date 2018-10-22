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


/// Generate a Latex figures for OC images
/// Author:  Wayne Schlei
/// Date:    1/15/2017


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

std::string me("");
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
    Convert << std::fixed << std::setprecision(6) << C; // Use some manipulators
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

//Write start of a new table
void writeFigureStart(FILE* f,const int tableCount)
{
    fprintf(f,"%% Figure number %d\n\n",tableCount);
    fprintf(f,"%%\\begin{sidewaysfigure}\n");
    fprintf(f,"\\begin{figure}[hp]\n");
    fprintf(f,"  \\centering\n");
}

//Write the end of a table (C is input based on range and steps)
void writeFigureEnd(FILE* f, std::string& domainStr, std::string& jcRangeStr,
                    const int sys, const int tableNumber)
{
    std::string sysStr = getSystemString(sys);
    fprintf(f,"  \\caption{Depiction of \\Poincare sections in ($x$,$\\dot{x}$) cooridnates over %s \n",domainStr.c_str() );  //domainStr = "the domain $D_{EM,A}$"
    fprintf(f,"             on $\\varSigma : y=0$ in the %s system with orbit convolution over %s.}\n",sysStr.c_str(),jcRangeStr.c_str() ); //jcRangeStr = "$C_{EM,A}$"
    
    
    //Table reference
    std::stringstream  os("");
    os << "   \\label{fig:ocSET" << sysStr << "_" << tableNumber << "}";
    std::string labelStr = os.str();
    fprintf(f,"%s\n",labelStr.c_str());
    fprintf(f,"\\end{figure}\n");
    fprintf(f,"%%\\end{sidewaysfigure}\n");
    fprintf(f,"\\pagebreak\n\n");
}

//Write a combo entry (L3 + Normal domain)
void writeComboEntry(
    FILE* f,
    std::string& folder1,
    std::string& folder2,
    std::string& filePrefix1,
    std::string& filePrefix2,
    const double& C)
{
    //With Underscore
    std::string jc_Str = getJCString(C);
    //Normal value entry
    std::stringstream jss("");
    jss << std::setprecision(6) << C;
    std::string jcStr = jss.str();
    
    std::stringstream pic1SS("");
    pic1SS << folder1.c_str() << "/" << filePrefix1.c_str() << jc_Str << ".png";
    std::stringstream pic2SS("");
    pic2SS << folder2.c_str() << "/" << filePrefix2.c_str() << jc_Str << ".png";
    
    //Subfigure Start
    fprintf(f,"  \\subfigure[$C = %s$]{\n",jcStr.c_str());
    //First figure (L3 side)
    fprintf(f,"     \\tikz[baseline=(a.north)]\\node[yscale=-1,inner sep=0,outer sep=0](a){\n"); //Flip image
    fprintf(f,"      \\includegraphics[width=0.49\\textwidth]{%s}};\n",pic1SS.str().c_str());
    //Second figure (main box)
    fprintf(f,"     \\tikz[baseline=(a.north)]\\node[yscale=-1,inner sep=0,outer sep=0](a){\n"); //Flip image
    fprintf(f,"      \\includegraphics[width=0.49\\textwidth]{%s}};}\n",pic2SS.str().c_str());
    
}

//Write a single image entry
void writeEntry(
    FILE* f,
    std::string& folder,
    std::string& filePrefix,
    const double& C)
{
    //With Underscore
    std::string jc_Str = getJCString(C);
    //Normal value entry
    std::stringstream jss("");
    jss << std::setprecision(6) << C;
    std::string jcStr = jss.str();
    
    std::stringstream picSS("");
    picSS << folder.c_str() << "/" << filePrefix.c_str() << jc_Str << ".png";
    
    //Subfigure Start
    fprintf(f,"  \\subfigure[$C = %s$]{\n",jcStr.c_str());
    fprintf(f,"     \\tikz[baseline=(a.north)]\\node[yscale=-1,inner sep=0,outer sep=0](a){\n"); //Flip image
    fprintf(f,"      \\includegraphics[width=0.46\\textwidth]{%s}};}\n",picSS.str().c_str());
    
}


void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Output figures of OC images for Latex\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -f  | --filePrefix <string>      *Filename prefix\n"
            << " -d  | --dir <string>             *Directory of image files\n"
            << " -f2 | --filePrefix2 <string>     Second Filename prefix\n"
            << " -d2 | --dir2 <string>            Second Directory of image files\n"
            << " -D  | --domainStr <string>       Domain string for figure caption\n"
            << " -cs | --jcRangeStr <string>      Jacobi range string for figure caption\n"
            << " -C0 | --jc0 <double>             Initial Jacobi value in range\n"
            << " -Cf | --jcf <double>             Final Jacobi value in range\n"
            << " -descend                         Print Jacobi in descending order\n"
            << " -nC | --numJC <int>              Number of Jacobi constant steps [C0,Cf]\n"
            << " -n  | --picsPerFigure <int>      Num images in one page [default = 8]\n"
            << " -s  | --system <int>             CR3BP system:\n"
            << "                                  0:Earth-Moon   | 1:Jupiter-Europa  |\n"
            << "                                  2:Saturn-Titan | 3:Saturn-Enceladus \n"
            << " -o  | --output <string>          Output LaTeX File\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    me = argv[0];
    int tableCounter = 0;
    int maxPicsPerPage = 8;
    int currentPageReturnCount = 0;
    int sysID = 0; //Default EM
    
    //Jacobi range values
    double jc0 = 2.91, jcf = 3.1932;
    int numC = 50;
    bool printDescending = false;
    
    //Default strings : EM Range A
    std::string domainStr("the domain $D_{EM,A}$");
    std::string jcRangeStr("$C_{EM,A}$");
    std::string folder("none");
    std::string filePrefix("none");
    bool writeCombo = false;
    std::string folder2("none");
    std::string filePrefix2("none");
    
    
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("");
        } else if (arg == "-f" || arg == "--filePrefix") {
            if (i == argc-1) {
                printUsageAndExit("missing file prefix");
            }
            filePrefix = argv[++i];
        } else if (arg == "-d" || arg == "--dir") {
            if (i == argc-1) {
                printUsageAndExit("missing directory string");
            }
            folder = argv[++i];
        } else if (arg == "-f2" || arg == "--filePrefix2") {
            if (i == argc-1) {
                printUsageAndExit("missing second file prefix");
            }
            filePrefix2 = argv[++i];
            writeCombo = true;
        } else if (arg == "-d2" || arg == "--dir2") {
            if (i == argc-1) {
                printUsageAndExit("missing second directory string");
            }
            folder2 = argv[++i];
            writeCombo = true;
        } else if (arg == "-D" || arg == "--domainStr") {
            if (i == argc-1) {
                printUsageAndExit("missing domain string");
            }
            domainStr = argv[++i];
        } else if (arg == "-cs" || arg == "--jcRangeStr") {
            if (i == argc-1) {
                printUsageAndExit("missing jc range string");
            }
            jcRangeStr = argv[++i];
        } else if (arg == "-C0" || arg == "--jc0") {
            if (i == argc-1) {
                printUsageAndExit("missing jc 0  value");
            }
            jc0 = atof(argv[++i]);
        } else if (arg == "-Cf" || arg == "--jcf") {
            if (i == argc-1) {
                printUsageAndExit("missing jc final  value");
            }
            jcf = atof(argv[++i]);
        } else if (arg == "-nC" || arg == "--numJC") {
            if (i == argc-1) {
                printUsageAndExit("missing num Jacobi values");
            }
            numC = atoi(argv[++i]);
        } else if (arg == "-descend") {
            printDescending = true;
        } else if (arg == "-n" || arg == "--picsPerFigure") {
            if (i == argc-1) {
                printUsageAndExit("missing num pics per figure");
            }
            maxPicsPerPage = atoi(argv[++i]);
        } else if (arg == "-s" || arg == "--system") {
            if (i == argc-1) {
                printUsageAndExit("missing system integer");
            }
            sysID = atoi(argv[++i]);
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit("missing output name");
            }
            texFile = argv[++i];
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if ( texFile=="none" || folder=="none" || filePrefix=="none" ) {
        printUsageAndExit("");
    }
    
    
    
    //Open a new file for writing data to LaTex figures
    FILE* f = fopen(texFile.c_str(), "w");
    
    int tableCount = 1, dataCounter = 0;
    
    //Start the first Figure
    writeFigureStart(f,tableCount);
    
    //Loop through all C values in range to write to file
    double deltaC = (jcf - jc0) / ((double) numC-1);
    double C = jc0;
    for (int i=0; i<numC; i++) {
        //Increment to new Jacobi value
        if (printDescending) {
            C = jcf - ((double) i) * deltaC;
        } else {
            C = jc0 + ((double) i) * deltaC;
        }
        
        //Note:  There is a maximum number of pics on a single page
        
        //Check if over, and reset if so
        if ( (dataCounter >= maxPicsPerPage) ) {
            //Close out the current figure
            writeFigureEnd(f,domainStr,jcRangeStr,sysID,tableCount);
            //Start a new figure
            tableCount++;
            writeFigureStart(f,tableCount);
            dataCounter = 0;
        }
        
        //Write Figure at C value
        if (writeCombo) {
            dataCounter += 2; //Add two images
            writeComboEntry(f,folder,folder2,filePrefix,filePrefix2,C);
        } else {
            writeEntry(f,folder,filePrefix,C);
            //Add some space between figs (if even counter)
            if(dataCounter == 0 || (dataCounter % 2) == 0) {
                fprintf(f,"   \\hspace{2mm}\n");
            } else {
                fprintf(f,"   \n");
            }
            dataCounter++;  //Add a single image
        }
        
    } //End for each C
    
    //Close out the current table when finished with data
    writeFigureEnd(f,domainStr,jcRangeStr,sysID,tableCount);
    
    
    //Close file
    fclose(f);
    
    return 0;
}
