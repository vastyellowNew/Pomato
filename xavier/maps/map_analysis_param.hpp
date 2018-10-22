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


#ifndef __XAVIER_MAP_ANALYSIS_PARAMS_HPP__
#define __XAVIER_MAP_ANALYSIS_PARAMS_HPP__


#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <boost/rational.hpp>
#include <boost/regex.hpp>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include "metric.hpp"

namespace xavier {

typedef metric<double, 2>           metric_type;
typedef grid<double, 2>             grid_type;
typedef grid_type::bounds_type      bounds_type;
typedef grid_type::vec_type         vec_type;
typedef grid_type::ivec_type        ivec_type;
typedef boost::rational<int>        rational_type;

struct map_analysis_param {
    typedef std::vector<nvis::vec2>           orbit_type;
    typedef std::pair<orbit_type, int>        tagged_orbit_type;
    typedef std::pair<nvis::vec2, int>        tagged_vector_type;
    
    map_analysis_param()
        : upsampling_factor(1), max_depth(1), min_period(1), max_period(12),
          samplingIntegTol(1.e-8), refinementIntegTol(1.e-12), refinementConvTol(1.e-8),
          refinementMaxIters(20), refinementPointsPerPeriod(5),
          max_angle(0.75*M_PI), maxDeltaTransSpeed(0.5), maxDeltaMapDisplacement(0.2,2.0),
          maxTimeChange(2.0), maxDeltaTimeFactor(0.15), min_fp_tol(1.e-4),
          nb_iterations(50), record(false), verbose(false), debug(false),
          linearMonodromy(true), lmin(3.e-4), ntEdgeDivisions(64), subcellResolution(6,6)
    {}
    
    nvis::ivec2       resolution;
    nvis::bbox2       bounds;
    metric_type       the_metric;
    int               upsampling_factor;
    int               max_depth;
    int               min_period;
    int               max_period;
    double            samplingIntegTol;
    double            refinementIntegTol;
    double            refinementConvTol;
    int               refinementMaxIters;
    int               refinementPointsPerPeriod; //For Multiple Shooting
    double            max_angle;
    double            maxDeltaTransSpeed;
    vec_type          maxDeltaMapDisplacement;
    double            maxTimeChange;
    double            maxDeltaTimeFactor;
    double            min_fp_tol;
    int               nb_iterations;
    std::vector<double>        winding_convexity_tols;
    std::vector<double>        winding_cell_maxDist;
    std::vector<nvis::vec2>                    edges;
    std::vector<tagged_orbit_type>             orbits;
    std::vector<std::pair<nvis::vec2, tagged_vector_type> >    vectors;
    bool              record;
    bool              verbose;
    bool              debug;
    bool              linearMonodromy;
    double            lmin;
    int               ntEdgeDivisions;
    nvis::ivec2       subcellResolution;
    std::map<double, rational_type>             valid_rationals;
    
    
    /// Read from file : NOTE, fixed format read!
    bool read(const char* filename)
    {
        FILE* f = fopen(filename, "r");
        if(!f) {
            std::cerr << " Bad Read for map_analysis_param...\n";
            return false;
        }
        //Ideally, we use boost::regex here, but I don't have the time...
        
        char buffer[200];
        int tempBool = 0;
        //Headers
        for(int i=0; i<4; i++) {
            fgets(buffer, 200, f);
        }
        //Parameters in order
        fscanf( f, "%s %d %d",buffer,&(resolution[0]),&(resolution[1]));
        fscanf( f, "%s %lf %lf %lf %lf",buffer,
                &(bounds.min()[0]),&(bounds.min()[1]),
                &(bounds.max()[0]),&(bounds.max()[1])       );
        fscanf( f, "%s %lf %lf %lf %lf",buffer,
                &(the_metric.bounds().min()[0]),&(the_metric.bounds().min()[1]),
                &(the_metric.bounds().max()[0]),&(the_metric.bounds().max()[1])  );
        fscanf( f, "%s %d",buffer,&tempBool);
        the_metric.periodic()[0] = (tempBool==1);
        fscanf( f, "%d",&tempBool);
        the_metric.periodic()[1] = (tempBool==1);
        fscanf( f, "%s %d",buffer,&min_period);
        fscanf( f, "%s %d",buffer,&max_period);
        fscanf( f, "%s %d",buffer,&nb_iterations);
        fscanf( f, "%s %lf",buffer,&samplingIntegTol);
        fscanf( f, "%s %lf",buffer,&refinementIntegTol);
        fscanf( f, "%s %lf",buffer,&refinementConvTol);
        fscanf( f, "%s %d",buffer,&refinementMaxIters);
        fscanf( f, "%s %d",buffer,&refinementPointsPerPeriod);
        fscanf( f, "%s %d",buffer,&tempBool);
        record = (tempBool==1);
        fscanf( f, "%s %d",buffer,&tempBool);
        verbose = (tempBool==1);
        fscanf( f, "%s %d",buffer,&tempBool);
        debug = (tempBool==1);
        fscanf( f, "%s %d",buffer,&tempBool);
        linearMonodromy = (tempBool==1);
        fscanf( f, "%s %lf",buffer,&lmin);
        fscanf( f, "%s %d",buffer,&ntEdgeDivisions);
        fscanf( f, "%s %d %d",buffer,&(subcellResolution[0]),&(subcellResolution[1]) );
        
        //Heuristics in order
        for(int i=0; i<3; i++) {
            fscanf(f, "%s", buffer);    //Read header lines
        }
        fscanf( f, "%s %d",buffer,&max_depth);
        fscanf( f, "%s %lf",buffer,&max_angle);
        fscanf( f, "%s %lf",buffer,&maxDeltaTransSpeed);
        fscanf( f, "%s %lf %lf",buffer,
                &(maxDeltaMapDisplacement[0]),&(maxDeltaMapDisplacement[1]));
        fscanf( f, "%s %lf",buffer,&maxTimeChange);
        fscanf( f, "%s %lf",buffer,&maxDeltaTimeFactor);
        fscanf( f, "%s %d",buffer,&tempBool);
        fscanf( f, "%s",buffer);
        winding_convexity_tols.clear();
        for(int i=0; i<tempBool; i++) {
            double value;
            fscanf( f, "%lf", &value);
            winding_convexity_tols.push_back(value);
        }
        fscanf( f,"%s",buffer);
        winding_cell_maxDist.clear();
        for(int i=0; i<tempBool; i++) {
            double value;
            fscanf( f, "%lf", &value);
            winding_cell_maxDist.push_back(value);
        }
        
        
        fclose(f);
        return true;
    }
    
    /// Write to file
    bool write(const char* filename) const
    {
        FILE* f = fopen(filename, "w"); //Open file
        if (!f) {
            std::cerr << " Bad filename for map_analysis_param...\n";
            return false;
        }
        
        //Header
        fprintf(f,"# Map Analysis Parameters File\n");
        
        //Domain
        fprintf(f,"-----------------------------------------\n");
        fprintf(f," PARAMETERS:\n");
        fprintf(f,"-----------------------------------------\n");
        fprintf(f,"RESOLUTION= %d %d\n",resolution[0],resolution[1]);
        fprintf(f,"BOUNDS= %.15f %.15f %.15f %.15f\n",
                bounds.min()[0],bounds.min()[1],
                bounds.max()[0],bounds.max()[1]);
        fprintf(f,"METRIC.BOUNDS= %.15f %.15f %.15f %.15f\n",
                the_metric.bounds().min()[0],the_metric.bounds().min()[1],
                the_metric.bounds().max()[0],the_metric.bounds().max()[1]);
        fprintf(f,"METRIC.PERIODIC= %d %d\n",
                (the_metric.periodic()[0])? 1 : 0,
                (the_metric.periodic()[1])? 1 : 0  );
        fprintf(f,"MIN_PERIOD= %d\n",min_period);
        fprintf(f,"MAX_PERIOD= %d\n",max_period);
        fprintf(f,"NUM_ITERS= %d\n",nb_iterations);
        fprintf(f,"SAMPLE_INTEG_TOL= %.15f\n",samplingIntegTol);
        fprintf(f,"REFINE_INTEG_TOL= %.15f\n",refinementIntegTol);
        fprintf(f,"REFINE_CONV_TOL= %.15f\n",refinementConvTol);
        fprintf(f,"REFINE_MAX_ITERS= %d\n",refinementMaxIters);
        fprintf(f,"REFINE_PTS_PER_P= %d\n",refinementPointsPerPeriod);
        fprintf(f,"RECORD= %d\n",(record)? 1 : 0);
        fprintf(f,"VERBOSE= %d\n",(verbose)? 1 : 0);
        fprintf(f,"DEBUG= %d\n",(debug)? 1 : 0);
        fprintf(f,"LINEAR_MONODROMY= %d\n",(linearMonodromy)? 1 : 0);
        fprintf(f,"LMIN= %.15f\n",lmin);
        fprintf(f,"NT_EDGE_DIV= %d\n",ntEdgeDivisions);
        fprintf(f,"SUBCELL_RES= %d %d\n",subcellResolution[0],subcellResolution[1]);
        
        
        fprintf(f,"\n\n");
        fprintf(f,"-----------------------------------------\n");
        fprintf(f," HEURISTICS:\n");
        fprintf(f,"-----------------------------------------\n");
        fprintf(f,"MAX_DEPTH= %d\n",max_depth);
        fprintf(f,"MAX_ANGLE= %.15f\n",max_angle);
        fprintf(f,"MAX_DELTA_TRANSSPEED= %.15f\n",maxDeltaTransSpeed);
        fprintf(f,"MAX_DELTA_DISPLACEMENT= %.15f %.15f\n",
                maxDeltaMapDisplacement[0],maxDeltaMapDisplacement[1]);
        fprintf(f,"MAX_DELTA_TIME= %.15f\n",maxTimeChange);
        fprintf(f,"MAX_DELTA_TIMEFACTOR= %.15f\n",maxDeltaTimeFactor);
        fprintf(f,"WINDING_SET_SIZE= %d\n",(int)winding_cell_maxDist.size());
        fprintf(f,"WINDING_CONVEXITY_TOLS= ");
        for(int i=0; i<(int) winding_convexity_tols.size(); i++) {
            fprintf(f,"%.15f ",winding_convexity_tols[i]);
        }
        fprintf(f,"\n");
        fprintf(f,"WINDING_MAX_DISTANCE= ");
        for(int i=0; i<(int) winding_cell_maxDist.size(); i++) {
            fprintf(f,"%.15f ",winding_cell_maxDist[i]);
        }
        fprintf(f,"\n");
        
        //Skipped:  valid_rationals, edges, orbits, vectors
        
        fclose(f);
        return true;
    }
};




} //end xavier




#endif  // __XAVIER_MAP_ANALYSIS_PARAMS_HPP__
