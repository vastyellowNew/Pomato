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


// EdgeRotation : Computes the rotation along a set of edges
// Author:  Wayne Schlei

#ifndef __EDGE_ROTATION_HPP__
#define __EDGE_ROTATION_HPP__

#include <iostream>
#include <iterator>
#include <vector>
#include <set>
#include <algorithm>
#include <memory>
#ifdef HX_HAS_STD
#include <boost/shared_ptr.hpp>
#endif
//xavierAPI
#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>
#include <data/edge.hpp>
#include <maps/definitions.hpp>
#include <maps/mapExceptions.hpp>
#include <maps/misc.hpp>
#include <maps/index.hpp>
#include <topology/SortableReturnData.hpp>
#include <topology/EdgeRotationFailure.hpp>
#include <topology/SectionTransversality.hpp>
#include <topology/EdgeRotationMapCalls.hpp>
#include <topology/AdaptiveEdgeRotation.hpp>


//OpenMP
#if _OPENMP
#include <omp.h>
#endif
//Avizo
#ifdef HX_HAS_STD
#include <QString>
#include <hxcore/HxMessage.h>
#include <hxcore/internal/HxWorkArea.h>
#endif

using namespace xavier;

namespace topology {

/// Macro for calling mapDisplacementUsingCache() in Edge_Rotation: throwing and storing EdgeRotationFailures
#define ERFMAPDISP(vOut,mapDispCall) \
    try { \
        erFail.reset(); \
        vOut = mapDispCall; \
    } catch (ERotFailure& eRotFail) { \
        switch(eRotFail.type) { \
            case ERotFailure::BACKWARD_MAP : \
                useSectionSeparationVersion = true; \
                break; \
            case ERotFailure::FORWARD_MAP : \
                useSectionSeparationVersion = true; \
                break; \
            case ERotFailure::BACKWARD_SECTION_SEP : \
                useSectionSeparationVersion = true; \
                break; \
            case ERotFailure::FORWARD_SECTION_SEP : \
                useSectionSeparationVersion = true; \
                break; \
            case ERotFailure::BACKWARD_SECTION_SEP_NODE : \
                useSectionSeparationVersion = true; \
                break; \
            case ERotFailure::FORWARD_SECTION_SEP_NODE : \
                useSectionSeparationVersion = true; \
                break; \
            case ERotFailure::BACKWARD_SINGULARITY : \
                useSectionSeparationVersion = true; \
                break; \
            case ERotFailure::FORWARD_SINGULARITY : \
                useSectionSeparationVersion = true; \
                break; \
            case ERotFailure::DOUBLE_PERIOD_FIXED_POINT : \
        {if (2*p <= theMapParams.max_period) { \
                    FP aGuess; \
                    aGuess.pos = eRotFail.failurePos; aGuess.K = p; \
                    cached_fps[thread_id].push_back( aGuess ); \
                } \
                if(verbose) std::cout << "FPSuspected : " << eRotFail.what() << " at " << eRotFail.where() << "\n"; \
                useSectionSeparationVersion = true; \
            }break; \
            case ERotFailure::FIXED_POINT_SUSPECTED : \
            {FP aGuess; \
                aGuess.pos = eRotFail.failurePos; aGuess.K = p; \
                cached_fps[thread_id].push_back( aGuess ); \
                if(verbose) std::cout << "FPSuspected : " << eRotFail.what() << " at " << eRotFail.where() << "\n"; \
                useSectionSeparationVersion = true; \
            }break; \
            default : \
                std::cout << "Unknown ERotFailure is thrown at Starting Point Calls.\n"; \
                std::cout << " ->ERotFailure : " << eRotFail.what() << " at " << eRotFail.where() << "\n"; \
                useSectionSeparationVersion = true; \
                break; \
        } \
    }
//Note on FATAL failure:  It's possible that a map state could evaluate as "INVALID"
//with "cr3bp: imaginary velocity" at this point due to numerical precision errors.
//                std::cout << " ->Moving to next period indicating rotation failure...\n";
//              e.rotFailures.push_back( eRotFail );
//              errorOccurred = true;
//              continue;

/// Edge Rotation function : Computes the rotation of the map-tangent vector along edges stored in a set [SEDGE => std::set< edge_type >]
template <class DATASET, class MAP, class PARAM, class SEDGE, class FP>
void Edge_Rotation(const DATASET& mapData, const MAP& theMap, const PARAM& theMapParams,
                   SEDGE& edges, std::vector< FP >& fpGuesses,
                   int nb_sub=4)
{
    typedef typename SEDGE::value_type                edge_type;
    typedef typename SEDGE::value_type::index_type    ivec_type;
    typedef typename MAP::lvec_type                   vec_type;
    typedef typename MAP::rhs_type                    RHStype;
    typedef typename DATASET::data_type               OrbitData;
    typedef EdgeRotationFailure<vec_type>             ERotFailure;
    
    int counter = 0;
    int nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    /** Note on lmin:
    * Minimum distance between consecutive sample points along
    * an edge to track the continuous rotation of a vector
    */
    double lmin = theMapParams.lmin;
    //int ntEdgeDivisions = theMapParams.ntEdgeDivisions;
    bool verbose = theMapParams.verbose;
    std::vector<edge_type>*  cached_edges = new std::vector<edge_type>[nbthreads];
    std::vector< FP >*  cached_fps = new std::vector< FP >[nbthreads];
    nvis::timer _timer;
    typename SEDGE::iterator edgeIt;
    static const int s = RHStype::numSingularities;
    typedef nvis::fixed_vector<double,s+1>                   ExtendedMapDataVec;
    typedef SortableReturnData<vec_type,ExtendedMapDataVec>  SortableData;
    
    //Copy to a std::vector
    std::vector<edge_type> copyOfEdges(edges.begin(),edges.end());
    size_t N = copyOfEdges.size();
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n=0; n<(int)N; n++) {
            int thread_id = 0;
            std::ostringstream os;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            //Create an edge copy
            edge_type e( copyOfEdges[n] );
            ivec_type i0 = e[0];
            ivec_type i1 = e[1];
            vec_type start0 = mapData.getVertex(i0);
            vec_type start1 = mapData.getVertex(i1);
            
            
            //Warn if we have no relevant map data (i.e., invalid initial sampling)
            if (!mapData(i0).hasSteps() || !mapData(i1).hasSteps() ) {
                std::cerr << " Edge_Rotation: WARNING!! Edge [" << i0 << "," << i1 << "] has no returns for corner node :\n";
                if (!mapData(i0).hasSteps()) {
                    std::cerr << "       Node 0 : " << i0 << " at " << start0 << "\n";
                }
                if (!mapData(i0).isICValid) {
                    std::cerr << "       Node 0 : is INFEASIBLE!\n";
                }
                if (!mapData(i1).hasSteps()) {
                    std::cerr << "       Node 1 : " << i1 << " at " << start1 << "\n";
                }
                if (!mapData(i1).isICValid) {
                    std::cerr << "       Node 1 : is INFEASIBLE!\n";
                }
            }
            //Do not store an edge outside feasible region!
            // -> Make sure to run removeInvalidLeaves() on Adaptive Grid BEFORE setting up cells.
            if (!mapData(i0).isICValid || !mapData(i1).isICValid) {
                continue;    //if(non-valid IC) skip this edge
            }
            
            //Compute the rotation amount occurring on an edge given a period:
            // Note:  We start with largest period first and store data to save on rerunning the map.
            std::set<int>::reverse_iterator rit;
            rit = e.periods.rbegin();
            //We want to store data in hopes that we can just look it up instead of recomputing
            std::set<SortableData> deltaStorageSet;  //Forward steps
            std::set<SortableData> betaStorageSet;   //Backward steps
            //Add the two end nodes for p => PROPAGATE, DON'T USE STORAGE
            //deltaStorageSet.insert( SortableData(start0,mapData,i0,*rit) );
            //deltaStorageSet.insert( SortableData(start1,mapData,i1,*rit) );
            
            //Debug:  Report what's being stored
            if (theMapParams.debug) {
                std::cout << "****************************************************\n";
                std::cout << " Edge : " << n << "\n";
                std::cout << "****************************************************\n";
                std::cout << " Original Delta states:\n";
                std::cout << "  p =" << (*rit) << "\n";
                std::cout << "      x0  = " << start0 << "     x1  = " << start1 << "\n";
                vec_type px0 = mapData(i0).steps[(*rit)-1];
                vec_type px1 = mapData(i1).steps[(*rit)-1];
                std::cout << "  P^p(x0) = " << px0 << " P^p(x1) = " << px1 << "\n";
                std::cout << "  delta0  = " << px0-start0 << " delta1  = " << px1 - start1 << "\n";
            }
            
            //Local copy per thread
#ifdef HX_HAS_STD
            boost::shared_ptr<MAP> pmap( theMap.clone() );
#else
            std::shared_ptr<MAP> pmap( theMap.clone() );
#endif
            
            //Step Direction for using interpolation to fix map errors
            vec_type edgeStep = start1 - start0;
            double edgeLength = nvis::norm(edgeStep);
            edgeStep /= edgeLength;
            
            
            //On Map Functions:
            //mapDisplacement(x,pmap,p,theMapParams);
            //   -> runs map and returns the vec_type v = displacement of pth iterate
            //map_complete(x,returns,p);
            //   -> returns = return_map_info using passive tracking
            //mapDisplacementUsingCache(x,pmap,p,theMapParams,forwardStorageSet,backwardStorageSet);
            //   -> Calls from storage set if available
            //   -> If this can't find it in storage, will run map and store new data
            
            
            //For each period, starting with largest, compute rotation
            for (; rit!=e.periods.rend(); ++rit) {
                int p = *rit;
                //Assume by default that transversality assumption is valid
                e.transverseSection.insert( std::pair<int,bool>(p,true) );
                
                //Flag for using separated section version of adaptiveRotationAngle()
                bool useSectionSeparationVersion = false;
                //Flag for errors
                bool errorOccurred = false;
                
                //Evaluate Tangents (not just Displacements) on corners
                vec_type delta0, delta1;
                ERotFailure erFail;
                //Mapping with built-in error catching
                // -> Note: if we get singularity on corner, we need ntAdaptiveRotation right away
                // -> Macro is merely for clarity, but it stores both fixed and unfixed failures yet throws unfixed
                ERFMAPDISP(delta0,mapDisplacementUsingCache(start0,*pmap,p,theMapParams,deltaStorageSet,erFail));
                ERFMAPDISP(delta1,mapDisplacementUsingCache(start1,*pmap,p,theMapParams,deltaStorageSet,erFail));
                
                //vec_type beta0, beta1;
                //ERFMAPDISP(beta0,mapDisplacementUsingCache(start0,*pmap,-p,theMapParams,betaStorageSet,erFail));
                //ERFMAPDISP(beta1,mapDisplacementUsingCache(start1,*pmap,-p,theMapParams,betaStorageSet,erFail));
                
                
                //vec_type v0 = (delta0-beta0)/2.0;
                //vec_type v1 = (delta1-beta1)/2.0;
                vec_type v0 = delta0, v1 = delta1;
                if (theMapParams.debug) {
                    std::cout << "===================================================\n";
                    std::cout << " Period : " << p << "\n";
                    std::cout << "===================================================\n";
                    std::cout << " Corner Node Data :\n";
                    std::cout << "  x0    = " << start0 << "    x1 = " << start1 << "\n";
                    std::cout << "  delta0= " << delta0 << " delta1= " << delta1 << "\n";
                    //std::cout << "  beta0 = " << beta0 << "  beta1 = " << beta1  << "\n";
                    //std::cout << "  eta0  = " << v0 <<    "  eta1  + " << v1 << "\n";
                    
                }
                
                
                
                //Rotation value
                double theta = 0;
                if (lmin > 0) {//-----------------------------------------------------------------------------
                    //EndPoints
                    vec_type x0 = start0;
                    vec_type x1 = start1;
                    //Automatically add some subdivisions (default is 4)
                    double du = 1./(double)(nb_sub+1);
                    vec_type* x = new vec_type[nb_sub+2];
                    vec_type* v = new vec_type[nb_sub+2];
                    double u = du;
                    //Loop through automatic subdivisions and propagate map
                    for (int i=1; i<=nb_sub ; ++i, u+=du) {
                        vec_type delta;
                        x[i] = (1.-u)*x0 + u*x1;
                        //Forward Mapping
                        ERFMAPDISP(delta,mapDisplacementUsingCache(x[i],*pmap,p,theMapParams,deltaStorageSet,erFail));
                        v[i] = delta;
                        //Backward Mapping
                        //vec_type beta;
                        //ERFMAPDISP(beta,mapDisplacementUsingCache(x[i],*pmap,-p,theMapParams,betaStorageSet,erFail));
                        //Map Tangent vector
                        //v[i] = (delta - beta) / 2.0;
                        
                        //We need to check for fixed points here
                        ERotFailure fpTestFail;
                        if (isFixedPointSuspected(x[i],v[i],p,theMapParams,fpTestFail)) {
                            //if (isFixedPointSuspectedTangent(x[i],delta,beta,v[i],p,theMapParams,fpTestFail)) {
                            //Send fixed point guess to cache
                            if(verbose) std::cout << "ERotFailure: " << fpTestFail.what()
                                                      << " at " << fpTestFail.where() << "\n";
                            if(fpTestFail.type == ERotFailure::DOUBLE_PERIOD_FIXED_POINT) {
                                if (2*p <= theMapParams.max_period) {
                                    //Add point to guesses
                                    FP aGuess;
                                    aGuess.pos = x[i];
                                    aGuess.K = p;
                                    cached_fps[thread_id].push_back( aGuess );
                                }
                            } else if(fpTestFail.type == ERotFailure::FIXED_POINT_SUSPECTED) {
                                //Add point to guesses
                                FP aGuess;
                                aGuess.pos = x[i];
                                aGuess.K = p;
                                cached_fps[thread_id].push_back( aGuess );
                            }
                            //Mark as separation conditions
                            useSectionSeparationVersion = true;
                            //if mark as failure
                            //errorOccurred = true;
                            //e.mapError.insert( std::pair<int,bool>(p,errorOccurred) );
                            //e.rotFailures.push_back( fpTestFail );
                        }//End FPSuspected test
                    } //End 4 subdivision loop
                    x[0] = x0;
                    x[nb_sub+1] = x1;
                    v[0] = v0;
                    v[nb_sub+1] = v1;
                    //If mapping problem occurred that couldn't be fixed, we have to stop this period
                    if (errorOccurred) {
                        e.rotationAngleMap.insert( std::pair<int,double>(p,0.0) );
                        continue; //Go to next period
                    }
                    
                    
                    
                    //Adaptively subdivide along edge segments with adaptiveRotationAngle()
                    std::vector<ERotFailure> stableManCrossFails;
                    try {
                        if (!useSectionSeparationVersion) {
                            //Default version:  Calls itself, but doesn't handle transversality assumption violations
                            for (int kk=0 ; kk<=nb_sub ; ++kk)
                                theta += adaptiveRotationAngle(x[kk], v[kk], x[kk+1], v[kk+1],
                                                               *pmap, p, theMapParams, deltaStorageSet, betaStorageSet,
                                                               stableManCrossFails);
                            //if (theMapParams.debug) std::cout << "total dtheta for edge = " << theta << std::endl;
                            //cached_angles[thread_id].push_back(pair_type(e, theta));
                        }
                    } catch (ERotFailure& eRotFail) {
                        //Store errors
                        switch(eRotFail.type) {
                            case ERotFailure::BACKWARD_SECTION_SEP :
                                //If this is a SectionSeparation-type failure,
                                // we need to call the separate algorithm
                                useSectionSeparationVersion = true;
                                break;
                            case ERotFailure::FORWARD_SECTION_SEP :
                                //If this is a SectionSeparation-type failure,
                                // we need to call the separate algorithm
                                useSectionSeparationVersion = true;
                                break;
                            case ERotFailure::BACKWARD_SINGULARITY :
                                //If this is a SectionSeparation-type failure,
                                // we need to call the separate algorithm
                                useSectionSeparationVersion = true;
                                break;
                            case ERotFailure::FORWARD_SINGULARITY :
                                //If this is a SectionSeparation-type failure,
                                // we need to call the separate algorithm
                                useSectionSeparationVersion = true;
                                break;
                            case ERotFailure::DOUBLE_PERIOD_FIXED_POINT : {
                                if (2*p <= theMapParams.max_period) {
                                    //Add point to guesses
                                    FP aGuess;
                                    aGuess.pos = eRotFail.failurePos;
                                    aGuess.K = p;
                                    cached_fps[thread_id].push_back( aGuess );
                                }
                                if(verbose) {
                                    std::cout << "FPSuspected : " << eRotFail.what() << " at " << eRotFail.where() << "\n";
                                }
                                //Treat this as a separation point
                                useSectionSeparationVersion = true;
                            }
                            break;
                            case ERotFailure::FIXED_POINT_SUSPECTED : {
                                //Add point to guesses
                                FP aGuess;
                                aGuess.pos = eRotFail.failurePos;
                                aGuess.K = p;
                                cached_fps[thread_id].push_back( aGuess );
                                if(verbose) {
                                    std::cout << "FPSuspected : " << eRotFail.what() << " at " << eRotFail.where() << "\n";
                                }
                            }
                            break;
                            default :
                                //Unfixed map errors here mean something is not being detected properly
                                if(verbose) {
                                    std::cout << "FATAL Edge_Rotation Error: Unknown ERotFailure is thrown at adaptiveRotationAngle()\n";
                                    std::cout << " ->ERotFailure : " << eRotFail.what() << " at " << eRotFail.where() << "\n";
                                    std::cout << " ->Moving to next period indicating rotation failure...\n";
                                }
                                e.rotFailures.push_back( eRotFail );
                                errorOccurred = true; //For map errors
                                e.mapError.insert( std::pair<int,bool>(p,errorOccurred) );
                                e.rotationAngleMap.insert( std::pair<int,double>(p,0.0) );
                                continue; //To next period
                                break;
                        }
                        
                        //Check for special stable manifold crossing case:
                        if( (int)stableManCrossFails.size() > 0 &&
                                stableManCrossFails[0].type == ERotFailure::STABLE_MANIFOLD_CROSS) {
                            for(int smc=0; smc<(int)stableManCrossFails.size(); smc++) {
                                //Check midpoints
                                vec_type deltaMid(50.0);
                                try {
                                    //Should always be available
                                    deltaMid = mapDisplacementUsingCache(stableManCrossFails[smc].failurePos,*pmap,p,
                                                                         theMapParams,deltaStorageSet);
                                } catch(...) {
                                }
                                //Add point to guesses if within error tolerances (10*e_fps)
                                if(nvis::norm(deltaMid) <= 10.0*theMapParams.min_fp_tol) {
                                    FP aGuess;
                                    aGuess.pos = stableManCrossFails[smc].failurePos;
                                    aGuess.K = p;
                                    cached_fps[thread_id].push_back( aGuess );
                                }
                                e.rotFailures.push_back( stableManCrossFails[smc] );
                                if(verbose) std::cout << "Special Unresolved : " << stableManCrossFails[smc].what()
                                                          << " at " << stableManCrossFails[smc].where() << "\n";
                            }
                            
                        }
                        
                    }
                    /*catch(...) {
                      //These are unknown errors
                      theta = 0.0;
                      //Store Unknown ERotFailure
                      if (theMapParams.debug) {
                          std::cout << " EdgeRot Failure:  An unknown error occurred in adaptive subdivision!\n"
                             << " SubDivID : " << kk << " Period : " << p << "\n"
                             << " Between Coords: x0 = " << x[kk] << ", x1 = " << x[kk+1] << std::endl;
                      }
                      ERotFailure notSure("Unknown failure in adaptive subdivision", (x[kk]+x[kk+1])/2.0,p)
                      e.rotFailures.push_back( notSure );
                      errorOccurred = true;
                    }*/
                    
                    //Employ alternative version of adaptiveRotationAngle() that computes the rotation
                    // assuming that section transversality is violated
                    if(useSectionSeparationVersion) {
                        //Store that edge had to useSectionSeparationVersion
                        e.transverseSection[p] = false;
                        std::vector<vec_type> xVec, vVec;
                        for(int i=0; i<nb_sub+2; i++) {
                            xVec.push_back( x[i] );
                            vVec.push_back( v[i] );
                        }
                        try {
                            theta = adaptiveRotationAngleNonTransverse(e, edgeLength, xVec, vVec,
                                    *pmap, p, theMapParams, cached_fps[thread_id], deltaStorageSet, betaStorageSet);
                        } catch(ERotFailure& eRotFail) {
                            errorOccurred = true;
                            //Store the error
                            if(verbose) {
                                std::cout << "FATAL Edge_Rotation Error: ERotFailure is thrown at adaptiveRotationAngleNonTransverse()\n";
                                std::cout << " ->ERotFailure : " << eRotFail.what() << " at " << eRotFail.where() << "\n";
                                std::cout << " ->Moving to next period indicating rotation failure...\n";
                            }
                            e.rotFailures.push_back( eRotFail );
                            e.mapError.insert( std::pair<int,bool>(p,errorOccurred) );
                            e.rotationAngleMap.insert( std::pair<int,double>(p,0.0) );
                            continue; //To next period
                            
                        }
                        /// Just throw any other type of error to create runtime error
                        /*catch(...) {
                          //Assume all other errors are runtime errors
                          errorOccurred = true;
                          ERotFailure uf;
                          uf.period = p;
                          uf.failurePos = (x[nb_sub+1] - x[0])/2.0; //midpoint
                          uf.setMessage( "Unknown failure during adaptiveRotationAngleNonTransverse() call" );
                          if(verbose) {
                            std::cout << "FATAL Edge_Rotation Error: Unknown RUNTIMEERROR is thrown at adaptiveRotationAngleNonTransverse()\n";
                            std::cout << " ->Failure : " << uf.what() << " at midpoint " << uf.where() << "\n";
                            std::cout << " ->Moving to next period indicating rotation failure...\n";
                          }
                          e.rotFailures.push_back( uf );
                          e.mapError.insert( std::pair<int,bool>(p,errorOccurred) );
                          e.rotationAngleMap.insert( std::pair<int,double>(p,0.0) );
                          continue; //To next period
                        }*/
                    }
                    
                    
                    //Store Information (Angle,MapError,RotFailures,xMinDelta,minDelta):
                    if (errorOccurred) {
                        theta = 0.0;    //Shouldn't be error at this point, but just for safety
                    }
                    e.mapError.insert( std::pair<int,bool>(p,errorOccurred) ); //made it with no fatal problems
                    e.rotationAngleMap.insert( std::pair<int,double>(p,theta) );
                    
                    //Loop through FORWARD cache to find minimum delta points on edge (useful for fp guesses)
                    typename std::set<SortableData>::iterator dataIT;
                    vec_type xMinDelta(0.0),minDelta(0.0);
                    double minDeltaNorm = 1000.0;
                    for (dataIT = deltaStorageSet.begin(); dataIT != deltaStorageSet.end(); ++dataIT) {
                        //Check for minimum delta (used in non-transverse cells for sampling fixed point guesses)
                        vec_type delta(1000,1000);
                        if( dataIT->isThere(p) ) {
                            delta = theMapParams.the_metric.displacement(dataIT->x0,dataIT->returns[std::abs(p)-1]);
                        } else {
                            continue;    //Likely a singularity point, so skip
                        }
                        if (nvis::norm(delta) < minDeltaNorm ) {
                            xMinDelta = dataIT->x0;
                            minDelta = delta;
                            minDeltaNorm = nvis::norm(delta);
                        }
                    }
                    //Store the minimum delta information for guess generation if non-transverse section
                    e.xMinDeltaMap.insert( std::pair<int,vec_type>(p,xMinDelta));
                    e.minDeltaMap.insert( std::pair<int,vec_type>(p,minDelta));
                    
                    //Done with x,v so free space
                    delete[] x;
                    delete[] v;
                    
                } else { //----------------------------------------------------------------------------
                    //Quick and simple: -> And thus, hardly ever works!!!
                    theta = signed_angle(v0, v1);
                    //if (theMapParams.debug) std::cout << "direct angle = " << theta << std::endl;
                    if (fabs(theta) < theMapParams.max_angle) {
                        e.rotationAngleMap.insert( std::pair<int,double>(p,theta) );
                        e.mapError.insert( std::pair<int,bool>(p,false) );
                        //cached_angles[thread_id].push_back(std::pair<p_edge_type, double>(edges[n], theta));
                        //if (theMapParams.debug) std::cout << "valid angle. done" << std::endl;
                    } else {
                        //if (theMapParams.debug) std::cout << "angle is too large." << std::endl;
                        //cached_edges[thread_id].push_back(edges[n]);
                        e.rotationAngleMap.insert( std::pair<int,double>(p,theta) );
                        e.mapError.insert( std::pair<int,bool>(p,true) );
                        //Build an edge rotation failure object
                        ERotFailure eRotFail(ERotFailure::UNKNOWN,(start0+start1)/2.0,p); //use unknown error
                        e.rotFailures.push_back( eRotFail );
                    }
                }//------------------------------------------------------------------------------------
                
            } //End per period loop
            
            //After each period is added, we add the edge to the cache
            cached_edges[thread_id].push_back(e);
            
            #pragma omp atomic
            ++counter;
            //Output
            //if (thread_id == 0) {
            std::cerr << "\rComputed rotation angles for " << counter << " / "
                      << N << " edges (" << 100*counter/N
                      << "%)          \r" << std::flush;
            //}
            //Avizo Output if available
#ifdef HX_HAS_STD
            if (thread_id == 0) {
                //Update progress info
                float progress = (float)counter/(float)N;
                QString infoText = "Edge Rotation Computation at ";
                QString progText;
                progText.setNum(100.*progress);
                infoText += progText + "\% complete";
                theWorkArea->setProgressInfo(infoText);
                //Update progress slider
                theWorkArea->setProgressValue(progress);
            }
            
#endif
        }//End parallel loop
    }//End Parallel Statement
    
    
    //After the per-Edge computation is done, collect the cached edges per thread
    edges.clear();
    for (int i=0 ; i<nbthreads ; ++i) {
        typename std::vector<edge_type>::iterator vit;
        for (vit = cached_edges[i].begin(); vit != cached_edges[i].end(); ++vit) {
            edges.insert(*vit);
        }
        typename std::vector<FP>::iterator fpit;
        for (fpit=cached_fps[i].begin(); fpit!=cached_fps[i].end(); ++fpit) {
            fpGuesses.push_back(*fpit);
        }
    }
    std::cerr << "PMATE : Processed " << edges.size() << " / " << N << " edges.\n";
    if (N != edges.size()) {
        std::cerr << "        WARNING!!! Some edges missing!!!!\n";
    }
    double dt = _timer.elapsed();
    std::cerr << "        Edge Rotation computation took "
              << dt << " s. (" << dt / (double)N << " s/edge) \n";
    //Free up data
    for (int i=0; i<nbthreads; ++i) {
        cached_edges[i].clear();
        cached_fps[i].clear();
    }
    delete[] cached_edges;
    delete[] cached_fps;
    
}

/** Compute the composite edge rotation after each of the TerminalEdges' rotation
 *  is computed (i.e., call AFTER a call to Edge_Rotation()).  This is necessary for
 *  Adaptive Grid edges in cells
 *  Author: Wayne Schlei
 */
template< typename SCOMPOSITE, typename SEDGE>
void compositeEdgeRotation(SCOMPOSITE& compositeEdgeSet, SEDGE& baseEdgeSet)
{
    typedef typename SEDGE::iterator         BSetIter;
    typedef typename SCOMPOSITE::iterator     CSetIter;
    typedef typename SEDGE::value_type        TerminalEdge;
    typedef typename SCOMPOSITE::value_type     CompositeEdge;
    typedef typename CompositeEdge::NodeSetType::iterator NodeSetIter;
    
    int nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    
    std::vector<CompositeEdge>* cached_edges = new std::vector<CompositeEdge>[nbthreads];
    
    nvis::timer _timer;
    
    //Copy to a std::vector so we can edit objects
    std::vector<CompositeEdge> copyOfEdges(compositeEdgeSet.begin(),compositeEdgeSet.end());
    size_t N = copyOfEdges.size();
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n=0; n<(int)N; n++) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            //Create an edge copy
            CompositeEdge e( copyOfEdges[n] );
            
            //Tally Quantities from each sub edge
            //int numSegments = (int)e.nodeSet.size() - 1;
            std::set<int>::iterator pit;
            //For each period for the mother edge
            for(pit=e.periods.begin(); pit!=e.periods.end(); ++pit) {
                int p = *pit;
                double theta = 0.0;
                bool mapErrorOccurred = false;
                bool sectionTransversality = true;
                NodeSetIter nit, nit2;
                nit = e.nodeSet.begin();
                nit2 = nit;
                ++nit2;
                double minDeltaNorm = 1000.0;
                nvis::vec2 xMinDelta(0.0), minDelta(0.0);
                for(; nit2 != e.nodeSet.end(); ++nit,++nit2 ) {
                    //Find the child edge
                    TerminalEdge tempEdge((*nit),(*nit2));
                    BSetIter subEdgeIt = baseEdgeSet.find( tempEdge );
                    //Compile information that is already computed
                    theta += subEdgeIt->rotationAngleMap.find(p)->second;
                    if (subEdgeIt->mapError.find(p)->second) {
                        mapErrorOccurred = true;
                    }
                    if (!(subEdgeIt->transverseSection.find(p)->second)) {
                        sectionTransversality = false;
                    }
                    //Find minimum point per period only if there is no map error
                    if (!(subEdgeIt->mapError.find(p)->second) ) {
                        if(nvis::norm(subEdgeIt->minDeltaMap.find(p)->second) < minDeltaNorm) {
                            xMinDelta = subEdgeIt->xMinDeltaMap.find(p)->second;
                            minDelta = subEdgeIt->minDeltaMap.find(p)->second;
                            minDeltaNorm = nvis::norm(minDelta);
                        }
                    }
                }
                //Assign rotation angle and map errors per period
                if(mapErrorOccurred) {
                    theta = 0.0;
                }
                e.rotationAngleMap.insert( std::pair<int,double>(p,theta) );
                e.mapError.insert( std::pair<int,bool>(p,mapErrorOccurred) );
                e.transverseSection.insert( std::pair<int,bool>(p,sectionTransversality) );
                e.xMinDeltaMap.insert( std::pair<int,nvis::vec2>(p,xMinDelta) );
                e.minDeltaMap.insert( std::pair<int,nvis::vec2>(p,minDelta) );
                //Ignore the assignment of failure locations (handled in base edges)
            }
            
            //Store to cache
            cached_edges[thread_id].push_back( e );
        }
        
    }//End parallel
    
    //After per-CompEdge compute, reassemble the edge set based on new information
    compositeEdgeSet.clear();
    for(int i=0; i<nbthreads; ++i) {
        typename std::vector<CompositeEdge>::iterator vit;
        for (vit=cached_edges[i].begin(); vit!=cached_edges[i].end(); ++vit) {
            compositeEdgeSet.insert( *vit );
        }
    }
    //Delete temp cache
    for(int i=0; i<nbthreads; ++i) {
        cached_edges[i].clear();
    }
    delete[] cached_edges;
}


} //end xavier

#endif
