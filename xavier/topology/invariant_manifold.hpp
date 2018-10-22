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


//////////////////////////////////////////////////////////////////////////////
// invariant_manifold.hpp
// Author:  Xavier Tricoche & Wayne Schlei
// Date:  2/21/2013
// Purpose:  Utilize the ManBVP algorithm for computing manifolds on a
// Poincare section (from Computing 1D Global Manifolds By Continuation by
// England et al.).
// Also now (3/15/2015) contains the call to the ManCurve algorithm which
// treats the manifold computation as a curve-refinement problem rather
// than a boundary value problem like ManBVP.  Should be more efficient
// for computationally expensive maps.
//
// Update notice:
// The original version of this has been modified to make algorithm
// parameters more tractable for CR3BP.  Also, the stopping criteria was
// modified to include leaving a region of interest.
//
// Future Work:
// 1)This propagates manifolds on a 2D surface of section; for work with higher
// dimensional sections, be sure to conform to a template that specifies dimensions.
// 2)This stops when you leave a region of interest, but it might be nice to
// check if a manifold will come back (or add ability to add segements that return).
///////////////////////////////////////////////////////////////////////////////

#ifndef INVARIANT_MANIFOLD_HPP
#define INVARIANT_MANIFOLD_HPP

//STL
#include <vector>
#include <list>
#include <iostream>
#include <sstream>
//NVIS
#include <math/fixed_vector.hpp>
#include <math/bezier_spline.hpp>
#include <math/bounding_box.hpp>
//Poincare API
#include <maps/metric.hpp>
#include <maps/fixpoints.hpp>
#include <maps/map_analysis.hpp>
#include <topology/ManifoldClasses.hpp>
#include <topology/EnglandManifold.hpp>
#include <topology/STHManifold.hpp>
//OpenMP
#if _OPENMP
#include <omp.h>
#endif

using namespace xavier;

namespace topology {



//Functions: -----------------------------------------------------------------------------
/// Computing minimum distance between chain points
template <class FPCHAIN>
double min_dist(const std::vector<FPCHAIN>& chains, const xavier::metric_type& _metric)
{
    std::vector<double> dist;
    std::vector<nvis::vec2> all_pos;
    
    for (int i=0 ; i<(int)chains.size() ; ++i) {
        for (int j=0 ; j<(int)chains[i].size() ; ++j) {
            all_pos.push_back(chains[i][j].pos);
        }
    }
    
    for (int i=0 ; i<(int)all_pos.size()-1 ; ++i) {
        for (int j=i+1 ; j<(int)all_pos.size() ; ++j) {
            dist.push_back(_metric.distance(all_pos[i], all_pos[j]));
        }
    }
    
    return *std::min_element(dist.begin(), dist.end());
}
/// Check if manifold segements were computed properly)
template <class SEPX>
bool sanity_check(const SEPX& sep, double alpha, const metric_type& _metric)
{

    //double huge = 0.05*_metric.diameter();
    double huge = 100; //Need a better number here (but _metric.diameter() = 0 for CR3BP)
    nvis::vec2 p0, p1, p2;
    const std::vector<nvis::vec2>& m = sep.manifold;
    if (m.size() < 2) {
        std::cerr << "\n\nnot enough points at " << m[0] << "\n\n\n" << std::endl;
        return false;
    }
    int c = 0;
    p1 = m[c++];
    p2 = m[c++];
    if (_metric.distance(p1, p2) > huge) {
        std::cerr << "segment too large (" << _metric.distance(p1, p2) << ") at p0 = " << p0 << "\n\n\n" << std::endl;
        return false;
    }
    while (c < (int)m.size()) {
        p0 = p1;
        p1 = p2;
        p2 = m[c++];
        if (_metric.distance(p1, p2) > huge) {
            std::cerr << "segment too large (" << _metric.distance(p1, p2) << ") at m[0] = " << m[0] << "\n\n\n" << std::endl;
            return false;
        }
    }
    return true;
}


/// Computing stable/unstable manifolds for every fp_chain - Used England Manifold
template< typename MAP >
bool compute_separatrix(std::vector<fp_chain>& all_p_chains,
                        std::vector<Separatrix>& separatrices,
                        unsigned int chain_id, const MAP& map,
                        const ManifoldSettings& settings, const map_analysis_param& param,
                        std::vector< std::vector<nvis::vec2> >& broken_manifolds)
{
    /* Upgrade thought:
     * Rework code to run as a massive openmp loop
     * This would require a data structure that stores everything each manifold (2 or 4) of each
     * periodic orbit (p on [1,pmax]) of each fixed-point chain (all_p_chains.size()) would
     * require as an isolated computation.
     */
    typedef typename MAP::section_type::lvec_type lvec_type;
    
    separatrices.clear();
    //param.verbose = true;
    
    //Copy the metric
    const metric_type _metric(param.the_metric);
    
    //Copy of section
    typename MAP::section_type theSection(map.section());
    
    //Make sure id is in range
    if (!all_p_chains.size() || chain_id >= all_p_chains.size()) {
        return true;
    }
    //Run for only one chain at present
    const fp_chain& chain = all_p_chains[chain_id];
    
    if (!chain.saddle()) {
        return true;    // that was easy!
    }
    
    int period = chain.period();
    
    if (param.verbose)
        std::cout << "\n\nprocessing saddle chain of period " << period
                  << " starting at " << chain[0].pos << std::endl;
                  
    //Use data structure for indicating settings
    /*
    toggles.delta_min = (period == 1 ? 1.e-4 : 1.e-3);
    double eps = 0.0001 * _metric.diameter();
    double sdelta_max = 0.005*_metric.diameter();
    double alpha_max = 0.1;
    double liberal_alpha_max = M_PI/2.;
    double delta_min = (period == 1 ? 1.0e-4 : 1.0e-3) * _metric.diameter();
    double delta_length = _metric.diameter();
    */
    manifold_type man1, man2, man3, man4;
    std::vector<int> bks1, bks2, bks3, bks4;
    std::vector<int>::iterator bkit;
    
    std::cerr << "compute_Separatrix called for period " << period << std::endl;
    
    ManifoldStopper<typename MAP::section_type> stop(all_p_chains, _metric, param.bounds, settings);
    Perturbation mType = UNSTABLE_PLUS;
    
    try {
        for (int i=0 ; i<(int)chain.size() ; ++i) {
            ManifoldProgress prog;
            std::cerr << "\n\nprocessing saddle " << i+1 << " from " << chain.size() << " at " << chain[i].pos << std::endl;
            
            std::cerr << "\nunstable manifold #0..." << std::endl;
            mType = UNSTABLE_PLUS;
            stop.start(fp_index_type(chain_id, i));
            prog.length = 0;
            ManBVP(man1, bks1, prog, map, stop, period, _metric, chain[i], mType, settings, param);
            separatrices.push_back(Separatrix());
            std::copy(man1.begin(), man1.end(), std::back_inserter(separatrices.back().manifold));
            for(bkit=bks1.begin(); bkit!=bks1.end(); ++bkit) {
                separatrices.back().breakIDs.insert(*bkit);
            }
            separatrices.back().start = fp_index_type(chain_id, i);
            separatrices.back().end = stop.connection();
            separatrices.back().forward = true;
            separatrices.back().length = prog.length;
            if (!sanity_check(separatrices.back(), settings.liberal_alpha_max, _metric)) {
            
                broken_manifolds.push_back(separatrices.back().manifold);
                
                std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
                separatrices.back().manifold.resize(2);
                lvec_type& x0 = separatrices.back().manifold[0];
                lvec_type& x1 = separatrices.back().manifold[1];
                x0 = chain[i].pos;
                x1 = x0 + 0.001 * _metric.displacement(x0, x1) / _metric.distance(x0, x1);
            }
            std::cerr << "length = " << prog.length << ", " << man1.size() << " vertices\n";
            
            //Stable Manifold if mirror theorem applies - note: also need symmetric orbit too!
            if ( theSection.isSymmetric() ) {
                //Symmetric section, so our stable manifolds are simply found through the Mirror Theorem
                std::cerr << "\n    Symmetric section: Applying Mirror Thm to get stable manifold #0..." << std::endl;
                manifold_type::iterator mit;
                std::vector<int>::iterator bkit;
                //Apply Mirror Thm to man1 to get man3
                man3.clear();
                bks3.clear();
                for (mit = man1.begin(); mit!=man1.end(); mit++) {
                    man3.push_back( theSection.mirror(*mit) );
                }
                for (bkit = bks1.begin(); bkit!=bks1.end(); bkit++) {
                    bks3.push_back( *bkit );
                }
                separatrices.push_back(Separatrix());
                std::copy(man3.begin(),man3.end(), std::back_inserter(separatrices.back().manifold));
                for(bkit=bks3.begin(); bkit!=bks3.end(); ++bkit) {
                    separatrices.back().breakIDs.insert(*bkit);
                }
                //Search for nearest fp since we actually changed to a new chain point by applying Mirror Thm
                //---------------------------------------------------------------------------------------------------------
                separatrices.back().start = fp_index_type(chain_id,
                                            (int) chain.closest(theSection.mirror(chain[i].pos)));
                separatrices.back().end = fp_index_type(chain_id,
                                                        (int) chain.closest(theSection.mirror(chain[stop.connection().second].pos)) );
                //---------------------------------------------------------------------------------------------------------
                separatrices.back().forward = false;
                separatrices.back().length = prog.length;
                if (!sanity_check(separatrices.back(), settings.liberal_alpha_max, _metric)) {
                
                    broken_manifolds.push_back(separatrices.back().manifold);
                    
                    std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
                    separatrices.back().manifold.resize(2);
                    nvis::vec2& x0 = separatrices.back().manifold[0];
                    lvec_type& x1 = separatrices.back().manifold[1];
                    x0 = chain[i].pos;
                    x1 = x0 + 0.001 * _metric.displacement(x0, x1) / _metric.distance(x0, x1);
                }
                std::cerr << "length = " << prog.length << ", " << man3.size() << " vertices\n";
            }
            
            std::cerr << "\nunstable manifold #1..." << std::endl;
            mType = UNSTABLE_MINUS;
            stop.start(fp_index_type(chain_id, i));
            prog.length = 0;
            ManBVP(man2, bks2, prog, map, stop, period, _metric, chain[i], mType, settings, param);
            separatrices.push_back(Separatrix());
            std::copy(man2.begin(), man2.end(), std::back_inserter(separatrices.back().manifold));
            for(bkit=bks2.begin(); bkit!=bks2.end(); ++bkit) {
                separatrices.back().breakIDs.insert(*bkit);
            }
            separatrices.back().start = fp_index_type(chain_id, i);
            separatrices.back().end = stop.connection();
            separatrices.back().forward = true;
            separatrices.back().length = prog.length;
            if (!sanity_check(separatrices.back(), settings.liberal_alpha_max, _metric)) {
            
                broken_manifolds.push_back(separatrices.back().manifold);
                
                std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
                separatrices.back().manifold.resize(2);
                lvec_type& x0 = separatrices.back().manifold[0];
                lvec_type& x1 = separatrices.back().manifold[1];
                x0 = chain[i].pos;
                x1 = x0 + 0.001 * _metric.displacement(x0, x1) / _metric.distance(x0, x1);
            }
            std::cerr << "length = " << prog.length << ", " << man2.size() << " vertices\n";
            
            //Stable Manifold if mirror theorem applies
            if ( theSection.isSymmetric() ) {
                //Symmetric section, so our stable manifolds are simply found through the Mirror Theorem
                std::cerr << "\n    Symmetric section: Applying Mirror Thm to get stable manifold #1..." << std::endl;
                manifold_type::iterator mit;
                std::vector<int>::iterator bkit;
                //Apply Mirror Thm to man1 to get man3
                man4.clear();
                bks4.clear();
                for (mit = man2.begin(); mit!=man2.end(); mit++) {
                    man4.push_back( theSection.mirror(*mit) );
                }
                for (bkit = bks2.begin(); bkit!=bks2.end(); bkit++) {
                    bks4.push_back( *bkit );
                }
                separatrices.push_back(Separatrix());
                std::copy(man4.begin(),man4.end(), std::back_inserter(separatrices.back().manifold));
                for(bkit=bks4.begin(); bkit!=bks4.end(); ++bkit) {
                    separatrices.back().breakIDs.insert(*bkit);
                }
                //Search for nearest fp since we actually changed to a new chain point by applying Mirror Thm
                //-----------------------------------------------------------------------------------------------
                separatrices.back().start = fp_index_type(chain_id, (int) chain.closest(theSection.mirror(chain[i].pos)));
                separatrices.back().end = fp_index_type(chain_id,
                                                        (int) chain.closest(theSection.mirror(chain[stop.connection().second].pos))
                                                       );
                //----------------------------------------------------------------------------------------------
                separatrices.back().forward = false;
                separatrices.back().length = prog.length;
                if (!sanity_check(separatrices.back(), settings.liberal_alpha_max, _metric)) {
                
                    broken_manifolds.push_back(separatrices.back().manifold);
                    
                    std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
                    separatrices.back().manifold.resize(2);
                    lvec_type& x0 = separatrices.back().manifold[0];
                    lvec_type& x1 = separatrices.back().manifold[1];
                    x0 = chain[i].pos;
                    x1 = x0 + 0.001 * _metric.displacement(x0, x1) / _metric.distance(x0, x1);
                }
                std::cerr << "length = " << prog.length << ", " << man3.size() << " vertices\n";
            }
            
            
            if ( !theSection.isSymmetric() ) {
                //Run ManBVP for both stable manifolds
                std::cerr << "\nstable manifold #0..." << std::endl;
                mType = STABLE_PLUS;
                stop.start(fp_index_type(chain_id, i));
                prog.length = 0;
                ManBVP(man3, bks3, prog, map, stop, -period, _metric, chain[i], mType, settings, param);
                separatrices.push_back(Separatrix());
                std::copy(man3.begin(), man3.end(), std::back_inserter(separatrices.back().manifold));
                for(bkit=bks3.begin(); bkit!=bks3.end(); ++bkit) {
                    separatrices.back().breakIDs.insert(*bkit);
                }
                separatrices.back().start = fp_index_type(chain_id, i);
                separatrices.back().end = stop.connection();
                separatrices.back().forward = false;
                separatrices.back().length = prog.length;
                if (!sanity_check(separatrices.back(), settings.liberal_alpha_max, _metric)) {
                
                    broken_manifolds.push_back(separatrices.back().manifold);
                    
                    std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
                    separatrices.back().manifold.resize(2);
                    lvec_type& x0 = separatrices.back().manifold[0];
                    lvec_type& x1 = separatrices.back().manifold[1];
                    x0 = chain[i].pos;
                    x1 = x0 + 0.001 * _metric.displacement(x0, x1) / _metric.distance(x0, x1);
                }
                std::cerr << "length = " << prog.length << ", " << man3.size() << " vertices\n";
                
                std::cerr << "\nstable manifold #1..." << std::endl;
                mType = STABLE_MINUS;
                stop.start(fp_index_type(chain_id, i));
                prog.length = 0;
                ManBVP(man4, bks4, prog, map, stop, -period, _metric, chain[i], mType, settings, param);
                separatrices.push_back(Separatrix());
                std::copy(man4.begin(), man4.end(), std::back_inserter(separatrices.back().manifold));
                for(bkit=bks4.begin(); bkit!=bks4.end(); ++bkit) {
                    separatrices.back().breakIDs.insert(*bkit);
                }
                separatrices.back().start = fp_index_type(chain_id, i);
                separatrices.back().end = stop.connection();
                separatrices.back().forward = false;
                separatrices.back().length = prog.length;
                if (!sanity_check(separatrices.back(), settings.liberal_alpha_max, _metric)) {
                
                    broken_manifolds.push_back(separatrices.back().manifold);
                    
                    std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
                    separatrices.back().manifold.resize(2);
                    lvec_type& x0 = separatrices.back().manifold[0];
                    lvec_type& x1 = separatrices.back().manifold[1];
                    x0 = chain[i].pos;
                    x1 = x0 + 0.001 * _metric.displacement(x0, x1) / _metric.distance(x0, x1);
                }
                std::cerr << "length = " << prog.length << ", " << man4.size() << " vertices\n";
            }
        }
    } catch(...) {
        std::cerr << "Caught an error while running ManBVP\n";
        return false;
    }
    
    
    return true;
}

/// Computing stable/unstable manifolds for a SINGLE orbit Using STH Manifold (curve-refinement)
template< typename MAP, typename PARAM >
bool compute_separatrix(std::vector<fp_chain>& all_p_chains,
                        std::vector<Separatrix>& separatrices,
                        unsigned int chain_id, const MAP& theMap,
                        const ManifoldSettings& settings, const PARAM& params
                       )
{
    typedef typename MAP::section_type::lvec_type lvec_type;
    typedef std::pair<int,Perturbation>           RunPair;
    
    //For printing perturbation
    std::map<Perturbation,std::string> pertNameMap;
    pertNameMap[UNSTABLE_PLUS] = std::string("UNSTABLE_PLUS");
    pertNameMap[UNSTABLE_MINUS] = std::string("UNSTABLE_MINUS");
    pertNameMap[STABLE_PLUS] = std::string("STABLE_PLUS");
    pertNameMap[STABLE_MINUS] = std::string("STABLE_MINUS");
    
    separatrices.clear();
    //params.verbose = true;
    
    //Copy the metric
    const metric_type& theMetric = params.the_metric;
    
    //Copy of section
    typename MAP::section_type theSection = theMap.section();
    
    //Make sure id is in range
    if (!all_p_chains.size() || chain_id >= all_p_chains.size()) {
        return true;
    }
    //Run for only one chain at present
    const fp_chain& chain = all_p_chains[chain_id];
    
    if (!chain.saddle()) {
        return true;
    }
    
    int period = chain.period();
    
    if (params.verbose)
        std::cout << "\n\nprocessing saddle chain of period " << period
                  << " starting at " << chain[0].pos << std::endl;
                  
    std::cerr << "compute_separatrix called for period " << period << std::endl;
    
    ManifoldStopper<typename MAP::section_type> stop(all_p_chains, theMetric, params.bounds, settings);
    
    //Create the execution packets per thread
    std::vector<RunPair> runData;
    for(int i=0; i<(int)chain.size(); ++i) {
        runData.push_back( RunPair(i,UNSTABLE_PLUS) );
        runData.push_back( RunPair(i,UNSTABLE_MINUS) );
        runData.push_back( RunPair(i,STABLE_PLUS) );
        runData.push_back( RunPair(i,STABLE_MINUS) );
    }
    int numManifolds = (int)runData.size();
    int counter = 0;
    int nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    std::vector<Separatrix>* cachedSEPs = new std::vector<Separatrix>[nbthreads];
    
    //Compute the Manifolds with ManCurve (STHManifold)
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n=0; n<numManifolds; n++) {
            int i = runData[n].first;
            Perturbation mType = runData[n].second;
            std::ostringstream os;
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            //Propagation direction (stable or unstable)
            bool fwd = (mType == UNSTABLE_PLUS || mType == UNSTABLE_MINUS) ? true : false;
            ManifoldProgress prog;
            //std::cerr << "\n\nprocessing saddle " << i+1 << " from " << chain.size() << " at " << chain[i].pos << std::endl;
            
            //Create a new Stopper
            ManifoldStopper<typename MAP::section_type> nStopper( stop );
            nStopper.start(fp_index_type(chain_id, i));
            prog.length = 0;
            //Compute Separatrix
            Separatrix tempSeparatrix;
            ManCurve(tempSeparatrix, prog, theMap, nStopper, params, chain[i], mType, settings);
            tempSeparatrix.start = fp_index_type(chain_id, i);
            tempSeparatrix.end = nStopper.connection();
            tempSeparatrix.forward = fwd;
            tempSeparatrix.length = prog.length;
            cachedSEPs[thread_id].push_back(tempSeparatrix);
            //std::cerr << "length = " << prog.length << ", " << tempSeparatrix.manifold.size() << " vertices\n";
            //Use Mirror Theorem IF theSection.isSymmetric() and chain[i].isSymmetric()
            
            #pragma omp atomic
            ++counter;
            
            //Update text
            os << "\rManifoldComplete: n=" << counter << " / " << numManifolds << " ("
               << 100*n/numManifolds << "%) : " << pertNameMap[mType]
               << " with Arclength = " << prog.length << "\r";
            std::cout << os.str() << std::flush;
            //Avizo?
        }
        /*}
        catch(...) {
            std::cerr << "Caught an error while running ManCurve\n";
            return false;
        }*/
    } //end parallel statement
    
    //Copy information per thread to output
    for(int i=0; i<nbthreads; i++) {
        typename std::vector<Separatrix>::iterator sit;
        for(sit=cachedSEPs[i].begin(); sit!=cachedSEPs[i].end(); ++sit) {
            separatrices.push_back( *sit );
        }
    }
    for(int i=0; i<nbthreads; i++) {
        cachedSEPs[i].clear();
    }
    delete[] cachedSEPs;
    //Worked correctly
    return true;
}

/// Extended run data for parallel runs
struct SepInput {
    SepInput(const int& chain_id, const int& fp_id, const Perturbation& type) :
        orbitID(chain_id), fpID(fp_id), mType(type) {}
        
    int orbitID;
    int fpID;
    Perturbation mType;
};

/// Computing stable/unstable manifolds for ALL input orbits Using STH Manifold (curve-refinement)
/// with ONE manifold per thread
template< typename MAP, typename PARAM >
bool compute_separatrix(std::vector<fp_chain>& all_p_chains,
                        std::vector<Separatrix>& separatrices,
                        const MAP& theMap,
                        const ManifoldSettings& settings, const PARAM& params
                       )
{
    typedef typename MAP::section_type::lvec_type lvec_type;
    
    //For printing perturbation
    std::map<Perturbation,std::string> pertNameMap;
    pertNameMap[UNSTABLE_PLUS] = std::string("UNSTABLE_PLUS");
    pertNameMap[UNSTABLE_MINUS] = std::string("UNSTABLE_MINUS");
    pertNameMap[STABLE_PLUS] = std::string("STABLE_PLUS");
    pertNameMap[STABLE_MINUS] = std::string("STABLE_MINUS");
    
    separatrices.clear();
    //params.verbose = true;
    
    //Copy the metric
    const metric_type& theMetric = params.the_metric;
    
    //Copy of section
    typename MAP::section_type theSection = theMap.section();
    
    //Make sure some exist
    if (!all_p_chains.size()) {
        return true;
    }
    
    //Stopper
    ManifoldStopper<typename MAP::section_type> stop(all_p_chains, theMetric, params.bounds, settings);
    
    //Setup run data for each fixed point
    std::vector<SepInput> runData;
    for(int chain_id=0; chain_id<(int)all_p_chains.size(); ++chain_id) {
        const fp_chain& chain = all_p_chains[chain_id];
        if (!chain.saddle()) {
            continue;
        }
        //int period = chain.period();
        //Create the execution packets per thread
        for(int i=0; i<(int)chain.size(); ++i) {
            runData.push_back( SepInput(chain_id,i,UNSTABLE_PLUS) );
            runData.push_back( SepInput(chain_id,i,UNSTABLE_MINUS) );
            runData.push_back( SepInput(chain_id,i,STABLE_PLUS) );
            runData.push_back( SepInput(chain_id,i,STABLE_MINUS) );
        }
    }
    
    
    int numManifolds = (int)runData.size();
    int counter = 0;
    int nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    std::vector<Separatrix>* cachedSEPs = new std::vector<Separatrix>[nbthreads];
    
    //Compute the Manifolds with ManCurve (STHManifold)
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n=0; n<numManifolds; n++) {
            int chain_id = runData[n].orbitID;
            int i = runData[n].fpID;
            Perturbation mType = runData[n].mType;
            const fp_chain& chain = all_p_chains[chain_id];
            std::ostringstream os;
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            //Propagation direction (stable or unstable)
            bool fwd = (mType == UNSTABLE_PLUS || mType == UNSTABLE_MINUS) ? true : false;
            ManifoldProgress prog;
            //std::cerr << "\n\nprocessing saddle " << i+1 << " from " << chain.size() << " at " << chain[i].pos << std::endl;
            
            //Create a new stopper
            ManifoldStopper<typename MAP::section_type> nStopper( stop );
            nStopper.start(fp_index_type(chain_id, i));
            prog.length = 0;
            Separatrix tempSeparatrix;
            ManCurve(tempSeparatrix, prog, theMap, nStopper, params, chain[i], mType, settings);
            tempSeparatrix.start = fp_index_type(chain_id, i);
            tempSeparatrix.end = nStopper.connection();
            tempSeparatrix.forward = fwd;
            tempSeparatrix.length = prog.length;
            cachedSEPs[thread_id].push_back(tempSeparatrix);
            //std::cerr << "length = " << prog.length << ", " << tempSeparatrix.manifold.size() << " vertices\n";
            //Use Mirror Theorem IF theSection.isSymmetric() and chain[i].isSymmetric()
            
            #pragma omp atomic
            ++counter;
            
            //Update text
            os << "\rManifoldComplete: n=" << counter << " / " << numManifolds << " ("
               << 100*n/numManifolds << "%) : " << pertNameMap[mType]
               << " with Arclength = " << prog.length << "\r";
            std::cout << os.str() << std::flush;
            //Avizo?
        }
        /*}
        catch(...) {
            std::cerr << "Caught an error while running ManCurve\n";
            return false;
        }*/
    } //end parallel statement
    
    //Copy information per thread to output
    for(int i=0; i<nbthreads; i++) {
        typename std::vector<Separatrix>::iterator sit;
        for(sit=cachedSEPs[i].begin(); sit!=cachedSEPs[i].end(); ++sit) {
            separatrices.push_back( *sit );
        }
    }
    
    //Worked correctly
    return true;
}


/// Using ManCurve Algorith to compute ALL manifolds using a Queue in OpenMP (tasks)
/*template< typename MAP, typename PARAM >
bool compute_separatrix_wTasks(std::vector<fp_chain>& all_p_chains,
                        std::vector<Separatrix>& separatrices,
                        const MAP& theMap,
                        const ManifoldSettings& settings, const PARAM& params
                        )
{
    typedef typename MAP::section_type::lvec_type lvec_type;


    //For printing perturbation
    std::map<Perturbation,std::string> pertNameMap;
    pertNameMap[UNSTABLE_PLUS] = std::string("UNSTABLE_PLUS");
    pertNameMap[UNSTABLE_MINUS] = std::string("UNSTABLE_MINUS");
    pertNameMap[STABLE_PLUS] = std::string("STABLE_PLUS");
    pertNameMap[STABLE_MINUS] = std::string("STABLE_MINUS");

    //Caches per thread
    int nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    std::vector<ManifoldMapData> mapDataPerThread[nbthreads];

    //Setup Initial Manifold Objects (4 per fixed point)
    bool allManifoldsDone = falfse;
    //Assign kids from same orbit (order is important)

    //Assign kids from Mirror Theorem

    //Main Processing Loop
    #pragma omp parallel
    {
      //Switch to single construct for scheduler
      #pragma omp single
      {
        //Anything here without "task" directive will execute on a single thread...
        while (!allManifoldsDone) {
            //For Each Manifold(Adult)
            for(int m=0;m<numManifolds;m++){
              //If(kid), do nothing since parents do the work

              //If (beginning) run start procedure

              //Evaluate downstream-curve criteria between points on current working segment

              //If curveOK, evaluate stopping criteria
                 //Continue to next working seg or stop

              //Subdivision if !curveOK

              //Run map calls within task directives
              #pragma omp task depend(out:a)
              {
                manifoldObjects[m].taskManifolds(nbthreads,mapDataPerThread);
              }
              //Wait unitl the map calls are processed
                  //#pragma omp taskwait //Note: Waits for ALL tasks to complete
              //Pull data from threads
              #pragma omp task depend(in:a)


              //Update kids

              //If kids are mature, need to grow up
            }
            //Evaluate total complete (SUM( max(arclength_i,max_arclength) / max_arclength)/numManifolds)
        }

      } //End single construct
    } //End Parallel loop

    //Indicate success
    return true;
}*/

} // topology

#endif
