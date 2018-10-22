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


#ifndef __MAP_ANALYSIS_HPP__
#define __MAP_ANALYSIS_HPP__

// std
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <iomanip> //C++ 11
#include <sstream>
#include <set>
#include <iterator>
#include <algorithm>
#include <queue>
#include <memory> //C++ 11
// boost
#include <boost/format.hpp>
#include <boost/limits.hpp>
#include <boost/rational.hpp>
//#include <boost/shared_ptr.hpp> //boost alternative?
// nvis
#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>
// xavier
#include <data/grid.hpp>
#include <data/raster_data.hpp>
#include <data/adaptive_grid.hpp>
#include <math/angle.hpp>
#include <math/sort.hpp>
#include <math/rational.hpp>
#include <maps/definitions.hpp>
#include <maps/mapExceptions.hpp>
#include <maps/index.hpp>
#include <maps/adaptiveGridMapInterface.hpp>
#include <maps/cellConvexityChecker.hpp>
// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/NumericalDiff>
#include <unsupported/Eigen/NonLinearOptimization>
#include <maps/models.hpp>
//OpenMP
#if _OPENMP
#include <omp.h>
#endif
//Avizo
#if HX_HAS_STD
#include <hxcore/HxMessage.h>
#include <hxcore/internal/HxWorkArea.h>
#endif

using namespace Eigen;

namespace {
template<typename T>
inline T _sign(T x)
{
    return (x<0) ? -1 : 1;
}
}

namespace xavier {

/// Compute the average distance between integration steps for period p
double average_distance(const std::vector<nvis::vec2>& steps, int p,
                        const metric_type& _metric);
                        
///Best period based on the average distance for integration steps ->Best for periodic domains
int best_period(const std::vector<vec_type>& steps, int maxp, const metric_type& m);

///Computing best approximate toroidal periods based on average distance of all points ->Best for periodic domains
void best_periods(std::vector<int>& periods, const std::vector<vec_type>& steps,
                  int maxp, int nper, const metric_type& m);
                  
/// Computes the set of possible periods based on the input set valid_rationals
template<class DATA, typename IVEC>
void period_range(std::vector<int>& ps,
                  const DATA& dataset,
                  const IVEC& cell_id,
                  const std::map<double, rational_type>& valid_rationals);
                  
/// Set the range of possible periods based on what's available in dataset
template<class DATA, typename IVEC>
void period_range(std::vector<int>& ps,
                  const DATA& dataset,
                  const IVEC& cell_id);
                  
/// Compute the set of possible periods based on the max_p value and a continued fraction algorithm
template<class DATA, typename IVEC>
void period_range(std::vector<int>& ps,
                  const DATA& dataset,
                  const IVEC& cell_id,
                  const int& maxp, const bool watchCell=false);
                  
/// A function to get the winding ratios at a given point in a dataset (the winding numbers are attached)
template<class DATASET,typename IVEC>
void getWindingRatios(const DATASET& data,const IVEC& id,const int p_max,
                      std::vector< boost::rational<int> >& alpha);
                      
/// A function to generate the plausible factors that might also be cell periods
template <typename VEC>
VEC GenerateFactors(int n)
{
    VEC factors;
    //factors.push_back(1);
    //factors.push_back(n);
    for(int i = 2; i * i <= n; ++i) {
        if(n % i == 0) {
            factors.push_back(i);
            if(i * i != n) {
                factors.push_back(n / i);
            }
        }
    }
    
    std::sort(factors.begin(), factors.end());
    return factors;
}


/// Assign data statement for a parallel implementation
template<class DATASET,class RTYPE>
void parallelAssignData(DATASET& data, const int i, const ivec_type& id, const vec_type& x0,
                        const std::vector<RTYPE>& steps, const std::vector<double>& wn);
                        
///Fill map iterates of a starting point
template<typename MAP>
void iterate(const nvis::vec2& x0, std::vector<nvis::vec2>& steps,
             const MAP& pmap, int niter);
/// Compute and Return the map displacement vector paramInput
template<typename MAP>
nvis::vec2 mapDisplacement(const nvis::vec2& x0, const MAP& pmap, int period,
                           const map_analysis_param& params);
/// Compute and Return the map displacement vector metricInput
template<typename MAP>
nvis::vec2 mapDisplacement(const nvis::vec2& x0, const MAP& pmap, int period,
                           const metric_type& _metric);
                           
/** NOT USED!!! Processing a cell's periods and Poincare indices for possible periods.
    */
template<typename MAP, class DATA, class GRID, typename IVEC>
void process_cell(std::vector<std::pair<int, int> >& pids,
                  const DATA& dataset,
                  const GRID& grid,
                  const IVEC& cell_id, const MAP& pmap, double dx,
                  map_analysis_param& params);
                  
/** Sampling from a Poincare map given a mesh, map parameters,
    and a Tracking system (for winding number or sf)*/
template<typename MAP, typename Tracker>
void sample_raster(dataset_type& dataset,
                   const grid_type& grid,
                   const MAP& pmap, map_analysis_param& params,
                   bool compute_sf = true);
                   
                   
                   
/** Sampling the map for seeds needed by the grid.  This maps points and
 *  assigns the returns and winding numbers to points within cells as
 *  part of the AdaptiveGrid class (in GRID)
 */
template<typename MAP, typename TRACKER, typename GRID, typename DATASET, typename IDTYPE, typename WINDING_VEC>
void sampleMapForAdaptiveGrid(GRID& grid, DATASET& dataSet,
                              const std::vector<IDTYPE>& seeds,
                              const MAP& pmap, const int& maxp,
                              const int& depth, bool subres = false);
                              
/** Adaptive Sampling from a Poincare map given a mesh, map parameters, convexity checking functor
    and a Tracking system (Uses a quad-tree-esque structure)*/
template<typename MAP, typename Tracker, typename GRID, typename DATASET, typename WINDING_VEC, typename CONVEXITY_FUNCTOR>
void adaptiveSampleRaster(GRID& grid, DATASET& dataSet,
                          CONVEXITY_FUNCTOR& convexFunc,
                          const MAP& pmap, map_analysis_param& params,
                          bool compute_sf = true);
                          
/** Get Data from the DATASET object (AdaptiveGridNodeData<>) for the adaptive grid
*/
template<typename MAP, typename GRID, typename DATASET, typename IDTYPE, typename WINDING_VEC>
void getDataForAdaptiveGrid(GRID& grid, DATASET& dataSet, const std::vector<IDTYPE>& seeds,
                            const MAP& pmap, const int& maxp, const int& depth, bool subres = false);
                            
                            
/** Adaptive Sampling from a Poincare map given a data file with return information, map parameters, convexity checking functor
    (Uses a quad-tree-esque structure).  The given data file is assumed to have "orbit_data" per node.
*/
template<typename MAP, typename GRID, typename DATASET, typename WINDING_VEC, typename CONVEXITY_FUNCTOR>
void adaptiveRasterFromFile(const char* inputDataFile, GRID& grid, DATASET& dataSet,
                            typename GRID::bounds_type& bounds, nvis::ivec2& baseDim, int& maxDepth,
                            CONVEXITY_FUNCTOR& convexFunc,
                            const MAP& pmap, double& ham, map_analysis_param& params,
                            bool compute_sf = true);
                            
                            
//In future, move these following functions to a ModelFitterClass:
/** Computing a saddle guess given map displacement data (forward and backward)
 *  using nonlinear model fitting.  Option is for a linear(false) or quadratic
 *  model (true).  Set x0 as cell center or your best guess
 */
template<typename MapPosType>
MapPosType saddleGuessFromModelFit(const MapPosType& x0, const std::vector<MapPosType>& xVals,
                                   const std::vector<MapPosType>& delta, const std::vector<MapPosType>& beta,
                                   VectorXd& parameterVector,
                                   bool fitOption = false, bool verbose = false);
                                   
/** Computing a saddle guess given map displacement data (forward and backward)
 *  using nonlinear model fitting WITH MORE INFORMATION.  Option is for a
 *  linear(false) or quadratic model (true). Set x0 as cell center or your best guess
 */
template<typename MapPosType>
MapPosType saddleGuessFromModelFit(const MapPosType& x0, const std::vector<MapPosType>& xVals,
                                   const std::vector<MapPosType>& delta, const std::vector<MapPosType>& beta,
                                   VectorXd& parameterVector,
                                   int info, int nfev, int njev, double fnorm, double covfac,  //outputs from levMar
                                   bool fitOption = false, bool verbose = false);
                                   
/** Computing a saddle guess given map displacement data (forward and backward)
 *  using nonlinear model fitting.  This uses the normalized map tangent for
 *  the model.  Set x0 as cell center or your best guess
 */
template<typename MapPosType>
MapPosType saddleGuessFromMapTangentNorm(
    const MapPosType& x0, const std::vector<MapPosType>& xVals,
    const std::vector<MapPosType>& delta, const std::vector<MapPosType>& beta,
    VectorXd& parameterVector,
    bool verbose=false);
    
} // namespace xavier

//---------------------------------------------------------------------------------------
//Compiler needs access to the source implementations for templates
//---------------------------------------------------------------------------------------

/** A function to get the winding ratios at a given point in a dataset
 * (raster_data - the winding numbers are attached)
 */
template<class DATASET,typename IVEC>
void xavier::getWindingRatios(const DATASET& data, const IVEC& id, const int p_max,
                              std::vector< boost::rational<int> >& alpha)
{
    //load the winding number vector
    std::vector<double> wns = data.getVertexValue(id);
    //Compute the rational approximation
    alpha.clear();
    for(int i=0; i<(int)wns.size(); i++) {
        alpha.push_back( xavier::rational_approx_CF<int>(wns[i],p_max) );
    }
}

/// Parallel data assignment
template<class DATASET,class RTYPE>
void xavier::parallelAssignData(DATASET& data, const int i, const ivec_type& id, const vec_type& x0,
                                const std::vector<RTYPE>& steps, const std::vector<double>& wn)
{
    //Used in an "ordered" assignment statement (each thread does this in order)
    // i = the point index in the run (ignored)
    
    //Set the steps
    data(id).steps.push_back(x0);
    for (int ii=0 ; ii<(int)steps.size() ; ++ii) {
        data(id).steps.push_back(steps[ii].x);
    }
    //Set the winding number
    std::copy(wn.begin(),wn.end(),std::back_inserter(data(id).wn));
}

///Fill map iterates of a starting point
template<typename MAP>
void xavier::iterate(const nvis::vec2& x0, std::vector<nvis::vec2>& steps,
                     const MAP& pmap, int niter)
{
    MAP* amap = pmap.clone();
    std::vector<nvis::vec2> tmp;
    steps.clear();
    steps.push_back(x0);
    try {
        amap->map(x0, tmp, niter);
    } catch(...) {
        return;
    }
    std::copy(tmp.begin(), tmp.end(), std::back_inserter(steps));
}
/// Compute and Return the map displacement vector
template<typename MAP>
nvis::vec2 xavier::mapDisplacement(const nvis::vec2& x0, const MAP& pmap, int period,
                                   const metric_type& _metric)
{
    MAP* amap = pmap.clone();
    std::vector<nvis::vec2> steps;
    try {
        amap->map(x0, steps, period);
    } catch(MapUndefined err) {
        throw err;
    } catch(FailedIntegration<typename MAP::lvec_type, typename MAP::gvec_type> err) {
        throw err;
    } catch(MapUnderflow err) {
        throw err;
    } catch(...) {
        MapUnknownError err("invalid result of map iteration");
        err.where = x0;
        throw err;
    }
    if ((int)steps.size() < std::abs(period)) {
        MapUnderflow err("invalid result of map iteration");
        err.where = x0;
        throw err;
    } else {
        nvis::vec2 v = _metric.displacement(x0, steps[std::abs(period)-1]);
        return v;
    }
}
/// Compute and Return the map displacement vector
template<typename MAP>
nvis::vec2 xavier::mapDisplacement(const nvis::vec2& x0, const MAP& pmap, int period,
                                   const map_analysis_param& params)
{
    MAP* amap = pmap.clone();
    std::vector<nvis::vec2> steps;
    try {
        amap->map(x0, steps, period);
    } catch(MapUndefined err) {
        throw err;
    } catch(FailedIntegration<typename MAP::lvec_type, typename MAP::gvec_type> err) {
        throw err;
    } catch(MapUnderflow err) {
        throw err;
    } catch(...) {
        MapUnknownError err("invalid result of map iteration");
        err.where = x0;
        throw err;
    }
    if ((int)steps.size() < std::abs(period)) {
        MapUnderflow err("invalid result of map iteration");
        err.where = x0;
        throw err;
    } else {
        nvis::vec2 v = params.the_metric.displacement(x0, steps[std::abs(period)-1]);
        ///Note: This recording step only works in the old StdMap format!
        /*if (params.record) {
            map_analysis_param::tagged_vector_type tv(v, period);
            params.vectors.push_back(std::pair<nvis::vec2,
                                     map_analysis_param::tagged_vector_type>(x0, tv));
        }*/
        return v;
    }
}


/// Processing a cell for a period-index pairing (StandardMap method - NOT USED IN CR3BP!!!)
template<typename MAP,class DATA, class GRID, typename IVEC>
void xavier::process_cell(std::vector<std::pair<int, int> >& pids,
                          const DATA& dataset,
                          const GRID& mesh,
                          const IVEC& cell_id, const MAP& pmap, double dx,
                          map_analysis_param& params)
{
    const vec_type& step                 = mesh.spacing();
    const bounds_type& bounds         = mesh.bounds();
    const ivec_type& resolution = mesh.dimensions();
    const metric_type& _metric         = params.the_metric;
    
    const std::pair<int, int> default_result(0, 0);
    
    const int ids[][2] = { {0,0}, {1,0}, {1,1}, {0,1} };
    
    pids.clear();
    
    //Compute list of possible periods
    std::vector<int> periods;
    period_range<DATA,IVEC>(periods, dataset, cell_id, params.valid_rationals);
    std::set<int> available_periods(periods.begin(), periods.end());
    
    std::cerr << "available periods: " << std::endl;
    std::copy(available_periods.begin(), available_periods.end(), std::ostream_iterator<int>(std::cerr, ", "));
    std::cerr << std::endl;
    
    //For each available period, evaluate Poincare index
    for (std::set<int>::const_iterator it=available_periods.begin();
            it!=available_periods.end() ; ++it) {
        if (params.verbose) {
            std::cerr << "considering period " << *it << std::endl;
        }
        
        // compute Poincare index of the map for current period
        int p = *it;
        if (p == 0) {
            return;
        }
        int pidx = poincare_index(mesh, cell_id, dataset, pmap, p, dx, params);
        if (pidx != 0) {
            pids.push_back(std::pair<int, int>(pidx, p));
        }
    }
} //NOT USED in CR3BP or current versions....

/// Sampling from a Poincare map given a mesh, map parameters, and a Tracking system (for winding number or sf)
template<typename MAP, typename TRACKER>
void xavier::sample_raster(dataset_type& dataset,
                           const grid_type& mesh,
                           const MAP& pmap,
                           map_analysis_param& params,
                           bool compute_sf)
{
    typedef typename MAP::return_type return_type;
    
    const vec_type& step         = mesh.spacing();
    const bounds_type& bounds         = mesh.bounds();
    const ivec_type& resolution = mesh.dimensions();
    
    int niter = params.nb_iterations;
    int npoints = resolution[0]*resolution[1];
    
    std::cerr << npoints << " points to sample\n";
    
    //size_t nbthreads = 1;
//#if _OPENMP
//    nbthreads = omp_get_max_threads();
//#endif

    size_t counter = 0;
    
    #pragma omp parallel
    {
        #pragma omp for ordered schedule(dynamic,1)
        for(int n = 0 ; n < npoints ; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            int j = n / resolution[0];
            int i = n % resolution[0];
            ivec_type id(i,j);
            
            vec_type x0 = bounds.min() + step * vec_type(i,j);
            //Checking the node
            //std::ostringstream osNode;
            //osNode << "Node " << id << " has x0 = " << x0 << " and is " << (dataset.inGrid(id) ? "in the grid" : " OUTSIDE GRID") << "\n";
            //std::cerr << osNode.str();
            
            //Create a map engine and a tracker per thread
            MAP* amap = pmap.clone();
            TRACKER trackingSystem(vec_type(0,0));
            // Computing Winding Number (or safety-factor) by tracking rotations about origin
            std::vector<return_type> tmp;
            try {
                amap->template map_and_track_complete<TRACKER>(x0, tmp, niter, trackingSystem);
            } catch(std::runtime_error& e) {
                //Map is not defined for the whole niters ->maybe add a flag?
                //                std::ostringstream os;
                //                os << "Map: exception caught: " << e.what() << "\n";
                //                os << "   Node " << id << " at " << x0 << " is invalid\n";
                //                std::cerr << os.str();
                continue;
            } catch(std::exception& e) {
                //                std::ostringstream os;
                //                os << "Map: exception caught: " << e.what () << "\n";
                //                os << "   Node " << id << " at " << x0 << " is invalid\n";
                //                std::cerr << os.str();
                continue;
            }
            #pragma omp atomic
            ++counter;
            
            if ((int)tmp.size() < niter) {
                std::ostringstream os;
                os << '\n' << "only " << tmp.size() << " vertices / " << niter << " computed\n";
                std::cerr << os.str();
                //continue;
            }
            
            //Winding number
            std::vector<double> wn;
            if (compute_sf) {
                wn = trackingSystem.getWindingFactors(tmp.size(),tmp.back().delta_theta);
                //Old Format
                //fabs(tmp.back().delta_theta[0]/(2.*M_PI*(tmp.size()-1)));
            }
            
            //Use an ordered statement to build dataset values
            //#pragma omp ordered
            //parallelAssignData<dataset_type,return_type>(dataset,n,id,x0,tmp,wn);
            //Store data to dataset
            dataset(id).steps.push_back(x0);
            for (int i=0 ; i<(int)tmp.size() ; ++i) {
                dataset(id).steps.push_back(tmp[i].x);
            }
            if (compute_sf) {
                // wn = q/p, whereby q= # poloidal rotations, p= # toroidal rotations
                dataset(id).wn = trackingSystem.getWindingFactors(tmp.size(),tmp.back().delta_theta);
                
                //Old Format
                //fabs(tmp.back().delta_theta[0]/(2.*M_PI*(tmp.size()-1)));
                //Call a function embedded in Trakcer instead
                //dataset(id).wn = tracker.getWindingNumbers(tmp.back(),tmp.size()-1);
            }
            
            //Lets double check it's being output correctly
            //if (true) {
            //std::ostringstream os;
            //os << " Node " << id << ": tmp.size() = " << tmp.size() << " Last element = " << tmp.back().x << "\n";
            //os << " DataSet " << id << ": valid() = " << ( (dataset(id).valid())? "true":"false") << " with "
            //   << " Last element = " << dataset(id).steps.back() << "\n";
            //std::cout << os.str();
            //}
            
            //Screen Output (forced)
            if (true) {
                std::ostringstream os;
                os << "\rcompleted " << counter << " / " << npoints << " ("
                   << 100.*(float)counter/(float)npoints << "%), "
                   << counter << " valid orbits             \r" << std::flush;
                std::cout << os.str();
            }
            
            //Avizo Output if applicable
#ifdef HX_HAS_STD
            if (thread_id == 0) {
                //Update progress info
                float progress = (float)counter/(float)npoints;
                QString infoText = "Initial Sampling at ";
                QString progText;
                progText.setNum(100.*progress);
                infoText += progText + "\% complete";
                theWorkArea->setProgressInfo(infoText);
                //Update progress slider
                theWorkArea->setProgressValue(progress);
            }
            
#endif
        }//End parallel loop
    }//End parallel statement
    
    std::cout << '\n';
}


/// Adaptive Sampling via a QuadTree structure with a priority_queue given a mesh, map parameters, and a Tracking system
template<typename MAP, typename TRACKER, typename GRID, typename DATASET, typename WINDING_VEC, typename CONVEXITY_FUNCTOR>
void xavier::adaptiveSampleRaster(
    GRID& theGrid,                                 //AdaptiveGrid structure - stores grid information
    DATASET& theData,                              //AdaptiveGridNodeData - stores "orbit_data" per node
    CONVEXITY_FUNCTOR& convexFunc,                 //The convexity functor determining if we should subdivide a cell
    const MAP& theMap,                             //The map
    map_analysis_param& params,
    bool compute_sf)
{
    typedef typename MAP::return_type               return_type;
    typedef typename GRID::Node                     node_type;
    typedef typename GRID::id_type                  id_type;
    typedef typename std::vector<id_type>::iterator CellVectorIterator;
    typedef CompareCellConvexity<CONVEXITY_FUNCTOR> ConvexCompare;
    typedef std::priority_queue<id_type, std::vector<id_type>, ConvexCompare> ConvexPQ;
    
    //const vec_type& step         = mesh.spacing();
    //const bounds_type& bounds    = mesh.bounds();
    //const ivec_type& resolution  = mesh.dimensions();
    //Adaptive Grid definitions
    //AdaptiveGrid< std::vector<double>, DefalutValueTraits > theGrid(mesh.bounds(),resolution[0],resolution[1]);
    
    std::vector<id_type> seed_ids;
    
    int niter = params.nb_iterations;
    
    
    size_t nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    
    //---------------------------------------------------------------------------------------
    //Fill Initial grid information by sampling nodes (d=0) and subpoints
    int depth = 0;
    int dIdx = 2;
    seed_ids.clear();
    theGrid.getUndefinedVertexes(seed_ids);
    int npoints = (int) seed_ids.size();
    std::cerr << " Initial Grid Sampling : " << npoints << " corners to sample\n";
#ifdef HX_HAS_STD
    QString updateStr(" Initial Grid Sampling "),sStr;
    sStr.setNum(npoints);
    updateStr += " : " + sStr + " corners to sample";
    theWorkArea->startWorking(updateStr);
#endif
    //Sample the given seeds at d=0 (runs parallel)
    sampleMapForAdaptiveGrid<MAP,TRACKER,GRID,DATASET,id_type,WINDING_VEC>
    (theGrid,theData,seed_ids,theMap,niter,-1,false);
#ifdef HX_HAS_STD
    QString fillStr = " Initial Grid Sampling : Filling empty cells with info ...";
    theWorkArea->setProgressInfo(fillStr);
#endif
    std::vector<id_type> leafIDs;
    theGrid.getLeaves(leafIDs);
    //Determine Empty Leaves
    std::vector<id_type> emptyCellIDs;
    std::vector<bool>    needNewSubPoints;
    for (int i=0; i<(int)leafIDs.size(); ++i) {
        //Indicate "To Sample" if the grid cell is empty
        if ( theGrid.getCellData(leafIDs[i]).empty()) {
            emptyCellIDs.push_back(leafIDs[i]);
            //Sample extra if it needs more points within
            needNewSubPoints.push_back( (depth <=1)? true : false );
        }
        // Note: Cells with only one or two returns need more information!
        //       Especially at the first couple depth levels.
        else if ( (depth < 1) ) {
            //Leaf needs additional sampling
            emptyCellIDs.push_back(leafIDs[i]);
            needNewSubPoints.push_back(true);
        }
        
    }
    
    //Create new seeds for the empty cells at center (and subcell centers - 5 total)
    std::vector<id_type> new_seeds(emptyCellIDs.size());
    std::vector<id_type> new_seed_down2Centers;
    for (int i=0; i<(int)emptyCellIDs.size(); i++) {
        //Center for d+1
        new_seeds[i] = theGrid.sub_id(emptyCellIDs[i]) + id_type(1,1,0);
        //Also need centers for d+2 cells For top level
        if (needNewSubPoints[i]) {
            // Only for top level -> Better resolution within cell (stores for later)
            id_type cellIDd2 = theGrid.sub_id( theGrid.sub_id(emptyCellIDs[i]) );
            new_seed_down2Centers.push_back( cellIDd2 + id_type(1,1,0) );
            new_seed_down2Centers.push_back( cellIDd2 + id_type(3,1,0) );
            new_seed_down2Centers.push_back( cellIDd2 + id_type(1,3,0) );
            new_seed_down2Centers.push_back( cellIDd2 + id_type(3,3,0) );
        }
    }
    std::copy(new_seed_down2Centers.begin(),new_seed_down2Centers.end(),
              std::back_inserter(new_seeds));
              
    //Compute map/winding numbers for new seeds
#ifdef HX_HAS_STD
    QString dZeroStr = " Initial Grid Sampling : Filling empty cells with information...";
    theWorkArea->setProgressInfo(dZeroStr);
#endif
    sampleMapForAdaptiveGrid<MAP,TRACKER,GRID,DATASET,id_type,WINDING_VEC>
    (theGrid,theData,new_seeds,theMap,niter,-1,true);
    
    //---------------------------------------------------------------------------------------
    
    
    //Utilize a priority queue to work through all cells for subdivision
    ConvexCompare theComparer(convexFunc);
    ConvexPQ  pQueue(theComparer);
    theGrid.getLeaves(leafIDs);
    std::vector<id_type> queueSeedIDs;
    //Fill priority_queue with ALL leaves
    for(int i=0; i<(int)leafIDs.size(); i++) {
        pQueue.push( leafIDs[i] );
        queueSeedIDs.push_back(leafIDs[i]);
    }
    //Debug
    /*std::cout << " Priority Queue Setup :\n";
    std::cout << "   Top = " << pQueue.top() << "\n";
    std::cout << "   Data->valid() = " << theData.valid(pQueue.top()) << "\n";
    std::cout << "   Grid->valid() = " << theGrid.valid(pQueue.top()) << "\n";
    std::cout << "   Grid->isCellValid() = " << theGrid.isCellValid(pQueue.top()) << "\n";
    std::cout << "   Grid->allCornersInvalid() = " << theGrid.allCornersInvalid(pQueue.top()) << "\n";
    std::cout << "   Grid->allCornersValid() = " << theGrid.allCornersValid(pQueue.top()) << "\n";
    std::cout << "   Convexity = " << convexFunc(pQueue.top()) << "\n";*/
    
    
    // Operate On Priority Queue until all elements are gone OR cells are convex
    //---------------------------------------------------------------------------------------
    int queueIter = 0;
    while( !pQueue.empty() && !convexFunc(pQueue.top()) ) {
        //Grab the first 64 for parallel processing
        std::vector<id_type> topCells;
        CellVectorIterator   cit;
        int numCells = std::min(64,(int)pQueue.size()); //Up to 64 cells for processing
        for(int count=0; count<numCells; count++) {
            topCells.push_back( pQueue.top() );
            pQueue.pop();
        }
        
        //Determine which cells must be refined (parallel)
        std::vector<id_type>* refineIDsPerThread = new std::vector<id_type>[nbthreads];
        int nLeaves = (int) topCells.size();
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic,1)
            for (int n=0; n<nLeaves; ++n) {
                int th_id = 0;
#ifdef _OPENMP
                th_id = omp_get_thread_num();
#endif
                //Make a clone of the convexity checker
                CONVEXITY_FUNCTOR* cfunc = convexFunc.clone();
                bool convexityCheck = (*cfunc)(topCells[n]);
                //Refine only non-convex cells that are valid or partially valid
                if ( !convexityCheck && !(theGrid.allCornersInvalid(topCells[n])) ) {
                    //Cell fail convexity check, need to refine
                    refineIDsPerThread[th_id].push_back( topCells[n] );
                }
            }
            
        }
        std::vector<id_type> toRefine;
        for (int k=0; k<(int)nbthreads; ++k) {
            std::copy(refineIDsPerThread[k].begin(),
                      refineIDsPerThread[k].end(),
                      std::back_inserter(toRefine) );
        }
        // Refine identified top-of-queue cells
        for(int i=0; i<(int)toRefine.size(); i++) {
            theGrid.refine(toRefine[i]);
        }
        std::cout << " Queue Batch " << queueIter << " : "
                  << toRefine.size() << " cells have been refined.\n";
#ifdef HX_HAS_STD
        QString updateStr2(" Queue Batch "),dStr;
        dStr.setNum(queueIter);
        sStr.setNum(npoints);
        updateStr2 += dStr + " : " + sStr + " node points to sample";
        theWorkArea->setProgressInfo(updateStr2);
#endif
        //Get the undefined vertexes with no data
        seed_ids.clear();
        theGrid.getUndefinedVertexes(seed_ids);
        int npoints = (int) seed_ids.size();
        std::cout << " Queue Batch " << queueIter << " : " << npoints << " nodes to sample\n";
        //Sample all seeds for grid & insert into data grid
        sampleMapForAdaptiveGrid<MAP,TRACKER,GRID,DATASET,id_type,WINDING_VEC>
        (theGrid,theData,seed_ids,theMap,niter,queueIter,false);
        //Locate empty cells/sample/insert
        leafIDs.clear();
        theGrid.getLeaves(leafIDs);
        emptyCellIDs.clear();
        needNewSubPoints.clear();
        for (int i=0; i<(int)leafIDs.size(); ++i) {
            depth = leafIDs[i][ dIdx ];
            //Indicate "To Sample" if the grid cell is empty
            if ( theGrid.getCellData(leafIDs[i]).empty()) {
                emptyCellIDs.push_back(leafIDs[i]);
                //Sample extra if it needs more points within
                needNewSubPoints.push_back( (depth <=1)? true : false );
            }
            // Note: Cells with only one or two returns need more information!
            //       Especially at the first couple depth levels.
            else if ( (depth <= 1) &&
                      ((int)theGrid.getCellData(leafIDs[i]).size() < 3) ) {
                //Leaf needs additional sampling
                emptyCellIDs.push_back(leafIDs[i]);
                needNewSubPoints.push_back(true);
            }
        }
        std::vector<id_type> newSeeds(emptyCellIDs.size());
        new_seed_down2Centers.clear();
        for(int i=0; i<(int)emptyCellIDs.size(); i++) {
            //Center for d+1
            newSeeds[i] = theGrid.sub_id(emptyCellIDs[i]) + id_type(1,1,0);
            //Also centers for d+2 in appropriate levels
            if (needNewSubPoints[i]) {
                //Only for top levels ->better resolution within cells
                id_type cellIDd2 = theGrid.sub_id( theGrid.sub_id(emptyCellIDs[i]) );
                new_seed_down2Centers.push_back( cellIDd2 + id_type(1,1,0) );
                new_seed_down2Centers.push_back( cellIDd2 + id_type(3,1,0) );
                new_seed_down2Centers.push_back( cellIDd2 + id_type(1,3,0) );
                new_seed_down2Centers.push_back( cellIDd2 + id_type(3,3,0) );
            }
        }
        std::copy(new_seed_down2Centers.begin(),new_seed_down2Centers.end(),
                  std::back_inserter(newSeeds));
                  
        //Compute map/winding numbers for new seeds
        npoints = (int) newSeeds.size();
        std::cout << " Queue Batch " << queueIter << " : " << npoints << " internal cell points to sample\n";
        sampleMapForAdaptiveGrid<MAP,TRACKER,GRID,DATASET,id_type,WINDING_VEC>
        (theGrid,theData,newSeeds,theMap,niter,queueIter,true);
        
        //Reassemble priority_queue with new information (convexity metric)
        while (!pQueue.empty()) {
            pQueue.pop();    //Clear the queue and reassemble
        }
        queueSeedIDs.clear();
        leafIDs.clear();
        theGrid.getLeaves(leafIDs);
        int numConvexLeaves = 0;
        for(int i=0; i<(int)leafIDs.size(); i++) {
            /*//Debug for (0,1,0) cell ->giving error
            if( nvis::all(leafIDs[i] == id_type(0,1,0)) ) {
             std::cout << " Queue Batch " << queueIter <<" : Debug testing (0,1,0)\n";
             std::cout << "   Data->valid() = " << theData.valid(leafIDs[i]) << "\n";
             std::cout << "   Grid->valid() = " << theGrid.valid(leafIDs[i]) << "\n";
             std::cout << "   Grid->isCellValid() = " << theGrid.isCellValid(leafIDs[i]) << "\n";
             std::cout << "   Grid->allCornersInvalid() = " << theGrid.allCornersInvalid(leafIDs[i]) << "\n";
             std::cout << "   Grid->allCornersValid() = " << theGrid.allCornersValid(leafIDs[i]) << "\n";
             std::cout << "   Convexity = " << convexFunc(leafIDs[i]) << "\n";
            }*/
            //Insert valid or partially valid cells that can be refined
            if ( (theGrid.isCellValid(leafIDs[i])) ||
                    !(theGrid.allCornersInvalid(leafIDs[i]))
               ) {
                //Only insert cells below max depth
                depth = leafIDs[i][dIdx];
                if (depth < params.max_depth) {
                    //Add leaf back into queue
                    pQueue.push(leafIDs[i]);
                    queueSeedIDs.push_back(leafIDs[i]);
                    if(convexFunc(leafIDs[i])) {
                        numConvexLeaves++;
                    }
                }
            }
        }
        
        //Debug : Output info about the top of the queue
        std::cout << " Queue Batch " << queueIter << " [DEBUG]\n";
        std::cout << "-----------------------------------------------------------------\n";
        std::cout << " Top Cell : " << pQueue.top() << "\n";
        std::cout << "    Convexity : " << convexFunc(pQueue.top()) << "\n";
        std::cout << "    Validity : " << theGrid.isCellValid(pQueue.top()) << "\n";
        std::cout << "    AllCornersInvalid : " << theGrid.allCornersInvalid(pQueue.top()) << "\n";
        std::cout << "-----------------------------------------------------------------\n";
        
        //Output
        float progress = (float) numConvexLeaves / (float) pQueue.size();
        std::cout << " Queue Batch " << queueIter << " : "
                  << " Approximate progress is " << numConvexLeaves << "/" << pQueue.size()
                  << " leaves are convex. (" << progress*100.0 << "%)\n";
                  
        //Avizo Output if applicable
#ifdef HX_HAS_STD
        //Update progress info
        theMsg->printf("Adaptive Grid: Queue Batch %d refined %d cells",queueIter,toRefine.size());
        QString infoText = "Adaptive Grid:  Queue Batch ";
        QString iterText, completeText;
        iterText.setNum(queueIter);
        completeText.setNum(progress*100.0);
        infoText += iterText + " : " + completeText + "% complete";
        theWorkArea->setProgressInfo(infoText);
        theWorkArea->setProgressValue(progress);
        //Update depth interval (done with this depth)
        //theWorkArea->stopWorking(); //End interval current depth
#endif
        
        queueIter++;
    }//End while loop
    
    //At the end, check if there are any additional conditions on cells
    // with convexFunc.validityAfterRefinement() function
    leafIDs.clear();
    theGrid.getLeaves(leafIDs);
    for(int i=0; i<(int)leafIDs.size(); i++) {
        convexFunc.validityAfterRefinement(leafIDs[i]);
        //Note:  This function adjusts some validity settings within theGrid.
    }
    
#ifdef HX_HAS_STD
    theWorkArea->setProgressInfo(QString("Sampling Reached Max Depth"));
    theWorkArea->stopWorking();
#endif
    
    std::cout << '\n';
    
}

///Sample the map for the adaptive grid
template<typename MAP, typename TRACKER, typename GRID, typename DATASET, typename IDTYPE, typename WINDING_VEC>
void xavier::sampleMapForAdaptiveGrid(
    GRID& grid, DATASET& dataSet, const std::vector<IDTYPE>& seeds,
    const MAP& pmap, const int& maxp, const int& iter, bool subres)
{
    typedef typename MAP::return_type  return_type;
    typedef typename MAP::lvec_type    MapPos_type;
    typedef typename MAP::gvec_type    FullState;
    typedef IDTYPE                     id_type;
    typedef std::vector<MapPos_type>   curve_type;
    typedef OrbitWindingPair<WINDING_VEC, MapPos_type, IDTYPE> OrbitInfo;
    typedef typename GRID::traits_type                         ValueTraits;
    typedef typename GRID::data_type                           data_pair_type;
    typedef typename DATASET::data_type                        OrbitDataType;
    
    
    int nseeds = (int) seeds.size();
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::vector<OrbitInfo>* orbits = new std::vector<OrbitInfo>[nthreads];
    ValueTraits tempTrait;
    WINDING_VEC invalid_value = tempTrait.invalid;
    nvis::timer timer;
    int started=0, completed=0, success=0;
    int cmp, st, suc;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (int n=0 ; n<nseeds ; ++n) {
        
            int th_id = 0;
#ifdef _OPENMP
            th_id = omp_get_thread_num();
#endif
            
            orbits[th_id].push_back(OrbitInfo(seeds[n]));
            OrbitInfo&  the_info  = orbits[th_id].back();
            curve_type& the_curve = the_info.orbit.second;
            
#pragma openmp atomic
            ++started;
            
            st = started;
            cmp = completed;
            suc = success;
            
            //Create a map engine and a tracker per thread
            MAP* amap = pmap.clone();
            TRACKER trackingSystem( MapPos_type(0.0) );  //About origin
            std::vector<return_type> tmp;
            typename std::vector<return_type>::iterator tmpIT;
            MapPos_type seedMapPos = grid.getVertex(seeds[n]);
            WINDING_VEC wn; //Winding number is actually a std::vector<double>
            bool isICValid = true;
            
            //Run the conversion check to see if the initial condition is in valid region
            FullState x0;
            try {
                x0 = amap->section().unproject(seedMapPos);
            } catch (...) {
                //IC is invalid if error is thrown
                isICValid = false;
            }
            
            if (isICValid) {
                //Compute Winding Number by running the map for maxp iterates
                try {
                    amap->template map_and_track_complete< TRACKER >(
                        seedMapPos,tmp,maxp,trackingSystem);
                } catch(...) {
                    //Mark as invalid
                    wn = GRID::traits_type::invalid;
                }
                //Check if reached the right number of returns
                if((int)tmp.size() < maxp) {
                    wn = trackingSystem.getDefaultWindingFactors();
                } else {
                    //Compute the winding numbers
                    wn = trackingSystem.getWindingFactors((int)tmp.size(),tmp.back().delta_theta);
                }
            } else {
                //IC is invalid, so set accordingly
                //Mark as invalid
                wn = GRID::traits_type::invalid;
            }
            
            std::ostringstream os;
            double dt = timer.elapsed();
            os << "Iter: " << iter << "\r" << cmp << "/" << nseeds << " complete (" << (float)cmp*100/nseeds << "%) in "
               << dt << "s. (" << (float)cmp/dt << "Hz) / " << st << " orbits in progress / "
               << suc << " (" << 100.*(float)suc/(float)cmp << "%) \r";
            std::cout << os.str() << std::flush;
            
            //Debug : (0,1,0) is valid when it shouldn't be?
            /*if( nvis::all(seeds[n] == IDTYPE(0,1,0)) ) {
              std::cout << " Propagating IDTYPE = " << seeds[n] << ":\n";
              std::cout << "  Map Pos = " << seedMapPos << " \n";
              std::cout << "  isICValid = " << isICValid << " \n";
              std::cout << "  Max p = " << maxp << "\n";
              std::cout << "  Number returns achieved = " << tmp.size() << "\n";
              nvis::vec3 wSet(wn[0],wn[1],wn[2]);
              std::cout << "  Set Winding Number = " << wSet << "\n";
            }*/
            
            //Store to data cache (winding number(s) <-> iterates (mapPos) pair)
            the_info.orbit.first = wn;
            for(tmpIT=tmp.begin(); tmpIT!=tmp.end(); tmpIT++) {
                the_curve.push_back(tmpIT->x);
            }
            the_info.isICValid = isICValid;
            
            #pragma omp atomic
            --started;
            #pragma omp atomic
            ++completed;
            
            bool isValueOK = true;
            for(int i=0; i<(int)the_info.orbit.first.size(); i++)
                if (the_info.orbit.first[i] == invalid_value[i]) {
                    isValueOK = false;
                }
                
            if(isValueOK)
#pragma openmp atomic
                ++success;
                
            st = started;
            cmp = completed;
            suc = success;
            
            os.clear();
            os.str("");
            dt = timer.elapsed();
            if (th_id == 0) {
                /*os << "(" << th_id << ")" << "Iter " << iter << ": " << cmp << "/" << nseeds
                        << " orbits completed (" << (float)cmp*100/nseeds << "%) in "
                        << dt << "s. (" << (float)cmp/dt << "Hz) / " << st << " orbits in progress / "
                        << suc << " (" << 100.*(float)suc/(float)cmp << "%) successful integrations         \n";
                std::cout << os.str() << std::flush;*/
                //Avizo Output if applicable
#ifdef HX_HAS_STD
                //Update progress info
                float progress = (float)cmp/(float)nseeds;
                QString infoText = "Iter ";
                QString progText,iterText;
                iterText.setNum(iter);
                progText.setNum(100.*progress);
                infoText += iterText + " : " + ((subres)? "EmptyCell" : "") + " Sampling at " + progText + "\% complete";
                //theWorkArea->setProgressInfo(infoText);
                //Update progress slider
                //theWorkArea->setProgressValue(progress);
#endif
            }//End thread0 print
        }
    } //End parallel
    double elapsed = timer.elapsed();
    std::cout << "\nComputation took " << elapsed << " s. (" << (float)nseeds/elapsed << " Hz)\n";
    std::cout << success << " orbits (" << 100.*(float)success/(float)nseeds << "%) were successfully integrated\n";
    
    timer.restart();
    for (int n=0 ; n<nthreads ; ++n) {
        for (int i=0 ; i<(int)orbits[n].size() ; ++i) {
            const OrbitInfo&  info   = orbits[n][i];
            const id_type&    id     = info.seed_id;
            const WINDING_VEC& w     = info.orbit.first;
            const curve_type& curve  = info.orbit.second;
            OrbitDataType orbitData;
            
            //Set the Vertex Values - Winding number
            grid.setVertexValue(id, w);
            orbitData.wn = w;
            if (subres) {
                grid.insertCellValue(grid.getVertex(id), w, id[2]);
            }
            //Insert downstream points into adaptive grid with same w
            for (int k=0 ; k<(int)curve.size() ; ++k) {
                //Safe call to insert into grid
                try {
                    grid.insertCellValue(curve[k], w, id[2]);
                } catch(...) {}
                //Also assign points into dataSet data
                orbitData.steps.push_back(curve[k]);
            }
            //Set whether or not the IC is in the valid region
            orbitData.isICValid = info.isICValid;
            
            //Add the id<->orbit_data pair to dataSet
            dataSet.addData(id, orbitData);
        }
    }
    //Update grid validity based on new data
    grid.updateCellValidity();
    
    //Done
    elapsed = timer.elapsed();
    //std::cout << "Data structure update took " << elapsed << " s. (" << (float)nseeds/elapsed << " Hz)\n";
}


/// Adaptive Sampling via a QuadTree-esque structure from a file containing node data (as if it were run from the map)
/// Note:  Works like adaptiveSampleRaster() but with dataset calls instead of map calls.
///        This is not particularly fast; it may be wiser to load data and establish
///        the grid at the same time.
template<typename MAP, typename GRID, typename DATASET, typename WINDING_VEC, typename CONVEXITY_FUNCTOR>
void xavier::adaptiveRasterFromFile(
    const char* inputDataFile,              //Input data file generated from AdaptiveGridNodeData<>.writeDataToFile()
    GRID& theGrid,                          //AdaptiveGrid structure - stores grid information
    DATASET& theData,                       //AdaptiveGridNodeData - stores "orbit_data" per node
    typename GRID::bounds_type& bounds,     //Bounds of the initial grid - read from file
    nvis::ivec2& baseDim,                   //Base Dims (depth=0) for adaptive grid
    int& maxDepth,                          //Maximum depth
    CONVEXITY_FUNCTOR& convexFunc,          //The convexity functor determining if we should subdivide a cell
    const MAP& theMap,                      //The map
    double& ham,                            //True Hamiltonian value loaded by file
    map_analysis_param& params,
    bool compute_sf)
{
    typedef typename MAP::return_type             return_type;
    typedef typename GRID::Node                   node_type;
    typedef typename GRID::id_type                id_type;
    typedef typename GRID::bounds_type            bounds_type;
    typedef typename GRID::pos_type               pos_type;
    typedef typename std::vector<id_type>::iterator CellVectorIterator;
    typedef CompareCellConvexity<CONVEXITY_FUNCTOR> ConvexCompare;
    typedef std::priority_queue<id_type, std::vector<id_type>, ConvexCompare> ConvexPQ;
    
    //Build the data set by reading from file
    int niter = params.nb_iterations;
    theData.readData(inputDataFile, theGrid, bounds, baseDim, ham, params);
    
    std::vector<id_type> seed_ids;
    
    size_t nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    
#ifdef HX_HAS_STD
    //Subdivide process by depth levels
    //theWorkArea->subdivide( params.max_depth+1 );
#endif
    
    
    //---------------------------------------------------------------------------------------
    //Fill Initial grid information by sampling nodes (d=0) and subpoints
    int depth = 0;
    int dIdx = 2;
    seed_ids.clear();
    theGrid.getUndefinedVertexes(seed_ids);
    int npoints = (int) seed_ids.size();
    std::cerr << " Initial Grid Sampling : " << npoints << " corners to sample\n";
#ifdef HX_HAS_STD
    QString updateStr(" Initial Grid Sampling "),sStr;
    sStr.setNum(npoints);
    updateStr += " : " + sStr + " corners to sample";
    theWorkArea->startWorking(updateStr);
#endif
    //Get the data from grid (in place of 'Sampling')
    getDataForAdaptiveGrid<MAP,GRID,DATASET,id_type,WINDING_VEC>
    (theGrid,theData,seed_ids,theMap,niter,-1,false);
#ifdef HX_HAS_STD
    QString fillStr = " Initial Grid Sampling : Filling empty cells with info ...";
    theWorkArea->setProgressInfo(fillStr);
#endif
    std::vector<id_type> leafIDs;
    theGrid.getLeaves(leafIDs);
    //Determine Empty Leaves
    std::vector<id_type> emptyCellIDs;
    std::vector<bool>    needNewSubPoints;
    for (int i=0; i<(int)leafIDs.size(); ++i) {
        //Indicate "To Sample" if the grid cell is empty
        if ( theGrid.getCellData(leafIDs[i]).empty()) {
            emptyCellIDs.push_back(leafIDs[i]);
            //Sample extra if it needs more points within
            needNewSubPoints.push_back( (depth <=1)? true : false );
        }
        // Note: Cells with only one or two returns need more information!
        //       Especially at the first couple depth levels.
        else if ( (depth < 1) ) {
            //Leaf needs additional sampling
            emptyCellIDs.push_back(leafIDs[i]);
            needNewSubPoints.push_back(true);
        }
        
    }
    
    //Create new seeds for the empty cells at center (and subcell centers - 5 total)
    std::vector<id_type> new_seeds(emptyCellIDs.size());
    std::vector<id_type> new_seed_down2Centers;
    for (int i=0; i<(int)emptyCellIDs.size(); i++) {
        //Center for d+1
        new_seeds[i] = theGrid.sub_id(emptyCellIDs[i]) + id_type(1,1,0);
        //Also need centers for d+2 cells For top level
        if (needNewSubPoints[i]) {
            // Only for top level -> Better resolution within cell (stores for later)
            id_type cellIDd2 = theGrid.sub_id( theGrid.sub_id(emptyCellIDs[i]) );
            new_seed_down2Centers.push_back( cellIDd2 + id_type(1,1,0) );
            new_seed_down2Centers.push_back( cellIDd2 + id_type(3,1,0) );
            new_seed_down2Centers.push_back( cellIDd2 + id_type(1,3,0) );
            new_seed_down2Centers.push_back( cellIDd2 + id_type(3,3,0) );
        }
    }
    std::copy(new_seed_down2Centers.begin(),new_seed_down2Centers.end(),
              std::back_inserter(new_seeds));
              
    //Compute map/winding numbers for new seeds
#ifdef HX_HAS_STD
    QString dZeroStr = " Initial Grid Sampling : Filling empty cells with information...";
    theWorkArea->setProgressInfo(dZeroStr);
#endif
    getDataForAdaptiveGrid<MAP,GRID,DATASET,id_type,WINDING_VEC>
    (theGrid,theData,new_seeds,theMap,niter,-1,true);
    
    //---------------------------------------------------------------------------------------
    
    
    //Utilize a priority queue to work through all cells for subdivision
    ConvexCompare theComparer(convexFunc);
    ConvexPQ  pQueue(theComparer);
    theGrid.getLeaves(leafIDs);
    std::vector<id_type> queueSeedIDs;
    //Fill priority_queue with all leaves
    for(int i=0; i<(int)leafIDs.size(); i++) {
        pQueue.push( leafIDs[i] );
        queueSeedIDs.push_back(leafIDs[i]);
    }
    //Debug
    /*std::cout << " Priority Queue Setup :\n";
    std::cout << "   Top = " << pQueue.top() << "\n";
    std::cout << "   Data->valid() = " << theData.valid(pQueue.top()) << "\n";
    std::cout << "   Grid->valid() = " << theGrid.valid(pQueue.top()) << "\n";
    std::cout << "   Grid->isCellValid() = " << theGrid.isCellValid(pQueue.top()) << "\n";
    std::cout << "   Grid->allCornersInvalid() = " << theGrid.allCornersInvalid(pQueue.top()) << "\n";
    std::cout << "   Grid->allCornersValid() = " << theGrid.allCornersValid(pQueue.top()) << "\n";
    std::cout << "   Convexity = " << convexFunc(pQueue.top()) << "\n";*/
    
    // Operate On Priority Queue until all elements are gone OR cells are convex
    //---------------------------------------------------------------------------------------
    int queueIter = 0;
    while( !pQueue.empty() && !convexFunc(pQueue.top()) ) {
        //Grab the first 64 for parallel processing
        std::vector<id_type> topCells;
        CellVectorIterator   cit;
        int numCells = std::min(64,(int)pQueue.size()); //Up to 64 cells for processing
        for(int count=0; count<numCells; count++) {
            topCells.push_back( pQueue.top() );
            pQueue.pop();
        }
        
        //Determine which cells must be refined (parallel)
        std::vector<id_type>* refineIDsPerThread = new std::vector<id_type>[nbthreads];
        int nLeaves = (int) topCells.size();
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic,1)
            for (int n=0; n<nLeaves; ++n) {
                int th_id = 0;
#ifdef _OPENMP
                th_id = omp_get_thread_num();
#endif
                //Make a clone of the convexity checker
                CONVEXITY_FUNCTOR* cfunc = convexFunc.clone();
                bool convexityCheck = (*cfunc)(topCells[n]);
                //Refine only non-convex cells that are valid or partially valid
                if ( !convexityCheck && !(theGrid.allCornersInvalid(topCells[n])) ) {
                    //Cell fail convexity check, need to refine
                    refineIDsPerThread[th_id].push_back( topCells[n] );
                }
            }
            
        }
        std::vector<id_type> toRefine;
        for (int k=0; k<(int)nbthreads; ++k) {
            std::copy(refineIDsPerThread[k].begin(),
                      refineIDsPerThread[k].end(),
                      std::back_inserter(toRefine) );
        }
        // Refine identified top-of-queue cells
        for(int i=0; i<(int)toRefine.size(); i++) {
            theGrid.refine(toRefine[i]);
        }
        std::cout << " Queue Batch " << queueIter << " : "
                  << toRefine.size() << " cells have been refined.\n";
#ifdef HX_HAS_STD
        QString updateStr2(" Queue Batch "),dStr;
        dStr.setNum(queueIter);
        sStr.setNum(npoints);
        updateStr2 += dStr + " : " + sStr + " node points to sample";
        theWorkArea->setProgressInfo(updateStr2);
#endif
        //Get the undefined vertexes with no data
        seed_ids.clear();
        theGrid.getUndefinedVertexes(seed_ids);
        int npoints = (int) seed_ids.size();
        std::cout << " Queue Batch " << queueIter << " : " << npoints << " nodes to sample\n";
        //Get data for seeds & insert into data grid
        getDataForAdaptiveGrid<MAP,GRID,DATASET,id_type,WINDING_VEC>
        (theGrid,theData,seed_ids,theMap,niter,queueIter,false);
        //Locate empty cells/sample/insert
        leafIDs.clear();
        theGrid.getLeaves(leafIDs);
        emptyCellIDs.clear();
        needNewSubPoints.clear();
        for (int i=0; i<(int)leafIDs.size(); ++i) {
            depth = leafIDs[i][ dIdx ];
            //Indicate "To Sample" if the grid cell is empty
            if ( theGrid.getCellData(leafIDs[i]).empty()) {
                emptyCellIDs.push_back(leafIDs[i]);
                //Sample extra if it needs more points within
                needNewSubPoints.push_back( (depth <=1)? true : false );
            }
            // Note: Cells with only one or two returns need more information!
            //       Especially at the first couple depth levels.
            else if ( (depth <= 1) &&
                      ((int)theGrid.getCellData(leafIDs[i]).size() < 3) ) {
                //Leaf needs additional sampling
                emptyCellIDs.push_back(leafIDs[i]);
                needNewSubPoints.push_back(true);
            }
        }
        std::vector<id_type> newSeeds(emptyCellIDs.size());
        new_seed_down2Centers.clear();
        for(int i=0; i<(int)emptyCellIDs.size(); i++) {
            //Center for d+1
            newSeeds[i] = theGrid.sub_id(emptyCellIDs[i]) + id_type(1,1,0);
            //Also centers for d+2 in appropriate levels
            if (needNewSubPoints[i]) {
                //Only for top levels ->better resolution within cells
                id_type cellIDd2 = theGrid.sub_id( theGrid.sub_id(emptyCellIDs[i]) );
                new_seed_down2Centers.push_back( cellIDd2 + id_type(1,1,0) );
                new_seed_down2Centers.push_back( cellIDd2 + id_type(3,1,0) );
                new_seed_down2Centers.push_back( cellIDd2 + id_type(1,3,0) );
                new_seed_down2Centers.push_back( cellIDd2 + id_type(3,3,0) );
            }
        }
        std::copy(new_seed_down2Centers.begin(),new_seed_down2Centers.end(),
                  std::back_inserter(newSeeds));
                  
        //Compute map/winding numbers for new seeds
        npoints = (int) newSeeds.size();
        std::cout << " Queue Batch " << queueIter << " : " << npoints << " internal cell points to sample\n";
        getDataForAdaptiveGrid<MAP,GRID,DATASET,id_type,WINDING_VEC>
        (theGrid,theData,newSeeds,theMap,niter,queueIter,true);
        
        //Reassemble priority_queue with new information (convexity metric)
        while (!pQueue.empty()) {
            pQueue.pop();    //Clear the queue and reassemble
        }
        queueSeedIDs.clear();
        leafIDs.clear();
        theGrid.getLeaves(leafIDs);
        int numConvexLeaves = 0;
        for(int i=0; i<(int)leafIDs.size(); i++) {
            /*//Debug for (0,1,0) cell ->giving error
                if( nvis::all(leafIDs[i] == id_type(0,1,0)) ) {
                 std::cout << " Queue Batch " << queueIter <<" : Debug testing (0,1,0)\n";
                 std::cout << "   Data->valid() = " << theData.valid(leafIDs[i]) << "\n";
                 std::cout << "   Grid->valid() = " << theGrid.valid(leafIDs[i]) << "\n";
                 std::cout << "   Grid->isCellValid() = " << theGrid.isCellValid(leafIDs[i]) << "\n";
                 std::cout << "   Grid->allCornersInvalid() = " << theGrid.allCornersInvalid(leafIDs[i]) << "\n";
                 std::cout << "   Grid->allCornersValid() = " << theGrid.allCornersValid(leafIDs[i]) << "\n";
                 std::cout << "   Convexity = " << convexFunc(leafIDs[i]) << "\n";
                }*/
            //Insert valid or partially valid cells that can be refined
            if ( (theGrid.isCellValid(leafIDs[i])) ||
                    !(theGrid.allCornersInvalid(leafIDs[i]))
               ) {
                //Insert cells below max depth
                depth = leafIDs[i][dIdx];
                if (depth < params.max_depth) {
                    //Add leaf back into queue
                    pQueue.push(leafIDs[i]);
                    queueSeedIDs.push_back(leafIDs[i]);
                    if(convexFunc(leafIDs[i])) {
                        numConvexLeaves++;
                    }
                }
            }
        }
        
        //Debug : Output info about the top of the queue
        std::cout << " Queue Batch " << queueIter << " [DEBUG]\n";
        std::cout << "-----------------------------------------------------------------\n";
        std::cout << " Top Cell : " << pQueue.top() << "\n";
        std::cout << "    Convexity : " << convexFunc(pQueue.top()) << "\n";
        std::cout << "    Validity : " << theGrid.isCellValid(pQueue.top()) << "\n";
        std::cout << "    AllCornersInvalid : " << theGrid.allCornersInvalid(pQueue.top()) << "\n";
        std::cout << "-----------------------------------------------------------------\n";
        
        //Output
        float progress = (float) numConvexLeaves / (float) pQueue.size();
        std::cout << " Queue Batch " << queueIter << " : "
                  << " Approximate progress is " << numConvexLeaves << "/" << pQueue.size()
                  << " leaves are convex. (" << progress*100.0 << "%)\n";
                  
        //Avizo Output if applicable
#ifdef HX_HAS_STD
        //Update progress info
        theMsg->printf("Adaptive Grid: Queue Batch %d refined %d cells",queueIter,toRefine.size());
        QString infoText = "Adaptive Grid:  Queue Batch ";
        QString iterText, completeText;
        iterText.setNum(queueIter);
        completeText.setNum(progress*100.0);
        infoText += iterText + " : " + completeText + "% complete";
        theWorkArea->setProgressInfo(infoText);
        theWorkArea->setProgressValue(progress);
        //Update depth interval (done with this depth)
        //theWorkArea->stopWorking(); //End interval current depth
#endif
        
        queueIter++;
    }//End while loop
    
    
    //At the end, check if there are any additional conditions on cells
    // with convexFunc.validityAfterRefinement() function
    leafIDs.clear();
    theGrid.getLeaves(leafIDs);
    for(int i=0; i<(int)leafIDs.size(); i++) {
        convexFunc.validityAfterRefinement(leafIDs[i]);
        //Note:  This function adjusts some validity settings within theGrid.
    }
    
#ifdef HX_HAS_STD
    theWorkArea->setProgressInfo(QString("Sampling Reached Max Depth"));
    theWorkArea->stopWorking();
#endif
    
    
    std::cout << '\n';
    
}

///Get the map from DATASET (AdaptiveGridNodeData) for the adaptive grid
template<typename MAP, typename GRID, typename DATASET, typename IDTYPE, typename WINDING_VEC>
void xavier::getDataForAdaptiveGrid(GRID& grid, DATASET& dataSet, const std::vector<IDTYPE>& seeds,
                                    const MAP& pmap, const int& maxp, const int& depth, bool subres)
{
    typedef typename MAP::return_type         return_type;
    typedef typename MAP::lvec_type           MapPos_type;
    typedef IDTYPE                            id_type;
    typedef std::vector<MapPos_type>          curve_type;
    typedef OrbitWindingPair<WINDING_VEC, MapPos_type, IDTYPE> OrbitInfo;
    typedef typename GRID::traits_type        ValueTraits;
    typedef typename GRID::data_type          data_pair_type;
    typedef typename DATASET::data_type       OrbitDataType;
    
    
    int nseeds = (int)seeds.size();
    
    nvis::timer timer;
    timer.restart();
    //Just gather the data from the pre-computed map
    for (int n=0 ; n<nseeds ; ++n) {
        const id_type&    id      = seeds[n];
        OrbitDataType temp;
        temp.wn = ValueTraits::invalid;
        OrbitDataType& data = temp;
        WINDING_VEC& w      = temp.wn;
        curve_type& curve   = temp.steps;
        bool& icValid   = temp.isICValid;
        try {
            data   = dataSet(id);
            w      = data.wn;
            curve  = data.steps;
            icValid = data.isICValid;
        } catch (...) {
            //If error is thrown, it's likely invalid point anyway
            std::cerr << "AdaptiveGridNodeData Access Violation:: point " << id << " returned as invalid\n";
        }
        
        //Set the Vertex Values - Winding number
        grid.setVertexValue(id, w);
        if (subres) {
            grid.insertCellValue(grid.getVertex(id), w, id[2]);
        }
        //Insert downstream points into adaptive grid with same w
        for (int k=0 ; k<(int)curve.size() ; ++k) {
            //Safe call to insert into grid
            try {
                grid.insertCellValue(curve[k], w, id[2]);
            } catch(...) {}
            
        }
    }
    //Update grid validity based on new data
    grid.updateCellValidity();
    
    //double elapsed = timer.elapsed();
    //std::cout << "Data structure update took " << elapsed << " s. (" << (float)nseeds/elapsed << " Hz)\n";
}







//In future, move these functions to a class... may be better for nD case
/** Computing a saddle guess given map displacement data (forward and backward)
 *  using nonlinear model fitting.  Option is for a linear(false) or quadratic
 *  model (true).  Set x0 as cell center or your best guess
 */
template<typename MapPosType>
MapPosType xavier::saddleGuessFromModelFit(const MapPosType& x0, const std::vector<MapPosType>& xVals,
        const std::vector<MapPosType>& delta, const std::vector<MapPosType>& beta,
        VectorXd& params,
        bool fitOption, bool verbose)
{
    int numDataPoints = (int) xVals.size();
    int spaceDim = (int) x0.size();
    //Form tangent vectors
    std::vector<MapPosType> eta(numDataPoints);
    for(int i=0; i<numDataPoints; i++) {
        eta[i] = (delta[i]-beta[i])/2.0;
    }
    //for(int i=0;i<numDataPoints;i++) eta[i] = (delta[i]-beta[i]);
    
    int numParams = spaceDim*(spaceDim+1); //Linear
    if (fitOption) {
        //If Quadratic
        numParams += spaceDim*spaceDim*spaceDim; //More Parameters
        //Use the Linear result to seed the nonlinear model,
        // otherwise, it will give an erroneous result!
        //Quadratic Tensor is assumed to be zero.
        
        //WILL USE WHAT'S IN 'params'!
        
    } else {
        //Linear Model - Sets input parameters for first pass
        params.setZero(numParams);
        for(int i=0; i<spaceDim; i++) {
            params[i] = x0[i];    //The guess for location
        }
        
        //*** NOT GENERAL ***//
        //Uses Lyapunov's Theorem for Reciprical Pair eigenvalues
        // ----- Canonical Saddle -------------
        // Unstable <-0-> on x axis, Stable at ->0<- on y axis
        //params[spaceDim] = 2.0; //A(0,0)=A[0]
        //params[spaceDim + spaceDim + spaceDim-1] = -0.5; //A(1,1)=A[2+1]
        // ----- Cross 'x' Saddle -------------
        // Note: Suited for CR3BP by saddle behavior.
        // Puts Unstable mode(lam=2) with vu=[1,1] and Stable mode(lam=-2) with vs=[-1,1]
        params[spaceDim]            = 0; //A(0,0)=A[0]
        params[spaceDim+1]          = 2; //A(1,0)=A[1]
        params[spaceDim+spaceDim]   = 2; //A(0,1)=A[2]
        params[spaceDim+spaceDim+1] = 0; //A(1,1)=A[3]
        //Note: 2D on map space
    }
    
    
    MapPosType guessPoint(0,0);
    int info;
    double fnorm;
    
    if (fitOption) {
        //Use Quadratic Fit to tangent vectors
        typedef QuadModel2DFunctor FunctorType;
        //Model functor - set data values
        FunctorType quadModel(numParams,numDataPoints);
        for(int k=0; k<numDataPoints; k++) {
            VectorXd xPoint(spaceDim), yPoint(spaceDim);
            for(int i=0; i<spaceDim; i++) {
                xPoint[i] = xVals[k][i];
                yPoint[i] = eta[k][i];
            }
            quadModel.x.push_back(xPoint);
            quadModel.y.push_back(yPoint);
        }
        //Solution algorithm
        //LevenbergMarquardt< FunctorType > lm(quadModel); //Analytical Derivatives
        NumericalDiff< FunctorType > nDiffer(quadModel); //Numerical Derivatives
        LevenbergMarquardt< NumericalDiff<FunctorType> > lm(nDiffer);
        lm.parameters.ftol = 1.e-12;
        lm.parameters.xtol = 1.e-12;
        info = lm.minimize(params);
        fnorm = lm.fvec.blueNorm();
    } else {
        //Use Linear Fit
        typedef LinearModelFunctor FunctorType;
        //Model functor - set data values
        FunctorType linModel(spaceDim, numParams, numDataPoints);
        for(int k=0; k<numDataPoints; k++) {
            VectorXd xPoint(spaceDim), yPoint(spaceDim);
            for(int i=0; i<spaceDim; i++) {
                xPoint[i] = xVals[k][i];
                yPoint[i] = eta[k][i];
            }
            linModel.x.push_back(xPoint);
            linModel.y.push_back(yPoint);
        }
        //Solution algorithm
        //LevenbergMarquardt< FunctorType > lm(linModel); //Analytical Derivatives
        NumericalDiff< FunctorType > nDiffer(linModel); //Numerical Derivatives
        LevenbergMarquardt< NumericalDiff<FunctorType> > lm(nDiffer);
        lm.parameters.ftol = 1.e-12;
        lm.parameters.xtol = 1.e-12;
        info = lm.minimize(params);
        fnorm = lm.fvec.blueNorm();
    }
    if (verbose) std::cerr << "LevenbergMarquardt Soln: " << ((fitOption)? "Quadratic" : "Linear")
                               << " Info = " << info << " Fnorm = " << fnorm << "\n";
    for(int i=0; i<spaceDim; i++) {
        guessPoint[i] = params[i];
    }
    
    return guessPoint;
}

/** Computing a saddle guess given map displacement data (forward and backward)
 *  using nonlinear model fitting WITH MORE OUTPUT INFORMATION.  Option is for a
 *  linear(false) or quadratic model (true). Set x0 as cell center or your best guess
 */
template<typename MapPosType>
MapPosType xavier::saddleGuessFromModelFit(const MapPosType& x0, const std::vector<MapPosType>& xVals,
        const std::vector<MapPosType>& delta, const std::vector<MapPosType>& beta,
        VectorXd& params,
        int info, int nfev, int njev, double fnorm, double covfac,  //outputs from levMar
        bool fitOption, bool verbose)
{
    int numDataPoints = xVals.size();
    int spaceDim = x0.size();
    //Form tangent vectors
    std::vector<MapPosType> eta(numDataPoints);
    for(int i=0; i<numDataPoints; i++) {
        eta[i] = (delta[i]-beta[i])/2.0;
    }
    //for(int i=0;i<numDataPoints;i++) eta[i] = (delta[i]-beta[i]);
    
    int numParams = spaceDim*(spaceDim+1); //Linear
    if (fitOption) {
        //If Quadratic
        numParams += spaceDim*spaceDim*spaceDim; //More Parameters
        //Use the Linear result to seed the nonlinear model,
        // otherwise, it will give an erroneous result!
        //Quadratic Tensor is assumed to be zero.
        
        //WILL USE WHAT'S IN 'params'!
        
    } else {
        //Linear Model - Sets input parameters for first pass
        params.setZero(numParams);
        for(int i=0; i<spaceDim; i++) {
            params[i] = x0[i];    //The guess for location
        }
        
        //*** NOT GENERAL ***//
        //Uses Lyapunov's Theorem for Reciprical Pair eigenvalues
        // ----- Canonical Saddle -------------
        // Unstable <-0-> on x axis, Stable at ->0<- on y axis
        //params[spaceDim] = 2.0; //A(0,0)=A[0]
        //params[spaceDim + spaceDim + spaceDim-1] = -0.5; //A(1,1)=A[2+1]
        // ----- Cross 'x' Saddle -------------
        // Note: Suited for CR3BP by saddle behavior.
        // Puts Unstable mode(lam=2) with vu=[1,1] and Stable mode(lam=-2) with vs=[-1,1]
        params[spaceDim]            = 0; //A(0,0)=A[0]
        params[spaceDim+1]          = 2; //A(1,0)=A[1]
        params[spaceDim+spaceDim]   = 2; //A(0,1)=A[2]
        params[spaceDim+spaceDim+1] = 0; //A(1,1)=A[3]
        //Note: 2D on map space
    }
    
    MapPosType guessPoint(0,0);
    
    
    if (fitOption) {
        //Use Quadratic Fit to tangent vectors
        typedef QuadModel2DFunctor FunctorType;
        //Model functor - set data values
        FunctorType quadModel(numParams,numDataPoints);
        for(int k=0; k<numDataPoints; k++) {
            VectorXd xPoint(spaceDim), yPoint(spaceDim);
            for(int i=0; i<spaceDim; i++) {
                xPoint[i] = xVals[k][i];
                yPoint[i] = eta[k][i];
            }
            quadModel.x.push_back(xPoint);
            quadModel.y.push_back(yPoint);
        }
        //Solution algorithm
        //LevenbergMarquardt< FunctorType > lm(quadModel);
        NumericalDiff< FunctorType > nDiffer(quadModel); //Numerical Derivatives
        LevenbergMarquardt< NumericalDiff<FunctorType> > lm(nDiffer);
        lm.parameters.ftol = 1.e-12;
        lm.parameters.xtol = 1.e-12;
        info = lm.minimize(params);
        fnorm = lm.fvec.blueNorm();
        nfev = lm.nfev;
        njev = lm.njev;
        covfac = fnorm*fnorm/(numDataPoints-numParams);
    } else {
        //Use Linear Fit
        typedef LinearModelFunctor FunctorType;
        //Model functor - set data values
        FunctorType linModel(spaceDim, numParams, numDataPoints);
        for(int k=0; k<numDataPoints; k++) {
            VectorXd xPoint(spaceDim), yPoint(spaceDim);
            for(int i=0; i<spaceDim; i++) {
                xPoint[i] = xVals[k][i];
                yPoint[i] = eta[k][i];
            }
            linModel.x.push_back(xPoint);
            linModel.y.push_back(yPoint);
        }
        //Solution algorithm
        //LevenbergMarquardt< FunctorType > lm(linModel); //Analytical
        NumericalDiff< FunctorType > nDiffer(linModel); //Numerical Derivatives
        LevenbergMarquardt< NumericalDiff<FunctorType> > lm(nDiffer);
        lm.parameters.ftol = 1.e-12;
        lm.parameters.xtol = 1.e-12;
        info = lm.minimize(params);
        fnorm = lm.fvec.blueNorm();
        nfev = lm.nfev;
        njev = lm.njev;
        covfac = fnorm*fnorm/(numDataPoints-numParams);
    }
    if (verbose) std::cerr << "LevenbergMarquardt Soln: " << ((fitOption)? "Quadratic" : "Linear")
                               << " Info = " << info << " Fnorm = " << fnorm << "\n";
    for(int i=0; i<spaceDim; i++) {
        guessPoint[i] = params[i];
    }
    
    return guessPoint;
    
}


/** Computing a saddle guess given map displacement data (forward and backward)
 *  using nonlinear model fitting.  This uses the normalized map tangent for
 *  the model.  Set x0 as cell center or your best guess
 */
template<typename MapPosType>
MapPosType xavier::saddleGuessFromMapTangentNorm(
    const MapPosType& x0, const std::vector<MapPosType>& xVals,
    const std::vector<MapPosType>& delta, const std::vector<MapPosType>& beta,
    VectorXd& params,
    bool verbose)
{
    int numDataPoints = (int) xVals.size();
    int spaceDim = (int) x0.size();
    
    //Form NORMALIZED map tangent vectors
    std::vector<MapPosType> eta(numDataPoints);
    for(int i=0; i<numDataPoints; i++) {
        MapPosType tmp = delta[i]-beta[i];
        eta[i] = tmp / nvis::norm(tmp);
    }
    
    int numParams = spaceDim*(spaceDim+1); //Linear
    //Linear Model - Sets input parameters for first pass
    params.setZero(numParams);
    for(int i=0; i<spaceDim; i++) {
        params[i] = x0[i];    //The guess for location
    }
    
    //*** NOT GENERAL ***//
    //Uses Lyapunov's Theorem for Reciprical Pair eigenvalues
    // ----- Canonical Saddle -------------
    // Unstable <-0-> on x axis, Stable at ->0<- on y axis
    //params[spaceDim] = 2.0; //A(0,0)=A[0]
    //params[spaceDim + spaceDim + spaceDim-1] = -0.5; //A(1,1)=A[2+1]
    // ----- Cross 'x' Saddle -------------
    // Note: Suited for CR3BP by saddle behavior.
    // Puts Unstable mode(lam=2) with vu=[1,1] and Stable mode(lam=-2) with vs=[-1,1]
    params[spaceDim]            = 0; //A(0,0)=A[0]
    params[spaceDim+1]          = 2; //A(1,0)=A[1]
    params[spaceDim+spaceDim]   = 2; //A(0,1)=A[2]
    params[spaceDim+spaceDim+1] = 0; //A(1,1)=A[3]
    //Note: 2D on map space
    
    
    MapPosType guessPoint(0,0);
    int info;
    double fnorm;
    
    //Use Linear Fit
    typedef LinearModelFunctor FunctorType;
    //Model functor - set data values
    FunctorType linModel(spaceDim, numParams, numDataPoints);
    for(int k=0; k<numDataPoints; k++) {
        VectorXd xPoint(spaceDim), yPoint(spaceDim);
        for(int i=0; i<spaceDim; i++) {
            xPoint[i] = xVals[k][i];
            yPoint[i] = eta[k][i];
        }
        linModel.x.push_back(xPoint);
        linModel.y.push_back(yPoint);
    }
    //Solution algorithm
    //LevenbergMarquardt< FunctorType > lm(linModel); //Analytical Derivatives
    NumericalDiff< FunctorType > nDiffer(linModel); //Numerical Derivatives
    LevenbergMarquardt< NumericalDiff<FunctorType> > lm(nDiffer);
    lm.parameters.ftol = 1.e-12;
    lm.parameters.xtol = 1.e-12;
    info = lm.minimize(params);
    fnorm = lm.fvec.blueNorm();
    if (verbose) std::cerr << "LevenbergMarquardt Soln: " << "Linear Fit for MapTangentNorm"
                               << " Info = " << info << " Fnorm = " << fnorm << "\n";
    for(int i=0; i<spaceDim; i++) {
        guessPoint[i] = params[i];
    }
    
    return guessPoint;
}

/// Get Period range by running rational approximations of the winding number
template<class DATA, typename IVEC>
void xavier::period_range(std::vector<int>& ps,
                          const DATA& dataset, const IVEC& cell_id, const int& maxp, const bool watchCell)
{

    std::vector<double> wns;
    std::vector<double>::const_iterator cit;
    
    //std::cout << " Calling dataset(" << cell_id << ") :\n";
    //std::cout << "     " << cell_id << " => inGrid " << dataset.inGrid(cell_id) << "\n";
    //std::cout << "     " << cell_id +IVEC(1,0) << " => inGrid " << dataset.inGrid(cell_id+IVEC(1,0)) << "\n";
    //std::cout << "     " << cell_id +IVEC(1,1) << " => inGrid " << dataset.inGrid(cell_id+IVEC(1,1)) << "\n";
    //std::cout << "     " << cell_id +IVEC(0,1) << " => inGrid " << dataset.inGrid(cell_id+IVEC(0,1)) << "\n";
    
    for (cit = dataset(cell_id).wn.begin();
            cit != dataset(cell_id).wn.end(); ++cit) {
        wns.push_back( *cit );
    }
    for (cit = dataset(cell_id + IVEC(1,0)).wn.begin();
            cit != dataset(cell_id + IVEC(1,0)).wn.end(); ++cit) {
        wns.push_back( *cit );
    }
    for (cit = dataset(cell_id + IVEC(1,1)).wn.begin();
            cit != dataset(cell_id + IVEC(1,1)).wn.end(); ++cit) {
        wns.push_back( *cit );
    }
    for (cit = dataset(cell_id + IVEC(0,1)).wn.begin();
            cit != dataset(cell_id + IVEC(0,1)).wn.end(); ++cit) {
        wns.push_back( *cit );
    }
    //Old format - scalar winding number per point
    //wns.push_back(dataset(cell_id).wn);
    //wns.push_back(dataset(cell_id + IVEC(1,0)).wn);
    //wns.push_back(dataset(cell_id + IVEC(1,1)).wn);
    //wns.push_back(dataset(cell_id + IVEC(0,1)).wn);
    
    // REMINDER:
    // winding number = # poloidal rotations (+ or -) / # toroidal rotations
    // period = # toroidal rotations necessary to complete an integral
    //          number of poloidal rotations
    // Remark: in a discrete map, toroidal rotations are discrete iterations
    //         of the map
    
    //Output for debugging
    std::ostringstream os;
    os << "Winding Numbers present in this cell are: ";
    std::vector<double>::iterator it;
    for (it = wns.begin(); it != wns.end(); ++it) {
        os << *it << " ";
    }
    os << "\n";
    if(watchCell) {
        std::cerr << os.str() << std::endl;
    }
    
    //Fill the list of periods by evaluating the integer ratio to pmax
    ps.clear();
    os.str("");
    os.clear();
    os << "Computed Fractions: ";
    boost::rational<int> alpha;
    
    for (it = wns.begin(); it!=wns.end(); ++it) {
        //compute rational number via continued fraction
        alpha = rational_approx_CF<int>(*it,maxp);
        os << alpha.numerator() << "/" << alpha.denominator() << " ";
        ps.push_back((int) fabs((double)alpha.denominator()) );
    }
    if (watchCell) {
        std::cerr << os.str() << std::endl;
    }
    os.str("");
    os.clear();
    
    //Gather the possible factors based on period set
    //std::vector<int> maybes = GenerateFactors<std::vector<int> >(ps.back());
    
    
    //Screen Output
    if (watchCell) {
        std::cout << "\rCell [" << cell_id << "] has periods: ";
        std::vector<int>::iterator pit;
        for(pit = ps.begin(); pit!=ps.end(); ++pit) {
            os << (*pit) << " ";
        }
        std::cerr << os.str() << "\n";
    }
    
}

/// Set the range of valid periods based on Winding numbers at each vertex of a cell (StandardMap)
template<class DATA, typename IVEC>
void xavier::period_range(std::vector<int>& ps,
                          const DATA& dataset,
                          const IVEC& cell_id,
                          const std::map<double, rational_type>& valid_rationals)
{

    std::vector<double> wns;
    std::vector<double>::const_iterator it;
    for (it = dataset(cell_id).wn.begin();
            it != dataset(cell_id).wn.end(); ++it) {
        wns.push_back( *it );
    }
    for (it = dataset(cell_id + IVEC(1,0)).wn.begin();
            it != dataset(cell_id + IVEC(1,0)).wn.end(); ++it) {
        wns.push_back( *it );
    }
    for (it = dataset(cell_id + IVEC(1,1)).wn.begin();
            it != dataset(cell_id + IVEC(1,1)).wn.end(); ++it) {
        wns.push_back( *it );
    }
    for (it = dataset(cell_id + IVEC(0,1)).wn.begin();
            it != dataset(cell_id + IVEC(0,1)).wn.end(); ++it) {
        wns.push_back( *it );
    }
    //Old format - scalar winding number per point
    //wns.push_back(dataset(cell_id).wn);
    //wns.push_back(dataset(cell_id + IVEC(1,0)).wn);
    //wns.push_back(dataset(cell_id + IVEC(1,1)).wn);
    //wns.push_back(dataset(cell_id + IVEC(0,1)).wn);
    
    // REMINDER:
    // winding number = # poloidal rotations (+ or -) / # toroidal rotations
    // period = # toroidal rotations necessary to complete an integral
    //          number of poloidal rotations
    // Remark: in a discrete map, toroidal rotations are discrete iterations
    //         of the map
    
    std::ostringstream os;
    os << "safety factors present in this cell are: ";
    std::copy(wns.begin(), wns.end(), std::ostream_iterator<double>(os, " "));
    os << "\n";
    std::cerr << os.str() << std::flush;
    
    double min = *std::min_element(wns.begin(), wns.end());
    double max = *std::max_element(wns.begin(), wns.end());
    
    typedef std::map<double, rational_type>         map_type;
    typedef map_type::const_iterator                iterator_type;
    
    // identify all valid rationals that are contained in the interval [min, max]
    iterator_type min_it = valid_rationals.lower_bound(min);
    iterator_type max_it = valid_rationals.upper_bound(max);
    // include boundaries on both end
    --min_it;
    ++max_it;
    ps.clear();
    // return corresponding denominators as potential periods
    for (iterator_type it=min_it ; it!=max_it ; ++it) {
        ps.push_back(it->second.denominator());
    }
    
    std::cout << "valid periods for this cell = ";
    std::copy(ps.begin(), ps.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << "\n";
}

/// Determine Cell's period range given that periods are set in the topology dataset
template<class DATA, class IVEC>
void xavier::period_range(std::vector<int>& ps, const DATA& dataset,
                          const IVEC& cell_id)
{
    typedef typename DATA::data_type OrbitData;
    ps.clear();
    std::vector<int> _inter, valid;
    IVEC ids[] = {0, IVEC(1,0), IVEC(0,1), IVEC(1,1)};
    std::ostringstream os;
    os << "period_range: ";
    for (int i=0 ; i<4; ++i) {
        const OrbitData& od  = dataset(cell_id + ids[i]);
        if (!od.valid()) {
            os << "[" << i << ": invalid]\n";
            std::cerr << os.str();
            return;
        }
        os << "[" << i << ": ";
        std::copy(od.periods.begin(), od.periods.end(), std::ostream_iterator<int>(os, " "));
        os << "], ";
        if (!i) {
            std::copy(od.periods.begin(), od.periods.end(), std::back_inserter(valid));
        } else {
            std::set_intersection(valid.begin(), valid.end(), od.periods.begin(), od.periods.end(), std::back_inserter(_inter));
            valid.swap(_inter);
        }
    }
    os << "\n";
    std::cerr << os.str();
    
    std::copy(valid.begin(), valid.end(), std::back_inserter(ps));
}


#endif
