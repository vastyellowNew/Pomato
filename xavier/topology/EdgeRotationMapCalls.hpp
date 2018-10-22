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


// EdgeRotation Poincare Map calls:
//  These are specially wrapped versions of map calls to handle a variety of exceptions
//  that occur while tracking the rotation of the "Map Displacement vector (P^{p}(x) - x)
//  Note:  This are also used during the computation of invariant manifolds as
//  section transveraslity violations occur frequently.
//Author : Wayne Schlei (Purdue University)

#ifndef __EDGE_ROTATION_MAPCALLS_HPP
#define __EDGE_ROTATION_MAPCALLS_HPP

#include <vector>
#include <map>
#include <maps/mapExceptions.hpp>
#include <topology/CATtracker.hpp>
#include <topology/SortableReturnData.hpp>
#include <topology/EdgeRotationFailure.hpp>
#include <topology/SectionTransversality.hpp>

using namespace xavier;

namespace topology {

/// Compute and Return the map displacement vector using a cache for data lookup
template<class MAP, class PARAM, class SRDATA>
typename MAP::lvec_type mapDisplacementUsingCache( const typename MAP::lvec_type& x0,
        const MAP& pmap, const int& period, const PARAM& params,
        std::set<SRDATA>& cache)
{
    typedef typename MAP::lvec_type         vec_type;
    typedef typename MAP::lmat_type         lmat_type;
    typedef typename MAP::state_type        state_type;
    typedef typename MAP::xstate_type       xstate_type;
    typedef typename MAP::rhs_type          RHStype;
    typedef typename SRDATA::MapDataType    MapDataType;
    typedef CATtracker<xstate_type,RHStype::numSingularities> TrackerType;
    static const int S = RHStype::numSingularities;
    
    //Cache variables
    typename std::set<SRDATA>::iterator cit;
    cit = cache.find( SRDATA(x0) );
    
    if ( cit != cache.end() && cit->isThere(period) ) {
        //If this x0 is in the cache, use it's returns:
        vec_type v = params.the_metric.displacement(x0, (*cit).returns[std::abs(period)-1]);
        return v;
    } else {
        //Point is not in cache, so let's compute it and add it to the cache
        std::vector<vec_type> steps;
        //Create a tracking system
        const xstate_type* singularityPtr = pmap.rhs().singularities();
        TrackerType tracker;
        for(int k=0; k<S; k++) {
            std::pair<vec_type,lmat_type> tmp = pmap.section().project(singularityPtr[k]);
            tracker.setSingularity(k,tmp.first);
        }
        try {
            //pmap.map(x0, steps, period);
            //Map with tracking close approach distances to singularities
            pmap.map_and_track_complete(x0, steps, period, tracker);
        } catch( MapSingularityError<state_type>& err) {
            throw err;
        } catch(MapUndefined& err) {
            throw err;
        } catch(FailedIntegration<typename MAP::lvec_type, typename MAP::gvec_type>& err) {
            throw err;
        } catch(MapUnderflow& err) {
            throw err;
        } catch(...) {
            throw std::runtime_error("mapDisplacementUsingCache: Unknown error!");
        }
        
        
        //Done with computation
        if((int)steps.size() < std::abs(period) ) {
            MapUndefined err("Unable to compute required steps in mapDisplacementUsingCache()");
            err.where = x0;
            throw err;
        } else {
            //Compute the displacement
            vec_type v = params.the_metric.displacement(x0, steps[std::abs(period)-1]);
            //Add this integration to cache
            std::vector<MapDataType> catData = tracker.getResult();
            cache.insert( SRDATA(x0,steps,catData) );
            return v;
        }
    }
}


/** Compute and Return the map displacement vector using a cache for data lookup
 *  => Specialized Edge-Rotation version that throws EdgeRotationFailure errors
 *  and checks NonTransverse Section conditions (i.e., Section Separation)
 */
template<class MAP, class PARAM, class SRDATA>
typename MAP::lvec_type mapDisplacementUsingCache(
    const typename MAP::lvec_type& x0, //Initial condition on map
    const MAP& pmap, const int& period, const PARAM& params,
    std::set<SRDATA>& cache,
    EdgeRotationFailure<typename MAP::lvec_type>& theFailure  //Edge Rotation Failure class
)
{
    typedef typename MAP::lvec_type         vec_type;
    typedef typename MAP::lmat_type         lmat_type;
    typedef typename MAP::state_type        state_type;
    typedef typename MAP::xstate_type       xstate_type;
    typedef typename MAP::rhs_type          RHStype;
    typedef typename SRDATA::MapDataType    MapDataType;
    typedef EdgeRotationFailure<vec_type>   Discont;
    typedef CATtracker<xstate_type,RHStype::numSingularities> TrackerType;
    static const int S = RHStype::numSingularities;
    
    //Temporary objects in case we encounter an error
    bool forward = (period < 0) ? false : true;
    theFailure.reset();
    theFailure.setMessage("Currently no failure");
    theFailure.failurePos = x0;
    theFailure.period = period;
    theFailure.setType( ((forward)? Discont::FORWARD_MAP : Discont::BACKWARD_MAP) );
    
    //Cache variables
    typename std::set<SRDATA>::iterator cit;
    cit = cache.find( SRDATA(x0) );
    
    if ( cit != cache.end() && cit->isThere(period) ) {
        //If this x0 is in the cache AND there is enough returns, use it's data:
        vec_type vOut = params.the_metric.displacement(x0, (*cit).returns[std::abs(period)-1]);
        return vOut;
    } else {
        //Point is not in cache, so let's compute it and add it to the cache
        std::vector<vec_type> steps;
        //Create a tracking system
        const xstate_type* singularityPtr = pmap.rhs().singularities();
        TrackerType tracker;
        for(int k=0; k<S; k++) {
            std::pair<vec_type,lmat_type> tmp = pmap.section().project(singularityPtr[k]);
            tracker.setSingularity(k,tmp.first);
        }
        bool isValid = true;
        try {
            //Original method maps the point
            //pmap.map(x0, steps, period);
            //Map with tracking close approach distances to singularities
            pmap.map_and_track_complete(x0, steps, period, tracker);
        } catch(MapSingularityError<state_type>& err) {
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            //We have to throw these just like Section Separation
            isValid = false;
        } catch(MapUndefined& err) {
            //Technically, this could be called for various reasons
            // - May need to check to be sure it's a singularity point
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            isValid = false;
        } catch(FailedIntegration<typename MAP::lvec_type, typename MAP::gvec_type>& err) {
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            isValid = false;
        } catch(MapUnderflow& err) {
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            isValid = false;
        } catch(...) {
            //These are also likely singularity errors but haven't found out why yet...
            theFailure.setType( (forward)? Discont::FORWARD_MAP : Discont::BACKWARD_MAP);
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setMessage( "Mapping failure:  unknown reason thrown by pmap.map()" );
            isValid = false;
        }
        
        //If not valid (usually due to a singularity intersection),
        // store partial data but throw error
        if(!isValid) {
            //Note: Data is filled out up to a certain period
            std::vector<MapDataType> partialCATData = tracker.getResult();
            cache.insert( SRDATA(x0,steps,partialCATData) );
            throw theFailure;
        }
        
        //Compute the displacement
        vec_type vOut = params.the_metric.displacement(x0, steps[std::abs(period)-1]);
        //Add this integration to cache
        std::vector<MapDataType> catData = tracker.getResult();
        cache.insert( SRDATA(x0,steps,catData) );
        return vOut;
    }
}


/** Compute and Return the map vector using a cache for data lookup
 *  => Specialized version that throws Discontinuity (EdgeRotationFailure) errors
 *  and checks NonTransverse Section conditions (i.e., Section Separation)
 */
template<class MAP, class PARAM, class SRDATA>
typename MAP::lvec_type mapUsingCache(
    const typename MAP::lvec_type& x0, //Initial condition on map
    const MAP& pmap, const int& period, const PARAM& params,
    std::set<SRDATA>& cache,
    EdgeRotationFailure<typename MAP::lvec_type>& theFailure  //Edge Rotation Failure class
)
{
    typedef typename MAP::lvec_type         vec_type;
    typedef typename MAP::lmat_type         lmat_type;
    typedef typename MAP::state_type        state_type;
    typedef typename MAP::xstate_type       xstate_type;
    typedef typename MAP::rhs_type          RHStype;
    typedef typename SRDATA::MapDataType    MapDataType;
    typedef EdgeRotationFailure<vec_type>   Discont;
    typedef CATtracker<xstate_type,RHStype::numSingularities> TrackerType;
    static const int S = RHStype::numSingularities;
    
    //Temporary objects in case we encounter an error
    bool forward = (period < 0) ? false : true;
    theFailure.reset();
    theFailure.setMessage("Currently no failure");
    theFailure.failurePos = x0;
    theFailure.period = period;
    theFailure.setType( ((forward)? Discont::FORWARD_MAP : Discont::BACKWARD_MAP) );
    
    //Cache variables
    typename std::set<SRDATA>::iterator cit;
    cit = cache.find( SRDATA(x0) );
    
    if ( cit != cache.end() && cit->isThere(period) ) {
        //If this x0 is in the cache AND there is enough returns, use it's data:
        vec_type vOut = (*cit).returns[std::abs(period)-1];
        return vOut;
    } else {
        //Point is not in cache, so let's compute it and add it to the cache
        std::vector<vec_type> steps;
        //Create a tracking system
        const xstate_type* singularityPtr = pmap.rhs().singularities();
        TrackerType tracker;
        for(int k=0; k<S; k++) {
            std::pair<vec_type,lmat_type> tmp = pmap.section().project(singularityPtr[k]);
            tracker.setSingularity(k,tmp.first);
        }
        bool isValid = true;
        try {
            //Original method maps the point
            //pmap.map(x0, steps, period);
            //Map with tracking close approach distances to singularities
            pmap.map_and_track_complete(x0, steps, period, tracker);
        } catch(MapSingularityError<state_type>& err) {
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            //We have to throw these just like Section Separation
            isValid = false;
        } catch(MapUndefined& err) {
            //Technically, this could be called for various reasons
            // - May need to check to be sure it's a singularity point
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            isValid = false;
        } catch(FailedIntegration<typename MAP::lvec_type, typename MAP::gvec_type>& err) {
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            isValid = false;
        } catch(MapUnderflow& err) {
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            isValid = false;
        } catch(...) {
            //These are also likely singularity errors but haven't found out why yet...
            theFailure.setType( (forward)? Discont::FORWARD_MAP : Discont::BACKWARD_MAP);
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setMessage( "Mapping failure:  unknown reason thrown by pmap.map()" );
            isValid = false;
        }
        
        //If not valid (usually due to a singularity intersection),
        // store partial data but throw error
        if(!isValid) {
            //Note: Data is filled out up to a certain period
            std::vector<MapDataType> partialCATData = tracker.getResult();
            cache.insert( SRDATA(x0,steps,partialCATData) );
            throw theFailure;
        }
        
        //Return the map
        vec_type vOut = steps[std::abs(period)-1];
        //Add this integration to cache
        std::vector<MapDataType> catData = tracker.getResult();
        cache.insert( SRDATA(x0,steps,catData) );
        return vOut;
    }
}

/** Compute the map and return full SortableReturnData entry using a cache for data lookup
 *  => Specialized version that throws Discontinuity (EdgeRotationFailure) errors
 *  and checks NonTransverse Section conditions (i.e., Section Separation)
 */
template<class MAP, class PARAM, class SRDATA>
const SRDATA&
mapDataUsingCache(
    const typename MAP::lvec_type& x0, //Initial condition on map
    const MAP& pmap, const int& period, const PARAM& params,
    std::set<SRDATA>& cache,
    EdgeRotationFailure<typename MAP::lvec_type>& theFailure  //Edge Rotation Failure class
)
{
    typedef typename MAP::lvec_type         vec_type;
    typedef typename MAP::lmat_type         lmat_type;
    typedef typename MAP::state_type        state_type;
    typedef typename MAP::xstate_type       xstate_type;
    typedef typename MAP::rhs_type          RHStype;
    typedef typename SRDATA::MapDataType    MapDataType;
    typedef EdgeRotationFailure<vec_type>   Discont;
    typedef CATtracker<xstate_type,RHStype::numSingularities> TrackerType;
    static const int S = RHStype::numSingularities;
    
    //Temporary objects in case we encounter an error
    bool forward = (period < 0) ? false : true;
    theFailure.reset();
    theFailure.setMessage("Currently no failure");
    theFailure.failurePos = x0;
    theFailure.period = period;
    theFailure.setType( ((forward)? Discont::FORWARD_MAP : Discont::BACKWARD_MAP) );
    
    //Cache variables
    typename std::set<SRDATA>::iterator cit;
    cit = cache.find( SRDATA(x0) );
    
    if ( cit != cache.end() && cit->isThere(period) ) {
        //If this x0 is in the cache AND there is enough returns, use it's data:
        return (*cit);
    } else {
        //Point is not in cache, so let's compute it and add it to the cache
        std::vector<vec_type> steps;
        //Create a tracking system
        const xstate_type* singularityPtr = pmap.rhs().singularities();
        TrackerType tracker;
        for(int k=0; k<S; k++) {
            std::pair<vec_type,lmat_type> tmp = pmap.section().project(singularityPtr[k]);
            tracker.setSingularity(k,tmp.first);
        }
        bool isValid = true;
        try {
            //Original method maps the point
            //pmap.map(x0, steps, period);
            //Map with tracking close approach distances to singularities
            pmap.map_and_track_complete(x0, steps, period, tracker);
        } catch(MapSingularityError<state_type>& err) {
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            //We have to throw these just like Section Separation
            isValid = false;
        } catch(MapUndefined& err) {
            //Technically, this could be called for various reasons
            // - May need to check to be sure it's a singularity point
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            isValid = false;
        } catch(FailedIntegration<typename MAP::lvec_type, typename MAP::gvec_type>& err) {
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            isValid = false;
        } catch(MapUnderflow& err) {
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            isValid = false;
        } catch(...) {
            //These are also likely singularity errors but haven't found out why yet...
            theFailure.setType( (forward)? Discont::FORWARD_MAP : Discont::BACKWARD_MAP);
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setMessage( "Mapping failure:  unknown reason thrown by pmap.map()" );
            isValid = false;
        }
        
        //If not valid (usually due to a singularity intersection),
        // store partial data and return info
        if(!isValid) {
            //Note: Data is filled out up to a certain period
            std::vector<MapDataType> partialCATData = tracker.getResult();
            cache.insert( SRDATA(x0,steps,partialCATData) );
            throw theFailure;
        }
        
        //Add this integration to cache
        std::vector<MapDataType> catData = tracker.getResult();
        cache.insert( SRDATA(x0,steps,catData) );
        //Find actual element and return
        cit = cache.find( SRDATA(x0) );
        return (*cit);
    }
}

/// Map output data from a manifold mapping task
template<class MANID, class SRDATA, class MAPDISCONT>
class MapTaskOutputData {
public :
    MapTaskOutputData(const MANID& mid, const SRDATA& data, const MAPDISCONT& discont, const bool& f) :
        theData(data),
        eFail(discont),
        failureDetected(f),
        id(mid)
    {}
    
    /// SortableReturnData (with extended information from CATtracker)
    SRDATA theData;
    /// EdgeRotationFailure
    MAPDISCONT eFail;
    /// bool indicating if failure occurred during propagation
    bool failureDetected;
    /// Manifold/Job id (Use a pair of ids for manifold<->edgeID job)
    MANID id;
};


/** Compute and Return the map vector for a specified period.
 *  => This version stores a (manID,Data) pair to a per-thread cache.
 *     The 'Data' value is a structure composed of SortableReturnData and a MapDiscont
 */
template<class MAP, class PARAM, class SRDATA, typename MANID>
typename MAP::lvec_type mapUsingTaskCache(
    const typename MAP::lvec_type& x0, //Initial condition on map
    const MAP& pmap, const int& period, const PARAM& params,
    std::set<SRDATA>& srDataSet,
    const MANID& id,
    std::vector<MapTaskOutputData<MANID,SRDATA,EdgeRotationFailure<vec_type> > >& perThreadCache
)
{
    typedef typename MAP::lvec_type         vec_type;
    typedef typename MAP::lmat_type         lmat_type;
    typedef typename MAP::state_type        state_type;
    typedef typename MAP::xstate_type       xstate_type;
    typedef typename MAP::rhs_type          RHStype;
    typedef typename SRDATA::MapDataType    MapDataType;
    typedef EdgeRotationFailure<vec_type>   Discont;
    typedef MapTaskOutputData<MANID,SRDATA,Discont>           MTOData;
    typedef CATtracker<xstate_type,RHStype::numSingularities> TrackerType;
    static const int S = RHStype::numSingularities;
    
    //Temporary objects in case we encounter an error
    bool forward = (period < 0) ? false : true;
    Discont theFailure;
    theFailure.reset();
    theFailure.setMessage("Currently no failure");
    theFailure.failurePos = x0;
    theFailure.period = period;
    theFailure.setType( ((forward)? Discont::FORWARD_MAP : Discont::BACKWARD_MAP) );
    
    //Sortable return data cache variables
    typename std::set<SRDATA>::iterator cit;
    cit = srDataSet.find( SRDATA(x0) );
    
    //Check if we got lucky and it's already in the SRDATA cache
    if ( cit != srDataSet.end() && cit->isThere(period) ) {
        //If this x0 is in the cache AND there is enough returns, use it's data:
        vec_type vOut = (*cit).returns[std::abs(period)-1];
        // We still have to store the lookup data in the per-thread cache
        perThreadCache.push_back( MTOData(id,*cit,theFailure,false) );
        return vOut;
    }
    //If NOT (the likely case), we propagate and store to perThreadCache
    else {
        //Point is not in cache, so let's compute it and add it to the cache
        std::vector<vec_type> steps;
        //Create a tracking system
        const xstate_type* singularityPtr = pmap.rhs().singularities();
        TrackerType tracker;
        for(int k=0; k<S; k++) {
            std::pair<vec_type,lmat_type> tmp = pmap.section().project(singularityPtr[k]);
            tracker.setSingularity(k,tmp.first);
        }
        bool isValid = true;
        try {
            //Original method maps the point
            //pmap.map(x0, steps, period);
            //Map with tracking close approach distances to singularities
            pmap.map_and_track_complete(x0, steps, period, tracker);
        } catch(MapSingularityError<state_type>& err) {
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            //We have to throw these just like Section Separation
            isValid = false;
        } catch(MapUndefined& err) {
            //Technically, this could be called for various reasons
            // - May need to check to be sure it's a singularity point
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            isValid = false;
        } catch(FailedIntegration<typename MAP::lvec_type, typename MAP::gvec_type>& err) {
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            isValid = false;
        } catch(MapUnderflow& err) {
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setType( (forward)? Discont::FORWARD_SINGULARITY : Discont::BACKWARD_SINGULARITY);
            isValid = false;
        } catch(...) {
            //These are also likely singularity errors but haven't found out why yet...
            theFailure.setType( (forward)? Discont::FORWARD_MAP : Discont::BACKWARD_MAP);
            theFailure.failIterate = 1 + (int) steps.size();
            theFailure.setMessage( "Mapping failure:  unknown reason thrown by pmap.map()" );
            isValid = false;
        }
        
        //Add this integration to cache
        std::vector<MapDataType> catData = tracker.getResult();
        bool isNotValid = (isValid) ? false : true;
        perThreadCache.push_back( MTOData(id,SRDATA(x0,steps,catData),theFailure,isNotValid) );
        //Return the map if valid
        if(!isValid) {
            return vec_type(50);
        } else {
            return steps[std::abs(period)-1];
        }
    }
}
} // end xavier

#endif
