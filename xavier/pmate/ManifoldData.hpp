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


#ifndef PMATE_MANIFOLD_DATA_HPP
#define PMATE_MANIFOLD_DATA_HPP


#include <vector>
#include <set>
#include <queue>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <util/wall_timer.hpp>
#include <math/fixed_vector.hpp>
#include <math/bezier_spline.hpp>
#include <math/bounding_box.hpp>
#include <math/angle.hpp>
#include <math/intersection.hpp>
#include <maps/metric.hpp>
#include <maps/fixpoints.hpp>
#include <maps/map_analysis.hpp>
#include <pmate/FixedPointData.hpp>
#include <pmate/ManifoldConnection.hpp>
#include <design/Trajectory.hpp>
#include <design/PeriodicOrbit.hpp>
#include <topology/CATtracker.hpp>
#include <topology/EdgeRotationFailure.hpp>
#include <topology/EdgeRotationMapCalls.hpp>
#include <topology/ManifoldClasses.hpp>
#include <topology/MapManifoldSegment.hpp>
#include <topology/MapManifold.hpp>
#include <topology/invariant_manifold.hpp>
#include <topology/STHManifold.hpp>

//OpenMP
#if _OPENMP
#include <omp.h>
#endif

using namespace xavier;
//using namespace topology;

namespace pmate {

///Manifold Data Class that contains ALL manifold data from an input set of fixed points (FixedPointData)
//  Performs manifold advection, stores all MapManifold objects (and segments)
template <class SRDATA, class FP>
class ManifoldData {
public :
    typedef typename SRDATA::VecType                                           VecType;
    typedef typename FP::StateType                                             StateType;
    typedef FixedPointData<FP>                                                 FPDataType;
    typedef typename FPDataType::FixedPointChain                               FixedPointChain;
    typedef typename FPDataType::FPChainsContainer                             FPChainsContainer;
    typedef pmateDesign::Trajectory<StateType>                                 TrajectoryType;
    typedef pmateDesign::PeriodicOrbit<StateType>                              PeriodicOrbitType;
    typedef std::vector<PeriodicOrbitType>                                     POrbitsVector;
    typedef EdgeRotationFailure<VecType>                                       MapDiscont;
    typedef topology::MapManifold<SRDATA,FP>                                   ManifoldType;
    typedef typename ManifoldType::ManifoldSeg                                 ManifoldSegType;
    //Manifold propagation data for tasks
    typedef topology::MapTaskOutputData<boost::uuids::uuid,SRDATA,MapDiscont>  ManMapStorageData;
    typedef std::vector<ManMapStorageData>                                     MMStorageVector;
    //Advection propagation data during subdivision procedures
    typedef topology::AdaptiveEdge<VecType,VecType>                            EdgeType;
    typedef typename EdgeType::IdType                                          EdgeIdType;
    typedef std::pair<int,EdgeIdType>                                          JobIDPair;
    typedef topology::MapTaskOutputData<JobIDPair,SRDATA,MapDiscont>           MapJobStorageData;
    typedef std::vector<MapJobStorageData>                                     MJStorageVector;
    //Intersection Types
    typedef LineSegmentIntersection2D<ManifoldSegType>                         MSegIntersector;
    typedef ManifoldConnection<VecType>                                        ManifoldCon;

    ///Constructor
    ManifoldData(FPDataType& fpData) :
        theFixedPointData(&fpData),
        enableMaxPeriodCutoff(false),
        nthreads(
#ifdef _OPENMP
			omp_get_max_threads()
#else
			1
#endif
			)
    {}

    ///Blank constructor
    ManifoldData() :
        theFixedPointData(0),
        enableMaxPeriodCutoff(false),
        nthreads(1)
    {}

    // Computing (or I/O) the manifold objects ----------------------------------------------------------------
    /** Read ManifoldData (.im file extension)
     *  Note: This is a bit tricky because you need to have fixed points stored
     *        in memory before calling the "read()" function for manifolds.
     *        Hence, why the constructor requires a FixedPointData object.
     */
    void read(const char* filename);

    /// Write ManifoldData (.im file extension)
    bool write(const char* filename, const char* fpDataFileName) const;

    ///Compute all the manifolds
    template<class PMAP, class PARAM>
    void compute(PMAP& amap, PARAM& params, ManifoldSettings& mSettings);

    /// Check the manifolds if the advection is complete in all objects
    bool isAdvectionComplete()
    {
        for(int i=0; (int)mapManifolds.size(); i++) {
            if(!(mapManifolds[i].isComplete)) {
                return false;
            }
        }
        return true;
    }

    /// The manifold objects
    std::map<int,ManifoldType> mapManifolds;
    ///Clear manifold objects (avoid memory leak)
    void clearManifolds()
    {
        mapManifolds.clear();
    }

    /// Poincare map data cache storing data in sortable order
    /// Note: these are only used during creation; neglected upon load operation.
    std::set<SRDATA>  cache, backCache;

    ///Set the number of computational threads
    void setNumThreads(int n)
    {
        nthreads = n;
    }
    /// Set the FixedPointData pointer
    void setFixedPointData(FPDataType& fpData)
    {
        theFixedPointData = &fpData;
    }
    /// Set the maximum period cutoff flag
    void setMaxCutoffFlag(const bool opt)
    {
        enableMaxPeriodCutoff = opt;
    }
    /// Get the maximum period cutoff flag
    bool getMaxCutoffFlag() const
    {
        return enableMaxPeriodCutoff;
    }
    // -----------------------------------------------------------------------------------------------------------

    // Autonomous Manifold Connection Computation ----------------------------------------------------------------
    /// Connection information
    std::vector<ManifoldCon> hhConnections;
    /// Map-based Delta-V assisted connections (same OR different C)
    std::vector<ManifoldCon> augConnections;
    /// Set a ManifoldConneciton given intersection information and two manifoldIDs from THIS object
    void setManifoldCon(const int uID, const int uSegID, const double ut,
                        const int sID, const int sSegID, const double st, ManifoldCon& mc);
    /// Compute homoclinic/heteroclinic connections between the s/u manifolds of this object
    void computeConnections();
    /// Compute augmented connections between s/u manifolds of this object
    void computeConnections(const VecType& maxMapDisp);
    /// Compute connections with another ManifoldData object
    //void computeConnections(ManifoldData& other);
    /** Write ManifoldConnection data to file (the .imc extension)
     *  This requires a known .im file that stores current ManifoldData, but can also
     *  store multiple ManifoldData names if you are computing connections through
     *  multiple Jacobi constants at once.  (It doesn't have to be one at a time.)
     */
    bool writeConnections(const char* filename, std::vector<std::string>& manDataFileName) const;
    /** Read ManifoldConnection data (.imc file extension)
     *  Note: This is also a bit tricky because the corresponding ManifoldData
     *  (and secondary ManifoldData files) are also required to be preloaded in order
     *  to interact with the connection data.
     */
    bool readConnections(const char* filename);
    // -----------------------------------------------------------------------------------------------------------


private:
    /// FixedPointData object pointer that holds the orbits relevant to these manifolds
    FPDataType* theFixedPointData;
    /// Enable a maximum cutoff on allowable periods to evaluate (value from MapParams.max_period)
    bool enableMaxPeriodCutoff;
    /// Number of threads to use in parallel computation
    int nthreads;
    /// Build all the MapManifold objects (assumes "omp single" directive)
    template<class PMAP, class PARAM>
    void buildManifolds(PMAP& amap, PARAM& params, POrbitsVector& orbits,
                        ManifoldSettings& mSettings, MMStorageVector* ptCache);
    /** Build/Simulate the periodic orbits to get a complete simulation history
     *  Note: removes fp's that fail a periodic error check ||xf - x0|| > 1e-6
     */
    template<class PMAP, class PARAM>
    void buildOrbits(PMAP& amap, PARAM& params, POrbitsVector& orbits);
    /** \breif Evaluate stopping criteria for this manifold
     *  Evaluate stopping criteria for this manifold just after advecting the segments
     *  (return estimated completion fraction [0,1])
     *
     *  Note: KAM manifold detection requires that we have access to all manifold objects.
     */
    template<class PMAP, class PARAM>
    double evaluateStoppingCriteria(const int manifoldID, const PMAP& theMap,
                                    const PARAM& params, const ManifoldSettings& settings);
    /*/// Evaluating the Progenitor States for new segments
    template<class PMAP, class PARAM>
    void evaluateProgenitorStates(PMAP& amap, PARAM& params);*/

};

///Compute all the manifolds
template <class SRDATA, class FP>
template <class PMAP, class PARAM>
void ManifoldData<SRDATA,FP>::
compute(PMAP& amap, PARAM& params, ManifoldSettings& mSettings)
{
    PMAP* theMap = &amap;
    PARAM* theMapParams = &params;
    ManifoldSettings* theManSettings = &mSettings;
#if _OPENMP
    // No dynamic scheduling
    omp_set_dynamic(0);
    omp_set_num_threads( nthreads );
#endif

    // Create per-thread cache for mapping information with only manifold uuid
    MMStorageVector* mapDataCache = new MMStorageVector[nthreads];
    // Create per-thread cache for mapping jobs during advection process
    MJStorageVector* jobDataCache = new MJStorageVector[nthreads];
    // Create per-thread cache for propagating full periodic orbits
    POrbitsVector orbits;


    // Setup the parallel region
    #pragma omp parallel
    {
        // Job (task) construction on a single thread, but mapping work is processed in pool
        #pragma omp single
        {
            nvis::timer timer;
            // Run the periodic orbits to collect all internal points (between iterates)
            std::cerr << "-----------------------------------------------------------------\n";
            std::cerr << "PMATE:  Examining validity of input fixed points...\n";
            std::cerr << "-----------------------------------------------------------------\n";
            buildOrbits(amap,params,orbits);
            double dt = timer.elapsed();
            std::cerr << "Periodic orbits simulated and analyzed. (" << dt << " s = " << dt/60.0 << "min)\n";

            // Build the Manifold Objects and their first segments
            std::cerr << "-----------------------------------------------------------------\n";
            std::cerr << "PMATE:  Constructing Manifold objects...\n";
            std::cerr << "-----------------------------------------------------------------\n";
            timer.restart();
            buildManifolds(amap,params,orbits,mSettings,mapDataCache); //Note: submits tasks to propagate first segments
            dt = timer.elapsed();
            std::cerr << "Manifold objects built. (" << dt << " s = " << dt/60.0 << "min)\n";


            //Start advection process
            std::cerr << "\n\n";
            std::cerr << "-----------------------------------------------------------------\n";
            std::cerr << "PMATE:  Manifold objects built.  Starting Advection Process...\n";
            std::cerr << "-----------------------------------------------------------------\n";
            //Some accounting variables
            int manifoldsCompleted = 0;
            int manifoldsInProgress = 0;
            int numParentManifolds = 0;
            double estManProg = 0.0;

            //Process until all manifolds have propagated far enough
            // or triggered stopping condition (manifold reached max depth/arclength or saddle-loop)
            bool isQueueEmpty = false;
            //Create the manifold queue (as a vector of ints)
            std::vector<int> manIntVec;
            std::vector<double> manProg;
            typename std::map<int,ManifoldType>::iterator mapit;
            for(mapit=mapManifolds.begin(); mapit!=mapManifolds.end(); mapit++) {
                if (mapit->second.isComplete) {
                    //These failed to start and so we are finished with them...
                    manifoldsCompleted++;
                    manProg.push_back( 1.0 ); //Done with these...
                } else if (!(mapit->second.hasParent()) ) {
                    manIntVec.push_back( mapit->first ); //Store int for parent manifolds
                    manProg.push_back( 0.0 ); //Estimated Progress for each manifold
                }
            }
            numParentManifolds = manifoldsInProgress = (int) manIntVec.size();


            //Process until the manifolds are done (meeting stopping criteria)
            int numPasses = 0;
            while (!isQueueEmpty) {

                //Loop through all manifolds in queue to submit propagation jobs
                std::vector<JobIDPair> jobs;
                std::vector<int>::iterator manIT, tempIT;
                int manifoldCounter = 0, numManifolds = (int) manIntVec.size();
                for (manIT=manIntVec.begin(); manIT!=manIntVec.end(); ++manIT) {
                    mapManifolds[(*manIT)].getEdgeNodeMapTasks((*manIT),jobs);
                    manifoldCounter++;
                    /*std::ostringstream os;
                    os << "\rGenerating propagation tasks: " << manifoldCounter << "/"
                       << numManifolds << " (%" << (float)manifoldCounter*100/numManifolds << ") \r";
                    std::cout << os.str() << std::flush;*/
                }

                //Utilize OpenMP Task-Scheduling to evaluate map calls
                int numProps = (int) jobs.size();
                if(theMapParams->verbose) {
                    std::cout << "   " << numProps << " propagation jobs queued for processing...\n";
                }
                timer.restart();
                int propsComplete = 0;
                for (int n =0; n<numProps; n++) {
                    #pragma omp task shared(amap,params,jobs,jobDataCache,propsComplete)
                    {
                        //Find the manifold object referring to this job
                        ManifoldType& theManifold = mapManifolds[jobs[n].first];

                        //Run the propagation
                        theManifold.propagateManifoldArcTask(jobs[n],amap,params,cache,backCache,jobDataCache);

                        #pragma omp atomic
                        propsComplete++;

                        dt = timer.elapsed();
                        std::ostringstream os;
                        os << "\rManifold Task Propagation: " << propsComplete << "/" << numProps
                        << " (" << (float)propsComplete*100/numProps << "%) in "
                        << dt << "s. (" << (float)propsComplete/dt << "Hz)       \r";
                        #pragma omp critical
                        std::cout << os.str() << std::flush;
                    }
                }

                //Extract information from map calls and save back into MapManifold objects
                #pragma omp taskwait
                {
                    std::cout << "ManifoldData:  Completed Propagation pass [" << numPasses << "]\n";
                    std::cout << "               Assigning data to manifold objects from threads...\n";
                } //End taskwait barrier
                //Search the thread caches for relevant jobs/data
                for (int k=0; k<nthreads; k++) {
                    typename MJStorageVector::iterator mjit, tempit;
                    for(mjit=jobDataCache[k].begin(); mjit!=jobDataCache[k].end(); ++mjit) {
                        ManifoldType& thisMan = mapManifolds[mjit->id.first];
                        thisMan.assignEdgeNodeData(
                            mjit->id.second,mjit->theData,
                            mjit->eFail,mjit->failureDetected,
                            amap, params, cache, backCache
                        );

                        //When done with data, lets pull it out of cache
                        tempit = mjit;
                        tempit--;
                        jobDataCache[k].erase(mjit);
                        mjit = tempit;
                    }//End cache loop
                }//End thread loop


                //Evaluate new information to see where refinement is needed
                std::queue<int> manifoldIDQueue;
                //Note: we could make this a priority_queue with manifolds that
                //have small advancement having highest priority and working in chunks,
                // but that may slow down the overall process as the task bucket is
                // emptied and refilled more frequently.
                std::cout << "ManifoldData:  Checking which segments to refine or advect...\n";
                std::vector<std::pair<int,int> > msIDsToEvalProgState;
                if(theMapParams->debug) {
                    std::cout << "ManifoldData: Update state for manifolds:\n";
                    std::cout << std::setw(6) << "ID"
                              << std::setw(12) << "pType"
                              << std::setw(8) << "Proc"
                              << std::setw(8) << "Advect"
                              << std::setw(8) << "WSeg"
                              << std::setw(8) << "WSegD"
                              << std::setw(8) << "NWSegD"
                              << std::setw(10) << "NWSegSDL"
                              << std::setw(12) << "NumSegs"
                              << std::setw(15) << "NumSegsD"
                              << std::setw(10) << "ProgAtD" << "\n";
                }
                int numRefineChecks = (int) manIntVec.size();
                int refineProcessed = 0;
                double refineTime = 0.0, advectTime = 0.0, evalStopTime = 0.0;
                timer.restart();
                for (manIT=manIntVec.begin(); manIT!=manIntVec.end(); ++manIT) {
                    //Inform the user about time break down (seem to be slowing up here??)
                    //dt = timer.elapsed();
                    //double t0 = dt;
                    //std::ostringstream refine_os;
                    //refine_os << "\rManifoldData: " << refineProcessed << "/" << numRefineChecks
                    //     << " (" << (float)refineProcessed*100/numRefineChecks << "%) in "
                    //     << dt << "s. (" << (float)refineProcessed/dt << "Hz)      \r";
                    //std::cout << refine_os.str() << std::flush;
                    //The Current manifold
                    ManifoldType& thisMan = mapManifolds[(*manIT)];
                    //Evaluate refinement criteria
                    if(thisMan.isForward()) {
                        thisMan.refineSegmentsTest(cache,*theMap,*theMapParams,*theManSettings);
                    } else {
                        thisMan.refineSegmentsTest(backCache,*theMap,*theMapParams,*theManSettings);
                    }

                    //Print the update entry
                    if(theMapParams->debug) {
                        thisMan.printUpdateEntry();
                    }
                    //refineTime += timer.elapsed() - t0;
                    //t0 = timer.elapsed();

                    //Test to see if this manifold is done with the current working segment
                    if (thisMan.emptyNodes < 1) {
                        //Record time in advect

                        /*refine_os << "\rManifoldData: " << refineProcessed << "/" << numRefineChecks
                           << " (" << (float)refineProcessed*100/numRefineChecks << "%) in "
                           << dt << "s. (" << (float)refineProcessed/dt << "Hz) [Advection]      \r";
                        std::cout << refine_os.str() << std::flush;*/
                        //Process the manifold segments and store
                        if(!theMapParams->verbose && theMapParams->debug) {
                            std::cout << "---------------------------------------------------------\n";
                            std::cout << "ManifoldData:  Manifold Ready for Advection:\n";
                            thisMan.info();
                            std::cout << "---------------------------------------------------------\n";
                        }
                        if(thisMan.isForward()) {
                            //try {
                            thisMan.advectSegment(*theMap,*theMapParams,*theManSettings,cache);
                            /*} catch(std::exception& e) {
                              //This manifold is done due to an error.
                              std::cerr << " Forcing Manifold:\n";
                              thisMan.info();
                              std::cerr << " to quit due to error : " << e.what() << "\n";
                              thisMan.isComplete = true;
                            }*/
                        } else {
                            //try {
                            thisMan.advectSegment(*theMap,*theMapParams,*theManSettings,backCache);
                            /*} catch(std::exception& e) {
                              //This manifold is done due to an error.
                              std::cerr << " Forcing Manifold:\n";
                              thisMan.info();
                              std::cerr << " to quit due to error : " << e.what() << "\n";
                              thisMan.isComplete = true;
                            }*/
                        }
                        //advectTime += timer.elapsed() - t0;
                        //t0 = timer.elapsed();

                        /*refine_os << "\rManifoldData: " << refineProcessed << "/" << numRefineChecks
                           << " (" << (float)refineProcessed*100/numRefineChecks << "%) in "
                           << dt << "s. (" << (float)refineProcessed/dt << "Hz) [Stop Criteria]      \r";
                        std::cout << refine_os.str() << std::flush;*/

                        //Evaluate stopping criteria (as we complete the segment):
                        //Also estimate total completion [but may not be great as it's an exponential growth]
                        manProg[(*manIT)] = evaluateStoppingCriteria(
                                                (*manIT),*theMap,*theMapParams,*theManSettings);

                        //Add new children to queue if encountered transversality violations
                        thisMan.childMaturity(manifoldIDQueue);
                        //evalStopTime += timer.elapsed() - t0;

                    }
                    //refineProcessed++;
                }
                dt = timer.elapsed();
                //double totalTime = refineTime + advectTime + evalStopTime;
                //std::cout << "ManifoldData:  Refine/Advect/Stop took " << dt
                //          << "s. (" << (float) numRefineChecks / dt << "Hz) \n";
                //std::cout << "    Total = " << totalTime << "s.\n";
                //std::cout << "    Refine = " << refineTime << "s. (" << refineTime*100.0/totalTime << "%)\n";
                //std::cout << "    Advect = " << advectTime << "s. (" << advectTime*100.0/totalTime << "%)\n";
                //std::cout << "    EvalStop = " << evalStopTime << "s. (" << evalStopTime*100.0/totalTime << "%)\n";



                //Perform the the check on stopping after progenitor states are evaluated:
                for (manIT=manIntVec.begin(); manIT!=manIntVec.end(); ++manIT) {
                    ManifoldType& thisMan = mapManifolds[(*manIT)];
                    //Fill queue if manifold advection is incomplete
                    if (!(thisMan.isComplete) ) {
                        manifoldIDQueue.push((*manIT));
                    } else {
                        manifoldsCompleted++;
                    }
                }

                //Refill the processing elements with what remains unfinished:
                manIntVec.clear();
                manifoldsInProgress = (int) manifoldIDQueue.size();
                double partialManProg = 0.0;
                while(!manifoldIDQueue.empty()) {
                    manIntVec.push_back( manifoldIDQueue.front() );
                    partialManProg += manProg[ manifoldIDQueue.front() ];
                    manifoldIDQueue.pop();
                }

                estManProg = (double) manifoldsCompleted + (double) partialManProg;
                estManProg /= (double) mapManifolds.size();

                //Update user:
                std::cout << "---------------------------------------------------------\n";
                std::cout << "ManifoldData: Computation progress: \n";
                std::cout << "    " << (int) mapManifolds.size() << " Total Manifolds\n";
                std::cout << "    " << manifoldsInProgress << " in progress \n";
                std::cout << "    " << manifoldsCompleted << " Manifolds Complete\n";
                //std::cout << "    " << numParentManifolds << " Total 'Adult' Manifolds\n";
                std::cout << "    " << estManProg*100.0 << "% complete\n";
                std::cout << "---------------------------------------------------------\n";

                //While() loop exit condition
                isQueueEmpty = ((int)manIntVec.size() > 0) ? false : true;
                if (isQueueEmpty && theMapParams->verbose) {
                    std::cerr << "ManifoldData: Map Advection Terminated 'Normally'\n";
                }
                numPasses++;
            } //End while loop

            //----------------------------------------------------------------------------------
            // PROGENITOR STATES:
            // When all manifolds are DONE, we can evaluate the Progenitor States (all at the same time)
            //----------------------------------------------------------------------------------
            /*std::cout << "---------------------------------------------------------\n";
            std::cout << "ManifoldData:  Processing Progenitor States\n";
            std::cout << "---------------------------------------------------------\n";
            evaluateProgenitorStates(*theMap,*theMapParams);*/

        } //End Single thread
    } //End Parallel statement

    //Clear job storage containers
    for(int k=0; k<nthreads; k++) {
        mapDataCache[k].clear();
        jobDataCache[k].clear();
    }
    delete[] mapDataCache;
    delete[] jobDataCache;


    //Sanity check
    //Evaluate all manifolds to see if complete
    //allComplete = isAdvectionComplete();
}

/**Build/Simulate the periodic orbits to get a complete simulation history
  *  Note: removes fp's that fail a periodic error check ||xf - x0|| > 1e-6
  */
template <class SRDATA, class FP>
template <class PMAP, class PARAM>
void ManifoldData<SRDATA,FP>::
buildOrbits(
    PMAP& amap,
    PARAM& params,
    POrbitsVector& orbits
)
{
    typedef typename PMAP::return_state    ReturnState;
    typedef typename PMAP::lvec_type       vec_type;
    typedef typename PMAP::lmat_type       lmat_type;
    typedef typename PMAP::xstate_type     xstate_type;
    typedef typename PMAP::rhs_type        RHStype;
    typedef typename SRDATA::MapDataType   MapDataType;
    typedef topology::CATtracker<xstate_type,RHStype::numSingularities> TrackerType;
    static const int S = RHStype::numSingularities;
    int maxThreads = 1;
#if _OPENMP
    maxThreads = omp_get_max_threads();
#endif
    int numOrbits = theFixedPointData->getNumOrbits();
    //Setup Orbit Data
    std::vector<FP> fpDataVec;
    for(int k=0; k<numOrbits; k++) {
        //Always first point in chain
        FP theFP = theFixedPointData->getFixedPoint(k,0);
        fpDataVec.push_back(theFP);
    }
    //Per-thread cache for propagation result
    POrbitsVector* poSimDataCache = new POrbitsVector[maxThreads];
    //Propagate all orbits by submitting task jobs
    int propsComplete = 0;
    int numInvalid = 0;
    for(int k=0; k<numOrbits; k++) {
        #pragma omp task shared(amap,params,poSimDataCache,propsComplete,numInvalid)
        {
            //Get the fixpoint object (and data like ic, type, and period)
            FP theFP( fpDataVec[k] );
            //Propagate each orbit forward by given period to get full set of states
            bool valid = true;
            std::vector<ReturnState> itersOut, internalOut;
            // Have to run with a CAT tracker to check for sub-surface passes:
            //Create a tracking system
            const xstate_type* singularityPtr = amap.rhs().singularities();
            TrackerType tracker;
            for(int s=0; s<S; s++)
            {
                std::pair<vec_type,lmat_type> tmp = amap.section().project(singularityPtr[s]);
                tracker.setSingularity(s,tmp.first);
            }
            // Call the map to get internal state information
            try {
                //amap.flow_map(theFP.pos,itersOut,internalOut,theFP.K);
                amap.map_and_track_complete(theFP.pos, itersOut, internalOut, theFP.K, tracker);
            } catch(...)
            {
                //if anything screwy happens, then this is not a valid fixed point
                valid = false;
                if(params.debug) {
                    std::cout << "Periodic Orbit: Orbit " << k
                    << " INVALID: via mapping failure!\n";
                }
                //if((int)itersOut.size() < theFP.K) valid = false;
            }

            //Construct the PeriodicOrbit object
            PeriodicOrbitType po;
            po.fpIndex = k;
            if(!valid)
            {
                po.setValidity(false); //Propagation didn't work so it's bad orbit
                #pragma omp atomic
                numInvalid++;
            } else {
                //Check the closest approach distances to see if within bodies
                //Note: we should have a parameter for this
                //if(params.excludePIntersectors) {
                std::vector<MapDataType> catData = tracker.getResult();
                for(int s=0; s<S; s++)
                {
                    //Safe distance based on inputs (settings)
                    double rSafeDist = amap.rhs().getSingularitySafeDistance(s);
                    double minCA = catData[std::abs((int)theFP.K)-1][s];
                    if (minCA <= rSafeDist) {
                        if(params.debug) {
                            std::cout << "Periodic Orbit: Orbit " << k
                            << " INVALID: has minCA (" << minCA << ") < rSafeDist (" << rSafeDist << ") for s=" << s << "\n";
                        }
                        valid = false;
                    }
                }

                if(valid)
                {
                    for(int i=0; i<(int) internalOut.size(); i++) {
                        po.states.push_back( internalOut[i].x );
                        po.times.push_back( internalOut[i].t );
                    }
                    //Check for periodic error
                    double perr = po.getPeriodicError();
                    if ( perr > 2.e-2 ) {
                        //Then this is not a valid periodic orbit for computing manifolds
                        if(params.debug) {
                            std::cout << "Periodic Orbit: Orbit " << k
                                      << " INVALID: has periodic error (" << perr << ") > " << 2.e-2 << "\n";
                        }
                        valid = false;
                    }
                }
                if(!valid)
                {
                    po.setValidity(false);
                    #pragma omp atomic
                    numInvalid++;
                }

            }
            //Store the periodic orbit data to a per-thread cache
            int tid = 0;
#if _OPENMP
            tid = omp_get_thread_num();
#endif
            poSimDataCache[tid].push_back( po );

            #pragma omp atomic
            propsComplete++;

            //Update the user
            std::ostringstream os;
            os << "\rPeriodic Orbit Propagation: " << propsComplete << "/" << numOrbits
               << " (" << (float)propsComplete*100/numOrbits <<"%)    \r";
            #pragma omp critical
            std::cout << os.str() << std::flush;
        }
    }

    //Extract and save simulation information
    #pragma omp taskwait
    {
        //Update user about valid orbits
        std::cout << "Periodic Orbit Propagation: Complete. " << numInvalid << " orbits invalid ("
                  << (float)numInvalid*100/numOrbits <<"%)\n";
        //Loop through each thread cache
        for(int n=0; n<maxThreads; n++) {
            for(int i=0; i<(int)poSimDataCache[n].size(); i++) {
                //Do NOT save orbits that fail
                if( poSimDataCache[n][i].isValid() ) {
                    orbits.push_back( poSimDataCache[n][i] );
                }
            }
        }
        //Done with the cache
        for(int n=0; n<maxThreads; n++) {
            poSimDataCache[n].clear();
        }
        delete[] poSimDataCache;
    }


}


///Compute all the manifolds
template <class SRDATA, class FP>
template <class PMAP, class PARAM>
void ManifoldData<SRDATA,FP>::
buildManifolds(
    PMAP& amap,
    PARAM& params,
    typename ManifoldData<SRDATA,FP>::POrbitsVector& orbits,
    ManifoldSettings& mSettings,
    MMStorageVector* ptCache
)
{
    PMAP* theMap = &amap;
    PARAM* theMapParams = &params;
    ManifoldSettings* theManSettings = &mSettings;
    //Note: This is assumed to be within the "omp single" thread scope

    // For each saddle-type fixed point
    //int numOrbits = theFixedPointData->getNumOrbits();
    int numOrbits = (int) orbits.size(); //Only stores valid orbits
    int numOrbitsSkipped = 0;
    std::vector< FixedPointChain > fpChains;
    theFixedPointData->getData(fpChains);
    Perturbation mTypes[4] = { UNSTABLE_MINUS, UNSTABLE_PLUS, STABLE_MINUS, STABLE_PLUS };
    std::vector<int> parentIDs;
    int manID = 0;
    for(int j=0; j<numOrbits; ++j) {
        //FixedPointData index
        int fpID = orbits[j].fpIndex;
        //int fpID = j; //Not all orbits are necessarily stored

        //Only create manifolds for saddles
        if (!(fpChains[fpID][0].saddle)) {
            if(theMapParams->debug) {
                std::cerr << "Build Orbits: Skipping " << fpID << " as it was not a saddle\n";
            }
            continue;
        }

        //Filter out orbits that are outside some maximum stability index
        //  Note: used to cut computation costs for large runs
        double si = fpChains[fpID][0].si;
        if (fabs(si) >= mSettings.max_si) {
            if(theMapParams->debug) {
                std::cerr << "Build Orbits: Skipping " << fpID << " as SI = " << si << " is too high\n";
            }
            numOrbitsSkipped++;
            continue;
        }

        //Filter out orbits that are not trusted (have nan as eigenvectors)
        if ( ISNANCALL( fpChains[fpID][0].evec[0][0] ) ||
                ISNANCALL( fpChains[fpID][0].evec[0][1] ) ) {
            if(theMapParams->debug) {
                std::cerr << "Build Orbits: Skipping " << fpID << " as it has NAN evecs\n";
            }
            numOrbitsSkipped++;
            continue;
        }

        //Filter out orbits that are outside the max_period cutoff:
        if ( enableMaxPeriodCutoff ) {
            //|thePeriod| > max_period are skipped
            if( ( fpChains[fpID][0].K > params.max_period ) ||
                    ( (fpChains[fpID][0].eval[1] < 0.0) && (2*fpChains[fpID][0].K > params.max_period) ) ) {
                if(theMapParams->debug) {
                    std::cerr << "Build Orbits: Skipping " << fpID << " as p=" << fpChains[fpID][0].K << " is too big (or Mobius) \n";
                }
                numOrbitsSkipped++;
                continue;
            }
        }

        //For each fixed point within a given chain
        int p = fpChains[fpID][0].K;
        for (int i=0; i<p; i++) {
            //Assume all ADULT objects with no children
            for(int k=0; k<4; k++) {
                mapManifolds.insert(
                    std::pair<int,ManifoldType>(manID,ManifoldType(mTypes[k],topology::ADULT))
                );
                ManifoldType& newManifold = mapManifolds[manID];
                newManifold.baseFixedPoint = theFixedPointData->getFixedPoint(fpID,i);
                newManifold.theOrbit = orbits[j];
                newManifold.fpdOrbitIdx = fpID;
                newManifold.fpdPointIdx = i;
                newManifold.setPeriod();
                newManifold.manifoldID = manID;
                const double& lambdaMax = fabs(newManifold.baseFixedPoint.eval[1]);
                //Manifold Strength parameter
                newManifold.weakStrength = (lambdaMax >= 50.0)? false : true;
                //Or in-plane stability index:
                //newManifold.weakStrength = (fabs(newManifold.baseFixedPoint.si) >= 25.0)? false : true;
                parentIDs.push_back( manID );

                //Increment manifold identifier
                manID++;

            }

            // PARENT-CHILD CONCEPT:  ...TODO...
            // Setup base 2 manifolds (or 4 if Mirror theorem is N/A)
            // Link children to parents ...
            // BIG NOTE:::
            //  This may not be accurate for non saddle loops!!!
            //  Evidence shows that sub returns form different "manifold" arcs emanating from
            //  the subsequent fixed points of a given orbit.  (kinda makes sense too!).
            //  Possible problem:  Apart from KAM manifolds, sub-iterates don't seem to fall
            //   on the "sister" manifold lines...  Needs more investigation!
            //  Mirror theorem would still apply though.

        }
    }


    // Initial step values for each manifold
    bool allManifoldsStarted = false;
    int numParents = (int) parentIDs.size();
    //Initial step points -> correlated in index with adults (pointers)
    std::vector<VecType> initialStepPoints(numParents,VecType(0));
    std::vector<VecType> initialStepPMap(numParents,VecType(0));
    std::vector<int> mIdxToSeed(numParents,0);
    for(int i=0; i<numParents; i++) {
        mIdxToSeed[i] = i;
    }
    int faultyParents = 0;
    while (!allManifoldsStarted) {
        //Propagate the first segment current steps with openMP tasks
        int numProps = (int) mIdxToSeed.size();
        int propsComplete = 0;
        for(int i=0; i<numProps; i++) {
            #pragma omp task shared(amap,params,mSettings,initialStepPoints,propsComplete)
            {
                int mID = mIdxToSeed[i];
                ManifoldType& thisManifold =  mapManifolds[ parentIDs[ mID ] ];
                VecType localStep = initialStepPoints[ mID ];
                //Generate the first step-off point (if first time through loop)
                if(thisManifold.firstStepCount == 0)
                {
                    localStep = thisManifold.firstSegmentStep(amap,params,mSettings);
                    /*#pragma omp critical
                    {
                      //std::cout << " Build=>thisManifold : orbitID = " << thisManifold.fpdOrbitIdx << "\n"
                      //          << "                       fpID =    " << thisManifold.fpdPointIdx << "\n"
                      std::cout << " Build=>baseFP Position = " << thisManifold.baseFixedPoint.pos << "\n"
                                << "             First Step = " << localStep << "\n";
                    }*/
                }
                //Propagate point
                thisManifold.propagateManifoldArcTask(localStep,amap,params,cache,backCache,ptCache);
                initialStepPoints[mID] = localStep;
                #pragma omp atomic
                propsComplete++;
                std::ostringstream os;
                os << "\rStarter Seeds:" << propsComplete << "/" << numProps
                << " (" << (float)propsComplete*100/numProps << "%)       \r";
                #pragma omp critical
                std::cout << os.str() << std::flush;
            }
        }
        //Extract / save the propagations / check if further steps are needed
        #pragma omp taskwait
        {
            for(int i=0; i<numProps; i++) {
                //Pull out information pertaining to this round of data
                int mID = mIdxToSeed[i];
                for(int k=0; k<nthreads; k++) {
                    typename std::vector<ManMapStorageData>::iterator mit, tempit;
                    for(mit=ptCache[k].begin(); mit!=ptCache[k].end(); ++mit) {
                        //Look for relevant data by matching id
                        ManifoldType& thisManifold = mapManifolds[ parentIDs[ mID ] ];
                        if ( thisManifold.isSameManifold( mit->id ) ) {
                            //Check for failure
                            if(mit->failureDetected) {
                                initialStepPMap[ mID ] = VecType(0);
                                //std::cerr << "ManifoldData:  Detected Failure during Initial Step Propagation\n";
                            } else {
                                //Everything is ok, insert to cache
                                if(thisManifold.isForward()) {
                                    cache.insert( mit->theData );  //Does nothing if already there...
                                } else {
                                    backCache.insert( mit->theData );
                                }
                                initialStepPMap[ mID ] = mit->theData.getReturn(fabs(thisManifold.getThePeriod()));
                            }
                            //Once done with it, pull it out of ptCache
                            tempit = mit;
                            tempit--;
                            ptCache[k].erase(mit);
                            mit = tempit;
                        } //Otherwise, do nothing
                    }
                }
            }

            //Check if step meets alignment criteria
            // ->Could be tasks/but not a lot of gain for the check [and significant recoding]
            std::vector<int>::iterator idIT, tmpidIT;
            for(idIT=mIdxToSeed.begin(); idIT!=mIdxToSeed.end(); ++idIT) {
                int mID = *idIT;
                ManifoldType& thisManifold = mapManifolds[ parentIDs[mID] ];
                bool testSuccess = false;
                if ((thisManifold.firstStepCount + 1) < 20) {
                    //First 20 iterates, just try (20-i)*eps as the step
                    testSuccess = thisManifold.firstSegmentCheck(
                                      initialStepPoints[mID],initialStepPMap[mID],
                                      *theMap, *theMapParams, *theManSettings );
                } else if((thisManifold.firstStepCount + 1) < 26) {
                    //Use halving to find a really small step size (limit up to 1e-11)
                    // So step = eps / (2^(i-20))
                    testSuccess = thisManifold.firstSegmentCheckDivisor(
                                      initialStepPoints[mID],initialStepPMap[mID],
                                      *theMap, *theMapParams, *theManSettings );
                } else {
                    //Likely that there is a section separation happening close to the fixed point.
                    //In this case, we accept the minimum linear-step (from autoInitStep)
                    thisManifold.firstSegmentForcedStep( initialStepPoints[mID],
                                                         *theMap, *theMapParams, *theManSettings);
                    testSuccess = true;
                }
                //If yes, store first segment
                if(testSuccess) {
                    try {
                        if (thisManifold.isForward()) {
                            thisManifold.storeFirstSegment(initialStepPoints[mID],cache,*theMap,*theMapParams,*theManSettings);
                        } else {
                            thisManifold.storeFirstSegment(initialStepPoints[mID],backCache,*theMap,*theMapParams,*theManSettings);
                        }
                    } catch(...) {
                        //We could have a map singularity issue in these calls:
                        testSuccess = false;
                        //With checks for singularity issues, this SHOULD NEVER HAPPEN!!!
                        std::cerr << "ManifoldData: Unable to store first segment for manID = " << mID
                                  << " for orbitID = " << thisManifold.fpdOrbitIdx << " ptID = " << thisManifold.fpdPointIdx << "\n";
                    }
                    if (testSuccess) {
                        //Erase from processing
                        tmpidIT = idIT;
                        tmpidIT--;
                        mIdxToSeed.erase(idIT);
                        idIT = tmpidIT;
                    }
                }
                if(!testSuccess) {
                    //if not, the step is modified by firstSegmentCheck, so try again
                    //Check if too many steps occurred [SHOULD NOT BE CALLED NOW!]
                    if( (thisManifold.firstStepCount + 1) >= 28 ) {
                        //This is done and we can't do anything
                        std::cerr << "ManifoldData:  Initial step process exceeds allowed iterations for fp = "
                                  << thisManifold.baseFixedPoint.pos << "\n    (p=" << thisManifold.baseFixedPoint.K
                                  << "). Using forced starting condition...\n";
                        //throw std::runtime_error("Initial step process failed");

                        //Mark the manifold as invalid?
                        //thisManifold.isComplete = true;

                        //Remove from processing
                        tmpidIT = idIT;
                        tmpidIT--;
                        mIdxToSeed.erase(idIT);
                        idIT = tmpidIT;
                    }
                }
            } //End criteria check loop

            //Not done if there are still parents to adjust
            allManifoldsStarted = ((int)mIdxToSeed.size() > 0) ? false : true;

        } //End taskwait section

    } //End while loop

    //Set the complementary IDs for ALL manifolds (For SaddleLoop/KAM exit conditions)
    for(int i=0; i<(int)mapManifolds.size(); i++) {
        bool isUnstable = mapManifolds[i].isForward();
        int fpID = mapManifolds[i].fpdOrbitIdx;
        std::vector<int> compIDs;
        //Examine all other manifolds for opposite (stable/unstable) ones from same orbit
        for(int j=0; j<(int)mapManifolds.size(); j++) {
            if(i==j) {
                continue;
            }
            if( (fpID == mapManifolds[j].fpdOrbitIdx) &&
                    (isUnstable != mapManifolds[j].isForward()) ) {
                //Store the complementary manifold index
                compIDs.push_back(j);
            }
        }
        //Store the vector of indexes
        mapManifolds[i].complementaryIDs = compIDs;
        //Also, tally faulty manifolds (ones that couldn't start)
        if(mapManifolds[i].isComplete) {
            faultyParents++;
        }
    }

    //Prompt user about orbit information
    if (numOrbitsSkipped > 0 ) {
        std::cout << "ManifoldData:  " << numOrbitsSkipped << "/" << numOrbits << " ("
                  << 100.0*(float)numOrbitsSkipped/(float)numOrbits
                  << "%) orbits are removed by either max stability index constraint (|si| < "
                  << mSettings.max_si << "), max_period constraint, or due to eigenvector computation failure.\n";
    }
    //Prompt user about faulty parents
    if (faultyParents > 0) {
        int numManifolds = (int)mapManifolds.size();
        std::cout << "ManifoldData:  " << faultyParents << "/" << numManifolds << " ("
                  << 100.0*(float)numOrbitsSkipped/(float)numManifolds
                  << "%) manifolds failed to start correctly and will be skipped.\n";

    }

    //TODO:
    //After initial steps, fill out kids with relevant information...
}


///Evaluate stopping criteria for this manifold just after advecting the segments
template <class SRDATA, class FP>
template <class MAP, class PARAM>
double ManifoldData<SRDATA,FP>::
evaluateStoppingCriteria(
    const int manID, const MAP& theMap, const PARAM& params, const ManifoldSettings& settings)
{
    //Evaluate stopping criteria:
    ManifoldType& theManifold = mapManifolds[manID];
    //Exit if already complete (happens on manifolds without a good start)
    if(theManifold.isComplete) {
        return 1.0;
    }

    //-----------------------------------------------------------
    //Set complete flag if condition triggered
    //-----------------------------------------------------------
    //0)For 'weak' manifolds, we should see if we have at least a minimum arc length of 1.0:
    bool beyondWeakMinThreshold = false;
    double minLength = 0.75;
    double epsSaddleLoop = 50.0*settings.delta_tau_min;
    static const int numSing = MAP::rhs_type::numSingularities;

    if( settings.enableMaxArclength ) {
        minLength = std::min(0.75,settings.max_arc_length);
    }
    if((!theManifold.weakStrength) ||
            (theManifold.weakStrength && theManifold.progressLength >= minLength)
      ) {
        beyondWeakMinThreshold = true;
    }

    //1)Achieved max progress arc length (disable if not beyond 'weak' threshold)
    if( beyondWeakMinThreshold && settings.enableMaxArclength ) {
        if (theManifold.progressLength >= settings.max_arc_length ) {
            if(params.debug)
                std::cout << " Stop Condition: Progress : length = " << theManifold.progressLength
                          << " > MaxArclength = " << settings.max_arc_length << "\n";
            theManifold.isComplete = true;
        }
    }
    //2)Number of seed segments (disable if not beyond 'weak' threshold)
    if( beyondWeakMinThreshold && settings.enableMaxSeg ) {
        if (theManifold.numSeedSegments >= settings.maxNumSeedSeg ) {
            if(params.debug)
                std::cout << " Stop Condition: NumSeeds : NumSeeds = " << theManifold.numSeedSegments
                          << " >= MaxNumSeeds = " << settings.maxNumSeedSeg << "\n";
            theManifold.isComplete = true; //indicate stoppage
        }
    }
    //2a) Number of segments (always, even if not beyond 'weak' threshold)
    if (settings.enableMaxSeg && theManifold.numSegments >= settings.maxNumSeg ) {
        if(params.debug)
            std::cout << " Stop Condition: Max Segments : NumSegs = " << theManifold.numSegments
                      << " >= MaxNumSegs = " << settings.maxNumSeg << "\n";
        theManifold.isComplete = true; //indicate stoppage
    }

    //3)Achieved max seed length -
    if ( settings.enableMaxArclength && (theManifold.seedLength >= settings.maxSeedArclength)) {
        if(params.debug)
            std::cout << " Stop Condition: SeedLength : length = " << theManifold.seedLength
                      << " > MaxSeedArclength = " << settings.maxSeedArclength << "\n";
        theManifold.isComplete = true;
    }
    //Note: The first 3 tests won't exit this function if true,
    //      so we will still run the other tests/computations.  Important to perform
    //      the saddle-loop test to evaluate and crop data for design.

    //Next Depth in Invariant Manifold Tree
    int nextSegDepth = 1 + theManifold.getDepth(theManifold.workingSegID);

    //Split decisions based on manifold strength:
    if (!theManifold.weakStrength) { //----------------------------------------------------------------------

        //5)Saddle-Loop stopping conditions (indicating we found a homoclinic loop or saddle-loop)
        //  We want to stop here to avoid the repeated (or nearly repeated) advection of manifolds.

        //Skip this check if we have already gathered a hefty portion of the manifold (unlikely a saddle loop)
        /*if (!theManifold.isComplete && theManifold.progressLength < 10.0) { //Only intended for KAM-like manifolds
          // Current depth layer segments
          std::vector<int> currentLayerIds;
          theManifold.getAllSegmentIDsAtDepth(currentLayerIds);
          //Information to compute distances from primary : Asymptote detection
          std::vector<double> safeDistances;
          for(int k=0;k<numSing;k++) safeDistances.push_back( theMap.rhs().getSingularitySafeDistance(k) );
          typedef typename MAP::lmat_type mat_type;
          typedef typename MAP::xstate_type xstate_type;
          typedef std::pair<vec_type,mat_type>  ProjPair;
          const xstate_type* singPtr = theMap.rhs().singularities();

          //Compare to complementary manifold segments at their current depth layer
          //Note: this should be another task construct and dont in parallel, but I'm lazy
          for(int j=0;j<(int)theManifold.complementaryIDs.size();j++) {
            //Get the complementary manifold
            int compID = theManifold.complementaryIDs[j];
            ManifoldType& compMan = mapManifolds[compID];
            //Skip if this is a bad manifold that didn't start
            if ((int)compMan.segments.size()<1) continue;
            //Skip if this manifold is already done
            if (compMan.isComplete) continue;

            //For each segment at the current depth level
            std::vector<int> cmLayerIDs;
            compMan.getAllSegmentIDsAtDepth(cmLayerIDs);
            std::map<int,int> ccmMap; //Assume it's a continiguous group
            //We must first check if something is partially within the primary or on an asymptote
            // - these are NOT saddle loops!
            //Compare to find the first grouping
            bool overlapDetected = false;
            for(int jj=0; jj<(int)cmLayerIDs.size(); jj++) {
              //Test for an asymptote by seeing if x coord is within safe distance of singularity
              int cmID = cmLayerIDs[jj];
              ManifoldSegType& compSeg = compMan.segments[cmID];
              //Test to see if we have a vertical & long segment
              if(compSeg.isNearlyVertical() && (compSeg.length(params)>settings.delta_max) ) continue;
              //Test if an asymptote is likely by singularity interesection
              bool compAsymptoteLikely = false;
              for(int k=0;k<numSing;k++) {
                ProjPair singPair = theMap.section().project( singPtr[k] );
                VecType s = singPair.first;
                VecType r0 = compSeg[0] - s;
                VecType r1 = compSeg[1] - s;
                if(r0[0] < safeDistances[k] || r1[0] < safeDistances[k])  //HARD-CODED to CR3BP position
                  compAsymptoteLikely = true;
              }
              //We don't need test this segment because it is likely within an asymptote
              if (compAsymptoteLikely) continue;

              for(int i=0; i<(int)currentLayerIds.size();i++) {
                  //Note:  both vectors are in reverse order
                  int tID = currentLayerIds[i];
                  ManifoldSegType& seg = theManifold.segments[tID];
                  //First test to see if we have a asymptote via vertical & long segment
                  if(seg.isNearlyVertical() && (seg.length(params)>settings.delta_max) ) continue;
                  //Test if an asymptote is likely by singularity interesection
                  bool currentOnAsymptote = false;
                  for(int k=0;k<numSing;k++) {
                    ProjPair singPair = theMap.section().project( singPtr[k] );
                    VecType s = singPair.first;
                    VecType r0 = seg[0] - s;
                    VecType r1 = seg[1] - s;
                    if(r0[0] < safeDistances[k] || r1[0] < safeDistances[k])  //HARD-CODED to CR3BP position
                      currentOnAsymptote = true;
                  }
                  //We don't need test this segment because it is likely within an asymptote
                  if (currentOnAsymptote) continue;

                  //The saddle-loop test for non-asymptotic segments
                  if (theManifold.segments[tID].isSaddleLoop(compMan.segments[cmID],epsSaddleLoop))
                  {

                      //Trigger that saddle-loop is found
                      overlapDetected = true;
                      //Store segID pairing in order of current segment
                      ccmMap.insert( std::pair<int,int>(cmID,tID) );
                  }
              }
            }

            //If we detected overlap, store connection and quit
            if (overlapDetected) {
              // Find mid point segments
              int numInGroup = (int) ccmMap.size();
              int idx = std::floor( numInGroup / 2 );
              std::map<int,int>::iterator it;
              int lSegID, clSegID, tally = 0;
              for(it=ccmMap.begin();it!=ccmMap.end();it++) {
                  if(tally == idx) {
                      clSegID = it->first;
                      lSegID = it->second;
                      //Exit the loop
                      break;
                  }
                  tally++;
              }
              //Crop manifolds at the midpoint segments
              mapManifolds[manID].crop(params,settings,lSegID);
              mapManifolds[compID].crop(params,settings,clSegID);

              //Signal that this manifold and it's complementary manifold are both complete
              mapManifolds[manID].isComplete = true;
              mapManifolds[manID].homoclinicPartnerID = compID;
              mapManifolds[compID].isComplete = true;
              mapManifolds[compID].homoclinicPartnerID = manID;
              //Prompt user
              if(params.debug) {
                    std::cout << " Stop Condition: Saddle-Loop Exit (strong) detected\n";
                    std::string iType = (theManifold.isForward()) ? "Unstable" : "Stable";
                    std::string cType = (mapManifolds[compID].isForward()) ? "Unstable" : "Stable";
                    std::cout << "   Saddle-Loop Intersection between :\n";
                    std::cout << "     Manifold " << manID << " (" << iType << ") | Segment "
                              << (int)theManifold.segments.size() - 1 << " [end]\n";
                    std::cout << "     Manifold " << compID << " ( " << cType <<") | Segment "
                              << (int)mapManifolds[compID].segments.size() - 1 << "\n";
              }
              //End compManifold loop
              break;
            }

          } //End loop through complementary manifolds
        } //End if length is still short enough
        */

        //6) Reached Maximum Depth: Stop when next workingSeg is beyond the maxTreeDepth
        if ( nextSegDepth > settings.maxTreeDepth ) {
            if(params.debug)
                std::cout << " Stop Condition: MaxTreeDepth : Next depth = " << nextSegDepth
                          << " > maxTreeDepth = " << settings.maxTreeDepth << "\n";
            theManifold.isComplete = true;
        }
        //6b) Major Subdivision has occurred and we have surpassed the msDepthLevel
        if ( theManifold.majorSubdivision ) {
            if ( nextSegDepth > theManifold.majorSubDepth ) {
                if(params.debug)
                    std::cout << " Stop Condition: Major Subdivision Depth Exit : Next depth = " << nextSegDepth
                              << " > majorSubDepth = " << theManifold.majorSubDepth << "\n";
                theManifold.isComplete = true;
            }
        }
        //6c) Force quit long orbits (most likely large manifold subdivision, cut computation time)
        if ( fabs(theManifold.getPeriod()) > 3 && nextSegDepth > 4 ) {
            if(params.debug)
                std::cout << " Stop Condition: Long Orbit Depth Exit : Next depth = " << nextSegDepth
                          << " > 4 AND Fixed point period = " << theManifold.getPeriod() << " > 3\n";
            theManifold.isComplete = true;

        }

    } else { //----------------------------------------------------------------------------------------------
        //The 'weak' manifolds may be KAM manifolds:
        //5)Saddle-Loop (or KAM) stopping conditions [Test the last 5 segments]
        //Skip this test if we have progressed too far already
        if (!theManifold.isComplete && theManifold.progressLength < 10.0) {
            //For each complementary manifold
            for(int j=0; j<(int)theManifold.complementaryIDs.size(); j++) {
                int compID = theManifold.complementaryIDs[j];
                ManifoldType& compMan = mapManifolds[compID];
                //Skip if complementary manifold is already done
                if(compMan.isComplete) {
                    continue;
                }
                //Last segment on the current manifold
                ManifoldSegType& lastSeg = theManifold.segments.back();

                //Check for an overlap of last segment with up-to the last 5 segments of complementary manifold:
                int segCount = 0;
                typename std::vector<ManifoldSegType>::reverse_iterator rit = compMan.segments.rbegin();
                for(; rit!=compMan.segments.rend() && segCount < 5; rit++) {
                    ManifoldSegType& compSeg = (*rit);
                    if( lastSeg.isSaddleLoop(compSeg,epsSaddleLoop) ) {
                        if(params.debug) {
                            std::cout << " Stop Condition: Saddle-Loop (or WeakMan KAM Manifold) crossing detected\n";
                            std::string iType = (theManifold.isForward()) ? "Unstable" : "Stable";
                            std::string cType = (compMan.isForward()) ? "Unstable" : "Stable";
                            std::cout << "   KAM Intersection between :\n";
                            std::cout << "     Manifold " << manID << " (" << iType << ") | Segment "
                                      << (int)theManifold.segments.size() - 1 << " [end]\n";
                            std::cout << "     Manifold " << compID << " ( " << cType <<") | Segment "
                                      << compSeg.segID << "\n";
                            //Run local test for debugging:
                            MSegIntersector lsx(lastSeg,compSeg);
                            std::cout << "     Seg 0 : " << lastSeg[0] << " , " << lastSeg[1] << "\n";
                            std::cout << "     Seg 1 : " << compSeg[0] << " , " << compSeg[1] << "\n";
                            std::cout << "     Test Result : " << lsx.getSolutionTypeString() << "\n";
                            std::cout << "          Angle = " << lsx.getAngle()*180.0/M_PI << " deg \n";
                        }
                        //Signal that this manifold and it's complementary manifold are both complete
                        mapManifolds[manID].isComplete = true;
                        mapManifolds[manID].homoclinicPartnerID = compID;
                        mapManifolds[compID].isComplete = true;
                        mapManifolds[compID].homoclinicPartnerID = manID;
                    }
                    if (theManifold.isComplete) {
                        break;
                    }
                    segCount++;
                }
                //End search through complementary manifolds if KAM crossing found.
                if (theManifold.isComplete) {
                    break;
                }
            }
        }

        //Also test for maximum tree depth but modified for "slow" progress.
        //The modification tracks the number of subdivision layers instead of just tree depth.
        //Note: This only applies if weak manifold has advected the minimum distance
        nextSegDepth = theManifold.getNumSubdivisionLayers(theManifold.workingSegID); //Counts itself as layer
        if ( beyondWeakMinThreshold && (nextSegDepth >= settings.maxTreeDepth) ) {
            if(params.debug) {
                std::cout << " Stop Condition: WeakMan MaxTreeDepth : Next depth = "
                          << nextSegDepth << " > maxTreeDepth = " << settings.maxTreeDepth << "\n";
            }
            theManifold.isComplete = true;
        }
    }

    //-----------------------------------------------------------
    //Estimate how close to completion
    //-----------------------------------------------------------
    if (theManifold.isComplete) {
        return 1.0;    //Done
    }
    double depthDone = std::min(1.0, ((double) theManifold.getDepth(theManifold.workingSegID) -1.0)/ (double) settings.maxTreeDepth);
    if (fabs(theManifold.getPeriod()) > 3) {
        depthDone = std::min(depthDone,((double) theManifold.getDepth(theManifold.workingSegID) -1.0) / 4.0);
    }
    double seedDone = std::min(1.0, theManifold.seedLength/settings.maxSeedArclength );
    double prog = std::max(depthDone,seedDone);
    if(beyondWeakMinThreshold && settings.enableMaxArclength) {
        double arcDone = std::min(1.0, theManifold.progressLength/settings.max_arc_length );
        prog = std::max(prog,arcDone);
    }
    if(beyondWeakMinThreshold && settings.enableMaxSeg) {
        double segsDone = std::min(1.0, (double)theManifold.numSegments/(double)settings.maxNumSeg);
        prog = std::max(prog,segsDone);
        //double seedSegsDone = std::min(1.0, (double)theManifold.numSeedSegments/(double)settings.maxNumSeedSeg);
        //prog = std::max(prog,seedSegsDone);
    }

    if(!beyondWeakMinThreshold) {
        //We have to get to the minimum arc length, so we can just scale the result:
        prog *= std::min(1.0, theManifold.progressLength / minLength);
    }

    //Return the estimated completion percent for this manifold
    return prog;
}


/*
///Compute Progenitor States of Manifolds that are advecting
template <class SRDATA, class FP>
template <class PMAP, class PARAM>
void ManifoldData<SRDATA,FP>::
evaluateProgenitorStates(
  PMAP& amap,
  PARAM& params
) {
    typedef orbital::ProgenitorCompData<PMAP>           ProgData;
    typedef orbital::ProgenitorIntersection<VecType>    ProgIntersect;
    typedef orbital::ProgenitorState<StateType,VecType> ProgStateType;
    typedef std::pair<int,TrajectoryType>               TrajJobPair;
    typedef std::vector<TrajJobPair>                    TrajJobStorageVector;
    typedef std::pair<int,ProgData>                     PDJobPair;
    typedef std::vector<PDJobPair>                      PDJobStorageVector;

    //Gather all segments that were generated during manifold computation:
    std::vector<std::pair<int,int> > msIDs;
    for(int manID=0; manID<(int)mapManifolds.size(); manID++ ) {
      const ManifoldType& theManifold = mapManifolds[manID];
      //Gather each segment except for the 0th segment on all manifolds:
      for(int segID=1;segID<(int)theManifold.segments.size();segID++) {
        msIDs.push_back( std::pair<int,int>(manID,segID) );
      }
    }

    //Over-arching data structures that sort computation possibilities
    std::vector<ProgData>  pData;
    std::map<nvis::ivec3,nvis::ivec3,nvis::lexicographical_order> duplicateProgComps;
    std::vector<nvis::ivec2> toInterpolate;

    //For each manifold that is requiring advection
    for (int i=0; i<(int)msIDs.size(); i++) {
        //Evaluate current depth level of working segment
        const int manID = msIDs[i].first;
        const int segID = msIDs[i].second;
        ManifoldType &thisMan = mapManifolds[manID];
        const int depth = thisMan.getDepth(segID);
        if ( depth > thisMan.depthPS ) { //Stop New Progenitor States after d=depthPS (per manifold)
            //Interpolate information & extract from manifolds :
            toInterpolate.push_back( nvis::ivec2(manID,segID) );

        } else { //Have to compute Progenitor State
            //Check if progenitor data for the first point already exists
            int lastSegID = thisMan.segments[segID].previous;
            if ( thisMan.previousConnection(segID)  //Check if this seg is connected to previous
                 && lastSegID != 0                  //Have to run the first point
               ) {
              //This point is already entered into a progenitor state computation
              //(OR already possesses progrenitor state information from last run),
              //so store it's duplication into a post-processing copy procedure:
              duplicateProgComps.insert(
                 std::pair<nvis::ivec3,nvis::ivec3>(
                         nvis::ivec3(manID,segID,0),     //The one to evaluate
                         nvis::ivec3(manID,lastSegID,1) )//The source of data
              ); //ToCopy->Source Pairing

              //Note: This data will be computed in computation, but we are trying
              //to avoid a duplicated computation

            }
            else { //New segment/point entirely
              //Store that we need a Progenitor Computation (and relevant info)
              ProgData dat0(thisMan.segments[segID][0],  //1st point of segment
                          manID, segID, 0);
              pData.push_back(dat0);
            }

            //Always add the second point of segment as Progenitor Computation:
            ProgData dat1(thisMan.segments[segID][1],  //2nd point of segment
                          manID, segID, 1);
            pData.push_back(dat1);
        }
    }
    //Update user on how many duplicate jobs are extracted
    //std::cerr << "\rPSGuess Propagation: Removed " << (int) duplicateProgComps.size() << " propagation jobs\n";

    //Submit progenitor propagation jobs
    int numPS = (int) pData.size();
    nvis::timer timer;
    //Per-thread cache
    TrajJobStorageVector *tjCache = new TrajJobStorageVector[nthreads];
    int propsComplete = 0;
    for (int j=0;j<numPS; j++) {
        #pragma omp task shared(amap,params,pData,tjCache,propsComplete)
        {
          //Propagate to assumed Progenitor State storing all internal points
          ManifoldType &thisMan = mapManifolds[pData[j].manID];
          TrajectoryType progArc;
          thisMan.propagateToProgenitorGuess(amap,params,j,pData[j].segID,pData[j].ptID,tjCache);
          #pragma omp atomic
          propsComplete++;
          double dt = timer.elapsed();
          std::ostringstream os;
          os << "\rPSGuess Propagation Task:" << propsComplete << "/" << numPS
             << " (" << (float)propsComplete*100/numPS << "%) in "
             << dt << "s. (" << (float)propsComplete/dt << "Hz)     \r";
          #pragma omp critical
          std::cout << os.str() << std::flush;

        }
    }
    //Extract progenitor propagation jobs from per-thread caches
    #pragma omp taskwait
    {
        //Store Trajectory object for each node
        for(int k=0;k<nthreads;k++) {
          for(int j=0;j<(int)tjCache[k].size();j++) {
            int jobID = tjCache[k][j].first;
            //Store the trajectory into ProgenitorCompData object
            pData[jobID].path = tjCache[k][j].second;
            //Note: if failure occurs, this could have blank data in trajectory.
          }
        }
    }
    //Done with trajectory storage
    for(int k=0;k<nthreads;k++) tjCache[k].clear();
    delete[] tjCache;

    //Have to make sure the Hamiltonian Check Stopping condition is off in Map Engine
    bool isHamCheck = amap.isCheckingHamiltonian();
    amap.stopAtHamiltonianError(false); //Need to disable for corrections

    //Localize the Actual Progenitor State by either
    // a)Finding the intersection or
    // b)Solving the TPBVP
    timer.restart();
    int psComplete = 0, numCorrections = 0, numSucCor = 0;
    for (int j=0;j<numPS; j++) {
        #pragma omp task shared(amap,params,pData,psComplete,numCorrections,numSucCor)
        {
            //Check if we have no data:
            if (pData[j].path.getNumStates() <=1) {
              //There is not enough data to evaluate the Progenitor State,
              //which is likely due to some form of propagation error

              //Indicate failed Progenitor State information
              pData[j].progState.found = false;

              //continue; //Can't do anything else with this PS job
            }
            else {
              //Compute intersections by comparing PeriodicOrbit and Trajectory objects
              std::vector<StateType> mSt, oSt;
              std::vector<double> mt, ot;
              const PeriodicOrbitType &theOrbit = mapManifolds[ pData[j].manID ].theOrbit;
              const bool isUnstable = mapManifolds[ pData[j].manID ].isForward();
              pData[j].path.findIntersection(theOrbit, mSt, oSt, mt, ot );

              //If an intersection exists,
              bool solveTPBVP = true;
              if ( ((int) mSt.size())>0 ) {
                  //Store the lowest delta-V intersection as the Progenitor State
                  std::list<ProgIntersect> iSectionList;
                  for(int ix=0;ix<(int)mSt.size();ix++) {
                      ProgIntersect psCandidate;
                      psCandidate.pos = VecType(mSt[ix][0],mSt[ix][1]); //HARD-CODED TO 2D!!!
                      psCandidate.mV = VecType(mSt[ix][3],mSt[ix][4]);  //Vx,Vy
                      psCandidate.oV = VecType(oSt[ix][3],oSt[ix][4]);  //Vx,Vy
                      psCandidate.tof = (isUnstable)? (pData[j].path.times.back()-mt[ix]) : mt[ix];
                      psCandidate.alpha = theOrbit.moduloAlpha(ot[ix] / theOrbit.getTimeOfFlight()); //On [0,1]
                      iSectionList.push_back( psCandidate );
                  }
                  //Sort the list : 1)delta-V, then 2)Time-of-flight
                  iSectionList.sort();
                  //Remove any intersection that has a poor velocity differential
                  typename std::list<ProgIntersect>::iterator lit;
                  for(lit=iSectionList.begin();lit!=iSectionList.end();lit++) {
                    bool feasible = true;
                    //Check that delta-v angle is feasible (<30deg)
                    if(lit->getDeltaVAngle() > 30.0*M_PI/180.0) feasible = false;
                    //Check that delta-v magnitude is feasible (<100m/s in EM | <10% nd velocity)
                    if(lit->getDeltaVMag() > 0.1) feasible = false;
                    //Remove if not feasible
                    if(!feasible) {
                      lit = iSectionList.erase(lit);
                      lit--;
                    }
                  }
                  if((int) iSectionList.size()>0) {
                    solveTPBVP = false;
                    //Store the first solution in the sort as the Progenitor State
                    //(should be thread-safe as new memory is not being allocated and pData is SHARED)
                    pData[j].progState.found = true;
                    ProgIntersect &psSol = (*(iSectionList.begin()));
                    pData[j].progState.alpha = psSol.alpha;
                    pData[j].progState.dv = psSol.getDeltaV();
                    StateType orbitState(0.0); //HARD-CODED to CR3BP planar problem.
                    orbitState[0] = psSol.pos[0]; orbitState[1] = psSol.pos[1];
                    orbitState[3] = psSol.oV[0]; orbitState[4] = psSol.oV[1];
                    pData[j].progState.orbitState = orbitState;
                    //Time of flight from PS to Manifold Point
                    pData[j].progState.tof = (!isUnstable)? -1.0*psSol.tof : psSol.tof;
                  }
              }

              if(solveTPBVP) {
                  #pragma omp atomic
                  numCorrections++;
                  //No intersections, so seed a guess with the closest neighboring points
                  StateType pG(0.0), oG(0.0);
                  double pdt, odt;
                  pData[j].path.findNearestNode(theOrbit,pG,pdt,oG,odt); //Gathers exact nodes
                  //Seed the guess
                  double alpha0 = theOrbit.moduloAlpha( odt / theOrbit.getTimeOfFlight() ); //On [0,1]
                  nvis::vec2 dv0(0.0); //HARD-CODED type for planar problem!
                  dv0[0] = pG[3] - oG[3]; //delta-xdot
                  dv0[1] = pG[4] - oG[4]; //delta-ydot
                  double tof0 = (isUnstable)? (pData[j].path.times.back()-pdt) : pdt;
                  ManifoldType &thisMan = mapManifolds[pData[j].manID];
                  VecType phif = thisMan.segments[pData[j].segID][pData[j].ptID];
                  //Solve with ProgenitorState corrections:
                  thisMan.solveProgenitorState(amap,alpha0,dv0,tof0,phif,pData[j].progState);
                  //Output to Progenitor Jobs (if it didn't work, indicate with flag)
                  if(pData[j].progState.found) {
                    #pragma omp atomic
                    numSucCor++;
                  }
              }
            } //End valid propagation check

            //Update user
            #pragma omp atomic
            psComplete++;
            double dt = timer.elapsed();
            std::ostringstream os;
            os << "\rPSComputation Task:" << psComplete << "/" << numPS
              << " (" << (float)psComplete*100/numPS << "%) in "
              << dt << "s (" << (float)psComplete/dt << "Hz) with "
              << numSucCor << "/" << numCorrections << " TPBVPs.  \r";
            #pragma omp critical
            std::cout << os.str() << std::flush;
        }

    }
    //Extract progenitor state information and write to MapManifoldSegments
    #pragma omp taskwait
    {
        for (int j=0;j<numPS; j++) {
            const int manID = pData[j].manID;
            const int segID = pData[j].segID;
            //ManifoldSegType& thisSeg = mapManifolds[manID].segments[segID];
            if(pData[j].ptID>0) {
                mapManifolds[manID].segments[segID].ps1 = pData[j].progState;
            } else {
                mapManifolds[manID].segments[segID].ps0 = pData[j].progState;
            }
        }

        //Fill out the duplicates now that information is already available
        std::map<nvis::ivec3,nvis::ivec3>::iterator mit;
        for(mit=duplicateProgComps.begin();mit!=duplicateProgComps.end();++mit) {
          nvis::ivec3 psIdx = mit->first;
          nvis::ivec3 source = mit->second;
          //Set the progenitor state of the FIRST point from a source computation
          ProgStateType sourceProgState = mapManifolds[source[0]].segments[source[1]].ps1; //A previous seg
          mapManifolds[psIdx[0]].segments[psIdx[1]].ps0 = sourceProgState; //Copy data into object
        }

        //Perform interpolation on remaining data
        for(int k=0;k<(int)toInterpolate.size();k++) {
            const int manID = toInterpolate[k][0];
            const int segID = toInterpolate[k][1];
            //Perform Interpolation from parent segment data:
            ManifoldType& thisMan = mapManifolds[manID];
            const int parentID = thisMan.segments[segID].parent;
            //For Point 0
            const double tau0 = thisMan.segments[segID].tau0;
            thisMan.segments[segID].ps0.found = thisMan.isPSFound(parentID,tau0);
            thisMan.segments[segID].ps0.alpha = thisMan.getPSAlpha(parentID,tau0);
            thisMan.segments[segID].ps0.dv = thisMan.getPSDeltaV(parentID,tau0);
            thisMan.segments[segID].ps0.orbitState = thisMan.getPSOrbitState(parentID,tau0);
            thisMan.segments[segID].ps0.tof = thisMan.getTimeOfFlight(segID, 0.0);
            //For Point 1
            const double tau1 = thisMan.segments[segID].tau1;
            thisMan.segments[segID].ps1.found = thisMan.isPSFound(parentID,tau1);
            thisMan.segments[segID].ps1.alpha = thisMan.getPSAlpha(parentID,tau1);
            thisMan.segments[segID].ps1.dv = thisMan.getPSDeltaV(parentID,tau1);
            thisMan.segments[segID].ps1.orbitState = thisMan.getPSOrbitState(parentID,tau1);
            thisMan.segments[segID].ps1.tof = thisMan.getTimeOfFlight(segID, 1.0);
        }
    }
    //Re-enable Hamiltonian check in Map Engine if necessary
    if (isHamCheck) amap.stopAtHamiltonianError(isHamCheck);
}*/


/// Read manifolds
template <class SRDATA, class FP>
void ManifoldData<SRDATA,FP>::
read(const char* filename)
{
    //Assumes you have the correct FixedPointData loaded (part of constructor)

    FILE* f = fopen(filename, "r");
    if(!f) {
        std::cerr << "ManifoldData:  Unable to read file = " << filename << "!\n";
        throw std::runtime_error("Bad file while reading ManifoldData");
    }

    //Header
    char buffer[200];
    for(int i=0; i<2; i++) {
        fgets(buffer, 200, f);
    }

    int numManifolds = 0;
    double ham = 0.0;
    fscanf(f,"%d %lf",&numManifolds,&ham);
    //Run a test to make sure Hamiltonian values match
    if (ham != theFixedPointData->getJacobiConstant() ) {
        std::cerr << "ManifoldData: Hamiltonian value from file (" << ham
                  << ") does NOT match FixedPointData input (" << theFixedPointData->getJacobiConstant() << ")\n";
        throw std::runtime_error("ManifoldData: Hamiltonian mismatch!");
    } //Could also check FixedPointData


    //Read in manifolds (and segments)
    for (int k=0; k<numManifolds; k++) {
        //Create a Manifold object with a read() command
        mapManifolds.insert( std::pair<int,ManifoldType>(k, ManifoldType()) );
        ManifoldType& loadedManifold = mapManifolds[k];
        //Utilize the read() command in MapManifold object to set quantities
        loadedManifold.read(f,theFixedPointData);
    }

    //Once all read, set parental pointers

    fclose(f);
}

/// Write manifolds
template <class SRDATA, class FP>
bool ManifoldData<SRDATA,FP>::
write(const char* filename, const char* fpDataFileName) const
{
    FILE* f = fopen(filename, "w");

    if(!f) {
        std::cerr << "ManifoldData:  Unable to write to file = " << filename << "!\n";
        return false;
    }

    // Write header
    fprintf(f,"# Manifold Data for STHManifold Algorithm \n");

    //Write the file indicating the FixedPointData
    fprintf(f,"FixedPointDataFile= %s\n",fpDataFileName);
    //Write Hamiltonian value to ensure solution agreement
    int numManifolds = (int) mapManifolds.size();
    fprintf(f,"%d %.15f\n\n",numManifolds,theFixedPointData->getJacobiConstant());

    //Write manifold data and segments utilizing member write() function
    typename std::map<int,ManifoldType>::const_iterator it;
    for (it=mapManifolds.begin(); it!=mapManifolds.end(); ++it) {
        it->second.write(f);
    }

    fclose(f);

    return true;
}

/// Compute homoclinic/heteroclinic connections between the s/u manifolds of this object
template <class SRDATA, class FP>
void ManifoldData<SRDATA,FP>::computeConnections()
{
    typedef std::pair<int,int>       IntPair;
    typedef std::vector<ManifoldCon> ManConVec;
    //Assuming the manifold data is already available
    std::vector<IntPair> usManPairs;
    //Construct all unstable-stable possibility pairings
    int numManifolds = (int) mapManifolds.size();
    for(int u=0; u<numManifolds; u++) {
        //For all unstable manifolds
        if ( mapManifolds[u].isForward() ) {
            //Loop through all STABLE manifolds
            for(int s=0; s<numManifolds; s++) {
                //Pair to an unstable for an intersection test
                if(!(mapManifolds[s].isForward())) {
                    usManPairs.push_back( IntPair(u,s) );
                }
            }
        }
    }

    //Test for manifold segment intersections on all pairings
    int numPairs = (int) usManPairs.size();
    ManConVec* mcPerThread = new ManConVec[nthreads];
    #pragma omp parallel for schedule(dynamic,1)
    for(int i=0; i<numPairs; i++) {
        int tid = 0;
#if _OPENMP
        tid = omp_get_thread_num();
#endif
        int uID = usManPairs[i].first;
        int sID = usManPairs[i].second;
        ManifoldType& uMan = mapManifolds[uID];
        ManifoldType& sMan = mapManifolds[sID];
        //Note: there could be multiples of the same transfer as we are just looking at
        //the crossings on these manifold objects.

        for(int uSegID = 0; uSegID<(int)uMan.segments.size(); uSegID++) {
            //For each stable manifold
            for (int sSegID=0; sSegID<(int)sMan.segments.size(); sSegID++) {
                double uTau=0.0,sTau=0.0;
                bool isf = false;
                isf = uMan.segments[uSegID].getIntersection(sMan.segments[sSegID],uTau,sTau);
                if(isf) {
                    //Store to cache
                    ManifoldCon mc;
                    setManifoldCon(uID,uSegID,uTau,sID,sSegID,sTau,mc);
                    mcPerThread[tid].push_back(mc);
                }
            }
        }
    }

    //Unroll the cache into Connection data
    hhConnections.clear();
    for(int k=0; k<nthreads; k++) {
        for(int i=0; i<(int)mcPerThread[k].size(); i++) {
            hhConnections.push_back( mcPerThread[k][i] );
        }
    }
    for(int k=0; k<nthreads; k++) {
        mcPerThread[k].clear();
    }
    delete[] mcPerThread;

}

/// Compute augmented connections between s/u manifolds of this object
template <class SRDATA, class FP>
void ManifoldData<SRDATA,FP>::computeConnections(const typename ManifoldData<SRDATA,FP>::VecType& maxMapDisp)
{}

/// Compute connections with another ManifoldData object
//void computeConnections(ManifoldData& other);

/** Write ManifoldConnection data to file (the .imc extension)
 *  This requires a known .im file that stores current ManifoldData, but can also
 *  store multiple ManifoldData names if you are computing connections through
 *  multiple Jacobi constants at once.  (It doesn't have to be one at a time.)
 */
template <class SRDATA, class FP>
bool ManifoldData<SRDATA,FP>::writeConnections(const char* filename, std::vector<std::string>& manDataFileName) const
{
    return false;
}

/** Read ManifoldConnection data (.imc file extension)
 *  Note: This is also a bit tricky because the corresponding ManifoldData
 *  (and secondary ManifoldData files) are also required to be preloaded in order
 *  to interact with the connection data.
 */
template <class SRDATA, class FP>
bool ManifoldData<SRDATA,FP>::readConnections(const char* filename)
{
    return false;
}

/// Set a ManifoldConneciton given intersection information and two manifoldIDs from THIS object
template <class SRDATA, class FP>
void ManifoldData<SRDATA,FP>::
setManifoldCon(const int uID, const int uSegID, const double ut,
               const int sID, const int sSegID, const double st,
               typename ManifoldData<SRDATA,FP>::ManifoldCon& mc)
{
    mc.manifoldDataIdx = std::pair<int,int>(0,0);
    //Make a different function for different ManifoldData (i.e., different Jacobi constants)

    mc.uMan = std::pair<int,int>(uID,uSegID);
    mc.uTau = ut;
    mc.sMan = std::pair<int,int>(sID,sSegID);
    mc.sTau = st;
    //Connection map position at intersection (there are likely more than just this one)
    mc.pos = mapManifolds[uID].segments[uSegID].getPoint(ut);

    //Look up the base orbit index to determine if this is a homoclinic connection:
    mc.isHomoclinic = (mapManifolds[uID].fpdOrbitIdx == mapManifolds[sID].fpdOrbitIdx);

    //First check if progenitor states are available (they should be)
    bool uPSFound = mapManifolds[uID].isPSFound(uSegID,ut);
    bool sPSFound = mapManifolds[sID].isPSFound(sSegID,st);

    if (uPSFound && sPSFound) {
        //Time of Flight (PSUnstable to PSStable)
        mc.tof = mapManifolds[uID].getTimeOfFlight(uSegID,ut); //Should be positive
        mc.tof += fabs( mapManifolds[sID].getTimeOfFlight(sSegID,st) ); //Likely negative
        //Departure Delta-V (PSUnstable)
        mc.depDV = mapManifolds[uID].getPSDeltaV(uSegID,ut);
        //Arrival Delta-V (PSStable)
        mc.arrDV = mapManifolds[sID].getPSDeltaV(sSegID,st);
    } else {
        mc.tof = -1.0;
        mc.depDV = VecType(0.0);
        mc.arrDV = VecType(0.0);
    }

    //Connection Delta-V (which is zero for heteroclinic and homoclinic connections)
    mc.deltaV = VecType(0.0);
    //Note: Make another version of this for augmented connections (delta-V assisted)
}

}// end pmate


#endif  // MANIFOLD_MASTER_CONTROL_HPP
