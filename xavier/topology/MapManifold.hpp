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


#ifndef MAP_MANIFOLD_HPP
#define MAP_MANIFOLD_HPP

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <queue>
// boost for universal unique identifiers
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
// PMATE API
#include <math/fixed_vector.hpp>
#include <math/bezier_spline.hpp>
#include <math/bounding_box.hpp>
#include <maps/metric.hpp>
#include <maps/fixpoints.hpp>
#include <maps/map_analysis.hpp>
#include <pmate/FixedPointData.hpp>
#include <design/Trajectory.hpp>
#include <design/PeriodicOrbit.hpp>
#include <orbital/controller.hpp>
#include <orbital/progenitor.hpp>
#include <topology/EdgeRotationFailure.hpp>
#include <topology/SectionTransversality.hpp>
#include <topology/ManifoldClasses.hpp>
#include <topology/MapManifoldSegment.hpp>
#include <topology/invariant_manifold.hpp>
#include <topology/STHManifold.hpp>

//OpenMP
#if _OPENMP
#include <omp.h>
#endif

using namespace xavier;

namespace topology {

/// Defining how to handle a MapManifold child - TODO: Useful in Future!
enum ManifoldChildType {
    ADULT = 0,  //Not a child of any manifold
    MINOR,      //Child defined by a sub-return of parent (n < fp.period)
    MIRROR,     //Child defined by mirror theorem:
    // - UNSTABLE_PLUS  => STABLE_PLUS
    // - UNSTABLE_MINUS => STABLE_MINUS
    ADULTMINOR, //Minor Child that grew up
    ADULTMIRROR //Mirror Child that grew up
};

/// Manifold Object containing all segments that make a particular manifold
template<class SRDATA,class FP>
class MapManifold {
public:
    typedef MapManifold<SRDATA,FP>                                   SelfType;
    typedef SRDATA                                                   SortableData;
    typedef typename SRDATA::VecType                                 VecType;
    typedef typename SRDATA::MapDataType                             MapDataType;
    typedef typename FP::StateType                                   StateType;
    typedef pmateDesign::Trajectory<StateType>                       TrajectoryType;
    typedef pmateDesign::PeriodicOrbit<StateType>                    PeriodicOrbit;
    typedef EdgeRotationFailure<VecType>                             MapDiscont;
    typedef MapManifoldSegment<VecType,StateType,MapDataType>        ManifoldSeg;
    typedef std::vector< ManifoldSeg >                               SegmentVector;
    //Adaptive Edge data structure
    typedef MapTaskOutputData<boost::uuids::uuid,SRDATA,MapDiscont>  ManMapStorageData;
    typedef std::vector<ManMapStorageData>                           MMStorageVector;
    //Storage for (uuid,seedID) pairs for jobs
    typedef AdaptiveEdge<VecType,VecType>                            EdgeType;
    typedef typename EdgeType::IdType                                IdType;
    typedef std::pair<int,IdType>                                    JobIDPair;
    typedef MapTaskOutputData<JobIDPair,SRDATA,MapDiscont>           MapJobStorageData;
    typedef std::vector<MapJobStorageData>                           MJStorageVector;
    //Storage for Progenitor State Objects
    typedef std::pair<int,TrajectoryType>                            TrajJobPair;
    typedef std::vector<TrajJobPair>                                 TrajJobStorageVector;
    typedef orbital::ProgenitorState<StateType,VecType>              ProgState;
    
    /// Constructor
    MapManifold();
    
    ///Constructor useful for a copy
    MapManifold(Perturbation dir, ManifoldChildType mct = ADULT);
    
    ///Constructor with Info (Not very good for containers
    MapManifold(const int orbitID, const int fpID, const FP& fp,
                Perturbation dir, ManifoldChildType mct = ADULT);
                
    /// Copy Constructor
    MapManifold(const SelfType& other);
    
    /// ID for working segment location
    int workingSegID;
    /// Is the data for the working segment initialized
    bool workingSegInitialized;
    /// ID for current downstream segment
    int lastSegID;
    /// Total accumulated Arc-length of all segments
    double totalArcLength;
    /// Accumulated progress length
    double progressLength;
    /// Seeding segment accumulated arc length
    double seedLength;
    /// Number of total segments computed for this manifold
    int numSegments;
    /// Number of seeding segments utilized for advection
    int numSeedSegments;
    /// Maximum depth level where Progenitor States are computed
    int depthPS;
    /// Flag indicating the first subdivision has occurred (weakStrength manifolds)
    bool firstSubdivision;
    /// Flag indicating a major subdivision (more than 1000 segments generated) at current depth
    bool majorSubdivision;
    /// Depth at major subdivision
    int majorSubDepth;
    /// Flag indicating that we need to retry depth<=2 (zeroth or first new depth level) with a fixed dTauMin
    bool noNewSegsAtNewDepth;
    /// Enable input delta_tau_min to be a fixed dTauMin metric
    bool useFixedTauMin;
    
    /// Flag tracking if the manifold satisfies advection stoppage criteria
    bool isComplete;
    
    // Mother-Daughter manifold data/references ------------------------------
    /// Pointer to Parent MapManifold with NULL representing no parent
    int parent;
    /// Current Manifold Index
    int manifoldID;
    /// Does this manifold have a parent
    bool hasParent() const
    {
        return (parent>=0)? true : false;
    }
    /// Container for pointers to children
    std::vector<int> children;
    /// Type of child
    ManifoldChildType childType;
    /// Complementary manifold indexes - the opposite manifolds (stable v. unstable) for KAM tests
    std::vector<int> complementaryIDs;
    /// Manifold ID of a homoclinic pairing (or saddle-loop on the section). [-1 for no saddle-loop]
    int homoclinicPartnerID;
    /// Storage for segment IDs at depth (for quick reference)
    std::map<int,std::vector<int> >  idsAtDepth;
    /// Next manifold correlating to a single mapping (downstream) which belongs to another return of periodic orbit
    int nextManifold;
    /// Previous manifold correlating to a single mapping (upstream) which belongs to another return of periodic orbit
    int prevManifold;
    //-------------------------------------------------------------------------
    
    /// Pointer to base-orbit fixed point for this manifold
    FP baseFixedPoint;
    /// Orbit Index in FixedPointData
    int fpdOrbitIdx;
    /// Point Index in fixed point chain in FixedPointData
    int fpdPointIdx;
    /// The full periodic orbit history of states
    PeriodicOrbit theOrbit;
    /// Type of perturbation, aka Manifold type
    Perturbation mType;
    /// Indicator of saddle (fixed point) strength : 'true' means lambda_max < 10
    bool weakStrength;
    /// Number of times the first step is tried (iteratively)
    int firstStepCount;
    
    ///Container for segments which are ordered in sequential progress
    SegmentVector  segments;
    
    ///AdaptiveEdge object for working segment subdivision
    EdgeType  adaptiveEdge;
    ///Number of subdivisions for the current working segment
    int numSubdivisions;
    ///Number of empty nodes after subdivision (helps indicate if subdivision remains)
    int emptyNodes;
    
    // Member Functions : --------------------------------------------------------
    /// Int to Perturbation
    enum Perturbation intToPert(const int& i) const;
    /// Perturbation to Int
    int pertToInt(const enum Perturbation& pert) const;
    /// PertInt to string
    std::string pIntToString(const int& i) const;
    /// Int to Child
    enum ManifoldChildType intToChild(const int& i) const;
    /// Child to int
    int childToInt(const enum ManifoldChildType& mct) const;
    /// ChildInt to string
    std::string cIntToString(const int& i) const;
    
    // First Segment processing during buildManifolds() ---------------------------------------
    /// Autonomous first step from function to compute step-off from fixed point (f(lambda_max))
    void autoInitStep(const double& mup, const double& lambda,
                      double& eps, double& epsMax, bool superAlignment);
                      
    /// Compute the coefficients to the auto init step function
    void getAutoStepCoeff(const double& mup, double& a, double& b, double& c, double& d) const;
    
    ///Create an initial step off of fixed point (for parents only)
    template<class MAP, class PARAM>
    VecType firstSegmentStep(const MAP& theMap, const PARAM& params,const ManifoldSettings& settings);
    
    ///Check that the current step is satisfactory (based on ManifoldSettings).
    //May overwrite 'p" with new step!
    template<class MAP, class PARAM>
    bool firstSegmentCheck(VecType& p, const VecType& map_of_p,
                           const MAP& theMap,const PARAM& params,const ManifoldSettings& settings);
                           
    /**Check that the current step is satisfactory (based on ManifoldSettings)
     *  - updates steps by dividing by 2 each attempt
     *  - increases the decrease rate of the initial step size
     *  - stops after 6 tries -> (1e-9/2^6 = 1.5625e-11)
     */
    //May overwrite 'p" with new step!
    template<class MAP, class PARAM>
    bool firstSegmentCheckDivisor(VecType& p, const VecType& map_of_p,
                                  const MAP& theMap,const PARAM& params,const ManifoldSettings& settings);
                                  
    /// Force the first segment step to be minimum linear step
    template<class MAP, class PARAM>
    void firstSegmentForcedStep(VecType& p, const MAP& theMap, const PARAM& params,
                                const ManifoldSettings& settings);
                                
    ///Store the first segment once it passes checks
    template<class MAP, class PARAM>
    void storeFirstSegment(const VecType& p, std::set<SRDATA>& cache, const MAP& theMap,
                           const PARAM& params,const ManifoldSettings& settings);
    //--------------------------------------------------------------------------------------------------------
    
    // Propagation functions ----------------------------------------------------------
    ///Propagation call for point on this manifold utilizing data storage (taskStorage & dataCache)
    template<class MAP, class PARAM>
    void propagateManifoldArcTask(const VecType& x0, const MAP& theMap, const PARAM& params,
                                  std::set<SRDATA>& cache, std::set<SRDATA>& backCache, MMStorageVector* mapDataCache);
    ///Propagation call for node point (IdType) on manifold working segment utilizing task data storage
    template<class MAP, class PARAM>
    void propagateManifoldArcTask(const JobIDPair& jobID, const MAP& theMap, const PARAM& params,
                                  std::set<SRDATA>& cache, std::set<SRDATA>& backCache, MJStorageVector* jobDataCache);
                                  
    ///Propagation call for point on this manifold (returns internal points)
    // ->Typically for interactive propagation and MissionDesign guess construction
    //void propagateManifoldArc(const VecType& x0,...);
    //void propagateManifoldArc(const int& segID, const double& tau, ...);
    //--------------------------------------------------------------------------------------------------------
    
    
    // Processing Manifold Segments utilizing STHManifold algorithm (curve-refinement)------------------------
    
    /// Get propagation tasks for the current adaptive edge
    void getEdgeNodeMapTasks( const int& manID, std::vector<JobIDPair>& jobs );
    
    /// Assign propagation data from a map call to an edge node
    template<class MAP, class PARAM>
    void assignEdgeNodeData(const IdType& nodeID, const SRDATA& result,
                            const MapDiscont& mapDis, bool fail,
                            const MAP& theMap, const PARAM& params,
                            std::set<SRDATA>& cache, std::set<SRDATA>& backCache);
                            
    /// Check and refine the current segments in the AdaptiveEdge
    template<class MAP, class PARAM>
    void refineSegmentsTest(std::set<SRDATA>& cache, const MAP& theMap,
                            const PARAM& params, const ManifoldSettings& settings);
                            
    /// Use STHManifold technique on current edge (using queue)
    template<class MAP, class PARAM>
    void advectSegment(const MAP& theMap, const PARAM& params,
                       const ManifoldSettings& settings, std::set<SRDATA>& cache);
                       
    /// Get a delta_tau_min parameter based on current working segment edge length and allowed map displacement u_min
    template<class PARAM>
    double getAdaptiveDeltaTauMin(const PARAM& params, const double& u_min) const;
    
    /// Add any children to parent processing queue if they grow up
    void childMaturity(std::queue<int>& processingQueue);
    
    //--------------------------------------------------------------------------------------------------------
    
    ///Write this manifold object to file
    void write(FILE* f) const;
    
    ///Read this manifold object from a file
    void read(FILE* f, pmate::FixedPointData<FP>* fpDataPtr);
    
    ///Operator to check if Manifold objects are the same
    bool operator==(SelfType const& rhs) const
    {
        return id == rhs.id;
    }
    ///Testing to see if thread data comes from this manifold
    bool isSameManifold(const boost::uuids::uuid& oid) const
    {
        return id == oid;
    }
    
    /// Set period values once fixed point is set
    void setPeriod()
    {
        period = baseFixedPoint.K;
        thePeriod = (fwd)? period : -period;
        //Note:  We need to double period for any orbit that has a negative lambda_max
        if (baseFixedPoint.eval[1] < 0.0) {
            thePeriod *= 2.0;
        }
        //Such manifolds switch to their opposite direction subspace (like a Mobius strip)
    }
    
    ///Get the period (number of crossings of hyperplane) for this source orbit
    int getPeriod() const
    {
        return period;
    }
    ///Get the period utilized to compute the manifold (can be negative or double getPeriod())
    int getThePeriod() const
    {
        return thePeriod;
    }
    
    /// Is this manifold using forward mappings
    bool isForward() const
    {
        return fwd;
    }
    
    /// Print base-level information about this manifold
    void info() const;
    
    /// Print an entry for an manifold generation update table (for debugging)
    void printUpdateEntry(); // const;
    
    /// Get a list of all Map Discontinuities from segments
    void getMapDiscontinuityList(std::list<MapDiscont>& mdList) const;
    
    /// Get the current working segment from list
    ManifoldSeg& getWorkingSegment();
    
    /// Test if a segment exists
    bool segmentExists(const int segID) const;
    /// Test if next segment from given segment ID is available
    bool nextSegmentExists(const int segID) const;
    /// Test if the previous segment exists
    bool previousSegmentExists(const int segID) const;
    /// Test if there is a connection between this segment and the next
    bool nextConnection(const int segID) const;
    /// Test if there is a connection between this segment and the previous
    bool previousConnection(const int segID) const;
    
    /// Crop segments AFTER a given segment id
    template<class PARAM>
    void crop(const PARAM& params, const ManifoldSettings& settings, const int segID);
    
    /// Clear segments and reset manifold like new
    void clearManifold();
    
    //----------------------------------------------------------------------------------------------
    //Segment-Tree Traversal functions
    /** Get the accumulated UPSTREAM sum of a data value (recursive).
     *  Typically useful for accumulated time-of-flight. (Note, doesn't fully work for tof!)
     *  Inputs:
     *    segID =  Segment Index of data request
     *    tau =    Linear parameter of data request of current segment on [0.0,1.0]
     *    dataID = Index of requested data from ExtendedMapDataVec
     *  Note, the perturbation type determines value (Stable manifolds have delta-t < 0).
     */
    double getDataSum(const int segID, const double tau, const int dataID) const;
    
    /** \brief Get the time of flight
     *  Get the time of flight for a manifold point.  Time of flight is assumed to be
     *  in ExtendedMapDataVec at dataID=2.
     *  Inputs:
     *    segID =  Segment Index of data request
     *    tau =    Linear parameter of current segment(segID) on [0.0,1.0]
     *  [OFF] The Progenitor State time of flight for
     *  depth levels < depthPS are utilized to approximate the total value.
     *  Note:  Stable manifolds have tof < 0.
     */
    double getTimeOfFlight(const int segID, const double tau) const; //Invokes progenitor State
    
    /// Get the minimum UPSTREAM data value (recursive)
    double getMinData(const int segID, const double tau, const int dataID) const;
    
    /// Get the depth of current segment or whole manifold tree (segment 0 is on depth = 1)
    int getDepth(const int segID=-1) const;
    
    /// Get the actual number of layers of subdivision (skipping single-segment layers)
    int getNumSubdivisionLayers(const int segID = -1); // const;
    
    /// Get the linear parameter for indicated UPSTREAM segment or initial segment (recursive)
    double getSourceLinearParameter(const int segID, const double tau, const int stopID=-1) const;
    
    /// Track the destinations of a source segment
    void trackSourceDownstream(const int sourceID, const int segID,
                               std::vector<int>& downstreamIDs,
                               std::vector< std::pair<double,double> >& tauAtSource) const;
                               
    /// Get UPSTREAM map points utilizing linear interpolation from manifold (in order 0-level -> currentDepth)
    void getUpstreamPoints(const int segID, const double tau, std::vector<VecType>& mapStates) const;
    /// Get DOWNSTREAM map points utilizing linear interpolation from manifold
    void getDownstreamPoints(const int segID, const double tau, std::vector<VecType>& mapStates) const;
    /// Get ALL map points for a selected arc (both upstream and downstream, in tree-depth order)
    void getAllMapPoints(const int segID, const double tau, std::vector<VecType>& mapStates) const;
    
    /// Gather all segments at the current depth level
    void getAllSegmentIDsAtDepth(std::vector<int>& ids, const int depth = -1); // const;
    /// Gather all segments at the depth level of given segment
    void getAllSegmentIDsOnSameDepth(std::vector<int>& ids, const int segID = -1); // const;
    
    /** Propagate UPSTREAM to top of manifold tree employing interpolated points (i.e., to depth=0)
     *  Note, output is actually in order from depth=0 to depth=currentDepth since propagation
     *  is all downstream to conform with method theory and assumptions.
     */
    template<class PMAP>
    void propagateUpstream(PMAP& theMap, const int segID, const double tau,
                           std::vector<typename PMAP::state_type>& states,
                           std::vector<double>& times ) const;
    /** Propagate UPSTREAM to top of manifold tree employing interpolated points
     *  at each p-iterate (i.e., to depth=0)
     *  This version also returns the iterates (p) and sub-iterates between each p-th iterate
     *  for reference (usually visualization).
     */
    template<class PMAP>
    void propagateUpstream(PMAP& theMap, const int segID, const double tau,
                           std::vector<typename PMAP::state_type>& states, std::vector<double>& times,
                           std::vector<typename PMAP::lvec_type>& iters,   std::vector<double>& itTimes,
                           std::vector<typename PMAP::lvec_type>& subIts,  std::vector<double>& subitTimes) const;
                           
    /// Propagate DOWNSTREAM to bottom of manifold tree employing interpolated points (i.e., to depth=d_max)
    template<class PMAP>
    void propagateDownstream(PMAP& theMap, const int segID, const double tau,
                             std::vector<typename PMAP::state_type>& states,
                             std::vector<double>& times ) const;
    /** Propagate DOWNSTREAM to bottom of manifold tree employing interpolated points
     *  at each p-iterate (i.e., to depth=d_max)
     *  This version also returns the iterates (p) and sub-iterates between each p-th iterate
     *  for reference (usually visualization).
     */
    template<class PMAP>
    void propagateDownstream(PMAP& theMap, const int segID, const double tau,
                             std::vector<typename PMAP::state_type>& states, std::vector<double>& times,
                             std::vector<typename PMAP::lvec_type>& iters,   std::vector<double>& itTimes,
                             std::vector<typename PMAP::lvec_type>& subIts,  std::vector<double>& subitTimes) const;
                             
    /** Propagate DOWNSTREAM to bottom of manifold tree and number of elements beyond
     *  employing interpolated points at each p-iterate (i.e., to depth=d_max).
     *  This version also returns the iterates (p) and sub-iterates between each p-th iterate
     *  for reference (usually visualization).
     */
    template<class PMAP>
    void propagateDownstream(PMAP& theMap, const int segID, const double tau, const int extendedIterates,
                             std::vector<typename PMAP::state_type>& states, std::vector<double>& times,
                             std::vector<typename PMAP::lvec_type>& iters,   std::vector<double>& itTimes,
                             std::vector<typename PMAP::lvec_type>& subIts,  std::vector<double>& subitTimes) const;
                             
                             
    //----------------------------------------------------------------------------------------------
    // Progenitor State Functions:  (THESE ARE CURRENTLY OFF!  NEEDS FUTURE WORK!!!!)
    //----------------------------------------------------------------------------------------------
    /// Is the requested progenitor state available (meaning viable solutions on both endpoints of segment)
    bool isPSFound(const int segID, const double tau) const;
    
    /// Get the Progenitor State DeltaV (Arrival or Departure)
    //VecType getPSDeltaV(const int segID, const double tau) const;
    /// Get the Progenitor State Alpha Parameter (location on orbit)
    //double getPSAlpha(const int segID, const double tau) const;
    /** \brief Get the Progenitor State Orbit State (full state on orbit)
     *  Linearly interpolate the orbit state at the progenitor location.
     *  Note: Using alpha and a pre-computed orbit is safer.
     */
    //StateType getPSOrbitState(const int segID, const double tau) const;
    
    /** \brief Use Shooting to solve the TPBVP for the progenitor state
     *  Use an implicit single shooting process to solve for (alpha,dV,TOF) values to
     *  indicate the progenitor state (manifold origin point on an orbit) given a
     *  DOWNSTREAM manifold state.
     *  Inputs:
     *    1) Map engine (theMap)
     *    2) Alpha (Guess)
     *    3) dV (Guess - 2D)
     *    4) TimeOfFlight (Guess) from Progenitor State to Manifold Point (+/-)
     *    5) Manifold Point (phif)
     *    6) Convergence Tolerance (default: 1e-8)
     *    7) Verbose flag (default: false)
     *    8) Max number of iterations (default: 20)
     *  Output:
     *    Progenitor State Object that holds result
     */
    /*template<class PMAP>
    bool solveProgenitorState(PMAP& theMap, const double& alpha0, const nvis::vec2& dv0,
                  const double& dt0, const VecType& phif, ProgState& output,
                  const double  eps = 1.e-10, const bool verbose = true, const int maxiter=20) const;
    
    /// Propagate to the estimated Progenitor State given a segment index and point index
    template<class PMAP,class PARAM>
    void propagateToProgenitorGuess(PMAP& theMap, PARAM& theMapParams, const int jobID,
                  const int segID, const int ptID, TrajJobStorageVector* tjCache);
    
    /// Propagate to the Prgenitor State (upstream) given a segment index (also returns iterates)
    template<class PMAP>
    void propagateToProgenitorState(PMAP& theMap, const int segID, const double tau,
                                    std::vector<typename PMAP::state_type>& states,
                                    std::vector<double>& times,
                                    std::vector<typename PMAP::return_type>& iterData) const;   */
    //----------------------------------------------------------------------------------------------
    
private:
    boost::uuids::uuid id; //Note sure I need this!
    ///Period of fixed point
    int period;
    //Propagation direction (stable or unstable)
    bool fwd;
    //The actual propagation period for map calls (takes fwd into account)
    int thePeriod;
    //Plus or minus perturbation flag
    bool p_or_m;
    
    
};

///Constructor
template<class SRDATA,class FP>
MapManifold<SRDATA,FP>::
MapManifold(const int orbitID, const int fpID, const FP& fp, Perturbation dir, ManifoldChildType mct) :
    workingSegID(-1),
    workingSegInitialized(false),
    lastSegID(-1),
    totalArcLength(0.0),
    progressLength(0.0),
    seedLength(0.0),
    numSegments(0),
    numSeedSegments(0),
    depthPS(4),
    firstSubdivision(false),
    majorSubdivision(false),
    majorSubDepth(-1),
    noNewSegsAtNewDepth(true),
    useFixedTauMin(false),
    isComplete(false),
    parent(-1),
    manifoldID(-1),
    childType(mct),
    homoclinicPartnerID(-1),
    nextManifold(-1),
    prevManifold(-1),
    baseFixedPoint(fp),
    fpdOrbitIdx(orbitID),
    fpdPointIdx(fpID),
    mType(dir),
    weakStrength(false),
    firstStepCount(0),
    adaptiveEdge(fp.pos,fp.pos),
    numSubdivisions(0),
    emptyNodes(0),
    id( boost::uuids::random_generator()() ),
    period(fp.K),
    fwd( (dir == UNSTABLE_PLUS || dir == UNSTABLE_MINUS) ? true : false  ),
    thePeriod(fp.K),
    p_or_m( (dir == STABLE_PLUS || dir == UNSTABLE_PLUS) ? true : false  )
{
    //Determine saddle strength
    if (std::abs(fp.si) <= 10.0) {
        weakStrength = true;
    }
    /*//Debug constructor
    std::cerr << " MapManifold Constructor :\n";
    std::cerr << "   id = " << id << "\n";
    std::cerr << "   fpdOrbitIdx = " << fpOrbitIdx << "\n";
    std::cerr << "   fpdPointIdx = " << fpPointIdx << "\n";
    std::cerr << "   fp.pos = " << fp.pos << "\n";
    std::cerr << " bfp->pos = " << baseFixedPoint.pos << "\n";
    std::cerr << "   period = " << period << "\n";
    std::cerr << "thePeriod = " << thePeriod << "\n";
    std::cerr << "      fwd = " << fwd << "\n";
    std::cerr << "   p_or_m = " << p_or_m << "\n";*/
    
    //Pre-allocate segment arrays
    segments.reserve(2500);
}



///Constructor useful for creating a blank with some data
template<class SRDATA,class FP>
MapManifold<SRDATA,FP>::
MapManifold(Perturbation dir, ManifoldChildType mct) :
    workingSegID(-1),
    workingSegInitialized(false),
    lastSegID(-1),
    totalArcLength(0.0),
    progressLength(0.0),
    seedLength(0.0),
    numSegments(0),
    numSeedSegments(0),
    depthPS(4),
    firstSubdivision(false),
    majorSubdivision(false),
    majorSubDepth(-1),
    noNewSegsAtNewDepth(true),
    useFixedTauMin(false),
    isComplete(false),
    parent(-1),
    manifoldID(-1),
    childType(mct),
    homoclinicPartnerID(-1),
    nextManifold(-1),
    prevManifold(-1),
    baseFixedPoint(),
    fpdOrbitIdx(-1),
    fpdPointIdx(-1),
    mType(dir),
    weakStrength(false),
    firstStepCount(0),
    adaptiveEdge(VecType(0),VecType(0)),
    numSubdivisions(0),
    emptyNodes(0),
    id( boost::uuids::random_generator()() ),
    period(0),
    fwd( (dir == UNSTABLE_PLUS || dir == UNSTABLE_MINUS) ? true : false  ),
    thePeriod(0),
    p_or_m( (dir == STABLE_PLUS || dir == UNSTABLE_PLUS) ? true : false  )
{
    //Still need to set pointers & fpData
    //Pre-allocate segment arrays
    segments.reserve(2500);
    
}

/// Constructor (Default/Blank)
template<class SRDATA,class FP>
MapManifold<SRDATA,FP>::
MapManifold() :
    workingSegID(-1),
    workingSegInitialized(false),
    lastSegID(-1),
    totalArcLength(0.0),
    progressLength(0.0),
    seedLength(0.0),
    numSegments(0),
    numSeedSegments(0),
    depthPS(4),
    firstSubdivision(false),
    majorSubdivision(false),
    majorSubDepth(-1),
    noNewSegsAtNewDepth(true),
    useFixedTauMin(false),
    isComplete(false),
    parent(-1),
    manifoldID(-1),
    childType(ADULT),
    homoclinicPartnerID(-1),
    nextManifold(-1),
    prevManifold(-1),
    baseFixedPoint(),
    fpdOrbitIdx(-1),
    fpdPointIdx(-1),
    mType(UNSTABLE_PLUS),
    weakStrength(false),
    firstStepCount(0),
    adaptiveEdge(VecType(0),VecType(0)),
    numSubdivisions(0),
    emptyNodes(0),
    id( boost::uuids::random_generator()() ),
    period(0),
    fwd(true),
    thePeriod(0),
    p_or_m(true)
{
    //Call the read() function to set values
    //read(f,fpData);
    //Pre-allocate segment arrays
    segments.reserve(2500);
    
}

/// Copy Constructor
template<class SRDATA,class FP>
MapManifold<SRDATA,FP>::
MapManifold(const MapManifold<SRDATA,FP>& other) :
    workingSegID( other.workingSegID ),
    workingSegInitialized( other.workingSegInitialized ),
    lastSegID( other.lastSegID ),
    totalArcLength(other.totalArcLength),
    progressLength(other.progressLength),
    seedLength(other.seedLength),
    numSegments(other.numSegments),
    numSeedSegments(other.numSeedSegments),
    depthPS(other.depthPS),
    firstSubdivision(other.firstSubdivision),
    majorSubdivision(other.majorSubdivision),
    majorSubDepth(other.majorSubDepth),
    noNewSegsAtNewDepth(other.noNewSegsAtNewDepth),
    useFixedTauMin(other.useFixedTauMin),
    isComplete(other.isComplete),
    parent(other.parent),
    manifoldID(other.manifoldID),
    childType(other.childType),
    homoclinicPartnerID(other.homoclinicPartnerID),
    nextManifold(other.nextManifold),
    prevManifold(other.prevManifold),
    baseFixedPoint(other.baseFixedPoint),
    fpdOrbitIdx(other.fpdOrbitIdx),
    fpdPointIdx(other.fpdPointIdx),
    mType(other.mType),
    weakStrength(other.weakStrength),
    firstStepCount(other.firstStepCount),
    segments(other.segments),
    adaptiveEdge(VecType(0),VecType(0)),
    numSubdivisions(other.numSubdivisions),
    emptyNodes(other.emptyNodes),
    id( other.id ),
    period( other.period ),
    fwd( other.fwd ),
    thePeriod( other.thePeriod ),
    p_or_m(other.p_or_m)
{
    //Have to set the adaptiveEdge -> force restart on this working seg
    workingSegInitialized = false;
    if( workingSegID>=0 && workingSegID<(int)segments.size() ) {
        VecType& x0 = (segments[workingSegID])[0];
        VecType& x1 = (segments[workingSegID])[1];
        //Setup an AdaptiveEdge structure to handle the subdivision process
        adaptiveEdge.reset(x0,x1);
    }
    
}


/// Print base-level information about this manifold
template<class SRDATA,class FP>
void MapManifold<SRDATA,FP>::
info() const
{
    std::cerr << " Manifold (" << id << ") :\n";
    std::cerr << "     manifoldID = " << manifoldID << "\n";
    std::cerr << "     (orbitID,fpID) = (" << fpdOrbitIdx << " , " << fpdPointIdx << " )\n";
    std::cerr << "     thePeriod = " << thePeriod << "\n";
    std::cerr << "     Perturbation = " << pIntToString(pertToInt(mType))
              << "  Type = " << cIntToString(childToInt(childType)) << "\n";
    std::cerr << "     WorkingSegID = " << workingSegID << "\n";
    std::cerr << "        LastSegID = " << lastSegID << "\n";
    std::cerr << "        ArcLength = " << totalArcLength << "\n";
    std::cerr << "  Progress Length = " << progressLength << "\n";
    std::cerr << "   First Subdiv.  =" << firstSubdivision << "\n";
    std::cerr << "    Num Seed Segs = " << numSeedSegments << "\n";
    std::cerr << "       emptyNodes = " << emptyNodes << "\n";
    std::cerr << "  Major Sub Depth = " << majorSubDepth << "\n";
}

///Print an entry to an update table
template<class SRDATA,class FP>
void MapManifold<SRDATA,FP>::
printUpdateEntry()
{
    int currentDepth = getDepth( workingSegID );
    int nextSegDepth = (workingSegID+1 >= (int)segments.size())? currentDepth+1 : getDepth( workingSegID + 1);
    //Gather segments at current depth for information
    std::vector<int> ids;
    getAllSegmentIDsAtDepth(ids,currentDepth);
    int numSegsAtDepth = (int) ids.size();
    int segNumOnDepth = workingSegID - ids[0];
    float progressAtDepth = 100.0 * (float) segNumOnDepth / (float) numSegsAtDepth;
    std::ostringstream os;
    os << thePeriod << ":(" << fpdOrbitIdx << "," << fpdPointIdx << ")";
    //Write entry
    std::cout << std::setw(6)  << manifoldID
              << std::setw(12) << os.str()
              << std::setw(8)  << (emptyNodes >= 1)
              << std::setw(8)  << (emptyNodes == 0)
              << std::setw(8)  << workingSegID
              << std::setw(8)  << currentDepth
              << std::setw(8)  << nextSegDepth
              << std::setw(10)  << getNumSubdivisionLayers( workingSegID )
              << std::setw(12) << (int) segments.size()
              << std::setw(15) << numSegsAtDepth
              << std::setw(10)  << progressAtDepth << "\n";
}

/// Manual step function for initial steps
template<class SRDATA,class FP>
void MapManifold<SRDATA,FP>::
autoInitStep(const double& mup, const double& lambda, double& eps, double& epsMax, bool superAlignment)
{
    //Modify the initial coefficients based on mup:
    double a=0.0,b=0.0,c=0.0,d=0.0;
    getAutoStepCoeff(mup,a,b,c,d);
    //std::cerr << " autoInitStep: mup = " << mup << " [" << a << ", " << b << ", " << c << ", " << d << "]\n";
    
    //Assume eps and epsMax are input
    if (lambda <= 150.0) {
        eps = a*(1.0-exp( 2.0*(150.0-lambda)/(1.0-lambda) )) + b;
    } else {
        eps = c*(1.0-exp( (150.0-lambda)/200.0 )) + d;
    }
    //if(superAlignment) eps = 1.e-7;
    epsMax = 20.0*eps;
    //std::cerr << " autoInitStep: eps = " << eps << " , epsMax = " << epsMax << "\n";
}

/// Compute auto step coefficients for autoInitStep()
template<class SRDATA,class FP>
void MapManifold<SRDATA,FP>::
getAutoStepCoeff(const double& mup, double& a, double& b, double& c, double& d) const
{
    double a0 = 5e-5 - 1e-7;
    double b0 = 1e-7;
    double c0 = 1e-9 - 1e-7;
    double d0 = 1e-7;
    
    double mEM   = 0.012150571430596;
    double mSEnc = 1.8984152807945e-07;
    double r = log(mSEnc) - log(mEM);
    double dm = log(mup) - log(mEM);
    
    double aSE = 2e-5 - 1e-8;
    double bSE = 1e-8;
    double cSE = 1e-10 - 1e-8;
    double dSE = 1e-8;
    double da = (aSE - a0)/r;
    double db = (bSE - b0)/r;
    double dc = (cSE - c0)/r;
    double dd = (dSE - d0)/r;
    
    a = da*dm + a0;
    b = db*dm + b0;
    c = dc*dm + c0;
    d = dd*dm + d0;
    //Debug:
    //std::cerr << " AutoStepCoeff: mup = " << mup << " [" << a << ", " << b << ", " << c << ", " << d << "]\n";
}
/** Create the first step off from a fixed point (note: parent-only function)
 *  - Compute the first step (assuming eigenvectors are in correct direction)
 *    utilizing a stepping function (if manualStep is disabled)
 *    The step distance will scale as a function of Lyapunov exponent magnitude.
 *  - This is an independent step per manifold, and thus, can be
 *    submitted to task-scheduling.
 */
template <class SRDATA, class FP>
template <class MAP, class PARAM>
typename SRDATA::VecType MapManifold<SRDATA,FP>::
firstSegmentStep(const MAP& theMap, const PARAM& params, const ManifoldSettings& settings)
{
    typedef typename MAP::rhs_type             RHStype;
    static const int s = RHStype::numSingularities;
    typedef nvis::fixed_vector<double,s+1>     ExtendedMapDataVec;
    typedef SRDATA                             SortableData;
    
    
    //Utilize under a task call as long as the First MapManifoldSegment is initialized
    const metric_type& theMetric = params.the_metric;
    
    
    //Find the initial step => Use function relating max eigenvalue and eps
    double eps = settings.eps; //The step off magnitude (non-dim) from fixed point
    double epsMax = settings.sdelta_max;
    const double& lambda = fabs(baseFixedPoint.eval[1]);
    //    ->     Unstable eigval(i.e., max, indicates saddle strength)
    
    //Check for super-alignment (eigenvectors are nearly pointing the same direction on the map!)
    double caEvec = nvis::inner( baseFixedPoint.evec[1] , baseFixedPoint.evec[0] );
    bool superAlignment = false;
    if (weakStrength && caEvec >= cos(M_PI/180.0)) {
        superAlignment = true; //Linear analysis is becoming untrustworthy!
        //Need really really tight step size even if we have weak manifolds!
    }
    //Piecewise function - Hard coded for CR3BP
    if (!settings.manualStep) {
        autoInitStep(theMap.rhs().getMu(),lambda,eps,epsMax,superAlignment);
    }
    VecType evec = (fwd ? baseFixedPoint.evec[1] : baseFixedPoint.evec[0]);
    
    
    
    if (params.verbose || params.debug) {
        #pragma omp critical
        {
            std::cerr << " MapManifold->firstSegment():  Computing phi1 (intial guess phi1=phi0+eps*v0)\n";
            std::cerr << " MapManifold->firstSegment():  lambdaMax = " << lambda << " \n";
            std::cerr << "        mu = " << theMap.rhs().getMu() << "\n";
            std::cerr << "        eps = " << eps << "\n";
            std::cerr << "        epsMax = " << epsMax << "\n";
            std::cerr << "        Perturbation: type = " << (fwd ? "unstable" : "stable")
            << ", period = " << thePeriod << "\n";
            std::cerr << "        fp.EigenValue = " << (fwd? baseFixedPoint.eval[1] : baseFixedPoint.eval[0] ) << "\n";
            std::cerr << "        fp.EigenVector = " << (fwd? baseFixedPoint.evec[1] : baseFixedPoint.evec[0]) << "\n";
            std::cerr << "        EigenVector = " << evec << "\n";
            std::cerr << "                    :  evec = " << evec << ",  step = " << epsMax
            << " with (" << (p_or_m ? "+":"-") << ", eps =" << eps << ")\n";
            std::cerr << "        Super-algined = " << superAlignment << "\n";
        }
    }
    
    //Perturb the state based on stored eigenvector - Start with largest possible step
    VecType p = baseFixedPoint.pos + epsMax * (p_or_m ? 1 : -1) * evec;
    
    //Increment counter to indicate first step was called
    firstStepCount++;
    
    return p;
}

///Check and update start segment (note: parent-only function)
template <class SRDATA, class FP>
template <class MAP, class PARAM>
bool MapManifold<SRDATA,FP>::
firstSegmentCheck(
    VecType& p,               //The current step (is overwritten with new value)
    const VecType& map_of_p,  //The Poincare map of the current step (input)
    const MAP& theMap,                        //Poincare map engine
    const PARAM& params,                      //Poincare map parameters
    const ManifoldSettings& settings          //Manifold settings
)
{
    typedef typename MAP::rhs_type                           RHStype;
    
    const metric_type& theMetric = params.the_metric;
    
    
    //Find the initial step => Use function relating max eigenvalue to eps,epsMax
    double eps = settings.eps; //The step off magnitude (non-dim) from fixed point
    double epsMax = settings.sdelta_max; //Largest possible step
    const double& lambda = fabs(baseFixedPoint.eval[1]);
    //    ->     Unstable eigval(i.e., max, indicates saddle strength)
    //Check for super-alignment (eigenvectors are nearly pointing the same direction on the map!)
    double caEvec = nvis::inner( baseFixedPoint.evec[1] , baseFixedPoint.evec[0] );
    bool superAlignment = false;
    if (weakStrength && caEvec >= cos(M_PI/180.0)) {
        superAlignment = true; //Linear analysis is becoming untrustworthy!
        //Need really really tight step size even if we have weak manifolds!
    }
    //Piecewise function - Hard coded for CR3BP
    if (!settings.manualStep) {
        autoInitStep(theMap.rhs().getMu(),lambda,eps,epsMax,superAlignment);
    }
    VecType evec = (fwd ? baseFixedPoint.evec[1] : baseFixedPoint.evec[0]);
    
    
    bool isAligned = false;
    //The Check to see if we have alignment :
    //     (eigenvector and first step vector are aligned and not too far away)
    if (params.debug)
        std::cerr << "  Step " << firstStepCount
                  << "  BaseFPID = " << fpdOrbitIdx
                  << " (step = " << epsMax - (firstStepCount+1)*eps << "):  p = " << p;
    if ((firstStepCount+1) >= 20) {
        return false;    //Reached max number of steps
    }
    
    if (params.debug) {
        std::cerr << " map(p) = " << map_of_p << "\n";
    }
    VecType disp = theMetric.displacement(p, map_of_p);
    // std::cerr << "     at [u1-u0] = " << disp << " norm = " << nvis::norm(disp)
    //     << " ( where delta_min = " << settings.delta_min << ") \n";
    // std::cerr << "     (test : u1-u0 = " << map_of_p-p
    //     << " and theMetric.displacement(u0,u1) = " << disp << ")\n";
    disp /= nvis::norm(disp);
    VecType p0u0 = theMetric.displacement(baseFixedPoint.pos,p);
    //std::cerr << "     (test : u0-p0 = " << p-saddle.pos
    //     << " and theMetric.displacement(p0,u0) = " << p0u0 << ")\n";
    p0u0 /= nvis::norm(p0u0);
    double cosalpha = nvis::inner( p0u0 , disp );
    //std::cerr << "     at cosAngle = "
    // << cosalpha << " so alpha = " << acos(cosalpha)
    // << " (alpha_max = " << settings.alpha_max << ") \n";
    if (!isAligned && cosalpha > cos(settings.alpha_max)) {
        isAligned = true;
    } else if(!isAligned) {
        //Update the next iterate with a
        p -= eps * (p_or_m ? 1 : -1) * evec; //Overwrites input
        //May want to update by halving this distance every time
        firstStepCount++;
    }
    
    return isAligned;
    
}

///Check and update start segment by halfing the minimum step (up to 5 times) (note: parent-only function)
template <class SRDATA, class FP>
template <class MAP, class PARAM>
bool MapManifold<SRDATA,FP>::
firstSegmentCheckDivisor(
    VecType& p,               //The current step (is overwritten with new value)
    const VecType& map_of_p,  //The Poincare map of the current step (input)
    const MAP& theMap,                        //Poincare map engine
    const PARAM& params,                      //Poincare map parameters
    const ManifoldSettings& settings          //Manifold settings
)
{
    typedef typename MAP::rhs_type                           RHStype;
    
    const metric_type& theMetric = params.the_metric;
    
    
    //Find the initial step => Use function relating max eigenvalue to eps,epsMax
    double eps = settings.eps; //The step off magnitude (non-dim) from fixed point
    double epsMax = settings.sdelta_max; //Largest possible step
    const double& lambda = fabs(baseFixedPoint.eval[1]);
    //    ->     Unstable eigval(i.e., max, indicates saddle strength)
    //Check for super-alignment (eigenvectors are nearly pointing the same direction on the map!)
    double caEvec = nvis::inner( baseFixedPoint.evec[1] , baseFixedPoint.evec[0] );
    bool superAlignment = false;
    if (weakStrength && caEvec >= cos(M_PI/180.0)) {
        superAlignment = true; //Linear analysis is becoming untrustworthy!
        //Need really really tight step size even if we have weak manifolds!
    }
    //Piecewise function - Hard coded for CR3BP
    if (!settings.manualStep) {
        autoInitStep(theMap.rhs().getMu(),lambda,eps,epsMax,superAlignment);
    }
    VecType evec = (fwd ? baseFixedPoint.evec[1] : baseFixedPoint.evec[0]);
    
    
    bool isAligned = false;
    //The actual step-off is now eps/2^f
    int factor = firstStepCount - 20 + 1;
    double divisor = 1.0;
    for(int f=0; f<factor; f++) {
        divisor*=2.0;
    }
    double step = eps / divisor;
    
    //The Check to see if we have alignment :
    //     (eigenvector and first step vector are aligned and not too far away)
    if (params.debug)
        std::cerr << "  Divisor Step " << firstStepCount - 20
                  << "  BaseFPID = " << fpdOrbitIdx
                  << " (step = " << step  << "):  f = " << factor << " p = " << p;
    if ((firstStepCount-20+1) >= 7) {
        return false;    //Reached max number of steps
    }
    
    
    if (params.debug) {
        std::cerr << " map(p) = " << map_of_p << "\n";
    }
    VecType disp = theMetric.displacement(p, map_of_p);
    /*std::cerr << "     at [u1-u0] = " << disp << " norm = " << nvis::norm(disp)
        << " ( where delta_min = " << settings.delta_min << ") \n";
    std::cerr << "     (test : u1-u0 = " << map_of_p-p
        << " and theMetric.displacement(u0,u1) = " << disp << ")\n";*/
    disp /= nvis::norm(disp);
    VecType p0u0 = theMetric.displacement(baseFixedPoint.pos,p);
    /*std::cerr << "     (test : u0-p0 = " << p-baseFixedPoint.pos
         << " and theMetric.displacement(p0,u0) = " << p0u0 << ")\n";*/
    p0u0 /= nvis::norm(p0u0);
    double cosalpha = nvis::inner( p0u0 , disp );
    /*std::cerr << "     at cosAngle = "
     << cosalpha << " so alpha = " << acos(cosalpha)
     << " (alpha_max = " << settings.alpha_max << ") \n";*/
    if (!isAligned && cosalpha >= cos(settings.alpha_max)) {
        isAligned = true;
    } else if(!isAligned) {
        //Reset the step (p)
        p = baseFixedPoint.pos + (step/2.0) * (p_or_m ? 1.0 : -1.0) * evec;
        //Indicate that we need another iteration
        firstStepCount++;
    }
    
    return isAligned;
    
}


/// Force the first segment step to be minimum linear step
template <class SRDATA, class FP>
template <class MAP, class PARAM>
void MapManifold<SRDATA,FP>::
firstSegmentForcedStep(
    VecType& p,
    const MAP& theMap,
    const PARAM& params,
    const ManifoldSettings& settings)
{
    typedef typename MAP::rhs_type                           RHStype;
    
    //Find the initial step => Use function relating max eigenvalue to eps,epsMax
    double eps = settings.eps; //The step off magnitude (non-dim) from fixed point
    double epsMax = settings.sdelta_max; //Largest possible step
    const double& lambda = fabs(baseFixedPoint.eval[1]);
    //    ->     Unstable eigval(i.e., max, indicates saddle strength)
    
    //Piecewise function - Hard coded for CR3BP
    bool superAlignment = false;
    if (!settings.manualStep) {
        autoInitStep(theMap.rhs().getMu(),lambda,eps,epsMax,superAlignment);
    }
    VecType evec = (fwd ? baseFixedPoint.evec[1] : baseFixedPoint.evec[0]);
    
    //Force the step to be a fixed value
    p = baseFixedPoint.pos + eps * (p_or_m ? 1 : -1) * evec;
}



///Check and update start segment (note: parent-only function)
template <class SRDATA, class FP>
template <class MAP, class PARAM>
void MapManifold<SRDATA,FP>::
storeFirstSegment(
    const VecType& p,                         //The current step point on the surface of section
    std::set<SRDATA>& cache,                  //Data cache for storing segment
    const MAP& theMap,                        //Poincare map engine
    const PARAM& params,                      //Poincare map parameters
    const ManifoldSettings& settings          //Manifold settings
)
{
    typedef typename MAP::rhs_type        RHStype;
    typedef typename MAP::lbox_type       BoundsType;
    static const int s = RHStype::numSingularities;
    //Set analysis bounds for not counting large growth near singularities
    BoundsType bbox = settings.bounds;
    
    //Need to lookup some data based on the call (won't run theMap)
    MapDiscont eFail;
    SortableData pData;
    MapDataType data1;
    IdType nodeID(1,0); //inserted first segment second node
    bool failDetected = false;
    try {
        //Lookup the map data for the first step x0 = x +/- eps*v
        pData = mapDataUsingCache(p,theMap,thePeriod,params,cache,eFail);
        //const VecType& map_of_p = pData.getReturn(std::abs(thePeriod));
        //double dt1 = pData.data[std::abs(thePeriod)-1][s]; //last element is TOF
        data1 = pData.getDataAtReturn(std::abs(thePeriod));
        // -> this could potentially throw an error in sensitive problems! (Map signularity)
    } catch(MapDiscont& e) {
        //Detected failure (likely Map Singularity) in cache
        failDetected = true;
        eFail = e;
    } catch(...) {
        failDetected = true;
        MapDiscont nodeTransFail(
            (fwd)? MapDiscont::FORWARD_SECTION_SEP_NODE : MapDiscont::BACKWARD_SECTION_SEP_NODE,
            p, thePeriod);
        nodeTransFail.edgeID = nodeID;
        eFail = nodeTransFail;
    }
    data1[s] = 0.0; //Reset TOF to correlate with tree traversal
    //Truthfully, you need a corrections scheme to figure out what the real TOF should be to this node.
    
    //Store to Manifold Segment [store previous 'data' in segment]
    ManifoldSeg theFirstSeg(baseFixedPoint.pos,p,data1);
    segments.push_back( theFirstSeg );
    //Set some indicators to start advection process
    workingSegID = 0;
    lastSegID = workingSegID;
    //Initialize the depth counter
    std::vector<int> startingDepthVec;
    startingDepthVec.push_back(workingSegID);
    idsAtDepth.insert( std::pair<int,std::vector<int> >(1,startingDepthVec) );
    
    //Initialize arc-length to current working segment edge length
    const double edgeLength = segments[workingSegID].length(params);
    totalArcLength = edgeLength;
    VecType& x0 = (segments[workingSegID])[0];
    VecType& x1 = (segments[workingSegID])[1];
    if(params.debug) {
        std::cout << "First Segment: \n";
        std::cout << "   x0 = " << x0 << "\n";
        std::cout << "   x1 = " << x1 << "\n";
        std::cout << "  edgeLength = " << edgeLength << "\n";
    }
    if(bbox.inside(x0) && bbox.inside(x1)) {
        progressLength += edgeLength;
    }
    seedLength = edgeLength;
    numSegments = 1;
    numSeedSegments = 1;
    workingSegInitialized = false;
    
    //Setup an AdaptiveEdge structure to handle the subdivision process
    adaptiveEdge.reset(x0,x1);
    adaptiveEdge.refine(IdType(0,0)); //Default one refine command
    VecType xm = (x0+x1) / 2.0;
    adaptiveEdge.setVertex(IdType(1,1),xm);
    
    //Setup failure on first segment if detected
    if(failDetected) {
        adaptiveEdge.insertDiscontinuityNode(nodeID,eFail);
    }
}

///Propagation call for point on this manifold utilizing data storage (taskStorage & dataCache)
template <class SRDATA, class FP>
template <class MAP, class PARAM>
void MapManifold<SRDATA,FP>::
propagateManifoldArcTask(const VecType& x0, const MAP& theMap, const PARAM& params,
                         std::set<SRDATA>& cache, std::set<SRDATA>& backCache, MMStorageVector* perThreadCache)
{
    typedef typename MAP::rhs_type        RHStype;
    typedef SRDATA                        SortableData;
    
    //Note:  Utilize under a task call
    
    int tid = 0;
#if _OPENMP
    tid = omp_get_thread_num();
#endif
    
    /// The per-thread call to the map
    VecType p(0);
    if(fwd) {
        p = mapUsingTaskCache(x0, theMap, thePeriod, params, cache, id, perThreadCache[tid]);
    } else {
        p = mapUsingTaskCache(x0, theMap, thePeriod, params, backCache, id, perThreadCache[tid]);
    }
    //Don't do anything with p
}


///Propagation call for point on the workingSeg of this manifold utilizing data storage (taskStorage & dataCache)
template <class SRDATA, class FP>
template <class MAP, class PARAM>
void MapManifold<SRDATA,FP>::
propagateManifoldArcTask(const typename MapManifold<SRDATA,FP>::JobIDPair& jobID,
                         const MAP& theMap, const PARAM& params,
                         std::set<SRDATA>& cache, std::set<SRDATA>& backCache, MJStorageVector* perThreadCache)
{
    typedef typename MAP::rhs_type        RHStype;
    typedef SRDATA                        SortableData;
    
    //Note:  Utilize under a task call
    
    int tid = 0;
#if _OPENMP
    tid = omp_get_thread_num();
#endif
    
    //Lookup the x value from the Working Segment adaptive edge structure
    IdType nodeID = jobID.second;
    VecType x0 = adaptiveEdge.getVertex(nodeID);
    
    /// The per-thread call to the map
    VecType p(0);
    if(fwd) {
        p = mapUsingTaskCache(x0, theMap, thePeriod, params, cache, jobID, perThreadCache[tid]);
    } else {
        p = mapUsingTaskCache(x0, theMap, thePeriod, params, backCache, jobID, perThreadCache[tid]);
    }
    //Don't do anything with p
}


/*
/// Propagate to the estimated Progenitor State given a segment index and point index
template<class SRDATA, class FP>
template<class PMAP,class PARAM>
void MapManifold<SRDATA,FP>::
propagateToProgenitorGuess(
        PMAP& theMap,
        PARAM& theMapParams,
        const int jobID,
        const int segID,
        const int ptID,
        typename MapManifold<SRDATA,FP>::TrajJobStorageVector* tjCache)
{
    typedef typename PMAP::rhs_type     RHStype;
    typedef typename PMAP::return_type  return_type;
    typedef typename PMAP::return_state return_state;


    //NOTE:  This is employed within a openMP task call

    int tid = 0;
    #if _OPENMP
      tid = omp_get_thread_num();
    #endif

    //Get the Manifold Point
    VecType x0 = segments[segID][ptID];
    //Compute the approximate number of return periods
    int period = getPeriod();
    int currentDepth = getDepth(segID);
    int p = (-std::abs(thePeriod))*(currentDepth+1); //Negative indicates upstream, go to depth -1 to pass through orbit
    std::vector<return_type> rTypeIters;
    std::vector<return_state> out_rstates;
    //Try to simulate
    bool trajOK = true;
    try {
        theMap.flow_map(x0,rTypeIters,out_rstates,p);
    } catch(...) {
        //Don't need to do anything but report the error
        trajOK = false;
    }
    //Store information to a Trajectory object
    TrajectoryType progArc;
    progArc.setValidity( trajOK );
    if(trajOK) {
        //Fill out the trajectory information (note, always stored in forward time for consistency)
        std::vector<double> times;
        std::vector<StateType> states;
        for(int i=0;i<(int)out_rstates.size();i++) {
            states.push_back( out_rstates[i].x );
            times.push_back( out_rstates[i].t );
        }
        bool bwd =((times.back() - times[0])<0.0)? true : false;
        progArc.setStates(states,times,bwd);
    }

    //Store the result to a per-thread cache
    tjCache[tid].push_back( TrajJobPair(jobID,progArc) );
}*/

/** Submitting propagation tasks from the current process
 *  - Submits a (uuid,edgeID) object for processing
 */
template <class SRDATA, class FP>
void MapManifold<SRDATA,FP>::
getEdgeNodeMapTasks(const int& manID, std::vector<typename MapManifold<SRDATA,FP>::JobIDPair>& jobs )
{
    //The first pass needs the end points and midpoint
    if (!workingSegInitialized) {
        jobs.push_back( JobIDPair(manID, IdType(0,1)) ); //left point
        jobs.push_back( JobIDPair(manID, IdType(1,1)) ); //mid point
        jobs.push_back( JobIDPair(manID, IdType(1,0)) ); //right point
        emptyNodes = 3;
        workingSegInitialized = true; //Will be initialized after first set of map calls
    } else {
        //Otherwise, directly call from AdaptiveEdge structure
        std::vector<IdType> seedIDs;
        adaptiveEdge.getEmptyLeaves(seedIDs);
        emptyNodes = (int)seedIDs.size();
        for(int j=0; j<emptyNodes; j++) {
            jobs.push_back( JobIDPair(manID,seedIDs[j]) );
        }
    }
    //Debugging
    //std::cout << " Manifold : " << id << "\n";
    //std::cout << "   Submits " << emptyNodes << " jobs.\n";
}

/** Check the current AdaptiveEdge data to see if we need to subdivision.
 *  Note:  This also performs the refinement process [refine()] and
 *  inserts MapDiscontinuities, if necessary.
 *  Thus, this function should NOT be called in a parallel implementation.
 */
template <class SRDATA, class FP>
template <class MAP, class PARAM>
void MapManifold<SRDATA,FP>::
refineSegmentsTest(
    std::set<SRDATA>& cache,
    const MAP& theMap,
    const PARAM& params,
    const ManifoldSettings& settings
)
{
    typedef typename MAP::lvec_type            vec_type;
    typedef typename MAP::lmat_type            mat_type;
    typedef typename MAP::xstate_type          xstate_type;
    typedef typename MAP::lbox_type            BoundsType;
    typedef std::pair<vec_type,mat_type>       ProjPair;
    //typedef EdgeRotationFailure<vec_type>      MapDiscont;
    //Adaptive Edge data structures
    //typedef AdaptiveEdge<vec_type,vec_type>    EdgeType;
    typedef typename EdgeType::IdType          IdType;
    typedef std::pair<IdType,IdType>           AESegment;
    typedef typename EdgeType::DataType        DataPairs;
    //Number of singularities
    static const int numSing = MAP::rhs_type::numSingularities;
    
    //Parameters
    bool verbose = params.verbose;
    bool debug = params.debug;
    //const double& dTauMin = settings.delta_tau_min; //Minimum allowed distance between points on segment
    //Note:  Actual delta_tau_min is dependent on the current segment and a u_min on allowed displacement on section.
    const double dTauMin = getAdaptiveDeltaTauMin(params,settings.delta_tau_min);
    
    //Info about the current working segment:
    const double edgeLength = segments[workingSegID].length(params);
    
    //For a special print debug case:
    //debug = (params.debug && manifoldID == 0 && workingSegID == 77)? true : false;
    //debug = (params.debug && manifoldID == 8 && workingSegID == 0)? true : false;
    
    
    //The current segments
    std::vector<AESegment> segs;
    adaptiveEdge.getAllSegments(segs);
    typename std::vector<AESegment>::iterator segIT, segITm1, segITp1;
    segIT = segITm1 = segs.begin();
    //Check each segment for subdivision    ----------------------------
    for(; segIT!=segs.end(); ++segIT) {
        segITm1 = segIT;
        --segITm1;
        segITp1 = segIT;
        ++segITp1;
        //Work on previous and current segment
        bool firstSeg = (segIT==segs.begin());
        bool lastSeg = nvis::all(segIT->second == IdType(1,0));
        
        //if(debug) std::cout << "Testing segment : ManID = " << manifoldID << " WorkingSegID = " << workingSegID << "\n";
        //if(debug) std::cout << "   Points: " << segIT->first << " , " << segIT->second << "\n";
        
        //Checking if we have to subdivide based on three nodes
        vec_type xA = adaptiveEdge.getVertex(segIT->first);
        vec_type xB = adaptiveEdge.getVertex(segIT->second);
        vec_type xC = adaptiveEdge.getVertex(
                          (firstSeg)? segITp1->second : segITm1->first);
        //Get parameters
        double tauA = adaptiveEdge.getParameter(segIT->first);
        double tauB = adaptiveEdge.getParameter(segIT->second);
        //double tauC = adaptiveEdge.getParameter(
        //               (firstSeg)? segITp1->second : segITm1->first);
        bool subdivide = false;
        //If testing segment is too far apart, subdivide
        double dTau = tauB - tauA;
        //if(dTau > dTauMax) subdivide = true; //Don't add any extra subdivision
        //Check for a discontinuity point, subdivide if still above dTauMin
        if(!subdivide && (dTau > dTauMin)) {
            //Are Node 0 and Node 1 both discontinuities? stop
            if(adaptiveEdge.isDiscontinuity(segIT->first) &&
                    adaptiveEdge.isDiscontinuity(segIT->second)) {
                //Stop subdivision as we likely reached a special case
                //if(debug) std::cout << "   Seg has dTau > dTauMin && discontinuity at both nodes\n";
                continue;
            }
            //Is Node 0 or Node 1 a discontinuity point? subdivide?
            else if(adaptiveEdge.isDiscontinuity(segIT->first) ||
                    adaptiveEdge.isDiscontinuity(segIT->second)) {
                //Subdivide the segments around the node
                subdivide = true;
            }
        }
        //Check the downstream manifold for curve-refinement parameters
        if(!subdivide) {
            //Get map data (downstream)
            vec_type pA = adaptiveEdge.getValue(segIT->first);
            vec_type pB = adaptiveEdge.getValue(segIT->second);
            vec_type pC = adaptiveEdge.getValue(
                              (firstSeg)? segITp1->second : segITm1->first);
            /*if(debug) {
              std::cout << " Data Points:\n";
              std::cout << "   Seg: x0    = " << xA << " , x1    = " << xB << "\n";
              std::cout << "   Seg: P(x0) = " << pA << " , P(x1) = " << pB << "\n";
              std::cout << "   Prev: xC = " << xC << " , P(xC) = " << pC << "\n";
            }*/
            
            
            //Run the downstream curve-refinement checks
            bool curveOK = false;
            //If primary segment (A-B) has separation, then curve is not ok
            if( !detectSectionSeparation(xA,xB,cache,theMap,std::abs(thePeriod),params) ) {
                if(firstSeg) {
                    //On first seg, we have to use the next segment's results to test curve
                    curveOK = isManifoldCurveOK(
                                  theMap,std::abs(thePeriod),params,xA,xB,xC,pA,pB,pC,settings,cache,true);
                    //otherwise, curve is not ok
                } else {
                    //Test the previous segment (C-A) for separation
                    bool prevSegHasSep =
                        detectSectionSeparation(xC,xA,cache,theMap,std::abs(thePeriod),params);
                    //Also, test if the previous segment (C-A) is too short to use safely (<delta_min)
                    bool prevSegIsTooShort =
                        isSegmentTooShort(theMap,params,pC,pA,settings);
                    if(prevSegHasSep || prevSegIsTooShort) {
                        if(!lastSeg) {
                            //Try to use the next segment to test segment
                            vec_type xD = adaptiveEdge.getVertex(segITp1->second);
                            vec_type pD = adaptiveEdge.getValue(segITp1->second);
                            curveOK = isManifoldCurveOK(
                                          theMap,std::abs(thePeriod),params,xA,xB,xD,pA,pB,pD,settings,cache,true);
                        }
                    } else {
                        //Standard segment test
                        curveOK = isManifoldCurveOK(
                                      theMap,std::abs(thePeriod),params,xC,xA,xB,pC,pA,pB,settings,cache,false);
                    }
                }
            }
            //std::cout << "   Curve-Refinement = " << curveOK << "\n";
            
            //Decide to subdivide segment based on manifold curve params
            if (curveOK) {
                continue;   //Done with segment
            } else if ( dTau <= dTauMin ) {
                //std::cout << "  Segment at dTauMin: " << dTau << " < dTauMin=" << dTauMin << "\n";
                MapDiscont rotErr; //Empty discontinuity for checking errors
                //First, check if nodes are discontinuities (any type),
                //  if yes, no further errors to report
                if ( adaptiveEdge.isDiscontinuity(segIT->first) ||
                        adaptiveEdge.isDiscontinuity(segIT->second) ) {
                    continue; //Done
                }
                //Check for fixed point conditions (of same period)
                else if(
                    isFixedPointSuspected(xA,pA-xA,xB,pB-xB,std::abs(thePeriod),params,rotErr)
                ) {
                    if (verbose) {
                        std::cerr << "Period (" << std::abs(thePeriod) << "): Fixed point is suspected on the manifold segment!\n";
                    }
                    rotErr.edgeID = adaptiveEdge.computeChildID(segIT->first,1);
                    adaptiveEdge.insertDiscontinuity(rotErr);
                    continue; //Done
                }
                /// If at lmin without successful rotation conditions nor a fixed point suspected,
                /// check if it's a section-separation condition (test both forward and backward):
                bool discontinuityFound = false;
                //Separation
                if (detectSectionSeparation(xA,xB,cache,theMap,std::abs(thePeriod),params) ) {
                    // If we are at lmin, store the forward separation point
                    vec_type intersection = xA+xB;
                    intersection /= 2.0;
                    MapDiscont sepErr((fwd)? MapDiscont::FORWARD_SECTION_SEP : MapDiscont::BACKWARD_SECTION_SEP,
                                      intersection,std::abs(thePeriod));
                    sepErr.edgeID = adaptiveEdge.computeChildID(segIT->first,1);
                    if (detectSingularityIntersection(xA,xB,cache,theMap,std::abs(thePeriod),params) ) {
                        sepErr.setType((fwd)? MapDiscont::FORWARD_SINGULARITY : MapDiscont::BACKWARD_SINGULARITY);
                        //Data output for debugging
                        if(debug) {
                            typename std::set<SRDATA>::iterator cit0,cit1;
                            cit0 = cache.find( SRDATA(xA) ); //Note:  Should be found
                            cit1 = cache.find( SRDATA(xB) ); //Note:  Should be found
                            vec_type p0 = (*cit0).returns[std::abs(thePeriod)-1]; //P^p(x0)
                            vec_type p1 = (*cit1).returns[std::abs(thePeriod)-1]; //P^p(x1)
                            double dt0 = (*cit0).data[std::abs(thePeriod)-1][numSing];
                            double dt1 = (*cit1).data[std::abs(thePeriod)-1][numSing];
                            const xstate_type* singPtr = theMap.rhs().singularities();
                            std::cout << " Singularity detected :\n";
                            std::cout << "  period = " << thePeriod << "\n";
                            std::cout << "  x0 = " << xA << " x1 = " << xB << "\n";
                            std::cout << "  p0 = " << p0 << " p1 = " << p1 << "\n";
                            std::cout << "  dt0 = " << dt0 << " dt1 = " << dt1 << "\n";
                            for(int k=0; k<numSing; k++) {
                                std::cout << "  Singularity " << k << ":\n";
                                double singularitySafeDistance = theMap.rhs().getSingularitySafeDistance(k);
                                double minCloseApproachDist0 = (*cit0).data[std::abs(thePeriod)-1][k];
                                double minCloseApproachDist1 = (*cit1).data[std::abs(thePeriod)-1][k];
                                //Check for sign change in position coordinate from singularity
                                ProjPair singPair = theMap.section().project( singPtr[k] );
                                vec_type s = singPair.first;
                                vec_type r0 = p0 - s;
                                vec_type r1 = p1 - s;
                                std::cout << "    r0 = " << r0 << "  r1 = " << r1 << "\n";
                                std::cout << "    (r0[0]*r1[0])<0) = " << ((r0[0]*r1[0])<0) << "\n";
                                std::cout << "    (p0[1]*p1[1])<0) = " << ((p0[1]*p1[1])<0) << "\n";
                                std::cout << "    SafeDist = " << singularitySafeDistance
                                          << " MinCA0 = " << minCloseApproachDist0
                                          << " MinCA1 = " << minCloseApproachDist1 << "\n";
                            }
                        }
                    }
                    adaptiveEdge.insertDiscontinuity( sepErr );
                    discontinuityFound = true;
                    if(debug) {
                        std::cout << " Discontinuity Found: " << sepErr.what() << "\n";
                    }
                }
                //At this point, we don't know what went wrong, so should throw 'unknown'
                // HOWEVER, it is likely that these are still separation conditions (likely singularity
                // intersections OR data assignment errors).  Pass out as unknown separation
                if (!discontinuityFound) {
                    vec_type intersection = xA+xB;
                    intersection /= 2.0;
                    MapDiscont singErr(MapDiscont::UNRESOLVED,intersection,thePeriod);
                    singErr.edgeID = adaptiveEdge.computeChildID(segIT->first,1);
                    adaptiveEdge.insertDiscontinuity( singErr );
                    if(debug) {
                        std::cout << " Unknown Curve Discontinuity Found: " << singErr.what() << "\n";
                        std::cout << " Data Points:\n";
                        std::cout << "   Seg: x0    = " << xA << " , x1    = " << xB << "\n";
                        std::cout << "   Seg: P(x0) = " << pA << " , P(x1) = " << pB << "\n";
                        std::cout << "   Prev: xC = " << xC << " , P(xC) = " << pC << "\n";
                        
                        /*typename std::set<SRDATA>::iterator cit0,cit1;
                        cit0 = cache.find( SRDATA(xA) ); //Note:  Should be found
                        cit1 = cache.find( SRDATA(xB) ); //Note:  Should be found
                        vec_type p0 = (*cit0).returns[std::abs(thePeriod)-1]; //P^p(x0)
                        vec_type p1 = (*cit1).returns[std::abs(thePeriod)-1]; //P^p(x1)
                        double dt0 = (*cit0).data[std::abs(thePeriod)-1][numSing];
                        double dt1 = (*cit1).data[std::abs(thePeriod)-1][numSing];
                        const xstate_type* singPtr = theMap.rhs().singularities();
                        std::cout << " Singularity detected :\n";
                        std::cout << "  period = " << thePeriod << "\n";
                        std::cout << "  x0 = " << xA << " x1 = " << xB << "\n";
                        std::cout << "  p0 = " << p0 << " p1 = " << p1 << "\n";
                        std::cout << "  dt0 = " << dt0 << " dt1 = " << dt1 << "\n";
                        for(int k=0;k<numSing;k++) {
                              std::cout << "  Singularity " << k << ":\n";
                              double singularitySafeDistance = theMap.rhs().getSingularitySafeDistance(k);
                              double minCloseApproachDist0 = (*cit0).data[std::abs(thePeriod)-1][k];
                              double minCloseApproachDist1 = (*cit1).data[std::abs(thePeriod)-1][k];
                              //Check for sign change in position coordinate from singularity
                              ProjPair singPair = theMap.section().project( singPtr[k] );
                              vec_type s = singPair.first;
                              vec_type r0 = p0 - s;
                              vec_type r1 = p1 - s;
                              std::cout << "    r0 = " << r0 << "  r1 = " << r1 << "\n";
                              std::cout << "    (r0[0]*r1[0])<0) = " << ((r0[0]*r1[0])<0) << "\n";
                              std::cout << "    (p0[1]*p1[1])<0) = " << ((p0[1]*p1[1])<0) << "\n";
                              std::cout << "    SafeDist = " << singularitySafeDistance
                                << " MinCA0 = " << minCloseApproachDist0
                                << " MinCA1 = " << minCloseApproachDist1 << "\n";
                        }*/
                    } //End 'debug' write
                }
                continue; //Force a completion of this segment since at dTauMin
                
            } else {
                subdivide = true;    //Curve needs refinement
            }
        }
        
        //Run subdivision step
        if (subdivide) {
            numSubdivisions++;
            if (!firstSubdivision) {
                firstSubdivision = true;
            }
            //Call the refine command on the left/bottom node
            adaptiveEdge.refine(segIT->first);
            //Set the new vertex (could be midpoint or something else like Secant method)
            vec_type xMidPoint = xA+xB;
            xMidPoint /= 2.0;
            IdType newID = adaptiveEdge.getChildID(segIT->first,1);
            adaptiveEdge.setVertex(newID, xMidPoint);
            //std::cout << " Subdivision values: \n";
            //std::cout << "   New Point : " << newID << " -> " << xMidPoint << "\n";
        }
    } //End for each segment loop
    
    //Evaluate if there is empty information at nodes (signifying that more mappings are need)
    emptyNodes = adaptiveEdge.getNumEmptyLeaves();
}


/// Assign propagation data from a map call to an edge node
template <class SRDATA, class FP>
template <class MAP, class PARAM>
void MapManifold<SRDATA,FP>::
assignEdgeNodeData(
    const typename MapManifold<SRDATA,FP>::IdType& nodeID,  //Point ID on AdaptiveEdge
    const SRDATA& result,                                     //Map Result
    const typename MapManifold<SRDATA,FP>::MapDiscont& mapDis,         //MapDiscontinuity (if occurred)
    bool fail,                                                //Flag indicating mapping failure
    const MAP& theMap,                                        //Poincare Map Engine
    const PARAM& params,                                      //Map Analysis Parameters
    std::set<SRDATA>& cache,                                  //Forward map data cache
    std::set<SRDATA>& backCache                               //Backward map data cache
)
{
    typedef typename MAP::lvec_type            vec_type;
    typedef typename MAP::lmat_type            mat_type;
    typedef typename MAP::xstate_type          xstate_type;
    typedef typename MAP::lbox_type            BoundsType;
    typedef std::pair<vec_type,mat_type>       ProjPair;
    //Adaptive Edge data structures
    //typedef AdaptiveEdge<vec_type,vec_type>    EdgeType;
    typedef std::pair<IdType,IdType>           AESegment;
    typedef typename EdgeType::DataType        DataPairs;
    
    VecType x = adaptiveEdge.getVertex(nodeID);
    //Insert data to cache (does nothing if already there)
    if (fwd) {
        cache.insert( result );
    } else {
        backCache.insert( result );
    }
    
    //Storage objects
    MapDiscont ef = mapDis;
    VecType px(0.0);
    bool detectedSeparation = false;
    if(fail) {
        //Singularity caught
        ef.edgeID = nodeID;
        adaptiveEdge.insertDiscontinuityNode(nodeID,ef);
        px = VecType(50);
        detectedSeparation = true;
        if(params.debug) {
            std::cout << "   Manifold Mapping Encountered Error = " << ef.what() << "\n";
        }
    }
    //Otherwise, it's ok
    else {
        px = result.getReturn(std::abs(thePeriod));
    }
    
    //Run test for transversality condition (dydot(P^p(x)) != 0)
    if( !detectedSeparation ) {
        if( (fwd)? detectSectionSeparation(x,cache,theMap,std::abs(thePeriod),params) :
                detectSectionSeparation(x,backCache,theMap,std::abs(thePeriod),params) ) {
            //Create node failure
            MapDiscont nodeTransFail(
                (fwd)? MapDiscont::FORWARD_SECTION_SEP_NODE : MapDiscont::BACKWARD_SECTION_SEP_NODE,
                x, thePeriod);
            nodeTransFail.edgeID = nodeID;
            adaptiveEdge.insertDiscontinuityNode(nodeID,nodeTransFail);
            detectedSeparation = true;
            if(params.debug) {
                std::cout << "   Section Separation Node = " << nodeTransFail.where() << "\n";
            }
        }
    }
    
    //Manually set the mapping coord of phi0 (i.e., theFixedPoint) to be phi1
    /*if ( workingSegID==0 && nvis::all(nodeID == IdType(0,1)) ) {
         px = adaptiveEdge.getVertex( IdType(1,0) );
    }*/
    
    //Store data into adaptive structure
    adaptiveEdge.setValue( nodeID, px );
}

/** Advecting to the next segment using STHManifold technique:
 *  -> This will store necessary data to the list of manifold segments
 *     and construct tree connectivity if the subdivision process is complete.
 *     Also, updates arc-length computation.
 *  Returns a list of new segmentIDs
 */
template <class SRDATA, class FP>
template <class MAP, class PARAM>
void MapManifold<SRDATA,FP>::
advectSegment(const MAP& theMap,
              const PARAM& params,
              const ManifoldSettings& settings,
              std::set<SRDATA>& cache
             )
{
    typedef typename MAP::lvec_type            vec_type;
    typedef typename MAP::lmat_type            mat_type;
    typedef typename MAP::xstate_type          xstate_type;
    typedef typename MAP::lbox_type            BoundsType;
    typedef std::pair<vec_type,mat_type>       ProjPair;
    //Adaptive Edge data structures
    //typedef AdaptiveEdge<vec_type,vec_type>    EdgeType;
    typedef std::pair<IdType,IdType>           AESegment;
    typedef typename EdgeType::DataType        DataPairs;
    
    static const int numSing = MAP::rhs_type::numSingularities;
    
    //Parameters
    bool verbose = params.verbose;
    bool debug = params.debug;
    
    //const double& dTauMin = settings.delta_tau_min; //Minimum allowed distance between points on segment
    //Note:  Actual delta_tau_min is dependent on the current segment and a u_min on allowed displacement on section.
    const double dTauMin = getAdaptiveDeltaTauMin(params,settings.delta_tau_min);
    
    //Set analysis bounds for not counting large growth near singularities
    BoundsType bbox = settings.bounds;
    
    //See if subdivision occurred on the current step
    //if (emptyNodes > 0) return; //Advection did not occur
    
    //For a special print debug case:
    //debug = (params.debug && manifoldID == 0 && workingSegID == 77)? true : false;
    //debug = (params.debug && manifoldID == 8 && workingSegID == 0)? true : false;
    
    
    
    
    //Otherwise, the current working segment is done, so store new info
    //-------------------------------------------------------------------------------------------
    //MapDiscontinuities
    adaptiveEdge.getDiscontinuityList( segments[workingSegID].discontinuityList );
    segments[workingSegID].discontinuityList.sort( mapDiscontOnEdgeCompare<MapDiscont> );
    int wSegDepth = getDepth(workingSegID);
    //Need to make these unique to avoid 2+ duplicates
    segments[workingSegID].discontinuityList.unique( );
    if((int)(segments[workingSegID].discontinuityList.size()) > 0) {
        segments[workingSegID].transverseSection = false;
    }
    
    //Make new segments based on type: weakStrength vs. normal:
    typename std::set<SRDATA>::iterator cit;
    //Check if this is a weak strength manifold and haven't found first subdivision
    if (weakStrength && !firstSubdivision) { //lambdaMax < 50.0
        //Enhancing slow advection:
        if(debug) {
            std::cerr << " advectSegment(): Utilizing fast-advance for Weak-Strength Manifold.\n";
        }
        //The final point and its mapping
        VecType x = adaptiveEdge.getVertex(IdType(1,0)); //End point of adaptiveEdge
        cit = cache.find( SRDATA(x) );
        if (cit == cache.end()) {
            throw std::runtime_error("Could not find requested map in SRDATA cache");
        }
        const SRDATA& xData = (*cit);
        VecType px = xData.getReturn(std::abs(thePeriod));
        MapDataType pxMapData = xData.getDataAtReturn(std::abs(thePeriod));
        //The initial point and its mapping
        VecType x0 = adaptiveEdge.getVertex(IdType(0,0)); //Start of working segment
        cit = cache.find( SRDATA(x0) );
        if (cit == cache.end()) {
            throw std::runtime_error("Could not find requested map in SRDATA cache");
        }
        const SRDATA& x0Data = (*cit);
        VecType px0 = x0Data.getReturn(std::abs(thePeriod));
        MapDataType px0MapData = xData.getDataAtReturn(std::abs(thePeriod));
        
        //Current last segment end node is the first point
        const ManifoldSeg& workSeg = segments[workingSegID];
        const ManifoldSeg& lastSeg = segments[lastSegID];
        VecType xf = px0; //Standard add just the one segment downstream
        //Advancing initial segments: Add only the segment formed by map of last point and its map
        if (workingSegID == 0) {
            xf = lastSeg[1];
            px0MapData = lastSeg.data1;
        }
        //Compute the displacement from the last segment
        double lastSegGap = nvis::norm( params.the_metric.displacement(px0,lastSeg[1]) );
        
        
        //Test if this segment is 'viable' := ok to add (no section sep and satisfies curve-refinement)
        bool viableSegment = true;
        //Checking if segment is viable with repsect to section transverality violations
        if (viableSegment) {
        
            //Check for a discontinuity point, Are Node 0 and Node 1 a discontinuity? (Unlikely without subdivision)
            if(adaptiveEdge.isDiscontinuity(IdType(0,0)) ||
                    adaptiveEdge.isDiscontinuity(IdType(1,0)) ) {
                std::cerr << "MapManifold: WeakManifold detected node discontinuity without subdivision!\n";
                viableSegment = false;
            }
            //Also, check if it has separation (unlikely unless without subdivision)
            else if ( detectSectionSeparation(x0,x,cache,theMap,std::abs(thePeriod),params) ) {
                std::cerr << "MapManifold: WeakManifold detected section separation without subdivision!\n";
                viableSegment = false;
            }
            //Check if we could not subdivide this segment based on u_min
            // -> i.e., does it satisfy curve-refinement criteria?
            else {
                //Compare current segment to last segment for curve refinement (skips section sep)
                bool curveOK = isManifoldCurveOK( theMap, params, lastSeg[0], xf, px, settings);
                if(!curveOK) {
                    viableSegment = false;
                }
            }
        }
        //It's not viable if we have a transversality violation
        //Note:  There is no check for "weakStrength" manifolds for failing
        //curve-refinement criteria with previous step.  We don't have any information
        //regarding the previous segment.
        
        //Debugging
        if(debug) {
            std::cout << "--------------------------------------------------------\n";
            std::cout << " Manifold " << manifoldID << ":\n";
            std::cout << "   Weak-manifold Point insertion...\n";
            std::cout << "   Working Segment [" << workingSegID << "] Points:\n";
            std::cout << "    prevManifold = " << prevManifold << "\n";
            std::cout << "    nextManifold = " << nextManifold << "\n";
            std::cout << "    x0 = " << workSeg[0] << "  tau0 = " << 0.0 << "\n";
            std::cout << "    x1 = " << workSeg[1] << "  tau1 = " << 1.0 << "\n";
            std::cout << "    depth = " << wSegDepth << "\n";
            std::cout << "    dTauMin = " << dTauMin << " | u_min = " << settings.delta_tau_min << "\n";
            std::cout << "    lastSegGap = " << lastSegGap << "\n";
            if (viableSegment) {
                std::cout << "   Inserting Segment [ " << (int)segments.size() << "]\n";
            } else {
                std::cout << "   Unviable New Segment [ " << (int)segments.size() << "]\n";
            }
            std::cout << "    p = " << period << " (prop = " <<  thePeriod << ")\n";
            std::cout << "    parentID = " << workSeg.segID << " \n";
            std::cout << "    phi0 = " << xf << "  dt0 = " << px0MapData[numSing] << "\n";
            std::cout << "    phi1 = P(x1)" << px << "  dt1 = " << pxMapData[numSing] << "\n";
            std::cout << "    firstSubdivision = " << firstSubdivision << "\n";
            std::cout << "    FixedTauMin = " << useFixedTauMin << "\n";
            std::cout << "--------------------------------------------------------\n";
        }
        
        //Create the new 'viable' segment and add to container
        if (viableSegment) {
            int nextSegID = (int) segments.size();
            ManifoldSeg newSeg(
                nextSegID,    //Segment Index
                xf,   px,     //Positions P^p(x0),P^p(x1)
                0.0, 1.0,     //Linear Parameters from source that make this segment
                px0MapData, pxMapData,   //MapDataType (CA,time) for relevant return
                workSeg, lastSeg );
            //Add to segments container
            segments.push_back( newSeg );
            //Add new segment to manifold tree (via working segment)
            segments[workingSegID].addChild( nextSegID );
            //Update the last segment with new information
            segments[lastSegID].next = nextSegID;
            lastSegID = nextSegID; //Advance last segment to the new segments
            numSegments++;
            //Reserve some more space for segments
            if (segments.size() >= segments.capacity()) {
                int sz = (int) segments.size() + 2500;
                segments.reserve(sz);
            }
            //Check if there is an option for the next depth
            std::map<int,std::vector<int> >::iterator it = idsAtDepth.find(wSegDepth+1);
            if(it != idsAtDepth.end() ) {
                //Add this segment to depth
                it->second.push_back( lastSegID );
            } else {
                //We have to create a new vector for a new depth
                std::vector<int> segAtNewDepth;
                segAtNewDepth.push_back(lastSegID);
                idsAtDepth.insert( std::pair<int,std::vector<int> >(wSegDepth+1,segAtNewDepth) );
            }
            //Update arc-length
            totalArcLength += segments[lastSegID].length(params);
            bool segInsideBounds = true;
            if (!bbox.inside(x) || !bbox.inside(px) ) {
                segInsideBounds = false;
            }
            if (segInsideBounds) {
                progressLength += segments[lastSegID].length(params);
            }
        }
        // else, it is not a viable segment and so we just progress to next working segment
        
        //Note: Set child manifold data (minors and mirrors)
        
    } else {
        //Otherwise add all computed and viable segments to container
        //if(debug) std::cerr << " advectSegment(): Adding viable segments...\n";
        std::vector<AESegment> aeSegs;
        adaptiveEdge.getAllSegments(aeSegs);
        typename std::vector<AESegment>::iterator segIT,segITm1, segITp1;
        segIT = segITm1 = aeSegs.begin();
        int numNewSegments = 0;
        //Check each potential segment generated from current AdaptiveEdge  ----------------------------------
        for(segIT=aeSegs.begin(); segIT!=aeSegs.end(); ++segIT) {
            //Other segments
            segITm1 = segIT;
            --segITm1;
            segITp1 = segIT;
            ++segITp1;
            //Need a special setting if its the first segment of current AdaptiveEdge
            bool firstSeg = (segIT==aeSegs.begin());
            bool lastSeg = nvis::all(segIT->second == IdType(1,0));
            //Get parameters
            double tauA = adaptiveEdge.getParameter(segIT->first);
            double tauB = adaptiveEdge.getParameter(segIT->second);
            //Get Point information
            vec_type xA = adaptiveEdge.getVertex(segIT->first);
            vec_type xB = adaptiveEdge.getVertex(segIT->second);
            vec_type xC = adaptiveEdge.getVertex(
                              (firstSeg)? segITp1->second : segITm1->first);
            //Linear Parameter distance
            double dTau = tauB - tauA;
            //Check if potential segment is viable (no transversality violations at nodes or between nodes)
            bool viableSegment = true;
            //If we detect separation between end nodes:
            if (viableSegment) { //all segments run this
                //Check for a discontinuity point, Are Node 0 and Node 1 a discontinuity?
                if(adaptiveEdge.isDiscontinuity(segIT->first) ||
                        adaptiveEdge.isDiscontinuity(segIT->second)) {
                    viableSegment = false;
                }
            }
            //Run the downstream curve-refinement checks
            if (viableSegment) {
                //Get map data (downstream)
                vec_type pA = adaptiveEdge.getValue(segIT->first);
                vec_type pB = adaptiveEdge.getValue(segIT->second);
                vec_type pC = adaptiveEdge.getValue(
                                  (firstSeg)? segITp1->second : segITm1->first);
                                  
                bool curveOK = false;
                //If primary segment (A-B) has separation, then curve is not ok
                if( !detectSectionSeparation(xA,xB,cache,theMap,std::abs(thePeriod),params) ) {
                    if(firstSeg) {
                        //On first seg, we have to use the next segment's results to test curve
                        curveOK = isManifoldCurveOK(
                                      theMap,std::abs(thePeriod),params,xA,xB,xC,pA,pB,pC,settings,cache,true);
                        //otherwise, curve is not ok
                    } else {
                        //Test the previous segment (C-A) for separation
                        bool prevSegHasSep =
                            detectSectionSeparation(xC,xA,cache,theMap,std::abs(thePeriod),params);
                        //Also, test if the previous segment (C-A) is too short to use safely (<delta_min)
                        bool prevSegIsTooShort =
                            isSegmentTooShort(theMap,params,pC,pA,settings);
                        if(prevSegHasSep || prevSegIsTooShort) {
                            if(!lastSeg) {
                                //Try to use the next segment to test segment
                                vec_type xD = adaptiveEdge.getVertex(segITp1->second);
                                vec_type pD = adaptiveEdge.getValue(segITp1->second);
                                curveOK = isManifoldCurveOK(
                                              theMap,std::abs(thePeriod),params,xA,xB,xD,pA,pB,pD,settings,cache,true);
                            }
                        } else {
                            //Standard segment test
                            curveOK = isManifoldCurveOK(
                                          theMap,std::abs(thePeriod),params,xC,xA,xB,pC,pA,pB,settings,cache,false);
                        }
                    }
                }
                //If passed curve tests (with section sep), then it is a viable segment to add:
                viableSegment = curveOK;
                
            }
            //IF viable segment, store information
            if (viableSegment) {
                //Get position data (on working segment)
                VecType xA = adaptiveEdge.getVertex(segIT->first);
                VecType xB = adaptiveEdge.getVertex(segIT->second);
                //Get the Poincare map point (downstream)
                VecType pA = adaptiveEdge.getValue(segIT->first);
                VecType pB = adaptiveEdge.getValue(segIT->second);
                //Get the Poincare map data
                cit = cache.find( SRDATA(xA) );
                if (cit == cache.end()) {
                    throw std::runtime_error("Could not find requested map in SRDATA cache");
                }
                const SRDATA& xAData = (*cit);
                MapDataType pxAMapData = xAData.getDataAtReturn(std::abs(thePeriod));
                cit = cache.find( SRDATA(xB) );
                if (cit == cache.end()) {
                    throw std::runtime_error("Could not find requested map in SRDATA cache");
                }
                const SRDATA& xBData = (*cit);
                MapDataType pxBMapData = xBData.getDataAtReturn(std::abs(thePeriod));
                
                //Values for setting the new segments
                const ManifoldSeg& workSeg = segments[workingSegID]; //read-only!
                const ManifoldSeg& lastSeg = segments[lastSegID]; //read-only!
                VecType p0 = pA;
                MapDataType p0MapData = pxAMapData;
                //In case we are on the first new segment of this Manifold,
                //we use the last available information
                if (firstSeg && workingSegID == 0) {
                    p0 = lastSeg[1];
                    p0MapData = lastSeg.data1;
                    //Bad to use beyond first constructed segment!
                    //If there is a transversality violations on the first point or between the
                    //last segment and the current segment, we create an erroneous segment!!!
                    // ->It rarely happens, but it is a possible event!
                }
                
                //Display information in debug mode
                if(debug) {
                    std::cout << "--------------------------------------------------------\n";
                    std::cout << " Manifold " << manifoldID << ":\n";
                    std::cout << "   Working Segment [" << workingSegID << "] Points:\n";
                    //std::cout << "    prevManifold = " << prevManifold << "\n";
                    //std::cout << "    nextManifold = " << nextManifold << "\n";
                    std::cout << "    x0 = " << xA << "  tau0 = " << tauA << "\n";
                    std::cout << "    x1 = " << xB << "  tau1 = " << tauB << "\n";
                    std::cout << "    depth = " << wSegDepth << "\n";
                    std::cout << "    dTauMin = " << dTauMin << " | u_min = " << settings.delta_tau_min << "\n";
                    std::cout << "   Inserting Segment [ " << lastSegID + 1 << "]\n";
                    std::cout << "    p = " << period << " (prop = " <<  thePeriod << ")\n";
                    std::cout << "    parentID = " << workSeg.segID << " \n";
                    std::cout << "    phi0 = " << p0 << "  dt0 = " << p0MapData[numSing] << "\n";
                    std::cout << "    phi1 = " << pB << "  dt1 = " << pxBMapData[numSing] << "\n";
                    std::cout << "--------------------------------------------------------\n";
                }
                
                //Create next segment object (With parent and previous inputs)
                int nextSegID = (int) segments.size();
                ManifoldSeg newSeg( nextSegID,
                                    p0, pB, tauA, tauB, p0MapData, pxBMapData,
                                    workSeg, lastSeg );
                segments.push_back( newSeg );
                //Update current last segment with new information
                segments[lastSegID].next = nextSegID;
                lastSegID = nextSegID; //Advance last segment
                numSegments++;
                numNewSegments++;
                //Reserve space for more segments if necessary
                if (segments.size() >= segments.capacity()) {
                    int sz = (int) segments.size() + 2500;
                    segments.reserve(sz);
                }
                //Add children to working segment
                segments[workingSegID].addChild( nextSegID );
                
                //Update progress length
                totalArcLength += segments[lastSegID].length(params);
                bool segInsideBounds = true;
                if (!bbox.inside(pA) || !bbox.inside(pB) ) {
                    segInsideBounds = false;
                }
                if (segInsideBounds) {
                    progressLength += segments[lastSegID].length(params);
                }
                
                //Note: Set child manifold data (minors and mirrors)
            }
            //Otherwise, move to the next segment and don't add this one to list!
            
        } // End each segment loop
        
        //Add segments to idsAtDepth based on children:
        //Check if there is an option for the next depth
        std::map<int,std::vector<int> >::iterator it = idsAtDepth.find(wSegDepth+1);
        if(it != idsAtDepth.end() ) {
            //Reserve space for new elements
            int sz = (int) it->second.size();
            it->second.reserve(sz + numNewSegments);
            //Add these child segments to depth
            for(int i=0; i<numNewSegments; i++) {
                it->second.push_back( segments[workingSegID].children[i] );
            }
        } else {
            //We have to create a new vector for a new depth
            std::vector<int> segsAtNewDepth = segments[workingSegID].children;
            idsAtDepth.insert( std::pair<int,std::vector<int> >(wSegDepth+1,segsAtNewDepth) );
        }
        
        //Check for a large number of new segments which indicates a major subdivision:
        // Notes:  There are multiple checks here to exit manifolds that do lots of subdivision.
        //         However, only the first is enabled because the others are semi-included
        //         inside standard ManifoldSettings now.
        //         This check preps the algoritm to quit at the end of current depth.
        //
        // 1) If over 1000, likely too many to process next depth
        // 2) If past first depth, stop on next depth if generating over 500 new segments
        // 3) Stop the process shortly after generating more than 5000 segments (10,000 segments
        //    for a single manifold can be overbearing for visualization)
        // Future Work:  Probably should put in a ManifoldSettings member to enable this feature
        // during the progress check in ManifoldData.hpp.
        if ( (!majorSubdivision && numNewSegments >= 1000) ) {
            //|| (!majorSubdivision && getDepth(workingSegID) > 1 && numNewSegments >= 500) ||
            //(!majorSubdivision && numSegments >= 5000) ) {
            //A major subdivision has occurred and we need to trigger the manifold to stop
            //once the current depth level is finished processing.
            majorSubdivision = true;
            majorSubDepth = getDepth(workingSegID);
        }
        
        //Check if there are no new segments added at 'new' depth (noNewSegsAtNewDepth = true)
        if (numNewSegments > 0) {
            noNewSegsAtNewDepth = false;
        } else {
            noNewSegsAtNewDepth = true;
        }
    } //end advection type conditional
    
    //--------------------------------------------------------------------
    
    //Move the Working Segment to the next in the linked list
    if (lastSegID == workingSegID ) {
        //If this occurs and we are either
        //    1) on weakStrength manifold before 1st subdivision OR
        //    2) below depth=2 OR workSegID <= 5
        //and no new segments are created
        if (!weakStrength &&
                (getDepth(workingSegID)<=2 || workingSegID<=3 ) && //Can happen more with troublesome starters
                noNewSegsAtNewDepth && !useFixedTauMin) {
            //Return back to the first segment of depth = 1 and try using dTauMin as u_min
            // i.e., disable the adaptive measure
            useFixedTauMin = true;
            
            //Restart working segment
            //workingSegID = segments[0].children[0];
            const ManifoldSeg& theNextWorkSeg = segments[workingSegID];
            //seedLength += segments[workingSegID].length(params);
            numSubdivisions = 0;
            workingSegInitialized = false;
            //numSeedSegments++;
            //AdapitveEdge reset
            VecType x0 = theNextWorkSeg[0];
            VecType x1 = theNextWorkSeg[1];
            adaptiveEdge.reset(x0,x1);
            adaptiveEdge.refine(IdType(0,0)); //Default one refine command
            VecType xm = (x0+x1) / 2.0;
            adaptiveEdge.setVertex(IdType(1,1),xm);
            
            //Prompt user
            std::cerr << " -----------------------------------------------------------------------------------------\n";
            std::cerr << " Warning:  Difficulting in subdiving Depth = " << getDepth(workingSegID) << " due to u_min.\n";
            std::cerr << " Attemping more subdivision to create more segments for manifold:\n";
            info();
            std::cerr << " -----------------------------------------------------------------------------------------\n";
            return;
            
        } else if (weakStrength && !firstSubdivision && !useFixedTauMin) {
            //This is the weakStrength Manifold version reset:
            //Try using a fixed dTauMin (i.e., disable the adaptive measure)
            useFixedTauMin = true;
            //Force the weakStrength Manifold to know a subdivision has occurred
            firstSubdivision = true;
            //Stay at the same working segment:
            const ManifoldSeg& theNextWorkSeg = segments[workingSegID];
            //seedLength += segments[workingSegID].length(params);
            numSubdivisions = 0;
            workingSegInitialized = false;
            //numSeedSegments++;
            //AdapitveEdge reset
            VecType x0 = theNextWorkSeg[0];
            VecType x1 = theNextWorkSeg[1];
            adaptiveEdge.reset(x0,x1);
            adaptiveEdge.refine(IdType(0,0)); //Default one refine command
            VecType xm = (x0+x1) / 2.0;
            adaptiveEdge.setVertex(IdType(1,1),xm);
            
            //Prompt user
            std::cerr << " -----------------------------------------------------------------------------------------\n";
            std::cerr << " Warning:  Difficulting in subdiving WorkingSegment on WEAK MANIFOLD due to u_min.\n";
            std::cerr << " Attemping more subdivision to create more segments for manifold:\n";
            info();
            std::cerr << " -----------------------------------------------------------------------------------------\n";
            return;
            
        } else {
            //Can't do anything else and we have to exit
            std::cerr << " -----------------------------------------------------------------------------------------\n";
            info();
            std::cerr << "    Unable to move to next working segment! workingSegID = lastSegID = " << lastSegID << "\n";
            std::cerr << "    You may need to decrease u_min parameter to encourage more subdivision!\n";
            std::cerr << " -----------------------------------------------------------------------------------------\n";
            //throw std::runtime_error("Unable to move to next working segment!");
            
            //If this occurs, we can't advance the manifold anymore
            isComplete = true;
            return;
            
        }
    }
    // ->Note: assumes that at least 1 segment is added!
    workingSegID = segments[workingSegID].next;
    const ManifoldSeg& theNextWorkSeg = segments[workingSegID];
    seedLength += segments[workingSegID].length(params);
    numSubdivisions = 0;
    workingSegInitialized = false;
    numSeedSegments++;
    //AdapitveEdge reset
    VecType x0 = theNextWorkSeg[0];
    VecType x1 = theNextWorkSeg[1];
    adaptiveEdge.reset(x0,x1);
    adaptiveEdge.refine(IdType(0,0)); //Default one refine command
    VecType xm = (x0+x1) / 2.0;
    adaptiveEdge.setVertex(IdType(1,1),xm);
    if((!weakStrength && getDepth(workingSegID)>2) || (weakStrength && firstSubdivision) ) {
        useFixedTauMin = false;
    }
    
    //Progressed to next working segment
}


/// Get a delta_tau_min parameter based on current working segment edge length and allowed map displacement u_min
template <class SRDATA, class FP>
template <class PARAM>
double MapManifold<SRDATA,FP>::
getAdaptiveDeltaTauMin(const PARAM& params, const double& u_min) const
{
    //Force usage of u_min as dTauMin to encourage subdivision
    //if (useFixedTauMin) return u_min;
    //if (useFixedTauMin) return 1.e-2;  //Utilize up to 100 steps
    if (useFixedTauMin) {
        return 1.e-4;    //Utilize up to 10000 steps
    }
    //An adaptive form based on map length that speeds up computation
    double workSegLength = segments[workingSegID].length(params);
    return (u_min / workSegLength);
}

///Add any children to parent processing queue if they grow up
template <class SRDATA, class FP>
void MapManifold<SRDATA,FP>::
childMaturity(std::queue<int>& processingQueue)
{

}


/// Write manifold data from curve-refinement advection for a vector container of MapManifolds
template <class SRDATA, class FP>
void MapManifold<SRDATA,FP>::write(FILE* f) const
{
    //Assume file is already open for writing
    
    //Write common elements for manifold
    // UUID - unique identifier
    std::ostringstream ss;
    ss << id;
    fprintf(f,"Manifold: %s\n",ss.str().c_str());
    // orbitID fpID mType childType numSegs hasParent weakStrength workingSegIdx totalArcLength progressLength seedLength isComplete
    int numSegs = (int) segments.size();
    int pertInt = pertToInt(mType);
    int cInt = childToInt(childType);
    int hpInt = (hasParent()) ? 1 : 0;
    int numCIDs = (int)complementaryIDs.size();
    fprintf(f,"%d %d %d %d %d %d %d %d %d %d %d %d %d %.15f %.15f %.15f \n",
            fpdOrbitIdx, fpdPointIdx, pertInt, cInt,
            numSegs, hpInt, prevManifold, nextManifold, depthPS,
            (weakStrength)?1:0, workingSegID, numCIDs, homoclinicPartnerID,
            totalArcLength, progressLength, seedLength);
    // parentID
    fprintf(f,"ParentManifoldIndex= %d\n",parent);
    
    //Write the complementaryIDs for companion manifold indexes
    for (int cid=0; cid<numCIDs; cid++) {
        fprintf(f,"%d ", complementaryIDs[cid]);
    }
    if(numCIDs>0) {
        fprintf(f,"\n");
    }
    
    //Write the segments
    for (int ll=0; ll<numSegs; ll++) {
        segments[ll].write(f);
    }
    
} //End writeMapManifolds()

/// Load manifold data from file into a vector container of MapManifolds
template<class SRDATA, class FP>
void MapManifold<SRDATA,FP>::
read(FILE* f, pmate::FixedPointData<FP>* fpDataPtr)
{
    //Assume file is already open
    char buffer[200];
    //For creating uuids with strings (from file)
    boost::uuids::string_generator string_gen;
    
    //Load unique id
    fscanf(f,"%s ",buffer); //'Manifold: '
    fscanf(f,"%s",buffer);  //The actual uuid string
    id = string_gen(buffer);
    
    //Load common data and set members
    int tempBool = 0;
    int numSegs = 0;
    int pertInt = 0, cInt=0, hpInt=0, numCIDs=0;
    fscanf(f,"%d %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf",
           &fpdOrbitIdx, &fpdPointIdx, &pertInt, &cInt,
           &numSegs, &hpInt, &prevManifold, &nextManifold, &depthPS,
           &tempBool, &workingSegID, &numCIDs, &homoclinicPartnerID,
           &totalArcLength, &progressLength, &seedLength);
    //Types from integers
    mType = intToPert(pertInt);
    childType = intToChild(cInt);
    //Set baseFixedPoint
    baseFixedPoint =  fpDataPtr->getFixedPoint(fpdOrbitIdx,fpdPointIdx) ;
    
    //Other values inferred from this information
    firstStepCount = 1;
    numSegments = numSegs;
    numSeedSegments = workingSegID+1;
    weakStrength = (tempBool==1);
    isComplete = true; //Reset for restart
    workingSegInitialized = false; //Assume for restart
    fwd = (mType == UNSTABLE_PLUS || mType == UNSTABLE_MINUS) ? true : false;
    p_or_m = (mType == STABLE_PLUS || mType == UNSTABLE_PLUS) ? true : false;
    setPeriod();
    
    //Initialize idsAtDepth map container
    std::vector<int> startingDepthVec;
    startingDepthVec.push_back(0); //Initial working seg id
    idsAtDepth.insert( std::pair<int,std::vector<int> >(1,startingDepthVec) );
    
    //Read parentID
    fscanf(f,"%s %d",buffer,&parent);
    
    //Read complementaryIDs
    for(int cid=0; cid<numCIDs; cid++) {
        int id = 0;
        fscanf(f,"%d",&id);
        complementaryIDs.push_back(id);
    }
    
    //Load Segments
    for(int ll=0; ll<numSegs; ll++) {
        //Make segment with read constructor
        ManifoldSeg newSeg(f); //Calls read commands
        
        if(ll>0) {
            //Last Segment Info
            segments.back().next = newSeg.segID;
            //Current Segment Info
            newSeg.previous = segments.back().segID;
        }
        
        //Add to container
        segments.push_back( newSeg );
        
        //Add segments to idsAtDepth based on children:
        int wSegDepth = getDepth(newSeg.segID);
        //Check if there is an option for the next depth
        std::map<int,std::vector<int> >::iterator it = idsAtDepth.find(wSegDepth+1);
        if(it != idsAtDepth.end() ) {
            //Reserve space for new elements
            int sz = (int) it->second.size();
            int numKids = (int) newSeg.children.size();
            it->second.reserve(sz + numKids);
            //Add these child segments to depth
            for(int i=0; i<numKids; i++) {
                it->second.push_back( newSeg.children[i] );
            }
        } else {
            //We have to create a new vector for a new depth
            std::vector<int> segsAtNewDepth = newSeg.children;
            idsAtDepth.insert( std::pair<int,std::vector<int> >(wSegDepth+1,segsAtNewDepth) );
        }
        
    }
    
    //Set Last Segment ID
    lastSegID = (int) segments.size() - 1;
    
    
} //End readMapManifolds()

/// Get a list of all Map Disconts from segments
template<class SRDATA, class FP>
void MapManifold<SRDATA,FP>::
getMapDiscontinuityList(std::list<typename MapManifold<SRDATA,FP>::MapDiscont>& mdList) const
{
    typename std::list<MapDiscont>::const_iterator mdIT;
    int numSegs = (int) segments.size();
    for (int ll=0; ll<numSegs; ll++) {
        for (mdIT = segments[ll].discontinuityList.begin();
                mdIT !=segments[ll].discontinuityList.end(); mdIT++) {
            mdList.push_back( (*mdIT) );
        }
    }
}

/// Get the current working segment from list
template<class SRDATA, class FP>
typename MapManifold<SRDATA,FP>::ManifoldSeg&
MapManifold<SRDATA,FP>::getWorkingSegment()
{
    return segments[workingSegID];
}

/// Test if a segment exists
template<class SRDATA, class FP>
bool MapManifold<SRDATA,FP>::segmentExists(const int segID) const
{
    if (segID >= 0 && segID < (int) segments.size() ) {
        return true;
    }
    //Otherwise no
    return false;
}

/// Test if the next segment exists
template<class SRDATA, class FP>
bool MapManifold<SRDATA,FP>::nextSegmentExists(const int segID) const
{
    //First check if segment exists
    if (!segmentExists(segID)) {
        return false;
    }
    //Get the next id
    int nextID = segments[segID].next;
    //Then return it's existence test
    return segmentExists(nextID);
}

/// Test if the previous segment exists
template<class SRDATA, class FP>
bool MapManifold<SRDATA,FP>::previousSegmentExists(const int segID) const
{
    //First check if segment exists
    if (!segmentExists(segID)) {
        return false;
    }
    //Get the previous id
    int prevID = segments[segID].previous;
    //Then return it's existence test
    return segmentExists(prevID);
}

/// Check if there is a connection between a segment and it's next segment
template<class SRDATA, class FP>
bool MapManifold<SRDATA,FP>::nextConnection(const int segID) const
{
    //First check if segment and next segment exists
    if (!nextSegmentExists(segID)) {
        return false;    //Nothing to do
    }
    //Test if Point 1 of segment (segID) and Point 0 of next are the same:
    int nextID = segments[segID].next;
    return nvis::all( segments[segID][1] == segments[nextID][0] );
}


/// Check if there is a connection between a segment and it's previous segment
template<class SRDATA, class FP>
bool MapManifold<SRDATA,FP>::previousConnection(const int segID) const
{
    //First check if segment and previous segment exists
    if (!previousSegmentExists(segID)) {
        return false;    //Nothing to do
    }
    //Test if Point 0 of segment (segID) and Point 1 of previous are the same:
    int prevID = segments[segID].previous;
    return nvis::all( segments[segID][0] == segments[prevID][1] );
}




/// Manifold Perturbation Type to integer
template<class SRDATA, class FP>
int MapManifold<SRDATA,FP>::pertToInt(const enum xavier::Perturbation& pert) const
{
    int out = 0;
    switch(pert) {
        case UNSTABLE_PLUS :
            out = 0;
            break;
        case UNSTABLE_MINUS :
            out = 1;
            break;
        case STABLE_PLUS :
            out = 2;
            break;
        case STABLE_MINUS :
            out = 3;
            break;
    }
    return out;
}

/// PertInt to string
template<class SRDATA, class FP>
std::string MapManifold<SRDATA,FP>::pIntToString(const int& i) const
{
    std::string out("");
    switch(i) {
        case 0 :
            out = std::string("UNSTABLE_PLUS");
            break;
        case 1 :
            out = std::string("UNSTABLE_MINUS");
            break;
        case 2 :
            out = std::string("STABLE_PLUS");
            break;
        case 3 :
            out =  std::string("STABLE_MINUS");
            break;
    }
    return out;
}

/// Integer to Manifold Perturbation Type
template<class SRDATA, class FP>
#ifdef _WIN32
xavier::Perturbation MapManifold<SRDATA,FP>::intToPert(const int& i) const
#else
enum xavier::Perturbation MapManifold<SRDATA,FP>::intToPert(const int& i) const
#endif
{
    Perturbation out = UNSTABLE_PLUS;
    switch(i) {
        case 0 :
            out = UNSTABLE_PLUS;
            break;
        case 1 :
            out = UNSTABLE_MINUS;
            break;
        case 2 :
            out = STABLE_PLUS;
            break;
        case 3 :
            out = STABLE_MINUS;
            break;
    }
    return out;
}

/// Manifold Child Type to integer
template<class SRDATA, class FP>
int MapManifold<SRDATA,FP>::childToInt(const enum topology::ManifoldChildType& mct) const
{
    int out = 0;
    switch(mct) {
        case ADULT :
            out = 0;
            break;
        case MINOR :
            out = 1;
            break;
        case MIRROR :
            out = 2;
            break;
        case ADULTMINOR :
            out = 3;
            break;
        case ADULTMIRROR :
            out = 4;
            break;
    }
    return out;
}

/// Integer to Manifold Child Type
template<class SRDATA, class FP>
#ifdef _WIN32
topology::ManifoldChildType MapManifold<SRDATA,FP>::intToChild(const int& i) const
#else
enum topology::ManifoldChildType MapManifold<SRDATA,FP>::intToChild(const int& i) const
#endif
{
    ManifoldChildType out = ADULT;
    switch(i) {
        case 0 :
            out = ADULT;
            break;
        case 1 :
            out = MINOR;
            break;
        case 2 :
            out = MIRROR;
            break;
        case 3 :
            out = ADULTMINOR;
            break;
        case 4 :
            out = ADULTMIRROR;
            break;
    }
    return out;
}

template<class SRDATA, class FP>
std::string MapManifold<SRDATA,FP>::cIntToString(const int& i) const
{
    std::string out("");
    switch(i) {
        case 0 :
            out = std::string("ADULT");
            break;
        case 1 :
            out = std::string("MINOR");
            break;
        case 2 :
            out = std::string("MIRROR");
            break;
        case 3 :
            out = std::string("ADULTMINOR");
            break;
        case 4 :
            out = std::string("ADULTMIRROR");
            break;
    }
    return out;
}


// TREE TRAVERSAL -----------------------------------------------------------------------------------
///Get the accumulated UPSTREAM sum of a data value (recursive, usually TOF)
template<class SRDATA, class FP>
double MapManifold<SRDATA,FP>::
getDataSum(const int segID, const double tau, const int dataID) const
{
    //Check inputs and gather segment
    assert(segID >=0 && segID < (int) segments.size());
    const ManifoldSeg& thisSeg = segments[segID];
    assert(dataID >= 0 && dataID < (int) thisSeg.data0.size());
    assert(tau >=0.0 && tau <= 1.0);
    
    double dataValue = (1.0-tau)*thisSeg.data0[dataID] + tau*thisSeg.data1[dataID];
    
    //If segment has no parents, return this value
    if(!(thisSeg.hasParent())) {
        return dataValue;
    } else {
        //If we are not at the parent seg, we have to go up the tree
        //and get the parent values (and add to sum).
        
        //Compute the parent tau
        double nu = (1.0-tau)*thisSeg.tau0 + tau*thisSeg.tau1;
        //Recursive call
        return dataValue + getDataSum(thisSeg.parent, nu, dataID);
        
    }
}

/// Get the time of flight (which is usually dataID=2, tau is on [0,1] for segments[segID])
template<class SRDATA, class FP>
double MapManifold<SRDATA,FP>::
getTimeOfFlight(const int segID, const double tau) const
{
    //Check inputs and gather segment
    assert( segmentExists(segID) );
    const ManifoldSeg& thisSeg = segments[segID];
    assert(tau >=0.0 && tau <= 1.0);
    int dataID = 2; //Assume Time-of-flight is data[2]
    
    double timeValue = (1.0-tau)*thisSeg.data0[dataID] + tau*thisSeg.data1[dataID];
    
    /*
    //Track time of flight upstream until reaching a known d_ps
    int currentDepth = getDepth(segID);
    if ((currentDepth <= depthPS) || !(thisSeg.hasParent())) {
        //Return the interpolated progenitor state time-of-flight
        return (1.0-tau)*(thisSeg.ps0.tof) + tau*(thisSeg.ps1.tof);*/
    if (!(thisSeg.hasParent())) {
        //Return the value of time-of-flight on 0th segment
        //(which is 0.0 with assumption that 0th segment IS on the fixed point)
        return 0.0;
    } else {
        //Recursively add interpolated time of flight to next level up:
        //Compute the parent tau
        double nu = (1.0-tau)*thisSeg.tau0 + tau*thisSeg.tau1;
        return timeValue + getTimeOfFlight(thisSeg.parent, nu);
    }
}

/// Get the minimum UPSTREAM data value (recursive)
template<class SRDATA, class FP>
double MapManifold<SRDATA,FP>::
getMinData(const int segID, const double tau, const int dataID) const
{
    //Check inputs and gather segment
    assert(segID >=0 && segID < (int) segments.size());
    const ManifoldSeg& thisSeg = segments[segID];
    assert(dataID < (int) thisSeg.data0.size());
    assert(tau >=0.0 && tau <= 1.0);
    
    double dataValue = (1.0-tau)*thisSeg.data0[dataID] + tau*thisSeg.data1[dataID];
    
    //If segment has no parents, return this value
    if(!(thisSeg.hasParent())) {
        return dataValue;
    } else {
        //If we are not at the parent seg, we have to go up the tree
        //and get the parent values (and add to sum).
        
        //Compute the parent tau
        double nu = (1.0-tau)*thisSeg.tau0 + tau*thisSeg.tau1;
        //Recursive call to test minimum
        return  std::min( dataValue , getMinData(thisSeg.parent, nu, dataID) );
        
    }
}

/// Get the linear parameter for indicated UPSTREAM segment or initial segment (recursive)
// Note: u is on [0,1] and is the parameter indicating the point location on segID
template<class SRDATA, class FP>
double MapManifold<SRDATA,FP>::
getSourceLinearParameter(const int segID, const double u, const int stopID) const
{
    //Note: this goes UP the tree!
    
    //Check inputs and gather segment
    assert(segID >=0 && segID < (int) segments.size());
    const ManifoldSeg& thisSeg = segments[segID];
    assert(u >=0.0 && u <= 1.0);
    
    //Return the linear parameter
    if(!(thisSeg.hasParent()) || segID==stopID ) {
        return u;
    } else {
        //Evaluate the linear parameter upstream (using current segment u values)
        double nu = (1.0-u)*thisSeg.tau0 + u*thisSeg.tau1;
        //Move up
        return getSourceLinearParameter( thisSeg.parent, nu, stopID );
    }
    
}


/// Track the destinations of a source segment (recursive) [FIX for efficiency, may be buggy!]
template<class SRDATA, class FP>
void MapManifold<SRDATA,FP>::
trackSourceDownstream(
    const int sourceID,
    const int segID,
    std::vector<int>& downstreamIDs,
    std::vector< std::pair<double,double> >& tauAtSource
) const
{
    //Check inputs
    assert(sourceID >= 0 && sourceID < (int) segments.size());
    assert(segID >= 0 && segID < (int) segments.size());
    const ManifoldSeg& thisSeg = segments[segID];
    
    //Only continue if children exist
    if( thisSeg.hasKids() ) {
        //Store children for this segment.  For each child
        for(int k=0; k<(int)thisSeg.children.size(); k++) {
            //Add the child value
            int kidID = thisSeg.children[k];
            downstreamIDs.push_back( kidID );
            
            //Lookup tau on source segment and store
            double tauA = getSourceLinearParameter(kidID,0.0,sourceID);
            double tauB = getSourceLinearParameter(kidID,1.0,sourceID);
            tauAtSource.push_back( std::pair<double,double>(tauA,tauB) );
        }
        
        //Run downstream elements (helps order data better as a separate loop)
        for(int k=0; k<(int)thisSeg.children.size(); k++) {
            int kidID = thisSeg.children[k];
            const ManifoldSeg& kidSeg = segments[kidID];
            
            //Call tracking again for child
            trackSourceDownstream( sourceID, kidID, downstreamIDs, tauAtSource );
        }
        
    }
}

/// Get the depth of current segment or whole manifold tree
template<class SRDATA, class FP>
int MapManifold<SRDATA,FP>::
getDepth(const int segID) const
{
    //Check inputs and gather segment
    int sID = segID;
    if (segID == -1) {
        sID = segments.back().segID;
    }
    assert(sID >=0 && sID < (int) segments.size());
    const ManifoldSeg& thisSeg = segments[sID];
    
    int depth = 1;
    
    //If segment has no parents, return this value
    if(!(thisSeg.hasParent())) {
        return depth;
    } else {
        //If we are not at the parent segment, we have to go up the tree
        //and get the depth value (and add to sum).
        
        //Recursive call
        return depth + getDepth(thisSeg.parent);
        
    }
}

/// Get the actual number of layers of subdivision (skipping single-segment layers)
template<class SRDATA, class FP>
int MapManifold<SRDATA,FP>::
getNumSubdivisionLayers(const int segID) //const
{
    //Check inputs and gather segment
    int sID = segID;
    if (segID == -1) {
        sID = segments.back().segID;
    }
    assert(sID >=0 && sID < (int) segments.size());
    const ManifoldSeg& thisSeg = segments[sID];
    
    int layer = 0;
    
    //We have to look at all segments at the current depth
    std::vector<int> idsAtDepth;
    getAllSegmentIDsAtDepth(idsAtDepth,getDepth(segID));
    int maxKids = 0;
    bool kidsOnDepth = false;
    for(int j=0; j<(int)idsAtDepth.size(); j++) {
        if(segments[idsAtDepth[j]].hasKids()) {
            kidsOnDepth = true;
            int numKids = (int)segments[idsAtDepth[j]].children.size();
            if ( numKids > maxKids) {
                maxKids = numKids;
            }
        }
    }
    
    
    //Only a viable layer if there is more than one child,
    // but counts for any segment with children at the current depth
    if (kidsOnDepth) {
        if ( maxKids >= 2 ) {
            layer = 1;
        }
        //Lots of subdivision happening (trying to shorten computation)
        if ( maxKids >= 10 ) {
            layer = 2;    //Double count on purpose
        }
        if ( maxKids >= 20 ) {
            //Rather large layer that encompasses lots of flowering
            layer = 1 + maxKids / 10; //Compute as multiples of 10
            //Hopefully, this forces a stop condition as layers this deep are difficult to process
            //in a reasonable time span (< 1 day)
        }
    }
    
    //If segment has no parents
    if (!(thisSeg.hasParent())) {
        return layer;
    } else {
        //Not at 0th segment, so go up tree and add other viable layers.
        //Recursive call
        return layer + getNumSubdivisionLayers(thisSeg.parent);
    }
    
}

/// Gather all segments at the current depth level
template<class SRDATA, class FP>
void MapManifold<SRDATA,FP>::
getAllSegmentIDsAtDepth(std::vector<int>& ids, const int depth) //const //For c++ 11 only
{
    ids.clear();
    int d = depth;
    if (d==-1) {
        //Get the current depth level
        d = getDepth( -1 );
    }
    
    //Just return the vector from the id map
    ids = idsAtDepth[d];
    
    /* Search is insanely in efficient!
    //Search from front of segment list since (more segments are in higher depths)
    typename SegmentVector::iterator rit;
    //for(auto rit = segments.cbegin(); rit != segments.cend(); ++rit) { //C++ 11
    for (rit = segments.begin(); rit != segments.end(); ++rit) {
        int thisDepth = getDepth( rit->segID );
        if (thisDepth == d) {
            ids.push_back( rit->segID );
        } else if (thisDepth > d)
            break; //Done with the loop as all remaining levels are lower
     }*/
    
}

/// Gather all segments at a depth level (with segID specified, hopefully faster)
template<class SRDATA, class FP>
void MapManifold<SRDATA,FP>::
getAllSegmentIDsOnSameDepth(std::vector<int>& ids, const int segID) //const //For c++ 11 only
{
    ids.clear();
    int d = getDepth(segID);
    
    //Just return the vector from the id map
    ids = idsAtDepth[d];
    
    //Search is insanely inefficient!
    /*int s0 = segID;
    if(s0 == -1) s0 = (int) segments.size() - 1;
    std::list<int> idList;
    idList.push_back(s0);
    int s = s0;
    //Go previous until at another depth
    while (previousSegmentExists(s) ) {
      int prevID = segments[s].previous;
      if (getDepth(prevID) == d) {
        idList.push_back( prevID );
        s = prevID;
      } else break; //Otherwise done going back
    }
    //Go next until at another depth (or DNE)
    s = s0;
    while (nextSegmentExists(s) ) {
      int nextID = segments[s].next;
      if (getDepth(nextID) == d) {
        idList.push_back( nextID );
        s = nextID;
      } else break; //Otherwise done going forward
    }
    //Sort list
    idList.sort();
    std::list<int>::iterator it;
    for(it=idList.begin();it!=idList.end();++it) ids.push_back( *it );*/
}


// PROGENITOR STATE FUNCTIONS ------------------------------------------------------------------------
/// Is the requested progenitor state available (meaning viable solutions on both endpoints of segment)
template<class SRDATA, class FP>
bool MapManifold<SRDATA,FP>::isPSFound(const int segID, const double tau) const
{
    //Check inputs and gather segment
    assert(segID >=0 && segID < (int) segments.size());
    const ManifoldSeg& thisSeg = segments[segID];
    assert(tau >=0.0 && tau <= 1.0);
    int currentDepth = getDepth(segID);
    
    //Return the boolean
    bool ok = false;
    if (tau == 0.0) {
        ok = thisSeg.ps0.found;
    } else if (tau == 1.0) {
        ok = thisSeg.ps1.found;
    } else if (thisSeg.ps0.found && thisSeg.ps1.found) {
        ok = true;
    }
    return ok;
}

/*
/// Get the Progenitor State DeltaV (Arrival or Departure)
template<class SRDATA, class FP>
typename MapManifold<SRDATA,FP>::VecType
MapManifold<SRDATA,FP>::getPSDeltaV(const int segID, const double tau) const
{
  //Check inputs and gather segment
  assert(segID >=0 && segID < (int) segments.size());
  const ManifoldSeg& thisSeg = segments[segID];
  assert(tau >=0.0 && tau <= 1.0);
  int currentDepth = getDepth(segID);
  //Return the interpolated progenitor state dv
  return (1.0-tau)*(thisSeg.ps0.dv) + tau*(thisSeg.ps1.dv);
}

/// Get the Progenitor State Alpha Parameter (location on orbit)
template<class SRDATA, class FP>
double MapManifold<SRDATA,FP>::getPSAlpha(const int segID, const double tau) const
{
  //Check inputs and gather segment
  assert(segID >=0 && segID < (int) segments.size());
  const ManifoldSeg& thisSeg = segments[segID];
  assert(tau >=0.0 && tau <= 1.0);
  //Return the interpolated progenitor state alpha
  return (1.0-tau)*(thisSeg.ps0.alpha) + tau*(thisSeg.ps1.alpha);
}

/// Get the Progenitor State Orbit State (full state on orbit)
template<class SRDATA, class FP>
typename MapManifold<SRDATA,FP>::StateType
MapManifold<SRDATA,FP>::getPSOrbitState(const int segID, const double tau) const
{
  //Check inputs and gather segment
  assert(segID >=0 && segID < (int) segments.size());
  const ManifoldSeg& thisSeg = segments[segID];
  assert(tau >=0.0 && tau <= 1.0);
  if (theOrbit.isInit() && theOrbit.isValid()) {
    //Sample from the PeriodicOrbit object
    double alpha = getPSAlpha(segID,tau);
    StateType orbitState(0.0);
    theOrbit.getStateFromAlpha(alpha,orbitState);
    return orbitState;
  } else {
    //Return the interpolated progenitor state alpha
    return (1.0-tau)*(thisSeg.ps0.orbitState) + tau*(thisSeg.ps1.orbitState);
  }
}


/// Solve Progenitor State TPBVP with Single Shooting
template<class SRDATA, class FP>
template<class PMAP>
bool MapManifold<SRDATA,FP>::
solveProgenitorState(
    PMAP& pmap,                //The Map engine
    const double& alpha0,     //Initial guess of orbit parameter
    const nvis::vec2& dv0,    //Initial guess of delta-V vector from orbit to ProgenitorState
    const double& dt0,        //Initial guess for time of flight ProgenitorState to Manifold Point
    const typename MapManifold<SRDATA,FP>::VecType& phif, //The manifold point
    typename MapManifold<SRDATA,FP>::ProgState& output, //The Progenitor State Object result
    const double eps,  //Convergence tolerance
    const bool verbose,
    const int maxiter
) const {

    double hamD = pmap.rhs().desired_hamiltonian();

    //Setup Progenitor State Corrector:
    orbital::ProgenitorCorrector<PMAP,SelfType> corrector(pmap,(*this),alpha0,dv0,dt0,phif,hamD,verbose);
    corrector.setMaxIters(maxiter);
    corrector.setTolerance(eps);


    int thread_id = 0;
#if _OPENMP
    thread_id = omp_get_thread_num();
#endif
    corrector.setJobIndex(thread_id);
    orbital::CorrectionResult theResult;
    theResult = corrector.correct();
    bool found = theResult.converged;
    if (!found) {
        //Let's hope this never happens!!!
        //std::runtime_error rt("Unable to locate Progenitor State");
        //throw rt;
        return false;
    }
    else {
       //Result must pass sanity checks
       //if (!checkOK) return false;
    }

    //Gather the result
    corrector.getSolution(output);

    //If no errors,
    return true;
}

/// Propagate to the Prgenitor State (upstream) given a segment index (also returns iterates)
template<class SRDATA, class FP>
template<class PMAP>
void MapManifold<SRDATA,FP>::
propagateToProgenitorState(PMAP& theMap, const int segID, const double tau,
                                std::vector<typename PMAP::state_type>& states,
                                std::vector<double>& times,
                                std::vector<typename PMAP::return_type>& iterData) const
{
  typedef typename PMAP::return_type   return_type;
  typedef typename PMAP::return_state  return_state;

  //Gather the progenitor state time of flight (PS to PhiF)
  double tof = getTimeOfFlight(segID,tau);
  //Get the initial manifold point
  VecType phif = segments[segID].getPoint(tau);
  //Run the propagation
  //bool integOK = true;
  std::vector<return_state> internalData;
  try {
    theMap.flow_map(phif,-tof,iterData,internalData); //Opposite ToF to go from PhiF to PS
  } catch(...) {
     // Will return with no data in states and times!
     // Check propagation by seeing if 'times' has data.
     return;
  }

  //Otherwise set the states and times:
  for(int i=0;i<(int)internalData.size();i++) {
    states.push_back( internalData[i].x );
    times.push_back( internalData[i].t );
  }
}  */


/// Get UPSTREAM map points utilizing linear interpolation from manifold (returns 0->d order)
// Note: tau is on [0,1] indicating the point on segments[segID]
template<class SRDATA, class FP>
void MapManifold<SRDATA,FP>::
getUpstreamPoints(const int segID, const double tau,
                  std::vector<typename MapManifold<SRDATA,FP>::VecType>& mapStates) const
{
    //Assertions
    assert( segmentExists(segID) );
    assert( (tau >= 0.0 && tau <= 1.0) );
    
    //First we need the depth of the current segment
    //int sDepth = getDepth(segID);
    
    //Get the current depth point
    VecType pt0 = segments[segID].getPoint(tau);
    //Start a list because we have to reverse the order to make it ascending
    std::list<VecType> tempList;
    tempList.push_back(pt0);
    //For each depth level up, lookup map point
    double tempTau = tau;
    int tempSegID = segID;
    while ( segments[tempSegID].hasParent() ) { //currentDepth to 0
        //Get the linear parameter of the segment one level up
        int parentID = segments[tempSegID].parent;
        double nu = getSourceLinearParameter(tempSegID,tempTau,parentID);
        // Get the point of the parent (loop stops at d=1, so always a parent)
        VecType pt = segments[parentID].getPoint(nu);
        tempList.push_back(pt);
        //Move up level
        tempSegID = parentID;
        tempTau = nu;
    }
    
    //The list is in opposite order of desired output, so reverse for output:
    typename std::list<VecType>::reverse_iterator rit;
    for(rit=tempList.rbegin(); rit!=tempList.rend(); rit++) {
        mapStates.push_back( (*rit) );
    }
}

/// Get DOWNSTREAM map points utilizing linear interpolation from manifold
template<class SRDATA, class FP>
void MapManifold<SRDATA,FP>::
getDownstreamPoints(const int segID, const double tau,
                    std::vector<typename MapManifold<SRDATA,FP>::VecType>& mapStates) const
{
    //Assertions
    assert( segmentExists(segID) );
    assert( (tau >= 0.0 && tau <= 1.0) );
    //Note: It's possible that there are no downstream points as the transversality violations
    //may have caused no segment to exist that has the desired tau value
    
    //Do NOT add the initial selection point here! (Compatibility with getAllMapPoints() )
    
    //Loop until no child segments are found
    int cSegID = segID;
    double tempTau = tau;
    bool downstreamExists = true;
    while( segments[cSegID].hasKids() && downstreamExists ) {
        //Loop through kids until found viable range, or exit if no viable segments
        int numKids = (int) segments[cSegID].children.size();
        int downstreamID = -1;
        for(int i=0; i<numKids; i++) {
            int childID = segments[cSegID].children[i];
            //If we find the viable segment, save point and progress
            if ( segments[childID].isOnSegment(tempTau) ) {
                downstreamID = childID;
                //Get the point based on parent segment tau value
                double cSegTau = segments[childID].getLinearParam(tempTau); //returns on [0,1]
                VecType pt = segments[childID].getPoint(cSegTau);
                mapStates.push_back(pt);
                //Advance info downstream
                tempTau = cSegTau;
                cSegID = downstreamID;
                //End this do loop
                break;
            }
        }
        //If no segment found, we exit
        if(downstreamID < 0) {
            downstreamExists = false;
        }
    }
}

/// Crop segments AFTER a given segment id
template<class SRDATA, class FP>
template<class PARAM>
void MapManifold<SRDATA,FP>::
crop(
    const PARAM& params,               //Poincare map parameters
    const ManifoldSettings& settings,  //Manifold Settings for viable bounds
    const int lsegID                   //The last segment ID to keep
)
{
    //This only runs if there are segments AFTER lsegID
    int numSegs = (int) segments.size();
    double progDeficit = 0.0, totalDeficit = 0.0;
    if (lsegID < numSegs-1 ) {
        //Gather information from each segment being cropped
        for(int i=lsegID+1; i<numSegs; i++) {
            //Remove from kids in parent segs
            if(segments[i].hasParent()) {
                int pSegID = segments[i].parent;
                std::vector<int>::iterator sit;
                for(sit=segments[pSegID].children.begin(); sit!=segments[pSegID].children.end(); ++sit) {
                    //Pull out this kid
                    if (*sit == i) {
                        segments[pSegID].children.erase(sit);
                        --sit; //jump back to move to next element correctly
                    }
                }
            } //End of kid removal loop
            
            //Total length and progress length deficits
            totalDeficit += segments[i].length(params);
            bool segInsideBounds = true;
            if (!settings.bounds.inside(segments[i][0]) ||
                    !settings.bounds.inside(segments[i][1]) ) {
                segInsideBounds = false;
            }
            if (segInsideBounds) {
                progDeficit += segments[i].length(params);
            }
        }
        
        //Reset Manifold parameters based on pulling out information
        if(workingSegID > lsegID) {
            workingSegID = lsegID;
        }
        lastSegID = lsegID;
        totalArcLength -= totalDeficit;
        progressLength -= progDeficit;
        numSegments = lsegID + 1;
        
        //Erase all segments after lsegID
        segments.erase(segments.begin() + lsegID + 1, segments.end());
        
    }
}

/// Clear all segments and reset manifold like new
template<class SRDATA, class FP>
void MapManifold<SRDATA,FP>::
clearManifold()
{
    //Wipe out all segments
    segments.clear();
    //Reset data parameters as if we are starting over
    workingSegID = -1;
    workingSegInitialized = false;
    lastSegID = -1;
    totalArcLength = 0.0;
    progressLength = 0.0;
    seedLength = 0.0;
    numSegments = 0;
    numSeedSegments = 0;
    firstSubdivision = false;
    majorSubdivision = false;
    majorSubDepth = -1;
    noNewSegsAtNewDepth = true;
    useFixedTauMin = false;
    isComplete = false;  //The most important reset
    //nextManifold = -1;
    //prevManifold = -1;
    firstStepCount = 0;
    numSubdivisions = 0;
    emptyNodes = 0;
    adaptiveEdge.reset(baseFixedPoint.pos,baseFixedPoint.pos);
}

/// Get ALL map points for a selected arc (both upstream and downstream, in tree-depth order)
template<class SRDATA, class FP>
void MapManifold<SRDATA,FP>::
getAllMapPoints(const int segID, const double tau,
                std::vector<typename MapManifold<SRDATA,FP>::VecType>& mapStates) const
{
    //Calls both the upstream and downstream entries and places in ascending order
    std::vector<VecType> ups, downs;
    //Get Upstream states
    getUpstreamPoints(segID,tau,ups);
    //Get Downstream states
    getDownstreamPoints(segID,tau,downs);
    //Concatenate vectors:
    mapStates = ups;
    for(int i=0; i<(int)downs.size(); i++) {
        mapStates.push_back(downs[i]);
    }
}

/// Propagate UPSTREAM to top of manifold tree employing interpolated points (i.e., to depth=0)
template<class SRDATA, class FP>
template<class PMAP>
void MapManifold<SRDATA,FP>::
propagateUpstream(PMAP& theMap, const int segID, const double tau,
                  std::vector<typename PMAP::state_type>& states,
                  std::vector<double>& times ) const
{
    //Assertions
    assert( segmentExists(segID) );
    assert( (tau >= 0.0 && tau <= 1.0) );
    //Type defs
    typedef typename PMAP::return_state  ReturnState;
    typedef typename PMAP::gvec_type     gvec_type;
    //Get the upstream points correlated with the selected tau (on [0,1]) for segments[segID]
    std::vector<VecType> ups;
    getUpstreamPoints(segID,tau,ups);
    states.clear();
    times.clear();
    double dtTotal = 0.0;
    //Store the first state (Never stored within the mapping function)
    gvec_type y0 = theMap.section().unproject(ups[0]);
    states.push_back( theMap.getState(y0) );
    times.push_back( 0.0 );
    
    
    //For all but last point, propagate upstream to the next p
    int numReturns = (int) ups.size();
    //for(int k=numReturns-1; k>0; k--) { //For upstream propagations (aren't very reliable)
    for(int k=0; k<(numReturns-1); k++) { //For propagating DOWNSTREAM (which conforms with setup more appropriately)
        //Output
        std::vector<ReturnState> itersOut, internalOut;
        bool valid = true;
        
        //Call the map with manifold information
        try {
            //Propagate upstream  (upstream -> -1*manifold_period)
            //theMap.flow_map(ups[k], itersOut, internalOut, -thePeriod );
            //Propagate downstream
            theMap.flow_map(ups[k], itersOut, internalOut, thePeriod );
            
        } catch(...) {
            //Mapping Error (shouldn't happen)
            valid = false;
        }
        
        
        
        //Skip storage if mapping error (but shouldn't happen
        if(!valid) {
            continue;
        }
        
        //Store to overall output
        for(int i=0; i<(int)internalOut.size(); i++) {
            //Skip the first state (which is map state) after 1st mapping
            //if(k<(numReturns-1) && i==0) continue; NOte: flow_map doesn't store first state
            //Store points
            states.push_back( internalOut[i].x );
            times.push_back( internalOut[i].t + dtTotal);
        }
        //Update the total time of flight occurred
        dtTotal += internalOut.back().t - internalOut[0].t;
    }
    
    //Note on propagation:
    // Upstream - Points are stored from initial tau to 0th level
    // Downstream - Points are stored from 0th level to current depth
}

/// Propagate UPSTREAM to top of manifold tree employing interpolated points (i.e., to depth=0)
template<class SRDATA, class FP>
template<class PMAP>
void MapManifold<SRDATA,FP>::
propagateUpstream(PMAP& theMap, const int segID, const double tau,
                  std::vector<typename PMAP::state_type>& states, std::vector<double>& times,
                  std::vector<typename PMAP::lvec_type>& iters,   std::vector<double>& itTimes,
                  std::vector<typename PMAP::lvec_type>& subIts,  std::vector<double>& subitTimes) const
{
    //Assertions
    assert( segmentExists(segID) );
    assert( (tau >= 0.0 && tau <= 1.0) );
    
    //Types
    typedef typename PMAP::return_state  ReturnState;
    typedef typename PMAP::gvec_type     gvec_type;
    typedef typename PMAP::lvec_type     lvec_type;
    typedef typename PMAP::lmat_type     lmat_type;
    
    //Get the upstream points correlated with the selected tau (on [0,1]) for segments[segID]
    getUpstreamPoints(segID,tau,iters);
    states.clear();
    times.clear();
    double dtTotal = 0.0;
    itTimes.push_back(0.0);
    //Store the first state (Never stored within the mapping function)
    gvec_type y0 = theMap.section().unproject(iters[0]);
    states.push_back( theMap.getState(y0) );
    times.push_back( 0.0 );
    
    
    //For all but last point, propagate upstream to the next p
    int numReturns = (int) iters.size();
    //for(int k=numReturns-1; k>0; k--) { //For upstream propagations (aren't very reliable)
    for(int k=0; k<(numReturns-1); k++) { //For propagating DOWNSTREAM (which conforms with setup more appropriately)
        //Output
        std::vector<ReturnState> itersOut, internalOut;
        bool valid = true;
        
        //Call the map with manifold information
        try {
            //Propagate upstream  (upstream -> -1*manifold_period)
            //theMap.flow_map(iters[k], itersOut, internalOut, -thePeriod );
            //Propagate downstream
            theMap.flow_map(iters[k], itersOut, internalOut, thePeriod );
        } catch(...) {
            //Mapping Error (shouldn't happen)
            valid = false;
        }
        
        //Skip storage if mapping error (but shouldn't happen
        if(!valid) {
            continue;
        }
        
        //Store the sub returns if the period is greater than 1
        int p = abs(thePeriod);
        if (p > 1) {
            //Grab up to pth return
            for(int i=0; i<(p-1); i++) {
                std::pair<lvec_type,lmat_type> tempPair = theMap.section().project(itersOut[i].getState());
                subIts.push_back( tempPair.first );
                subitTimes.push_back( dtTotal + itersOut[i].t );
            }
        }
        
        //Store to overall output
        for(int i=0; i<(int)internalOut.size(); i++) {
            //Skip the first state (which is map state) after 1st mapping
            //if(k<(numReturns-1) && i==0) continue; NOte: flow_map doesn't store first state
            //Store points
            states.push_back( internalOut[i].x );
            times.push_back( internalOut[i].t + dtTotal);
        }
        //Update the total time of flight occurred
        dtTotal += internalOut.back().t - internalOut[0].t;
        itTimes.push_back( dtTotal );
    }
    
    //Note on propagation:
    // Upstream - Points are stored from initial tau to 0th level
    // Downstream - Points are stored from 0th level to current depth
}

/// Propagate DOWNSTREAM to bottom of manifold tree employing interpolated points (i.e., to depth=0)
template<class SRDATA, class FP>
template<class PMAP>
void MapManifold<SRDATA,FP>::
propagateDownstream(PMAP& theMap, const int segID, const double tau,
                    std::vector<typename PMAP::state_type>& states,
                    std::vector<double>& times ) const
{
    //Assertions
    assert( segmentExists(segID) );
    assert( (tau >= 0.0 && tau <= 1.0) );
    //Typedefs
    typedef typename PMAP::return_state  ReturnState;
    typedef typename PMAP::gvec_type     gvec_type;
    //Get the current point
    std::vector<VecType> downs;
    VecType pt0 = segments[segID].getPoint(tau);
    downs.push_back(pt0);
    
    //Gather the remaining points downstream
    getDownstreamPoints(segID,tau,downs);
    states.clear();
    times.clear();
    double dtTotal = 0.0;
    gvec_type y0 = theMap.section().unproject(downs[0]);
    states.push_back( theMap.getState(y0) );
    times.push_back( 0.0 );
    
    //For all but last point, propagate downstream to the next p
    for(int k=0; k<((int)downs.size()-1); k++) {
        //Output
        std::vector<ReturnState> itersOut, internalOut;
        bool valid = true;
        
        //Call the map with manifold information
        try {
            theMap.flow_map(downs[k], itersOut, internalOut, thePeriod );
        } catch(...) {
            //Mapping Error (shouldn't happen)
            valid = false;
        }
        
        //Skip storage if mapping error (but shouldn't happen
        if(!valid) {
            continue;
        }
        
        //Store to overall output
        for(int i=0; i<(int)internalOut.size(); i++) {
            //Skip the first state (which is map state) after 1st mapping
            if(k>0 && i==0) {
                continue;
            }
            //Store points
            states.push_back( internalOut[i].x );
            times.push_back( internalOut[i].t + dtTotal);
        }
        //Update the total time of flight occurred
        dtTotal += internalOut.back().t - internalOut[0].t;
    }
    
}

/// Propagate DOWNSTREAM to bottom of manifold tree employing interpolated points (i.e., to depth=0), also gather iterates
template<class SRDATA, class FP>
template<class PMAP>
void MapManifold<SRDATA,FP>::
propagateDownstream(PMAP& theMap, const int segID, const double tau,
                    std::vector<typename PMAP::state_type>& states, std::vector<double>& times,
                    std::vector<typename PMAP::lvec_type>& iters,   std::vector<double>& itTimes,
                    std::vector<typename PMAP::lvec_type>& subIts,  std::vector<double>& subitTimes) const
{
    //Assertions
    assert( segmentExists(segID) );
    assert( (tau >= 0.0 && tau <= 1.0) );
    //Type defs
    typedef typename PMAP::return_state  ReturnState;
    typedef typename PMAP::gvec_type     gvec_type;
    typedef typename PMAP::lvec_type     lvec_type;
    typedef typename PMAP::lmat_type     lmat_type;
    //Get the current point
    VecType pt0 = segments[segID].getPoint(tau);
    iters.push_back(pt0);
    itTimes.push_back(0.0);
    
    //Gather the remaining points downstream
    getDownstreamPoints(segID,tau,iters);
    states.clear();
    times.clear();
    double dtTotal = 0.0;
    gvec_type y0 = theMap.section().unproject(iters[0]);
    states.push_back( theMap.getState(y0) );
    times.push_back( 0.0 );
    
    //For all but last point, propagate downstream to the next p
    for(int k=0; k<((int)iters.size()-1); k++) {
        //Output
        std::vector<ReturnState> itersOut, internalOut;
        bool valid = true;
        
        //Call the map with manifold information
        try {
            theMap.flow_map(iters[k], itersOut, internalOut, thePeriod );
        } catch(...) {
            //Mapping Error (shouldn't happen)
            valid = false;
        }
        
        //Skip storage if mapping error (but shouldn't happen
        if(!valid) {
            continue;
        }
        
        //Store the sub returns if the period is greater than 1
        int p = abs(thePeriod);
        if (p > 1) {
            //Grab up to pth return
            for(int i=0; i<(p-1); i++) {
                std::pair<lvec_type,lmat_type> tempPair = theMap.section().project(itersOut[i].getState());
                subIts.push_back( tempPair.first );
                subitTimes.push_back( dtTotal + itersOut[i].t );
            }
        }
        
        //Store to overall output
        for(int i=0; i<(int)internalOut.size(); i++) {
            //Skip the first state (which is map state) after 1st mapping
            if(k>0 && i==0) {
                continue;
            }
            //Store points
            states.push_back( internalOut[i].x );
            times.push_back( internalOut[i].t + dtTotal);
        }
        //Update the total time of flight occurred
        dtTotal += internalOut.back().t - internalOut[0].t;
        itTimes.push_back( dtTotal );
    }
    
}

/// Propagate DOWNSTREAM to bottom of manifold tree employing interpolated points (i.e., to depth=0), also gather iterates
template<class SRDATA, class FP>
template<class PMAP>
void MapManifold<SRDATA,FP>::
propagateDownstream(PMAP& theMap, const int segID, const double tau, const int extendedIterates,
                    std::vector<typename PMAP::state_type>& states, std::vector<double>& times,
                    std::vector<typename PMAP::lvec_type>& iters,   std::vector<double>& itTimes,
                    std::vector<typename PMAP::lvec_type>& subIts,  std::vector<double>& subitTimes) const
{
    //Assertions
    assert( segmentExists(segID) );
    assert( (tau >= 0.0 && tau <= 1.0) );
    //Type defs
    typedef typename PMAP::return_state  ReturnState;
    typedef typename PMAP::gvec_type     gvec_type;
    typedef typename PMAP::lvec_type     lvec_type;
    typedef typename PMAP::lmat_type     lmat_type;
    //Get the current point
    VecType pt0 = segments[segID].getPoint(tau);
    iters.push_back(pt0);
    itTimes.push_back(0.0);
    
    //Gather the remaining points downstream
    getDownstreamPoints(segID,tau,iters);
    states.clear();
    times.clear();
    double dtTotal = 0.0;
    gvec_type y0 = theMap.section().unproject(iters[0]);
    states.push_back( theMap.getState(y0) );
    times.push_back( 0.0 );
    
    //For all but last point, propagate downstream to the next p
    for(int k=0; k<((int)iters.size()-1); k++) {
        //Output
        std::vector<ReturnState> itersOut, internalOut;
        bool valid = true;
        
        //Call the map with manifold information
        try {
            theMap.flow_map(iters[k], itersOut, internalOut, thePeriod );
        } catch(...) {
            //Mapping Error (shouldn't happen)
            valid = false;
        }
        
        //Skip storage if mapping error (but shouldn't happen
        if(!valid) {
            continue;    //May need to break loop here!!!!! See other functions too!!!!!
        }
        
        //Store the sub returns if the period is greater than 1
        int p = abs(thePeriod);
        if (p > 1) {
            //Grab up to pth return
            for(int i=0; i<(p-1); i++) {
                std::pair<lvec_type,lmat_type> tempPair = theMap.section().project(itersOut[i].getState());
                subIts.push_back( tempPair.first );
                subitTimes.push_back( dtTotal + itersOut[i].t );
            }
        }
        
        //Store to overall output
        for(int i=0; i<(int)internalOut.size(); i++) {
            //Skip the first state (which is map state) after 1st mapping
            if(k>0 && i==0) {
                continue;
            }
            //Store points
            states.push_back( internalOut[i].x );
            times.push_back( internalOut[i].t + dtTotal);
        }
        //Update the total time of flight occurred
        dtTotal += internalOut.back().t - internalOut[0].t;
        itTimes.push_back( dtTotal );
    }
    
    //Extended iterates (if greater than 0)
    if (extendedIterates > 0) {
        for(int k=0; k<extendedIterates; k++) {
            nvis::vec2 seedPt = iters.back();
            //Output
            std::vector<ReturnState> itersOut, internalOut;
            bool valid = true;
            
            //Call the map with manifold information
            try {
                theMap.flow_map(seedPt, itersOut, internalOut, thePeriod );
            } catch(...) {
                //Mapping Error (shouldn't happen)
                valid = false;
            }
            
            //Skip storage if mapping error (but shouldn't happen
            if(!valid) {
                continue;    //May need to break loop here!!!!! See other functions too!!!!!
            }
            
            //Store the sub returns if the period is greater than 1
            int p = abs(thePeriod);
            if (p > 1) {
                //Grab up to pth return
                for(int i=0; i<(p-1); i++) {
                    std::pair<lvec_type,lmat_type> tempPair = theMap.section().project(itersOut[i].getState());
                    subIts.push_back( tempPair.first );
                    subitTimes.push_back( dtTotal + itersOut[i].t );
                }
            }
            
            //Store to overall output
            for(int i=0; i<(int)internalOut.size(); i++) {
                //Skip the first state (which is map state) after 1st mapping
                if(k>0 && i==0) {
                    continue;
                }
                //Store points
                states.push_back( internalOut[i].x );
                times.push_back( internalOut[i].t + dtTotal);
            }
            //Add the last point to the iterates
            std::pair<lvec_type,lmat_type> tempPair = theMap.section().project(itersOut.back().getState());
            iters.push_back( tempPair.first );
            //Update the total time of flight occurred
            dtTotal += internalOut.back().t - internalOut[0].t;
            itTimes.push_back( dtTotal );
        }
    }
    
}


} // end topology

#endif   // MAP_MANIFOLD_HPP
