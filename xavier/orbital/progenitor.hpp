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


/////////////////////////////////////////////////////////////////
//  progenitor.hpp header file
//  Author:  Wayne Schlei
//  Date: 5/28/2016
//
//  Purpose:  Use differential corrections to solve for the
//  progenitor state (initial location on a periodic orbit)
//  for a given manifold state.
//
//  Returns: alpha (linear param), dV (dxdot,dydot) and TOF
/////////////////////////////////////////////////////////////////
#ifndef PROGENITOR_HPP
#define PROGENITOR_HPP

//Useful elements
#include <vector>
#include <iostream>
#include <memory>
#include <complex>
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>

//EIGEN API
#include <Eigen/Core>
#include <Eigen/Dense>
//#include <Eigen/Sparse>
//#include <Eigen/UmfPackSupport>

//OpenMP
#if _OPENMP
#include <omp.h>
#endif

//API for Maps
#include <math/angle.hpp>
#include <maps/poincare_map.hpp>
#include <maps/section.hpp>
#include <maps/fixpoints.hpp>
#include <orbital/controller.hpp>
#include <design/Trajectory.hpp>
#include <design/PeriodicOrbit.hpp>


#ifdef _WIN32
#include <boost/math/special_functions/fpclassify.hpp> // isnan for Windows
#endif

using namespace Eigen;
using namespace xavier;

namespace orbital {

/** \brief Progenitor State Data Object
 *  Data locating the progenitor state (or source point on periodic orbit) for
 *  an invariant manifold point on a surface of section.
 */
template <class STATE, class VEC>
class ProgenitorState {
public :
    /// Constructor
    ProgenitorState() : alpha(0.0), dv(0.0), tof(0.0), orbitState(0.0), found(false)
    {}
    
    /// The Progenitor State parameter on orbit
    double alpha;
    /// The departure/arrival delta-v vector at the Progenitor State from the orbit to manifold arc
    VEC dv;
    /// The time of flight (non-dimensional) from Progenitor State to Manifold State (negative for StableManifold)
    double tof;
    /// The departure/arrival state on the periodic orbit (to save from re-propagating on load for lookups)
    STATE orbitState;
    /// Boolean indicating if the Progenitor State computation succeeded
    bool found;
};


/** \brief Data storage for a Progenitor State computation
 *  Storage of Manifold State, manID, segID, pointID, and trajectory history
 *  relevant to the computation of a Progenitor State.
 */
template<class MAP>
class ProgenitorCompData {
public :
    typedef typename MAP::value_type                value_type;
    typedef typename MAP::rhs_type                  rhs_type;
    typedef typename rhs_type::state_type           state_type;
    typedef typename MAP::xstate_type               xstate_type;
    typedef typename MAP::section_type              section_type;
    typedef typename section_type::lvec_type        lvec_type;
    typedef typename section_type::lmat_type        lmat_type;
    typedef typename section_type::gvec_type        gvec_type;
    typedef typename section_type::gmat_type        gmat_type;
    typedef typename section_type::lgMatrix         lgmat_type;
    typedef pmateDesign::Trajectory<state_type>     Trajectory;
    
    ///Constructor
    ProgenitorCompData(const lvec_type& phif, const int mID, const int sID, const int pID) :
        path(), manID(mID), segID(sID), ptID(pID), phi(phif),
        progState()
    {}
    
    
    /// Trajectory object for storing estimated path 'back' to periodic orbit
    Trajectory  path;
    /// Manifold ID
    int manID;
    /// Segment ID
    int segID;
    /// Point ID on segment (either 0 or 1)
    int ptID;
    /// Manifold Point on Map
    lvec_type phi;
    
    ///Results of Progenitor Computation:
    ProgenitorState<state_type,lvec_type> progState;
    
};

/** \brief Class for sorting intersections in progenitor state generation.
 *  Data class that stores information from intersection tests between
 *  manifold arcs and the sire periodic orbit.  This also comes with
 *  an less-than ordering function that ranks intersections by
 *  1) The tangentiality (angle between velocity vectors)
 *  2) The time of flight (magnitude, always >=0.0)
 *  Note: this is to sift through POSSIBLE progenitor states that are
 *  generated through intersection testing. (Used ONLY in ManifoldData::evaluateProgenitorState())
 */
template<class VEC>
class ProgenitorIntersection {
public :
    typedef ProgenitorIntersection<VEC>  SelfType;
    
    /// Constructor
    ProgenitorIntersection() :
        pos(0.0), mV(0.0), oV(0.0), tof(0.0), alpha(0.0) {}
    ProgenitorIntersection(const VEC& x, const VEC& manVec, const VEC& orbitVec, const double& dt, const double& alf) :
        pos(x),
        mV(manVec),
        oV(orbitVec),
        tof(fabs(dt)),
        alpha(alf)
    {}
    
    /// The location of the intersection (position space)
    VEC pos;
    /// The velocity vector from the manifold arc
    VEC mV;
    /// The velocity vector from the periodic orbit
    VEC oV;
    /// The time of flight (magnitude, always positive) from manifold point to the progenitor state
    double tof;
    /// The linear parameter indicating the progenitor state location on the periodic orbit (time factor on [0,1])
    double alpha;
    
    /// Less than operator for working with sorting - smallest values have the best delta-V's and time of flights
    bool operator<(const SelfType& other) const
    {
        //First order by deltaV
        double dv = getDeltaVMag();
        double odv = other.getDeltaVMag();
        if (dv != odv ) {
            return (dv<odv);
        }
        //Then order by time of flight
        else if (tof != other.tof) {
            return (tof<other.tof);
        }
        //Then just order with position
        else {
            nvis::lexicographical_order compare;
            return compare(pos,other.pos);
        }
    }
    
    ///Compute the delta-V between velocity vectors (Always Orbit-To-ManifoldArc)
    VEC getDeltaV() const
    {
        return mV - oV;
    }
    ///Compute the Delta-V magnitude
    double getDeltaVMag() const
    {
        return nvis::norm( getDeltaV() );
    }
    
    ///Compute the angle differential between velocity vectors (angular separation in rads)
    double getDeltaVAngle() const
    {
        return unsigned_angle<2>(mV,oV);
    }
    
};

/** \brief Corrections Class for finding a Progenitor State for a manifold arc
 *  Corrections Class that locates a Progenitor (or origin) State
 *  on a periodic orbit for a given invariant manifold arc and a guess.
 *  Process is an implicit mapping single-shooting method that solves
 *  M=2 equations given N=3 unknowns.
 */
template<class MAP, class MANIFOLD>
class ProgenitorCorrector {
public :
    /// Type defs
    typedef typename MAP::value_type                value_type;
    typedef typename MAP::rhs_type                  rhs_type;
    typedef typename rhs_type::state_type           state_type;
    typedef typename MAP::xstate_type               xstate_type;
    typedef typename MAP::section_type              section_type;
    typedef typename section_type::lvec_type        lvec_type;
    typedef typename section_type::lmat_type        lmat_type;
    typedef typename section_type::gvec_type        gvec_type;
    typedef typename section_type::gmat_type        gmat_type;
    typedef typename section_type::lgMatrix         lgMatrix;
    typedef typename MAP::return_type               return_type;
    typedef typename MAP::return_state              return_state;
    typedef typename return_state::mat_type         return_mat;
    typedef typename MAP::solver_type               solver_type;
    typedef typename solver_type::step              step_type;
    typedef CorrectionsRegulator                    reg_type;
    //typedef Eigen::Triplet<double>                  Triple;
    //typedef Eigen::SparseMatrix<double>             SpMat;
    
    /// Constructor
    ProgenitorCorrector(
        MAP& pmap,                 //Map engine
        const MANIFOLD& mman,      //MapManifold object
        const value_type& alpha0,  //Initial linear paramter along orbit to indicate
        const lvec_type& dv0,      //Initial delta-V for depature onto manifold at progenitor state
        const value_type& dt,      //Projected time of flight FROM progenitor state TO manifold point
        const lvec_type& phiT,     //Manifold for which we are computing the progenitor state
        const double& Hval,        //Desired Hamiltonian value for this manifold arc
        bool verbose=true
    ) :
        _verbose(verbose), outputF(false), outputDF(false), isSolutionFound(false),
        alpha(alpha0), dv(dv0), tof(dt),
        _maxIters(20), zeroTimeTol(1.e-6), _tol(1.e-10), _prec(1.e-12),
        theMap(&pmap), theManifold(&mman), desired_H(Hval), jobID(0),
        X(VectorXd::Zero(1)), F(VectorXd::Zero(1)), dX(VectorXd::Zero(1)),
        Xstore(VectorXd::Zero(1))
    {
        setPrecision(_prec);
        //Have to call "setGuess()" to setup _x0/_period correctly
    }
    
    //Parameters:
    /// Maximum iterations
    void setMaxIters(int maxIters)
    {
        _maxIters = maxIters;
    }
    /// Convergence Tolerance
    void setTolerance(double tol)
    {
        _tol = tol;
    }
    /// Time-Period (if you know it)
    void setPeriod(double T)
    {
        tof = T;
    }
    /// Precision of integration
    void setPrecision(double& prec)
    {
        _prec = prec;
        theMap->setPrecision(_prec);
    }
    /// Job index (To tell which run is being displayed)
    void setJobIndex(const int id)
    {
        jobID = id;
    }
    
    /// Set an initial guess
    void setGuess(const double& alpha0, const lvec_type& dv0, const double& dt)
    {
        //Set parameters
        alpha = alpha0;
        dv = dv0;
        tof = dt;
    }
    
    /// Run corrections process
    CorrectionResult correct();
    
    
    /// Get solution (useful for 1-line set of variables)
    void getSolution(ProgenitorState<state_type,lvec_type>& ps) const
    {
        ps.found = isSolutionFound;
        ps.alpha = X(0);
        int dim = (int)dv.size();
        for(int i=0; i<dim; i++) {
            ps.dv[i] = X(i+1);
        }
        ps.tof = X.tail(1)[0];
        //Set the orbit state by calling Periodic Orbit object:
        theManifold->theOrbit.getStateFromAlpha(ps.alpha, ps.orbitState);
    }
    
    /// Printing options
    void setVerbose(const bool v, const bool outF, const bool outDF)
    {
        _verbose = v;
        outputF = outF;
        outputDF = outDF;
    }
    
private :
    /// Initialize:  Sampling the Initial State to fill out X
    void initialize()
    {
    
        const int n = (int)phiT.size() + 1;
        
        //Single shooting - forward or backward
        //Initialize F(X)
        F.resize(n);
        //Initialize F vector
        F.setZero();
        
        //Fill X values based on booleans
        X.resize(n+1);
        dX.resize(n+1);
        dX.setZero();
        //Fill Out initial values
        X(0) = alpha;
        for(int i=0; i<(int)dv.size(); i++) {
            X(i+1) = dv[i];
        }
        X(n) = tof;
        
    } //End initialize()
    
    
    /// Compute F and DF given X
    void computeF_and_DF(MatrixXd& DF);
    
    
    //Member variables
    bool _verbose,outputF,outputDF;
    bool isSolutionFound;
    double alpha; //Initial guess parameter on orbit
    lvec_type dv; //initial guess delta-V vector (dxdot,dydot)
    lvec_type phiT; //Manifold Point to hit from progenitor state
    double desired_H; //Hamiltonian value to maintain in this problem
    int _maxIters; //Max Iterations in corrections process
    double zeroTimeTol; //Tolerance on time approaching zero (failure condition)
    double _tol;  //Convergence Tolerance
    double _prec; //Integration precision
    MAP* theMap;   //The Poincare map (with RHS, ODESolver, and Section)
    const MANIFOLD* theManifold; //MapManifold object for this problem
    double tof; //Time of flight
    int jobID; //An index when multiple instances are running at once
    //Parameters for multi-dim problem
    VectorXd X,F,dX; //Free-variables, Constraints, Newton-Step
    VectorXd Xstore;  //Storage when needed
};


///The correction function for ProgenitorCorrector - main algorithm loop
template<class MAP, class MANIFOLD>
CorrectionResult ProgenitorCorrector<MAP,MANIFOLD>::correct()
{
    bool allDone = false;
    //Initialize some variables
    double normF = 1000.0;
    double normF0 = 1000.0;
    int count = 0;
    bool succeed = true;
    isSolutionFound = false;
    //Lets check if the guess is valid
    /*if (!isValid(_x0)) {
    allDone = true;
    return CorrectionResult(false,0,-1.0,-1.0);
    }*/
    
    
    //Initialize the run by computing X,F,DF for initial guess
    bool initOK = true;
    try {
        initialize();
    } catch(...) {
        //Note:  initialize() has calls to "theMap->integrate_state()",
        //       which can throw Singularity Intersection errors from the EOMs class.
        initOK = false;
    }
    MatrixXd DF = MatrixXd::Zero(F.size(),X.size());//(n*numPatchPts+1,n*numPatchPts+1) With all zeros
    if (_verbose) {
        //std::cerr << "Entering Periodic Orbit Corrections:\n";
        std::cerr << "(" << jobID << ")  Initial Guess a0 = " << alpha
                  << " dvMag = " <<  nvis::norm(dv) << " ToF= " << tof << "\n";
        if (!initOK) {
            std::cerr << "(" << jobID << ")  initialize() caught exception! Singularity crossing likely!\n";
        }
        //std::cerr << "(" << jobID << ")  DF.size() = ( " << n*numPatchPts+1 << " , " << n*numPatchPts+1 <<" )\n";
    }
    
    //Reset values
    normF = 1000.0; //Current error norm
    normF0 = 1000.0;//Error norm of first guess
    count = 0; //Current number of iterations
    succeed = true; //Success boolean
    
    if(initOK) {
    
        //Normal Loop for Iteratively solving system of nonlinear equations
        while (normF > _tol) {
            //Compute F and DF
            try {
                computeF_and_DF(DF);
                normF = F.norm();
            } catch(...) {
                //Error occurrs in mapping
                normF = 1000.0; //Force a quit with failure.
            }
            
            //Store normF0 if we are on the first pass
            if(count==0) {
                normF0=normF;
            }
            
            if (_verbose) {
                std::cerr << "(" << jobID << ") Iter " << count << " : ||F(X)|| = " << normF << "\n";
            }
            
            //Surpassed Maximum Iterations
            if (count >= _maxIters) {
                succeed = false;
                break;
            }
            //F(X) grew too big -> Likely to fail or jump to undesired behavior
            if (normF > 2.0) {
                succeed = false;
                break;
            }
            //Is time-of-flight approaching zero
            if (fabs(X.tail(1)[0]) <= zeroTimeTol) {
                succeed = false;
                break;
            }
            //Double-check |F(X)| isn't a nan
#ifdef _WIN32
            if (  boost::math::isnan(normF)  ) {
#else
            if (std::isnan(normF)) {
#endif
                succeed = false;
                break;
            }
            
            count++;
            //Solve for dX (Uses minNorm or just Newton-Raphson based on DF size)
            dX = DF.colPivHouseholderQr().solve(-F); //Non-sparse implementation
            
            
            //Update if we are not done
            if (normF > _tol) {
                X += dX;
            }
        } //End while loop
    } //End QN mode switch
    else {
        //Initialize guess failed:
        //  - Likely a really bad/sensitive guess that intersects a singularity
        //  - Probably not a useful guess anyway...
        succeed = false;
        normF = 1e4;
    } //End !initOK check
    
    //Temporary result information
    CorrectionResult tempResult(succeed,count,normF0,normF);
    
    
    //Store result (which satisfies ALL constriants)
    if (succeed) {
        if (_verbose) {
            std::cerr << " Progenitor State Correction Result: SUCCESS!\n";
            std::cerr << "----------------------------------------------------------------\n";
        }
        //STORE RESULT
        alpha = X(0);
        for(int i=0; i<(int)dv.size(); i++) {
            dv[i] = X(i+1);
        }
        tof =
            isSolutionFound = true;
    } else {
        //Bummer
        if (_verbose) {
            std::cerr << " Progenitor State Correction Result: FAILURE!\n";
            std::cerr << "----------------------------------------------------------------\n";
        }
        isSolutionFound = false;
    }
    //Return the overall result
    return tempResult;
}

/// Compute F and DF given X
template<class MAP, class MANIFOLD>
void ProgenitorCorrector<MAP,MANIFOLD>::computeF_and_DF(MatrixXd& DF)
{
    const int m = (int)phiT.size() + 1; //F(X) dimension
    static const int n = rhs_type::dimension; //state dimension
    static const int lDim = section_type::local_dimension; //Section dimension
    //outputF = _verbose;
    //outputDF = _verbose;
    
    double dt = X(m); //Last entry is time
    
    
    //Find state given alpha and delta-V vector
    const rhs_type& rhs = theMap->rhs();
    xstate_type y0(0.0);
    state_type icState(0.0);
    theManifold->theOrbit.getStateFromAlpha(X(0),icState); //State from PeriodicOrbit
    for (int i=0; i<n; i++) {
        y0[i] = icState[i];
    }
    xstate_type dxdt0 = rhs(0.0,y0); //State derivative at alpha
    double orbitPeriod = theManifold->baseFixedPoint.timePeriod; //NonDim time
    for (int i=0; i<n; i++) {
        y0[n+n*i+i] = 1.0;    //Identity for STM0
    }
    
    //Add the delta-v elements
    for(int i=0; i<(int)dv.size(); i++) {
        y0[n/2+i] += X(i+1);
    }
    
    /*if (_verbose) {
        std::cout << " : y0 = " << y0 << '\n';
    }*/
    //Run the state through dt time
    return_state finalInfo;
    try {
        finalInfo = theMap->integrate_state(y0,dt);
    } catch(...) {
        if(_verbose) {
            std::cerr << "Map error occurred in corrections process. Exiting.\n";
        }
        throw std::runtime_error("computeF_and_DF(): Error in Mapping...");
        return;
    }
    xstate_type yf = finalInfo.getState();
    /*if (_verbose) {
        std::cout << " : y0 = " << y0 << '\n';
        std::cout << " : yf = " << yf << '\n';
    }*/
    
    //State derivative
    xstate_type dxdt = rhs(0.0,yf);
    
    //Hitting the Map Point constraint
    std::pair<lvec_type, lmat_type> yfPair = theMap->section().project(yf);
    for(int i=0; i<(int)phiT.size(); i++) {
        F(i) = yfPair.first[i] - phiT[i];
    }
    //Hamiltonian Constant constraint
    F(m-1) = theMap->rhs().hamiltonian(0.0,y0) - desired_H;
    
    if(outputF) {
        for(int i=0; i<m; i++) {
            std::cerr << "     F(" << i << ") = " << F(i) << "\n";
        }
    }
    
    //DF(X) - alpha components:
    lgMatrix phiLocal = theMap->section().localProjection( finalInfo.getState() );
    for (int j = 0; j<lDim; j++) {
        double sum = 0.0;
        for (int i=0; i<n; i++) {
            sum += phiLocal[j][i]*dxdt0[i];
        }
        DF(j,0) = sum*orbitPeriod;
    }
    //DF(X) - delta-xdot & delta-ydot components:
    for (int j = 0; j<lDim; j++) {
        for (int i=0; i<(int)dv.size(); i++) {
            DF(j,1+i) = phiLocal[j][n/2+i];
        }
    }
    //DF(X) - time-of-flight components (state derivative):
    std::pair<lvec_type, lmat_type> dxdtProj = theMap->section().project(dxdt);
    for(int j=0; j<lDim; j++) {
        DF(j,1+lDim) = dxdtProj.first[j];
    }
    
    //Set the Hamiltonian constraint elements
    state_type hPartials = theMap->rhs().hamiltonian_first_partials(0.0,y0);
    for(int i=0; i<(int)dv.size(); i++) {
        DF(m-1,1+i) = hPartials[n/2+i];
    }
} //End computeF_and_DF();



} //namespace orbital

#endif //corrections.hpp
