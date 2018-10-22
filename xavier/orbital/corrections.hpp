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
//  Corrections header file
//  Author:  Wayne Schlei
//  Date: 12/7/2012
//
//  Purpose:  For differential corrections procedures using
//  multiple shooting.  (Fixed-point computation)
//
//  Modification:
//  1/13/2014-WRS:  Attempting to attach a control mechanism on
//  the corrections process that autonomously adjust constraints
//  and resolves to attain a solution.
/////////////////////////////////////////////////////////////////
#ifndef CORRECTIONS_HPP
#define CORRECTIONS_HPP

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
#include <maps/poincare_map.hpp>
#include <maps/section.hpp>
#include <maps/fixpoints.hpp>
#include <orbital/monodromy.hpp>
#include <orbital/controller.hpp>
#include <orbital/qnSingleShoot.hpp>
#include <orbital/qnMultiShoot.hpp>
#include <orbital/quasiNewtonRootFind.hpp>

#ifdef _WIN32
#include <boost/math/special_functions/fpclassify.hpp> // isnan for Windows
#endif

using namespace Eigen;
using namespace xavier;

namespace orbital {


/** Corrections Class for a Periodic Orbit
 *  Process is a variable-time multiple-shooting method that includes
 *  Hamiltonian, Dist-to-Section, and Periodicity constraints.
 *  This class is monitored and modified with a CorrectionsRegulator
 *  class that will turn certain constraints ON/OFF or enable a Quasi-Newton
 *  method adaptively to assist convergence.
 */
template<typename MAP>
class PeriodicOrbitCorrector {
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
    typedef typename MAP::return_type               return_type;
    typedef typename MAP::return_state              return_state;
    typedef typename return_state::mat_type         return_mat;
    typedef typename MAP::solver_type               solver_type;
    typedef typename solver_type::step              step_type;
    typedef CorrectionsRegulator                    reg_type;
    //typedef Eigen::Triplet<double>                  Triple;
    //typedef Eigen::SparseMatrix<double>             SpMat;
    
    /// Constructors
    PeriodicOrbitCorrector(MAP& pmap, const lvec_type& x0, const value_type& Hval,
                           const int& numCross, const int ptsPerCross=8, bool verbose=false) :
        _verbose(verbose), outputF(false), outputDF(false),
        _x0(x0), _p(numCross), _ptsPerp(ptsPerCross),
        _maxIters(20), zeroTimeTol(1.e-5), _tol(1.e-8), _prec(1.e-12),
        theMap(&pmap), _desired_H(Hval), _period(0.0), numPatchPts(ptsPerCross*numCross), jobID(0),
        X(VectorXd::Zero(1)), F(VectorXd::Zero(1)), dX(VectorXd::Zero(1)),
        Xstore(VectorXd::Zero(1)), dGydotdX(MatrixXd::Zero(1,1)), dGdlambda(VectorXd::Zero(1)),
        _Hlambda(Hval), gState(0.0), finalState(0.0), _homotopyMonitor()
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
        _period = T;
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
    void setGuess(const lvec_type& x0, const value_type Hdesired, int numCross)
    {
        //Set parameters
        _x0 = x0;
        _desired_H = Hdesired;
        _Hlambda = Hdesired;
        theMap->setHamiltonian( Hdesired );
        _p = numCross;
        _period = 0.0; //Reset
        numPatchPts = _ptsPerp*_p;
        
        //Run the map to find the period(time) of the orbit
        //->Map x0 for p iterates
        return_type r;
        try {
            r = theMap->map_complete(_x0, _p);
            _period += r.t;
        } catch (xavier::MapUndefined& m_u ) {
            std::cerr << "PeriodicOrbitCorrector:: Caught Map Undefined while propagating initial state\n";
            std::cerr << "  x0 = " << _x0 << " for p = " << _p << " iters;  (t = " << _period << ")\n";
            throw m_u;
        } catch (...) {
            std::cerr << "PeriodicOrbitCorrector:: Caught Error while propagating initial state\n";
            std::cerr << "  x0 = " << _x0 << " for p = " << _p << " iters;  (t = " << _period << ")\n";
            xavier::MapUndefined m_u;
            m_u.where = _x0;
            throw m_u;
        }
    }
    
    /// Unique orbit index for a fixed point
    static void generateOrbitID()
    {
#if _OPENMP
        #pragma omp atomic
        _orbitIndex++; //Static member
#else
        _orbitIndex++; //Static member
#endif
    }
    static void resetOrbitID()
    {
#if _OPENMP
        #pragma omp atomic
        _orbitIndex -= _orbitIndex; //Static member
#else
        _orbitIndex = 0; //Static member
#endif
    }
    static int getOrbitID()
    {
        return _orbitIndex;
    }
    
    
    /// Run corrections process
    CorrectionResult correct();
    
    /// Set fixed-point information
    //  - call parts of fixedpoints_wSTM.hpp to do the linear analysis
    bool setFixedPoints(fixpoint& fp, std::vector<fixpoint>& iterates, bool linearSTM=false)
    {
        //Assume the fixed point is found (fp is initial fixed point)
        generateOrbitID();
        fp.pos = _x0;
        fp.K = _p;
        fp.orbitID = _orbitIndex;
        
        //Call monodromy analysis to find fixed-point chain and eigen vals/vecs
        bool ok = orbital::monodromy_analysis<MAP>(*theMap, fp, iterates, linearSTM, _verbose);
        
        //Set the xavier::fixpoint data from iterates since it's reset
        if (ok) {
            fp = iterates[0];
        }
        
        //Return if orbit passed sanity check (and therefore if data is valid)
        return ok;
    }
    
    /// Get the correction mode from the homotopy regulator
    reg_type::OperatingMode getMode()
    {
        return _homotopyMonitor.mode;
    }
    /// Get the correction mode as a string
    std::string getModeName()
    {
        return _homotopyMonitor.getModeName();
    }
    
    /// Set the stopping mode (last method to try) for the Corrections Regulator object
    void setStopMode(reg_type::OperatingMode mode)
    {
        _homotopyMonitor.setStoppingMode(mode);
    }
    /// Set the operating mode
    void setOpMode(reg_type::OperatingMode mode)
    {
        _homotopyMonitor.setMode(mode);
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
    
        static const int n = rhs_type::dimension;
        reg_type* monitor = &_homotopyMonitor;
        
        //Switch based on OperatingMode
        /** Try to use outputs from regulator to determine what to do:
            - Skip initializing if still in loop
            - Call stored value
            - Various samples techniques
        */
        
        //Multiple Shooting by default
        bool singleShooting = false; //Indicating single shooting
        bool mixedSample = false; //Split sampling
        bool forwardSampling = true; //forward or mixed sampling
        bool useStorage = false; //Use stored distributed error or don't init
        
        
        switch (monitor->mode) {
            case reg_type::MULTIPLE_SHOOT :
                numPatchPts = _p * _ptsPerp;
                //Initialize F(X)
                F.resize(n*numPatchPts+1);
                break;
            case reg_type::RESAMPLE :
                //ptsPerCross is incremented by 1
                _ptsPerp++;
                numPatchPts = _p * _ptsPerp;
                //Use FORWARD sampling - initialized values
                //Initialize F(X)
                F.resize(n*numPatchPts+1);
                //forwardSampling = true;
                break;
            case reg_type::DISTRIBUTED_ERROR :
                //Use given value of numPatchPts and modified sampling
                forwardSampling = false;
                mixedSample = true; //Split backward&forward sampling
                numPatchPts = _p * _ptsPerp;
                //Initialize F(X)
                F.resize(n*numPatchPts+1);
                break;
            case reg_type::HAMSEC_OFF :
                numPatchPts = _p * _ptsPerp;
                if (!monitor->processingSteps) {
                    //Initialize F(X)
                    F.resize(n*numPatchPts-1);  //No Hamiltonian or section constraints
                    forwardSampling = false;
                    useStorage = true;
                } else {
                    //Don't initialize X, just use the current value
                    useStorage = false;
                    forwardSampling = false;
                    F.resize(n*numPatchPts+1); //Adding constraints back in
                }
                break;
            case reg_type::PERIODICITY_HOMOTOPY :
                numPatchPts = _p * _ptsPerp;
                if (!monitor->processingSteps) {
                    //Use forward sampling when initializing
                    forwardSampling = true;
                    useStorage = false;
                    //Initialize F(X)
                    F.resize(n*numPatchPts); //n*k (missing Hamiltonian)
                } else {
                    //Processing steps
                    forwardSampling = false;
                    useStorage = false;
                    if (monitor->lambda >= 1.0) {
                        //Add Hamiltonian constraint back for continuation
                        F.resize(n*numPatchPts+1);
                    }
                }
                break;
            case reg_type::QUASI_NEWTON :
                //Try on single Shooting
                singleShooting = true;
                numPatchPts = 1;
                useStorage = false;
                F.resize(n+1);
                //MultipleShooter
                //singleShooting = false;
                //useStorage = true;
                //F.resize(n*numPatchPts+1);
                forwardSampling = false;
                mixedSample = false;
                break;
            default :
                //Single shooting - forward or backward
                singleShooting = true;
                numPatchPts = 1;
                //Initialize F(X)
                F.resize(n+1);
        }
        //Initialize F vector
        F.setZero();
        
        
        
        //Fill X values based on booleans
        if (singleShooting) { //-------------------------------------------------
            X.resize(n+1);
            dX.resize(n+1);
            dX.setZero();
            //Fill Out initial values
            xstate_type y0 = theMap->section().unproject(_x0);
            for(int j=0; j<n; j++) {
                X(j) = y0[j];
            }
            X(n) = _period;
            if (monitor->mode == reg_type::RSINGLE_SHOOT) {
                X(n)=-_period;
            }
        } else if (forwardSampling) { //------------------------------------------
            //Initialize size
            X.resize(n*numPatchPts+1);
            dX.resize(n*numPatchPts+1);
            
            //Run integration to fill out X
            xstate_type y0 = theMap->section().unproject(_x0);
            return_state stateInfo(y0,0.0);
            //Fill initial value
            for (int j=0; j<n; j++) {
                X(j) = y0[j];
            }
            
            //Sample states
            //Option 1) Integrate FORWARD to get the remaining states
            //(endpoints of k-1 integrations for period/k time)
            for (int i=0; i<(numPatchPts-1); i++) {
                xstate_type y = stateInfo.getState();
                stateInfo = theMap->integrate_state(y,_period/(double)numPatchPts);
                for(int j=0; j<n; j++) {
                    X(n*(i+1)+j) = stateInfo.x[j];
                }
                //Debugging
                /*std::ostringstream os("");
                os << "\rInitialize(): MS Propagating point " << i << "/" << (numPatchPts-1)
                   << " for dt = " << _period/(double)numPatchPts << " \r" << std::flush;
                std::cout << os.str();*/
            }
            //Last entry is the period guess
            X(n*numPatchPts) = _period;
            
            //Store the end point of integration for use in homotopy method
            xstate_type y = stateInfo.getState();
            stateInfo = theMap->integrate_state(y,_period/(double)numPatchPts);
            for(int j=0; j<n; j++) {
                gState[j] = stateInfo.x[j];
            }
            
        } else if (mixedSample) { //---------------------------------------------
            //Sample with error focused away from the periodicity constraints
            //Initialize sizes
            X.resize(n*numPatchPts+1);
            dX.resize(n*numPatchPts+1);
            dX.setZero();
            
            //Run integration to fill out X
            //std::cerr << " Mixed Sampling:\n";
            xstate_type y0 = theMap->section().unproject(_x0);
            return_state stateInfo(y0,0.0);
            //Fill initial value
            for (int j=0; j<n; j++) {
                X(j) = y0[j];
            }
            
            //Option 2) Integrate FORWARD AND BACKWARD to get the states
            //starting from the initial state (reduces initial error in
            //periodicity constraints)
            
            //Integrate and fill initial guess based on k being odd or even
            if (numPatchPts % 2 != 0) { //Odd k
                //Integrate forward & backward halves
                //std::cerr << " numPatchPts is ODD so integrating halves...\n";
                return_state backStateInfo(y0,0.0);
                for(int i=0; i<(int)floor((double)numPatchPts/2); i++) {
                    //Forward
                    xstate_type y = stateInfo.getState();
                    stateInfo = theMap->integrate_state(y,_period/(double)numPatchPts);
                    for(int j=0; j<n; j++) {
                        X(n*(i+1)+j) = stateInfo.x[j];
                    }
                    //Debugging:
                    //std::cerr << "   Forward state " << i << "\n";
                    //std::cerr << "    y = " << stateInfo.x << " at Patch Point" << i+1 << "\n";
                    /* std::ostringstream os("");
                     os << "\rInitialize(): MixOddFOR Propagated point " << i << "/" << (int)floor((double)numPatchPts/2)
                        << " for dt = " << _period/(double)numPatchPts << " \r" << std::flush;
                     std::cout << os.str();*/
                    
                    //Backward Time
                    xstate_type yb = backStateInfo.getState();
                    backStateInfo = theMap->integrate_state(yb, -_period/(double)numPatchPts);
                    //Reverse Fill
                    for(int j=0; j<n; j++) {
                        X(n*(numPatchPts-1-i) + j) = backStateInfo.x[j];
                    }
                    //Debugging:
                    //std::cerr << "   Backward state " << i << "\n";
                    //std::cerr << "    y = " << backStateInfo.x << " at Patch Point" << (numPatchPts-1)-i << "\n";
                    /*std::ostringstream os2("");
                     os2 << "\rInitialize(): MixOddBACK Propagating point " << i << "/" << (int)floor((double)numPatchPts/2)
                        << " for dt = " << -_period/(double)numPatchPts << " \r" << std::flush;
                     std::cout << os2.str();*/
                }
            } else {    //Even k
                //Integrate Forward and backward to 1 pt prior end point
                //std::cerr << " numPatchPts is EVEN, so integrating halves -1 pt prior...\n";
                return_state backStateInfo(y0,0.0);
                for(int i=0; i<(numPatchPts/2 - 1); i++) {
                    //Forward
                    xstate_type y = stateInfo.getState();
                    stateInfo = theMap->integrate_state(y,_period/(double)numPatchPts);
                    for(int j=0; j<n; j++) {
                        X(n*(i+1)+j) = stateInfo.x[j];
                    }
                    //Debugging:
                    //std::cerr << "   Forward state " << i << "\n";
                    //std::cerr << "    y = " << stateInfo.x << " at Patch Point" << i+1 << "\n";
                    /*std::ostringstream os("");
                     os << "\rInitialize(): MixEvenFOR Propagating point " << i << "/" << (numPatchPts/2 - 1)
                        << " for dt = " << _period/(double)numPatchPts << " \r" << std::flush;
                     std::cout << os.str();*/
                    
                    //Backward Time
                    xstate_type yb = backStateInfo.getState();
                    backStateInfo = theMap->integrate_state(yb, -_period/(double)numPatchPts);
                    //Reverse Fill
                    for(int j=0; j<n; j++) {
                        X(n*(numPatchPts-1-i) + j) = backStateInfo.x[j];
                    }
                    //Debugging:
                    //std::cerr << "   Backward state " << i << "\n";
                    //std::cerr << "    y = " << backStateInfo.x << " at Patch Point" << (numPatchPts-1)-i << "\n";
                    /*std::ostringstream os2("");
                     os2 << "\rInitialize(): MixEvenBACK Propagating point " << i << "/" << (numPatchPts/2 - 1)
                        << " for dt = " << -_period/(double)numPatchPts << " \r" << std::flush;
                     std::cout << os2.str();*/
                }
                //Integrate last states and average values
                xstate_type y = stateInfo.getState();
                stateInfo = theMap->integrate_state(y,_period/(double)numPatchPts);
                xstate_type yb = backStateInfo.getState();
                backStateInfo = theMap->integrate_state(yb, -_period/(double)numPatchPts);
                //Fill with average
                int iMiddle = numPatchPts/2;
                for(int j=0; j<n; j++) {
                    X(n*(iMiddle)+j) = (backStateInfo.x[j] + stateInfo.x[j])/2.0;
                }
                //std::cerr << "   Averaged Midpoint " << iMiddle << "\n";
                //std::cerr << "   y = " << 0.5*(backStateInfo.x + stateInfo.x) << " at Patch Point " << iMiddle << "\n";
            }
            //Last entry is the period guess
            X(n*numPatchPts) = _period;
            //Force to storage : Assumes we run DISTRIBUTED_ERROR before HAMSEC_OFF
            Xstore.resize(n*numPatchPts+1); //Note, this won't apply if we have map errors
            Xstore = X;
            
        } else if (useStorage) { //------------------------------------------------
            //Use the stored value (from DISTRIBUTED_ERROR step) to fill
            //the design variable vector
            X.resize(n*numPatchPts+1);
            dX.resize(n*numPatchPts+1);
            dX.setZero();
            X.setZero();
            if(Xstore.size() > 1) {
                X = Xstore;    //Otherwise, this is overwritten.
            }
        }
        //else -> Use the current X as we are continuing to another step
    }//End initialize()
    
    /// Set the constant coefficients
    void setConstantDFCoeffs(MatrixXd& DF)
    {
        static const int n = rhs_type::dimension;
        reg_type* monitor = &_homotopyMonitor;
        if (outputDF) {
            std::cerr << " Fixed coefficients of DF = \n";
        }
        
        //-I's for continuity
        for(int kk=0; kk<(numPatchPts-1); kk++) {
            for (int jj=0; jj<n; jj++) {
                DF(jj+kk*n,jj+n*(kk+1)) = -1.0;
                if (outputDF) {
                    std::cerr << "  DF(" << jj+kk* n << " , " << jj+n*(kk+1)
                              << ") = " << DF(jj+kk*n,jj+n*(kk+1)) << "\n";
                }
            }
        }
        
        /** This assumes the (n-2)th state component is a skipped element.
          * Note: the equations are assumed to have a constant Hamiltonian,
          * i.e., at least ONE integral of motion.
          * In the future, this may have to be generalized or specified by
          * the user.  n-2 may not always be a good choice.  (Usually use
          * a perpendicular velocity component?)
        */
        //-1's for periodicity
        for(int kk=0; kk<(n-2); kk++) {
            if (monitor->mode == reg_type::PERIODICITY_HOMOTOPY) {
                DF(n*(numPatchPts-1)+kk,kk) = -(monitor->lambda);
            } else {
                DF(n*(numPatchPts-1)+kk,kk) = -1.0;
            }
            
            if (outputDF) {
                std::cerr << "  DF(" << n*(numPatchPts-1)+kk << " , " << kk
                          << ") = " << DF(n*(numPatchPts-1)+kk,kk) << "\n";
            }
            
        }
        dGydotdX(n-2) = -(monitor->lambda);//Homotopy
        if (monitor->mode == reg_type::PERIODICITY_HOMOTOPY) {
            DF(n*(numPatchPts-1)+(n-2),n-1) = -(monitor->lambda);
        } else {
            DF(n*(numPatchPts-1)+(n-2),n-1) = -1.0;    //last (n-1)th state equation
        }
        //Verbose output
        if (outputDF) {
            std::cerr << "  DF(" << n*(numPatchPts-1)+n-2 << " , " << n-1
                      << ") = " << DF(n*(numPatchPts-1)+n-2,n-1) << "\n";
        }
        
        //The rest of DF(X) is set with computeF_and_DF();
    }
    
    // A sparse representation for use with SparseMatrix in Eigen
    /*void setConstantDFCoeffs_Sparse(std::vector<Triple>& tuples) {
        static const int n = rhs_type::dimension;
    reg_type *monitor = &_homotopyMonitor;
    
        //Set the size of coefficients
        tuples.clear();
        tuples.reserve(n*(numPatchPts-1) + n-1 + n*n*(numPatchPts-1) + n*(n-1) + n*(numPatchPts-1) + n-1 + n + n);
    
       //-I's for continuity
        for(int kk=0; kk<(numPatchPts-1); kk++) {
        for (int jj=0;jj<n;jj++) {
                tuples.push_back(Triple(jj+kk*n,jj+n*(kk+1),-1.0) );
        }
        }
        //-1's for periodicity (n-1)
        for(int kk=0; kk<(n-2); kk++) {
            tuples.push_back( Triple(n*(numPatchPts-1)+kk,kk,-1.0) );
        }
        tuples.push_back( Triple(n*(numPatchPts-1)+(n-2),n-1,-1.0) );
        if (_verbose) {
            std::cerr << " Fixed coefficients of DF = \n";
            int idx1 = n*(numPatchPts-1)+(n-1);
            for (int i=0; i<idx1; i++) std::cerr << "   (" << i
                      << ") = (" << tuples[i].row() << ", "
                      << tuples[i].col() << ") "
                      << " = " << tuples[i].value() << "\n";
        }
        //The rest of the triples are set with computeF_and_DF_Sparse();
    }*/
    
    /// Compute F and DF given X
    void computeF_and_DF(MatrixXd& DF);
    //void computeF_and_DF_Sparse(std::vector<Triple>& tuples);
    
    /// Compute the change in X given an update (lambda) in periodicity homotopy method
    bool computePeriodicityHomotopyPrediction(MatrixXd& DF);
    
    /// Compute the change in X given an update (gamma) in H-homotopy method
    bool computeHamiltonianHomotopyPrediction(MatrixXd& DF);
    
    
    
    //Member variables
    static int _orbitIndex;
    bool _verbose,outputF,outputDF;
    lvec_type _x0; //initial guess
    int _p; //Number of crossings
    int _ptsPerp; //Patch points per crossing
    int _maxIters;
    double zeroTimeTol; //Tolerance on time going to zero
    double _tol;  //Convergence Tolerance
    double _prec; //Integration precision
    MAP* theMap;   //The Poincare map (with RHS, ODESolver, and Section)
    double _desired_H; //The desired Energy value
    double _period; //Time
    int numPatchPts; //Actual number of patch pts being used
    int jobID; //An index when multiple instances are running at once
    //Parameters for multi-dim problem
    VectorXd X,F,dX; //Free-variables, Constraints, Newton-Step
    VectorXd Xstore;  //Storage when needed
    MatrixXd dGydotdX;
    VectorXd dGdlambda; //Homotopic method vectors
    double _Hlambda; //Stored value of Hamiltonian after lambda = 1
    state_type gState, finalState; //Storage of end states for homotopic method
    reg_type _homotopyMonitor; //Corrections Regulator for using homotopic method
};

/// Initialize orbit index
template<typename MAP>
int PeriodicOrbitCorrector<MAP>::_orbitIndex = 0;

///The correction function for PeriodicOrbitCorrector - main algorithm loop
template<typename MAP>
CorrectionResult PeriodicOrbitCorrector<MAP>::correct()
{
    static const int n = rhs_type::dimension;
    reg_type* monitor = &_homotopyMonitor;
    monitor->setJobID(jobID);
    monitor->setVerbose(_verbose);
    bool allDone = false;
    //Initialize some variables
    double normF = 1000.0;
    double normF0 = 1000.0;
    int count = 0;
    bool succeed = true;
    //Lets check if the guess is valid
    /*if (!isValid(_x0)) {
    allDone = true;
    return CorrectionResult(false,0,-1.0,-1.0);
    }*/
    
    
    //Recipe and Homotopic Method CONTROL Loop
    while (!allDone && !monitor->modesExplored) {
    
        //Prompt the user what is going on
        if (_verbose) {
            monitor->verbose();
        }
        
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
        //SparseMatrix<double> DF; //Initialize, column-major order
        //DF.resize(F.size(),X.size());
        dGydotdX.resize(1,X.size()); //Row vec
        dGdlambda.resize(n*numPatchPts+1); //One more than F_HomotopyLambda(X)
        dGydotdX.setZero();
        dGdlambda.setZero();
        
        if (_verbose) {
            //std::cerr << "Entering Periodic Orbit Corrections:\n";
            std::cerr << "(" << jobID << ")  Initial Guess x0 = " << _x0 << "\n";
            gvec_type tempState(0.0);
            if (initOK) {
                for(int i=0; i<n; i++) {
                    tempState[i] = X[i];
                }
                std::pair<lvec_type, lmat_type> x0Pair = theMap->section().project(tempState);
                lvec_type x0Now = x0Pair.first;
                std::cerr << "(" << jobID << ")  Now: x0 ~= " << x0Now << "\n";
            } else {
                std::cerr << "(" << jobID << ")  initialize() caught exception! Singularity crossing likely!\n";
            }
            std::cerr << "(" << jobID << ")  Period = " << _p << "\n";
            //std::cerr << "(" << jobID << ")  DF.size() = ( " << n*numPatchPts+1 << " , " << n*numPatchPts+1 <<" )\n";
            std::cerr << "(" << jobID << ")  Num Patch Points = " << numPatchPts << "\n";
            //Monitor status
            std::cerr << "(" << jobID << ")  Monitor->hsConstraintsOn = "
                      << std::string((monitor->hsConstraintsOn)? "true" : "false") << "\n";
            std::cerr << "(" << jobID << ")  Monitor->processingSteps = "
                      << std::string((monitor->processingSteps)? "true" : "false") << "\n";
        }
        
        //Reset values
        normF = 1000.0; //Current error norm
        normF0 = 1000.0;//Error norm of first guess
        count = 0; //Current number of iterations
        succeed = true; //Success boolean
        
        //Split behavior if Quasi-Newton is enabled
        if (initOK && monitor->lineSearchBacktrack) {
            /*//Build RHS (Multiple Shooting)
            MultipleShooter<MAP> fpMSFunc(theMap,numPatchPts,_desired_H);
            //Build QuasiNewton root-finder object (from Num.Rec.)
            QuasiNewton< MultipleShooter<MAP> > qnSolver(fpMSFunc,_verbose);*/
            
            //Build RHS (Single Shooting)
            SingleShooter<MAP> fpSSFunc(*theMap,_desired_H);
            //Build QuasiNewton root-finder object
            QuasiNewton< SingleShooter<MAP> > qnSolver(fpSSFunc,_verbose);
            
            qnSolver.setJobID(jobID);
            qnSolver.setMaxIts(2*_maxIters);
            qnSolver.setTol(_tol);
            bool check = false;
            std::vector<double> xVec;
            for (int i=0; i<(n*numPatchPts+1); i++) {
                xVec.push_back( X(i) );
            }
            try {
                qnSolver.solve(xVec,normF,normF0,check);
            } catch(...) {
                if(_verbose) {
                    std::cerr << "(" << jobID << ") QN solver error.\n";
                }
                succeed = false;
            }
            //Output result
            normF = sqrt(2.0*normF);  //Actually f=0.5*F(X)*F(X) in QuasiNewton
            if (normF > _tol) {
                succeed = false;
            }
            for (int i=0; i<(n*numPatchPts+1); i++) {
                X(i) = xVec[i];
            }
            count = qnSolver.getNumIters();
            
        } else if(initOK) {
        
            //Normal Loop for Iteratively solving system of nonlinear equations
            while (normF > _tol) {
                //Compute F and DF
                bool integrationFailure = false;
                try {
                    computeF_and_DF(DF);
                } catch(...) {
                    integrationFailure = true;
                }
                //Sparse Implementation
                //std::vector<Triple> elements;
                //computeF_and_DF_Sparse(elements); //Will change F and DF given current X
                //DF.setFromTriplets(elements.begin(),elements.end());
                //if (_verbose) {
                //    std::cout << "  Sparse DF(X) Matrix: \n";
                //    std::cout << "  DF = " << DF << '\n';
                //}
                
                
                //Can't actually evaluate F(X) and DF(X) [usually singularity in EOMs]
                if (integrationFailure) {
                    if(_verbose) {
                        std::cerr << "(" << jobID << ") Iter " << count << " : F/DF evaluation failure! (period = " << _p << ")\n";
                    }
                    succeed = false;
                    break;
                }
                
                //Set ||F(X)||
                normF = F.norm();
                
                //Store normF0 if we are on the first pass
                if(count==0) {
                    normF0=normF;
                }
                
                if (_verbose) {
                    std::cerr << "(" << jobID << ") Iter " << count << " : ||F(X)|| = " << normF << " (period = " << _p << ")\n";
                }
                
                //Surpassed Maximum Iterations
                if (count >= _maxIters) {
                    succeed = false;
                    break;
                }
                //F(X) grew too big -> Likely to fail or jump to undesired behavior
                if (normF > 1.5) {
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
                //Time is going to zero -> A sad cheat by the corrector
                double dt = X(n*numPatchPts)/(double) numPatchPts; //Last entry is time
                if (std::fabs(dt) < zeroTimeTol) {
                    if(_verbose) {
                        std::cerr << "(" << jobID << ") Iter " << count << " : Propagation time went to zero.\n";
                    }
                    succeed = false;
                    break;
                }
                
                count++;
                //Solve for dX
                dX = DF.colPivHouseholderQr().solve(-F); //Non-sparse implementation
                
                //Sparse Matrix Solving
                //BiCGSTAB<SpMat> solver(DF);
                //UmfPackLU<SpMat> solver(DF);
                //dX = solver.solve(-F);
                //if (solver.info() != Success ) {
                //    std::cerr << " WARNING:  Solution to DF * dX = -F has failed!\n";
                //    std::cerr << "           Solver returned : " << solver.info() << "\n";
                //}
                
                
                //Update if we are not done
                if (normF > _tol) {
                    X += dX;
                }
            }
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
        
        //Update Regulator -----------------------------------------------------------
        allDone = monitor->update(tempResult);
        if (_verbose) {
            std::cerr << "(" << jobID << ")  AllDone = " << std::string((allDone) ? "true": "false") << "\n";
        }
        //After trying RESAMPLE, we need to store the active number of patch points
        if (monitor->mode == reg_type::DISTRIBUTED_ERROR) {
            //Go back to original choice and use for sampling
            if (!(monitor->useExtraPointPerRev)) {
                _ptsPerp--;
                numPatchPts = _p*_ptsPerp;
            }
            //Keep numPatchPts the same otherwise
        } else if (monitor->mode == reg_type::PERIODICITY_HOMOTOPY &&
                   monitor->processingSteps) { //Homotopy except first step
            double& lambda = monitor->lambda;
            double& gamma = monitor->gamma;
            if (_verbose) {
                std::cerr << "(" << jobID << ") Periodicity Homotopy Params: lambda = "
                          << lambda << " gamma = " << gamma << "\n";
            }
            //Compute a prediction of X based on the next step and apply
            if (lambda < 1) {
                //Periodicity Homotopy Method
                bool predictOK = true;
                try {
                    predictOK = computePeriodicityHomotopyPrediction(DF);
                } catch(...) {
                    predictOK = false;
                }
                //Break the loop if no nullspace or inability to solve
                if (!predictOK) {
                    allDone = true;
                    succeed = false;
                }
                if (_verbose) {
                    std::cerr << "(" << jobID << ") Periodicity Homotopy Params: deltaLambda = "
                              << monitor->dlam << ", Next lambda = " << monitor->lambda << "\n";
                }
            } else if (lambda>=1.0 && gamma < 1.0) {
                bool predictOK = true;
                if (gamma == 0.0) {
                    //Periodicity homotopy is finished so store the resulting Hval
                    xstate_type startPoint(0.0);
                    for (int i=0; i<n; i++) {
                        startPoint[i] = X(i);
                    }
                    try {
                        _Hlambda = theMap->rhs().hamiltonian(0.0,startPoint);//H at lambda = 1
                    } catch(...) {
                        predictOK = false;
                    }
                }
                //Continuation in Hamiltonian constraint
                try {
                    predictOK = computeHamiltonianHomotopyPrediction(DF);
                } catch(...) {
                    predictOK = false;
                }
                if (!predictOK) {
                    allDone = true;
                    succeed = false;
                    std::cerr << " Homotopy Update FAILED :  Check program...\n";
                }
                if (_verbose) {
                    std::cerr << "(" << jobID << ") Periodicity Homotopy Params: deltaGamma = "
                              << monitor->dgam << ", Next gamma = " << monitor->gamma << "\n";
                }
            }
        }//---------------------------------------------------------------------------
        
    } //End CONTROL Loop
    
    //Store result (which satisfies ALL constriants)
    if (monitor->solutionFound()) {
        if (_verbose) {
            //std::cerr << " Approach (" << monitor->getModeName() << ") SUCCESS!\n";
            std::cerr << "----------------------------------------------------------------\n";
        }
        _period = X(n*numPatchPts);
        if(monitor->mode == reg_type::RSINGLE_SHOOT) {
            _period = -_period;
        }
        gvec_type x0(0);
        for (int i=0; i<n; i++) {
            x0[i] = X(i);
        }
        std::pair<lvec_type, lmat_type> x0pair = theMap->section().project(x0);
        _x0 = x0pair.first;
    } else {
        //Bummer
        if (_verbose) {
            std::cerr << " Reached Last Option - Stopping...\n";
            std::cerr << "----------------------------------------------------------------\n";
        }
    }
    //Return the overall result
    return CorrectionResult(monitor->solutionFound(),count,normF0,normF);
}

/// Compute F and DF given X
template<typename MAP>
void PeriodicOrbitCorrector<MAP>::computeF_and_DF(MatrixXd& DF)
{
    static const int n = rhs_type::dimension;
    reg_type* monitor = &_homotopyMonitor;
    
    //outputF = _verbose;
    //outputDF = _verbose;
    
    //Integrate Each state
    double dt = X(n*numPatchPts)/(double) numPatchPts; //Last entry is time
    
    //Call the function to set the constants
    setConstantDFCoeffs(DF);
    int setCount = n*(numPatchPts-1) + n-1;
    
    //For each patch point -> We could multi-thread this if we aren't solving fp's per thread
    for(int kk=0; kk<numPatchPts; kk++) {
        xstate_type y0(0.0);
        for (int i=0; i<n; i++) {
            y0[i] = X(n*kk+i);
        }
        for (int i=0; i<n; i++) {
            y0[n+n*i+i] = 1.0;    //Identity for STM0
        }
        /*if (_verbose) {
            std::cout << "      Point " << kk;
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
            return;
        }
        xstate_type yf = finalInfo.getState();
        /*if (_verbose) {
            std::cout << "      After Integ - Point " << kk;
            std::cout << " : y0 = " << y0 << '\n';
        
            std::cout << "      After Integ - Point " << kk;
            std::cout << " : yf = " << yf << '\n';
        }*/
        
        //STM = finalInfo.J
        //State derivative
        const rhs_type& rhs = theMap->rhs();
        xstate_type dxdt = rhs(0.0,yf);
        
        //Partials/constraints relevant to first point
        if (kk == 0) {
            if (monitor->mode == reg_type::PERIODICITY_HOMOTOPY) {//-------------------
                //Homotopy forms
                //Section Distance constraint
                F(n*numPatchPts - 1) = theMap->section().distance(y0);
                if (outputF) {
                    std::cerr << "      Distance Constraint: F(" << n* numPatchPts-1
                              << ") = " << F(n*numPatchPts-1) << '\n';
                }
                xstate_type sectPartials = theMap->section().distance_first_partials(y0);
                for (int i=0; i<n; i++) {
                    DF(n*numPatchPts - 1,i) = sectPartials[i];
                    setCount++;
                }
                if (monitor->lambda >=1.0) {
                    //Hamiltonian Constraint
                    F(n*numPatchPts) = theMap->rhs().hamiltonian(0.0,y0)
                                       -((1.0-monitor->gamma)*_Hlambda + (monitor->gamma)*_desired_H);
                    if (outputF) {
                        std::cerr << "      Current JC = " << theMap->rhs().hamiltonian(0.0,y0)
                                  << "  F(" << n* numPatchPts << ") = " << F(n*numPatchPts) << '\n';
                    }
                    state_type hPartials = theMap->rhs().hamiltonian_first_partials(0.0,y0);
                    for (int i=0; i<n; i++) {
                        DF(n*numPatchPts,i) = hPartials[i];
                        setCount++;
                    }
                }
            } else if (monitor->hsConstraintsOn) { //----------------------------------
                //Section Distance constraint
                F(n*numPatchPts - 1) = theMap->section().distance(y0);
                if (outputF) {
                    std::cerr << "      Distance Constraint: F(" << n* numPatchPts-1
                              << ") = " << F(n*numPatchPts-1) << '\n';
                }
                xstate_type sectPartials = theMap->section().distance_first_partials(y0);
                for (int i=0; i<n; i++) {
                    DF(n*numPatchPts - 1,i) = sectPartials[i];
                    setCount++;
                }
                //Hamiltonian Constraint
                F(n*numPatchPts) = theMap->rhs().hamiltonian(0.0,y0) - _desired_H;
                if (outputF) {
                    std::cerr << "      Current JC = " << theMap->rhs().hamiltonian(0.0,y0)
                              << "  F(" << n* numPatchPts << ") = " << F(n*numPatchPts) << '\n';
                }
                state_type hPartials = theMap->rhs().hamiltonian_first_partials(0.0,y0);
                for (int i=0; i<n; i++) {
                    DF(n*numPatchPts,i) = hPartials[i];
                    setCount++;
                }
            }//------------------------------------------------------------------------
        } //kk==0
        
        if (kk < numPatchPts-1) {
            //Continuity Constraints
            for (int i=0; i<n; i++) {
                F(n*kk+i) = yf[i] - X(n*(kk+1)+i);
            }
            if (outputF) {
                std::cerr << "    Continuity Constraint for Point " << kk << " :\n";
                for (int i=0; i<n; i++) {
                    std::cerr << "      F(" << n* kk+i << ") = " << F(n*kk+i) << "\n";
                }
            }
            //DF entries
            for (int i = 0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    DF(n*kk+j,n*kk+i) = finalInfo.J[j][i];
                    setCount++;
                }
            }
            //State derivative
            for (int i=0; i<n; i++) {
                DF(n*kk+i,n*numPatchPts) = dxdt[i] / (double) numPatchPts;
                setCount++;
            }
        } else {
            //Periodicity for kk==numPatchPts-1
            double& lambda = monitor->lambda;
            //Fp(X)
            for (int i=0; i<n-2; i++) {
                F(n*kk+i) = yf[i] - ((1.0-lambda)*gState[i] + lambda*X(i));
                //Homotopy Update vector elements
                dGdlambda(n*kk+i) = gState[i] - X(i);
            }
            F(n*kk+n-2) = yf[n-1] - ((1.0-lambda)*gState[n-1] + lambda*X(n-1));
            dGdlambda(n*kk+n-2) = gState[n-1] - X(n-1);
            dGdlambda(n*(kk+1)) = gState[n-2] - X(n-2); //ydot eqn
            if (outputF) {
                std::cerr << "    Periodicity Constraints:\n";
                for (int i=0; i<n-1; i++) {
                    std::cerr << "       F(" << n* kk+i << ") = " << F(n*kk+i) << "\n";
                }
            }
            //DF entries
            if (numPatchPts == 1) { //Just add result to DF since it's phi-I
                for (int i=0; i<n; i++) {
                    for (int j=0; j<n-2; j++) {
                        DF(n*kk+j,n*kk+i) += finalInfo.J[j][i];
                        setCount++;
                    }
                }
                for (int i=0; i<n; i++) { //Last STM row
                    DF(n*kk+n-2,n*kk+i) += finalInfo.J[n-1][i];
                    setCount++;
                }
            } else { //Normal setting of STM elements ------------------------------------
                for (int i=0; i<n; i++) {
                    for (int j=0; j<n-2; j++) {
                        DF(n*kk+j,n*kk+i) = finalInfo.J[j][i];
                        setCount++;
                    }
                }
                for (int i=0; i<n; i++) { //Last STM row
                    DF(n*kk+n-2,n*kk+i) = finalInfo.J[n-1][i];
                    setCount++;
                }
            } //End numPatchPts setting if statement
            for (int i=0; i<n; i++) {
                dGydotdX(0,n*kk+i) = finalInfo.J[n-2][i];    //Homotopy
            }
            //State Derivs
            for (int i=0; i<n-2; i++) {
                DF(n*kk+i,n*numPatchPts) = dxdt[i] / (double) numPatchPts;
                setCount++;
            }
            dGydotdX(0,n*numPatchPts) = dxdt[n-2] / (double) numPatchPts; //Homotopy
            DF(n*kk+n-2,n*numPatchPts) = dxdt[n-1] / (double) numPatchPts;
            setCount++;
        }
    } //End for loop
}//End computeF_and_DF();

/*/// Compute F and DF given X - For a sparse representation
template<typename MAP>
void PeriodicOrbitCorrector<MAP>::computeF_and_DF_Sparse(std::vector<Triple>& tuples) {
    static const int n = rhs_type::dimension;

    //Integrate Each state
    double dt = X(n*numPatchPts)/numPatchPts; //Last entry is time

    //Call the function to set the constants
    setConstantDFCoeffs(tuples);
    int setCount = n*(numPatchPts-1) + n-1;
    int lastSetCount = setCount;

    //For each patch point -> We could multi-thread this if we aren't solving fp's per thread
    for(int kk=0; kk<numPatchPts; kk++) {
        xstate_type y0(0.0);
        for (int i=0; i<n; i++) y0[i] = X(n*kk+i);
        for (int i=0; i<n; i++) y0[n+n*i+i] = 1.0; //Identity for STM0
    if (_verbose) {
        std::cout << "      Point " << kk;
        std::cout << " : y0 = " << y0 << '\n';
    }
        //Run the state through dt time
        return_state finalInfo = theMap->integrate_state(y0,dt);
        xstate_type yf = finalInfo.getState();
    if (_verbose) {
        std::cout << "      After Integ - Point " << kk;
        std::cout << " : y0 = " << y0 << '\n';

        std::cout << "      After Integ - Point " << kk;
        std::cout << " : yf = " << yf << '\n';
    }

        //STM = finalInfo.J
        //State derivative
        const rhs_type& rhs = theMap->rhs();
        xstate_type dxdt = rhs(0.0,yf);

        //Partials/constraints relevant to first point
        if (kk == 0) {
            int idx = n*(numPatchPts-1)+(n-1) + n*n*(numPatchPts-1) + n*(n-1) + n*(numPatchPts-1) + (n-1); //Idx5
            //Section Distance constraint
            F(n*numPatchPts - 1) = theMap->section().distance(y0);
        if (outputF) {
        std::cerr << "      Distance Constraint: F(" << n*numPatchPts-1
        << ") = " << F(n*numPatchPts-1) << '\n';
        }
            xstate_type sectPartials = theMap->section().distance_first_partials(y0);
            for (int i=0; i<n; i++) {
                tuples.push_back( Triple(n*numPatchPts - 1,i,sectPartials[i]) );
                idx++;
        setCount++;
            }
            //Hamiltonian Constraint
            F(n*numPatchPts) = theMap->rhs().hamiltonian(0.0,y0) - _desired_H;
        if (outputF) {
            std::cerr << "      Current JC = " << theMap->rhs().hamiltonian(0.0,y0)
        << "  F(" << n*numPatchPts << ") = " << F(n*numPatchPts) << '\n';
        }
            state_type hPartials = theMap->rhs().hamiltonian_first_partials(0.0,y0);
            for (int i=0; i<n; i++) {
                tuples.push_back( Triple(n*numPatchPts,i,hPartials[i]) );
                idx++;
        setCount++;
            }
        }

        if (kk < numPatchPts-1) {
            //Continuity Constraints
            for (int i=0; i<n; i++) F(n*kk+i) = yf[i] - X(n*(kk+1)+i);
        if (outputF) {
        std::cerr << "    Continutity Constraint for Point " << kk << " :\n";
        for (int i=0;i<n;i++) {
            std::cerr << "      F(" << n*kk+i << ") = " << F(n*kk+i) << "\n";
        }
        }
            //DF entries
            int idx = n*(numPatchPts-1)+(n-1); //Idx1
            for (int i = 0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    tuples.push_back( Triple(n*kk+j,n*kk+i,  finalInfo.J[j][i]) );
                    idx++;
            setCount++;
                }
            }
            //State derivative
            idx = n*(numPatchPts-1)+(n-1) + n*n*(numPatchPts-1) + n*(n-1); //Idx3
            for (int i=0; i<n; i++) {
                tuples.push_back( Triple(n*kk+i,n*numPatchPts, dxdt[i] / numPatchPts) );
                idx++;
        setCount++;
            }
        } else {
            //Periodicity for kk==numPatchPts-1
            for (int i=0; i<n-2; i++) F(n*kk+i) = yf[i] - X(i);
            F(n*kk+n-2) = yf[n-1] - X(n-1);
        if (outputF) {
        std::cerr << "    Periodicity Constraints:\n";
        for (int i=0;i<n-1;i++) {
            std::cerr << "       F(" << n*kk+i << ") = " << F(n*kk+i) << "\n";
        }
        }
            //DF entries
            int idx = n*(numPatchPts-1)+(n-1) + n*n*(numPatchPts-1); //Idx2
            for (int i=0; i<n; i++) {
                for (int j=0; j<n-2; j++) {
                    tuples.push_back( Triple(n*kk+j,n*kk+i, finalInfo.J[j][i]) );
                    idx++;
            setCount++;
                }
            }
            for (int i=0; i<n; i++) { //Last STM row
                tuples.push_back( Triple(n*kk+n-2,n*kk+i, finalInfo.J[n-1][i]) );
                idx++;
        setCount++;
            }
            //State Derivs
            idx = n*(numPatchPts-1)+(n-1) + n*n*(numPatchPts-1) + n*(n-1) + n*(numPatchPts-1); //Idx4
            for (int i=0; i<n-2; i++) {
                tuples.push_back( Triple(n*kk+i,n*numPatchPts, dxdt[i] / numPatchPts) );
                idx++;
        setCount++;
            }
            tuples.push_back( Triple(n*kk+n-2,n*numPatchPts,dxdt[n-1]) );
        setCount++;
        }

    if(outputDF) {
        std::cerr << "      Values: " << lastSetCount << " : " << setCount-1 << "\n";
        //Print the values set since the last patch point
        for (int i = lastSetCount; i< setCount; i++) {
            int row, col; double value;
            row = tuples[i].row();
            col = tuples[i].col();
            value = tuples[i].value();
            std::cerr << "      DF( " << row << " , " << col << " ) = " << value << '\n';
        }
    }
    lastSetCount = setCount;

    }
    //Projected size comparison
    if (_verbose) {
        std::cerr << " Setting coefficients:  size() = " << tuples.size() << "\n";
    int idxProjected = n*(numPatchPts-1) + n-1 + n + n*n*(numPatchPts-1) + (n-1)*n + n +n*(numPatchPts-1) + (n-1);
    std::cerr << "   Projected size = " << idxProjected
        << " (offset: " << idxProjected - tuples.size() << " )\n";
    }

}//end computeF_and_DF();
*/

/// Compute the change in X given an update (lambda) in periodicity homotopy method
template<typename MAP>
bool PeriodicOrbitCorrector<MAP>::computePeriodicityHomotopyPrediction(MatrixXd& DF)
{
    static const int n = rhs_type::dimension;
    reg_type* monitor = &_homotopyMonitor;
    //Fill Vectors
    VectorXd dY;
    dY.resize(n*numPatchPts+2);
    
    //Fill DG(Y) matrix
    MatrixXd DG(n*numPatchPts+1,n*numPatchPts+2);
    DG.setZero();
    DG.topLeftCorner(n*numPatchPts,n*numPatchPts+1) << DF;
    DG.row(n*numPatchPts).head(n*numPatchPts+1) << dGydotdX; //Row
    DG.col(n*numPatchPts+1) << dGdlambda; //Column
    
    //Compute Nullspace of DG to get direction
    MatrixXd nullDG = DG.fullPivLu().kernel(); //Note: not normalized
    //Check for 1D nullspace
    if (nullDG.cols()>1 || nullDG.cols()==0) {
        return false;    //check failure
    }
    //Set update vector
    dY << nullDG.col(0);
    dY.normalize(); //Normalize in place
    if (dY.tail(1)[0] == 0.0) {
        return false;    //Nowhere to step
    }
    if (dY.tail(1)[0] < 0.0) {
        dY *= -1.0;
    }
    //Scale to get the appropriate step
    double& ds = monitor->dlam;
    double& lam = monitor->lambda;
    double deltaLambda = 0.0;
    //Compute change in lambda and update vector
    if (lam == 1.0) {
        //All done - don't update
        dX.setZero();
    } else if ((dY.tail(1)[0]*ds + lam > 1.0) || (lam > 0.99) ) {
        //Steps past 1 so enforce a step to lam = 1.0
        deltaLambda = 1.0-lam;
        dX = deltaLambda * dY.head(n*numPatchPts+1);
    } else if (dY.tail(1)[0] >= 0.9) {
        //We can take a full step based on monitor->dlam
        double sf = ds/dY.tail(1)[0];
        deltaLambda = ds;
        dX = sf*dY.head(n*numPatchPts+1);
    } else {
        //Sensitive areas so take scaled steps
        deltaLambda = dY.tail(1)[0]*ds;
        dX = ds*dY.head(n*numPatchPts+1);
    }
    
    //Apply update
    monitor->incrementLambda(deltaLambda);
    X += dX;
    return true;
    
}

/// Compute the change in X given an update (gamma) in H-homotopy method
template<typename MAP>
bool PeriodicOrbitCorrector<MAP>::computeHamiltonianHomotopyPrediction(MatrixXd& DF)
{
    static const int n = rhs_type::dimension;
    reg_type* monitor = &_homotopyMonitor;
    
    //Fill Vectors
    VectorXd dW;
    dW.resize(n*numPatchPts+2);
    //Fill DE(W) matrix
    MatrixXd DE(n*numPatchPts+1,n*numPatchPts+2);
    DE.setZero();
    DE.topLeftCorner(n*numPatchPts+1,n*numPatchPts+1) << DF;
    DE(n*numPatchPts,n*numPatchPts+1) = _Hlambda - _desired_H; //Only non-zero element in last column (DE/dgamma)
    //Compute Nullspace
    MatrixXd nullDE = DE.fullPivLu().kernel(); //not normalized
    if (nullDE.cols()>1 || nullDE.cols() == 0) {
        return false;    //nullspace issues
    }
    //Set nullspace vector
    dW << nullDE.col(0);
    dW.normalize();
    if (dW.tail(1)[0] == 0.0) {
        return false;    //Nowhere to step
    }
    if (dW.tail(1)[0] < 0.0) {
        dW *= -1.0;
    }
    //Scale to get the appropriate step
    double& ds = monitor->dgam;
    double& gam = monitor->gamma;
    double deltaGamma = 0.0;
    //Compute change in gamma and update vector
    if (gam == 1.0) {
        //All done - don't update
        dX.setZero();
    } else if ((dW.tail(1)[0]*ds + gam > 1.0) || gam > 0.99) {
        //Steps past 1 so enforce a step to gam = 1.0
        deltaGamma = 1.0 - gam;
        dX = deltaGamma * dW.head(n*numPatchPts+1);
    } else if (dW.tail(1)[0] >= 0.9) {
        //We can take full step
        deltaGamma = ds;
        double sf = ds/dW.tail(1)[0];
        dX = sf*dW.head(n*numPatchPts+1);
    } else {
        //Sensitive step, scale down
        deltaGamma = dW.tail(1)[0]*ds;
        dX = ds*dW.head(n*numPatchPts+1);
    }
    
    //Apply update
    monitor->incrementGamma(deltaGamma);
    X += dX;
    return true;
    
}


/// The replacement for meta_newton: - The Function we call to compute fixed points
///   Use Multiple Shooting to solve for Fixed points in higher dimensional or sensitive systems
template<typename MAP>
bool solveFixedPoints(MAP& pmap, const metric<double,2>& metric, const nvis::bbox2& bounds,
                      const nvis::vec2& first_guess, const double& Hval,
                      const int& depth, const int& period,
                      fixpoint& fp, std::vector<fixpoint>& iterates,
                      const double& eps, bool verbose = false, const int maxiter=20,
                      const int numPtsPerPeriod=8, const bool linearSTM = false,
                      CorrectionsRegulator::OperatingMode stopMode = CorrectionsRegulator::THE_END)
{
    typedef rhs_only_wrapper<MAP,2> simple_residual;
    
    simple_residual err(pmap, metric, period);
    nvis::vec2 x = first_guess, ferr;
    try {
        ferr = err(x);
    } catch(...) {
        return false;
    }
    
    //See if we can improve the seed point - Not sure if this should be on or not...
    /*if (nvis::norm(ferr) > eps)
    {
        if (verbose) {
            std::cerr << "  norm at seed point (" << first_guess << ") = " << nvis::norm(ferr) << " is too large\n";
        }
        //Run the seed check
        bool ok = find_seed(err, metric, bounds, x, depth);
        if (!ok) {
            if (verbose) {
                std::cerr << "Unable to find seed\n";
            }
            return false;
        }
        else if (verbose) {
            std::cerr << "  improved seed found:\n"
                      << "   norm at new seed (" << x << ") is " << nvis::norm(err(x)) << "\n";
        }
    }*/
    
    //Disable the Hamiltonian check for corrections process
    bool isHamCheck = pmap.isCheckingHamiltonian();
    pmap.stopAtHamiltonianError(false);
    
    //Run the corrections algorithm
    PeriodicOrbitCorrector<MAP> corrector(pmap, x, Hval, period, numPtsPerPeriod, verbose);
    corrector.setMaxIters(maxiter);
    corrector.setTolerance(eps);
    corrector.setStopMode(stopMode);
    //Debugging - manually set modes
    //corrector.setOpMode(CorrectionsRegulator::MULTIPLE_SHOOT);
    //corrector.setOpMode(CorrectionsRegulator::QUASI_NEWTON);
    //corrector.setStopMode(CorrectionsRegulator::HAMSEC_OFF);
    //corrector.setMaxIters(2);
    
    //Setup the guess
    bool guessFailed = false;
    try {
        corrector.setGuess(x,Hval,period);
    } catch(...) {
        guessFailed = true;
    }
    
    int thread_id = 0;
#if _OPENMP
    thread_id = omp_get_thread_num();
#endif
    corrector.setJobIndex(thread_id);
    //Stop here if guess failed
    if(guessFailed) {
        return false;
    }
    
    CorrectionResult theResult;
    try {
        theResult = corrector.correct();
    } catch(...) {
        //Some form of error is thrown: shouldn't happen, but seems to be!
        //Sometimes, this is an imaginary velocity ???
        return false;
    }
    bool found = theResult.converged;
    if (!found) {
        return false;
    } else {
        //Compute Map iterates and gather fixed-point information
        iterates.clear();
        bool checkOK = true;
        try {
            checkOK = corrector.setFixedPoints(fp,iterates,linearSTM);
        } catch(...) {
            checkOK = false;
        }
        //Result must pass sanity checks
        if (!checkOK) {
            return false;
        }
    }
    
    //if (!bounds.inside(x)) {
    //Do something if we travel outside bounds
    //} else {
    //}
    
    //Disable Hamiltonian check (if off beforehand)
    if(!isHamCheck) {
        pmap.stopAtHamiltonianError(false);
    }
    
    //If no errors,
    return true;
}

}//namespace orbital

#endif //corrections.hpp
