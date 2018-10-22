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


////////////////////////////////////////////////////////////////////
//  Corrections Control header file
//  Author:  Wayne Schlei
//  Date: 1/23/2014
//
//  Purpose:  For controlling/regulating differential corrections
//  procedures using multiple shooting for Fixed-point computation.
////////////////////////////////////////////////////////////////////
#ifndef CORRECTIONS_CONTROL_HPP
#define CORRECTIONS_CONTROL_HPP

#include <iostream>
#include <map>
#include <string>

namespace orbital {

///Corrections Data Structure
struct CorrectionResult {
    CorrectionResult(bool success, int numIters, double initError, double finalError) :
        converged(success), iterations(numIters), normF0(initError), normF(initError) {}
    CorrectionResult() : converged(false), iterations(0), normF0(100.0), normF(100.0) {}
    
    bool      converged;
    int       iterations;
    double    normF0, normF;
    //EnumConstraint culprit; //Try to identify the problem?
};


/// Class to monitor and modify the corrections progress
class CorrectionsRegulator {
public :
    typedef CorrectionsRegulator     self_type;
    /**Operating mode for the corrections process
     * corrections.hpp::Corrector class will work down this list, in order, until
     * the mode indicated by 'stoppingMode' or QUASI_NEWTON.
     */
    enum OperatingMode {
        SINGLE_SHOOT = 0,        //Single shooting from first point forward
        RSINGLE_SHOOT,        //Single shooting in backward time
        MULTIPLE_SHOOT,        //Multiple shooting (all forward)
        RESAMPLE,         //Add a point/rev & resample forward
        DISTRIBUTED_ERROR,    //Sample forward and backward to reduce error in periodicity constraints
        HAMSEC_OFF,        //Hamiltonian and Section constraints set to OFF
        PERIODICITY_HOMOTOPY,    //Transistion Periodicity constraint to ON (Then, H)
        QUASI_NEWTON,        //Try full problem with a linesearch for steps with backtracking
        THE_END            //Artifical last entry for stopping
    };
    
    ///Constructor: Default so mode is 0
    CorrectionsRegulator() :
        mode(self_type::SINGLE_SHOOT),
        modesExplored(false),
        useExtraPointPerRev(false),
        hsConstraintsOn(true),
        lineSearchBacktrack(false),
        processingSteps(false),
        stoppingMode(self_type::THE_END),
        found(false),
        lambda(1.0),
        dlam(0.1),
        gamma(1.0),
        dgam(0.2),
        jobID(0),
        _verbose(false),
        ppr0Err(1000.0), ppr1Err(1001.0)
    {
        //Build name map for enums
        modeName.insert( std::pair<OperatingMode,std::string>(SINGLE_SHOOT,"SINGLE_SHOOT") );
        modeName.insert( std::pair<OperatingMode,std::string>(RSINGLE_SHOOT,"RSINGLE_SHOOT") );
        modeName.insert( std::pair<OperatingMode,std::string>(MULTIPLE_SHOOT,"MULTIPLE_SHOOT") );
        modeName.insert( std::pair<OperatingMode,std::string>(RESAMPLE,"RESAMPLE") );
        modeName.insert( std::pair<OperatingMode,std::string>(DISTRIBUTED_ERROR,"DISTRIBUTED_ERROR") );
        modeName.insert( std::pair<OperatingMode,std::string>(HAMSEC_OFF,"HAMSEC_OFF") );
        modeName.insert( std::pair<OperatingMode,std::string>(PERIODICITY_HOMOTOPY,"PERIODICITY_HOMOTOPY") );
        modeName.insert( std::pair<OperatingMode,std::string>(QUASI_NEWTON,"QUASI_NEWTON") );
        modeName.insert( std::pair<OperatingMode,std::string>(THE_END,"THE_END") );
        /*modeName[SINGLE_SHOOT] = std::string("SINGLE_SHOOT");
        modeName[RSINGLE_SHOOT] = std::string("RSINGLE_SHOOT");
        modeName[MULTIPLE_SHOOT] = std::string("MULTIPLE_SHOOT");
        modeName[RESAMPLE] = std::string("RESAMPLE");
        modeName[DISTRIBUTED_ERROR] = std::string("DISTRIBUTED_ERROR");
        modeName[HAMSEC_OFF] = std::string("HAMSEC_OFF");
        modeName[PERIODICITY_HOMOTOPY] = std::string("PERIODICITY_HOMOTOPY");
        modeName[QUASI_NEWTON] = std::string("QUASI_NEWTON");
        modeName[THE_END] = std::string("THE_END");*/
    }
    /*///Constructor: Set mode
    CorrectionsRegulator(const int oppMode, const double lam, const double deltaLam ) :
          mode(oppMode),
          modesExplored(true),
          processingSteps(false),
          lambda(lam),
          dlam(deltaLam),
          gamma(lam),
          dgam(deltaLam)
    {}*/
    
    ///Set the job index
    void setJobID(const int& jID)
    {
        jobID = jID;
    }
    void setVerbose(const bool& v)
    {
        _verbose = v;
    }
    
    ///Update to compute new changes to module - Returns boolean stating if more work is needed
    bool update(const CorrectionResult& result)
    {
        if (mode == self_type::MULTIPLE_SHOOT) {
            //Store the initial normF result
            ppr0Err = result.normF0;
        } else if (mode == self_type::RESAMPLE) {
            //Store the initial normF
            ppr1Err = result.normF0;
            //Pick which sampling is better
            if (ppr0Err>ppr1Err) {
                useExtraPointPerRev = true;
            }
        }
        
        //The result of THIS STEP determines what to do next
        bool allDone = modeModifier(result.converged);
        //Check if this was the last mode
        bool lastMode = (mode == stoppingMode) ? true : false;
        
        //Check if we have reached/passed the stopping mode still with failure
        if(lastMode && !allDone) {
            allDone = true;
            modesExplored = true;
        }
        //Return whether or not there is still more to do:
        return allDone;
    }
    
    ///Did the method work and found a periodic solution
    bool solutionFound()
    {
        return found;
    }
    
    ///Increment lambda externally
    void incrementLambda(const double& delta)
    {
        lambda += delta;
    }
    ///Increment gamma externally
    void incrementGamma(const double& delta)
    {
        gamma += delta;
    }
    
    
    /// Function for outputing current configuration information
    void verbose() const
    {
        if (!_verbose) {
            return;
        }
        //Print the mode and the current step info
        std::cerr << "------------------------------------------------------------------------------------\n";
        switch (mode) {
            case SINGLE_SHOOT :
                std::cerr << "(" << jobID << ") Periodic Orbit Corrector:  1) SINGLE SHOOTING METHOD\n";
                //std::cerr << "    -Computing differential corrections based on full-period, forward propagation.\n";
                break;
            case RSINGLE_SHOOT :
                std::cerr << "(" << jobID << ") Periodic Orbit Corrector:  2) SINGLE SHOOTING METHOD with Reverse Time\n";
                //std::cerr << "    -Computing differential corrections based on full-period, backward propagation.\n";
                break;
            case MULTIPLE_SHOOT :
                std::cerr << "(" << jobID << ") Periodic Orbit Corrector: 3) MULTIPLE SHOOTING METHOD\n";
                /*std::cerr << "    -Differential correction with multiple 'patch' points with the indicated \n";
                std::cerr << "    number of points per revolution (which are actually distributed throughout\n";
                std::cerr << "    the forward propagation in time).\n";*/
                break;
            case RESAMPLE :
                std::cerr << "(" << jobID << ") Periodic Orbit Corrector: 4) MULTIPLE SHOOTING METHOD +1\n";
                /*std::cerr << "    -Differential correction with multiple 'patch' points with the indicated \n";
                std::cerr << "    number of points per revolution PLUS an extra patch point (which are \n";
                std::cerr << "    actually distributed throughout the forward propagation in time).\n";*/
                break;
            case DISTRIBUTED_ERROR :
                std::cerr << "(" << jobID << ") Periodic Orbit Corrector: 5) MULTIPLE SHOOTING METHOD with F/B Split Sampling\n";
                /*std::cerr << "    -Differential correction with multiple 'patch' points with sampling split\n";
                std::cerr << "    half-way between forward and backward time propagation.\n";*/
                break;
            case HAMSEC_OFF :
                std::cerr << "(" << jobID << ") Periodic Orbit Corrector: 6) MULTIPLE SHOOTING METHOD with H&Section Switch\n";
                /*std::cerr << "    -Differential correction with multiple 'patch' points with split sampling\n";
                std::cerr << "    half-way between forward and backward time propagation. The Hamiltonian and\n";
                std::cerr << "    Section constraints are OFF for the first run, then ON for the second.\n";*/
                break;
            case PERIODICITY_HOMOTOPY :
                std::cerr << "(" << jobID << ") Periodic Orbit Corrector: 7) MULTIPLE SHOOTING METHOD with PERIODICITY HOMOTOPY\n";
                /*std::cerr << "    -Differential correction with multiple 'patch' points with FORWARD sampling\n";
                std::cerr << "    and the application of homotopy to gradually enable the periodicity and \n";
                std::cerr << "    Hamiltonian constraints.\n";*/
                break;
            case QUASI_NEWTON :
                std::cerr << "(" << jobID << ") Periodic Orbit Corrector: 8) QUASI-NEWTON METHOD with Single Shooting\n";
                /*std::cerr << "    -Differential correction with multiple 'patch' points with split sampling\n";
                std::cerr << "    and the application of a linesearch with backtracking on the update step.\n";*/
                break;
            default :
                std::cerr << "(" << jobID << ") CorrectionsRegulator:  Modes explored.\n";
        }
    }
    
    /// Set the stopping mode - The mode after the last one you want to use.
    void setStoppingMode(OperatingMode& stopper)
    {
        stoppingMode = stopper;
    }
    /// Get the stopping mode
    OperatingMode getStoppingMode()
    {
        return stoppingMode;
    }
    /// Printing a method name to a string
    std::string getModeName()
    {
        return modeName[mode];
    }
    /// Set the operating mode (Important: use this to set it manually and not follow the recipe!)
    void setMode(OperatingMode& theMode)
    {
        switch (theMode) {
            default :
                //SINGLE_SHOOT, RSINGLE_SHOOT, MULTIPLE_SHOOT, DISTRIBUTED_ERROR
                processingSteps = false;
                hsConstraintsOn = true;
                lambda = gamma = 1.0;
                useExtraPointPerRev = false; //in update() this depends on other stuff
                break;
            case RESAMPLE :
                processingSteps = false;
                hsConstraintsOn = true;
                lambda = gamma = 1.0;
                useExtraPointPerRev = true;
                break;
            case HAMSEC_OFF :
                processingSteps = false;
                hsConstraintsOn = false;
                processingSteps = false;
                lambda = gamma = 1.0;
                useExtraPointPerRev = false;
                break;
            case PERIODICITY_HOMOTOPY :
                hsConstraintsOn = false;
                processingSteps = false;
                useExtraPointPerRev = false;
                //Set homotopy parameters to zero for processing
                lambda = 0.0;
                gamma = 0.0;
                break;
            case QUASI_NEWTON :
                hsConstraintsOn = true;
                processingSteps = false;
                useExtraPointPerRev = false;
                lineSearchBacktrack = true;
                lambda = gamma = 1.0;
                break;
        }
        mode = theMode;
    }
    
    //Member Variables: --------
    /// The Operating mode (what is currently being done)
    OperatingMode mode;
    /// Bool to see if all modes are explored
    bool modesExplored;
    /// Bool for whether or not to use k or k+1 pts per rev (distributed in time)
    bool useExtraPointPerRev;
    /// Bool for turning both Hamiltonian and section constraints ON/OFF
    bool hsConstraintsOn;
    /// Bool for using a line search with backtracking in Quasi-Newton
    bool lineSearchBacktrack;
    /// Bool indicating that the method is reusing the last solution as a guess for next step
    bool processingSteps;
    /// Indicator for what mode to stop on
    OperatingMode stoppingMode;
    /// Indicator for solution found
    bool found;
    ///Homotopy parameters - Periodicity
    double lambda, dlam;
    ///Homotopy parameters - Jacobi
    double gamma, dgam;
    ///Job ID
    int jobID;
    ///Verbose
    bool _verbose;
    ///Error values to compute what sampling to use
    double ppr0Err, ppr1Err;
    
private :
    std::map<OperatingMode,std::string> modeName;
    void stepFailedVerbose()
    {
        if (_verbose) {
            std::cerr << "(" << jobID << ") Approach " << getModeName() << " FAILED\n";
        }
    }
    void stepProcessingVerbose()
    {
        if (_verbose) {
            std::cerr << "(" << jobID << ") Approach " << getModeName() << " PROCESSING ...\n";
        }
    }
    void stepSuccessVerbose()
    {
        if (_verbose) {
            std::cerr << "(" << jobID << ") Approach " << getModeName() << " SUCCESS!\n";
        }
    }
    ///Modify the mode according to what happened
    bool modeModifier(const bool& converged)
    {
        bool allDone = false;
        switch(mode) {
            case SINGLE_SHOOT :
                if(converged) {
                    allDone = true;
                    found = true;
                } else {
                    stepFailedVerbose();
                    //Try backward sampling
                    mode = self_type::RSINGLE_SHOOT;
                    processingSteps = false;
                }
                break;
            case RSINGLE_SHOOT :
                if(converged) {
                    allDone = true;
                    found = true;
                } else {
                    stepFailedVerbose();
                    //Try multiple shooting
                    mode = self_type::MULTIPLE_SHOOT;
                    processingSteps = false;
                }
                break;
            case MULTIPLE_SHOOT :
                if(converged) {
                    allDone = true;
                    found = true;
                } else {
                    stepFailedVerbose();
                    //Go to RESAMPLE mode and try again
                    mode = self_type::RESAMPLE;
                    processingSteps = false;
                    useExtraPointPerRev = true;
                }
                break;
            case RESAMPLE :
                if(converged) {
                    allDone = true; //It worked!
                    found = true;
                } else {
                    stepFailedVerbose();
                    //Try the problem with different sampling
                    //that reduces peridocity-constraint error
                    mode = self_type::DISTRIBUTED_ERROR;
                    processingSteps = false;
                }
                break;
            case DISTRIBUTED_ERROR :
                if(converged) {
                    allDone = true; //It worked!
                    found = true;
                } else {
                    stepFailedVerbose();
                    //Try the problem with H&Sec constraints off
                    mode = self_type::HAMSEC_OFF;
                    hsConstraintsOn = false;
                    processingSteps = false;
                }
                break;
            case HAMSEC_OFF :
                if (converged && hsConstraintsOn) {
                    allDone = true; //It worked!
                    found = true;
                } else if(converged && !hsConstraintsOn) {
                    //Turn on constraints and reconverge
                    hsConstraintsOn = true;
                    processingSteps = true;
                    stepProcessingVerbose();
                } else {
                    stepFailedVerbose();
                    //Skipping periodicity homotopy as it is not the best for periodicity constraints
                    mode = self_type::QUASI_NEWTON;
                    hsConstraintsOn = true;
                    lineSearchBacktrack = true;
                    processingSteps = false;
                    
                    /*//Go to periodicity homotopy
                    mode = self_type::PERIODICITY_HOMOTOPY;
                    hsConstraintsOn = false;
                    processingSteps = false;
                    //Set homotopy parameters to zero for processing
                    lambda = 0.0;
                    gamma = 0.0;*/
                }
                break;
            case PERIODICITY_HOMOTOPY : //Temporarily skipping homotopy.......
                if (converged && lambda>=1.0 && gamma >=1.0 ) {
                    allDone = true; //It worked
                    found = true;
                } else if (converged && (lambda >= 1.0) ) {
                    //Resolve true orbit by switching to homotopy with JC constraint
                    processingSteps = true;
                    hsConstraintsOn = true;
                    stepProcessingVerbose();
                } else if (converged && (lambda<1.0)) {
                    //Step performed in corrections.hpp::correct() as part of prediction
                    processingSteps = true;
                    stepProcessingVerbose();
                } else {
                    //Couldn't solve
                    stepFailedVerbose();
                    //Let's try Quasi-Newton
                    mode = self_type::QUASI_NEWTON;
                    hsConstraintsOn = true;
                    lineSearchBacktrack = true;
                    lambda = 1.0;
                    gamma = 1.0;  //Turn constraints fully on
                }
                break;//...........................................................
            case QUASI_NEWTON :
                //Last option : converged or not, the buck stops here.
                modesExplored = true;
                allDone = true;
                if (converged) {
                    found = true;
                } else {
                    stepFailedVerbose();
                    mode = self_type::THE_END;
                }
                break;
            default :
                //Nothing to do cause it's THE_END!!!
                break;
        }
        if (found) {
            stepSuccessVerbose();
        }
        return allDone;
    } //End modeModifier()
};

}//Namespace orbital


#endif //CORRECTIONS_CONTROL_HPP
