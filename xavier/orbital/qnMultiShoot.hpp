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


/** Functional for Multiple Shooting Periodicity Targetting problem
 *  - Computing the multiple shooting periodicity constraints F(X)
 *  and Jacobian DF(X) matrix.
 *
 * Author - Wayne Schlei (Purdue University)
 * Note:  This is for use with Map Topology extraction code
 */

#ifndef MS_PERIODICITY_HPP
#define MS_PERIODICITY_HPP

#include <vector>
#include <limits>
//Eigen Lib
#include <Eigen/Core>
#include <Eigen/Dense>


using namespace std;

namespace orbital {

/// Functor for Quasi-Newton root-finding method
template<class MAP>
struct MultipleShooter {
    MAP* theMap; //Poincare_map class : from xavier nm
    MatrixXd DF;
    int dim, numPatchPts;
    double _desired_H;
    int rows, cols;
    
    ///Constructor
    MultipleShooter(MAP& amap,const int k,const double Hval) :
        theMap(&amap),
        dim(MAP::rhs_type::dimension),
        numPatchPts(k),
        _desired_H(Hval),
        rows(dim*numPatchPts+1),
        cols(dim*numPatchPts+1)
    {
        //Fill out DF that remains unchanged
        DF.resize(rows,cols);
        DF.setZero(); //Init to zero
        
        //Set constant DF elements -------------------------------
        int n=dim;
        for(int kk=0; kk<(numPatchPts-1); kk++) {
            //Continuity for all but last point
            for(int jj=0; jj<n; jj++) {
                DF(jj+kk*n,jj+n*(kk+1)) = -1.0;
            }
        }
        //Periodicity - Assume skips n-2 state
        for(int kk=0; kk<(n-2); kk++) {
            DF(n*(numPatchPts-1)+kk,kk) = -1.0;
        }
        DF(n*(numPatchPts-1)+(n-2),n-1) = -1.0;
        
    }
    
    ///Boolean to show if derivatives are supplied
    bool derivsSupplied()
    {
        return true;
    }
    
    ///Function call
    void eval(vector<double>& x,vector<double>& fvec)
    {
        typedef typename MAP::rhs_type     rhs_type;
        typedef typename MAP::state_type state_type;
        typedef typename MAP::xstate_type xstate_type;
        typedef typename MAP::return_state return_state;
        
        //Compute F and DF based on X
        fvec.clear();
        fvec.resize(rows,0.0);
        
        //Run integration for each arc ----------------------------------
        int n = dim;
        double dt = x[n*numPatchPts] / (double) numPatchPts; //Last entry is time-period
        for (int kk=0; kk<numPatchPts; kk++ ) {
            xstate_type y0(0.0);
            for (int i=0; i<n; i++) {
                y0[i] = x[n*kk+i];
            }
            for (int i=0; i<n; i++) {
                y0[n+n*i+i] = 1.0;    //STM0
            }
            return_state finalInfo;
            try {
                finalInfo = theMap->integrate_state(y0,dt);
            } catch(...) {
                std::cerr << "Map error occured during integration in Quasi-Netwon process.\n";
                throw("Integration Error");
            }
            xstate_type yf = finalInfo.getState();
            //Final state derivative
            const rhs_type& rhs = theMap->rhs();
            xstate_type dxdt = rhs(0.0,yf);
            
            //Constraints relative to first point
            if (kk==0) {
                //Section constraint
                fvec[n*numPatchPts - 1] =  theMap->section().distance(y0);
                xstate_type sectPartials = theMap->section().distance_first_partials(y0);
                
                //Hamiltonian constraint
                fvec[n*numPatchPts] = theMap->rhs().hamiltonian(0.0,y0) - _desired_H;
                state_type hPartials = theMap->rhs().hamiltonian_first_partials(0.0,y0);
                
                //DF for section/hamiltonian constraints
                for (int i=0; i<n; i++) {
                    DF(n*numPatchPts - 1, i ) = sectPartials[i];
                    DF(n*numPatchPts, i ) = hPartials[i];
                }
            }//End first point constraints
            
            //Continuity & Derivatives
            if (kk < numPatchPts-1) {
                for(int i=0; i<n; i++) {
                    fvec[n*kk+i] = yf[i] - x[n*(kk+1)+i];
                }
                for (int i=0; i<n; i++) {
                    for (int j=0; j<n; j++) {
                        DF(n*kk+j,n*kk+i) = finalInfo.J[j][i];
                    }
                }
                for (int i=0; i<n; i++) {
                    DF(n*kk+i,n*numPatchPts) = dxdt[i] / (double) numPatchPts;
                }
                
            } else {
                //Periodicity constraints
                for(int i=0; i<n-2; i++) {
                    fvec[n*kk+i] = yf[i] - x[i];
                }
                fvec[n*kk+n-2] = yf[n-1] - x[n-1];
                //Jacobian elements
                for(int i=0; i<n; i++) {
                    for (int j=0; j<n-2; j++) {
                        DF(n*kk+j,n*kk+i) = finalInfo.J[j][i];
                    }
                    DF(n*kk+n-2,n*kk+i) = finalInfo.J[n-1][i];
                }
                //State derivatives
                for (int i=0; i<n-2; i++) {
                    DF(n*kk+i,n*numPatchPts) = dxdt[i] / (double) numPatchPts;
                }
                DF(n*kk+n-2,n*numPatchPts) = dxdt[n-1] / (double) numPatchPts;
            }
            
        }
    }
    
    ///Jacobian Matrix : First-order derivatives
    void df(vector<double>& x, MatrixXd& df)
    {
        df = DF;
    }
};

} //orbital

#endif
