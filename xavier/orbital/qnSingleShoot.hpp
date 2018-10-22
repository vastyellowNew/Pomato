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


/** Functional for Single Shooting Periodicity Targeting problem
 *  - Computing the single shooting periodicity constraints F(X)
 *  and Jacobian DF(X) matrix.
 *
 * Author - Wayne Schlei (Purdue University)
 * Note:  This is for use with Map Topology extraction code
 */

#ifndef SS_PERIODICITY_HPP
#define SS_PERIODICITY_HPP

#include <vector>
#include <limits>
//Eigen Lib
#include <Eigen/Core>
#include <Eigen/Dense>


using namespace std;

namespace orbital {

template<class MAP>
struct SingleShooter { //Derivatives
    MAP* theMap; //Poincare_map class : from xavier nm
    MatrixXd DF;
    int dim;
    double _desired_H;
    int rows, cols;
    
    ///Constructor
    SingleShooter(MAP& amap,const double Hval) :
        theMap(&amap),
        dim(MAP::rhs_type::dimension),
        _desired_H(Hval),
        rows(dim+1),
        cols(dim+1)
    {
        DF.resize(rows,cols);
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
        fvec.reserve(rows);
        DF.setZero(); //Init to zero
        //Set constant DF elements
        for(int i=0; i<dim-2; i++) {
            DF(i,i) = -1.0;
        }
        DF(dim-2,dim-1) = -1.0; //Periodicity
        
        //Check integration time is valid (has a tendancy to go to zero for bad guesses)
        double dt = x.back(); //Last entry is time-period
        if (dt < 0.05) {
            std::cerr << "QN:  Single-shooting process contains a period that is impractically small or negative\n";
            throw("Invalid integration time for periodic orbit");
        }
        
        //Run integration
        xstate_type y0(0.0);
        for (int i=0; i<dim; i++) {
            y0[i] = x[i];
        }
        for (int i=0; i<dim; i++) {
            y0[dim+dim*i+i] = 1.0;    //STM0
        }
        return_state finalInfo;
        try {
            finalInfo = theMap->integrate_state(y0,dt);
        } catch(...) {
            std::cerr << "Map error occured during integration in Quasi-Netwon process.\n";
            throw("Integration Error");
        }
        xstate_type yf = finalInfo.getState();
        
        //Peridoicity constraints - dim-1 equations
        for (int i=0; i<dim-2; i++) {
            fvec.push_back(yf[i] - y0[i]);
        }
        fvec.push_back(yf[dim-1] - y0[dim-1]);
        
        //Final state derivative
        const rhs_type& rhs = theMap->rhs();
        xstate_type dxdt = rhs(0.0,yf);
        
        //Periodicity DF
        for (int j=0; j<dim-2; j++) {
            DF(j,dim) = dxdt[j];
            for (int i=0; i<dim; i++) {
                DF(j,i) += finalInfo.J[j][i];
            }
        }
        for (int i=0; i<dim; i++) {
            DF(dim-2,i) += finalInfo.J[dim-1][i];
        }
        DF(dim-2,dim) = dxdt[dim-1];
        
        //Section constraint
        fvec.push_back( theMap->section().distance(y0) );
        xstate_type sectPartials = theMap->section().distance_first_partials(y0);
        
        //Hamiltonian constraint
        fvec.push_back( theMap->rhs().hamiltonian(0.0,y0) - _desired_H );
        state_type hPartials = theMap->rhs().hamiltonian_first_partials(0.0,y0);
        
        //DF for section/hamiltonian constraints
        for (int i=0; i<dim; i++) {
            DF(dim - 1, i ) = sectPartials[i];
            DF(dim, i ) = hPartials[i];
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