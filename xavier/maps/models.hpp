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


/** Functor to represent the linear model around a saddle point.
    Used for Levengberg-Marquardt method in Eigen
*/
#ifndef __LINEAR_MODEL_HPP
#define __LINEAR_MODEL_HPP

#include <Eigen/Core>
#include <vector>

using namespace Eigen;

namespace xavier {

// Generic functor  (same as the provided functor in NonLinearOptimization.cpp)
template<typename _Scalar, int NX=Dynamic, int NY=Dynamic>
struct VectorFunctor {
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Matrix<Scalar,NY,InputsAtCompileTime> JacobianType;
    
    const int m_inputs, m_values;
    
    VectorFunctor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    VectorFunctor(int inputs, int values) : m_inputs(inputs), m_values(values) {}
    
    int inputs() const
    {
        return m_inputs;
    }
    int values() const
    {
        return m_values;
    }
    
};

/// Linear model functor
struct LinearModelFunctor : public VectorFunctor<double> {
public:
    ///Constructor
    LinearModelFunctor(const int& spaceDim,const int numParams, const int numData) :
        VectorFunctor<double>(numParams,numData),
        n(spaceDim), nParam(numParams), nData(numData)
    {}
    
    //Dimension of map space
    int n; //NX = n*(n+1)
    int nParam; //number of parameters
    int nData; //number of available data points
    //Members that store data
    std::vector<VectorXd> x; //x vectors
    std::vector<VectorXd> y; //y vectors
    /**When setting data points (x-input value, y-output value)
        be sure to write like this example
        \code
        LinearModelFunctor<numParams,numDataPoints> linModel(spaceDim);
        for (int k=0; k<numDataPoints; k++) {
            VectorXd xPoint(spaceDim), yPoint(spaceDim);
            //... Set xPoint, yPoint ...//
            xPoint << 1, 0, ... val[spaceDim-1];
            yPoint << 23, -20, ... val[spaceDim-1];
            linModel.x.push_back(xPoint);
            linModel.y.push_back(yPoint);
        }
        \endcode
    */
    //Required call for NonLinearOptimization routines
    int operator()(const VectorXd& beta, VectorXd& fvec) const
    {
        //make sure sizes are right
        assert(beta.size()==nParam);
        assert(fvec.size()==nData);
        //beta is parameter vector, consisting of (xm,A.row(0),A.row(1),...)
        VectorXd xs(n), xdot(n), diff(n);
        MatrixXd A(n,n);
        xs = beta.head(n);
        for (int j=0; j<n; j++) {
            A.row(j) = beta.segment(n+j*n,n);
        }
        //For each value provided, compute fvec value (=chi^2)
        for(int i=0; i<nData; i++) {
            //The Linear Model
            xdot = A*(x[i] - xs);
            //Chi^2 or error
            diff = xdot - y[i];
            fvec[i] = diff.dot(diff);
        }
        //success
        return 0;
    }
    
    //Analytical Derivatives:  NOTE: There is an error in these somewhere!
    /*int df(const VectorXd &beta, MatrixXd &fjac)
    {
        assert(beta.size() == nParam);
        assert(fjac.rows() == nData);
        assert(fjac.cols() == nParam);
    
        //Basically re-evaluate fvec given beta and data
        VectorXd xs(n), xdot(n), diff(n), dx(n);
        MatrixXd A(n,n);
        xs = beta.head(n);
        for (int j=0;j<n;j++) A.row(j) = beta.segment(n+j*n,n);
        //Compute each row => each data value
        for (int row=0; row<nData; row++) {
            MatrixXd paramJacobian(n,nParam);
            paramJacobian.setZero();
            //The Linear Model
            dx = x[row] - xs;
            xdot = A*dx;
            //Chi^2 or error
            diff = xdot - y[row];
            //Fill entries of parameter Jacobian (dxdot/dbeta)
            paramJacobian.leftCols(n) = -1.0*A; //dxdot/dxs
            for (int j=0; j<n; j++) {
                for (int i=0;i<n;i++) paramJacobian(j,j*n+i) = dx[i];
            }
    
            //Fill out Jacobian row
            fjac.row(row) = 2.0*diff.transpose()*paramJacobian;
    
        }
        //success
        return 0;
    }*/
    //Can do Numerical Differentiation using Eigen
    
    ///Evaluate based on parameters
    VectorXd eval(const VectorXd& beta, const VectorXd& xVals)
    {
        assert(beta.size() == nParam);
        assert(xVals.size() == n);
        VectorXd xs(n), xdot(n), diff(n);
        MatrixXd A(n,n);
        xs = beta.head(n);
        for (int j=0; j<n; j++) {
            A.row(j) = beta.segment(n+j*n,n);
        }
        xdot = A*(xVals - xs);
        return xdot;
    }
};


/// Quadratic model functor - force in 2D case for now
//  Unfotunately, this is hard to generalize (derivatives).  But might be
//  possible with complex-step derivatives
struct QuadModel2DFunctor : public VectorFunctor<double> {
public:
    ///Constructor
    QuadModel2DFunctor(const int numParams, const int numData) :
        VectorFunctor<double>(numParams,numData),
        n(2), nParam(numParams), nData(numData)
    {}
    
    //Dimension of map space
    int n; //NX = [n*(n+1) + n*n*n] parameters
    int nParam; //number of parameters
    int nData; //number of available data points
    
    //Members that store data
    std::vector<VectorXd> x; //x vectors
    std::vector<VectorXd> y; //y vectors
    /**When setting data points (x-input value, y-output value)
        be sure to write like this example
        \code
        QuadModelFunctor<numParams,numDataPoints> quadModel(spaceDim);
        for (int k=0; k<numDataPoints; k++) {
            VectorXd xPoint(spaceDim), yPoint(spaceDim);
            //... Set xPoint, yPoint ...//
            xPoint << 1, 0, ... val[spaceDim-1];
            yPoint << 23, -20, ... val[spaceDim-1];
            quadModel.x.push_back(xPoint);
            quadModel.y.push_back(yPoint);
        }
        \endcode
    */
    //Required call for NonLinearOptimization routines
    int operator()(const VectorXd& beta, VectorXd& fvec) const
    {
        //make sure sizes are right
        assert(beta.size()==nParam);
        assert(fvec.size()==nData);
        //beta is parameter vector, consisting of (xm,A.row(0),A.row(1),...)
        VectorXd xs(n), xdot(n), diff(n), dx(n);
        MatrixXd A(n,n), partialQdx(n,n);
        //Q->Tensor or n nxn matrices
        MatrixXd Jx(n,n), Jy(n,n);
        partialQdx.setZero(n,n);
        Jx.setZero();
        Jy.setZero();
        xs = beta.head(n);
        for (int j=0; j<n; j++) {
            A.row(j) =  beta.segment(n+j*n          , n);
            Jx.row(j) = beta.segment(n+n*n +j*n     , n);
            Jy.row(j) = beta.segment(n+n*n +n*n +j*n, n);
        }
        //The Quad Model
        //For each value provided, compute fvec value (=chi^2)
        for(int i=0; i<nData; i++) {
            assert(beta.size() == nParam);
            assert(fvec.size() == nData);
            dx = x[i] - xs;
            partialQdx = Jx*dx[0] + Jy*dx[1];
            xdot = A*dx + 0.5*partialQdx.transpose()*dx;
            //Chi^2 or error
            diff = xdot - y[i];
            fvec[i] = diff.dot(diff);
        }
        //success
        return 0;
    }
    
    //Analytical Derivatives : (Hard to generalize!)
    //Note: There's an error in here somewhere!
    /*int df(const VectorXd &beta, MatrixXd &fjac)
    {
        assert(beta.size() == nParam);
        assert(fjac.rows() == nData);
        assert(fjac.cols() == nParam);
        //Basically re-evaluate fvec given beta and data
        VectorXd xs(n), xdot(n), diff(n), dx(n);
        MatrixXd A(n,n), partialQdx(n,n);
        //Q->Tensor or n nxn matrices
        std::vector<MatrixXd> Q(n,MatrixXd(n,n));
        partialQdx.setZero(n,n);
        xs = beta.head(n);
        for (int j=0;j<n;j++) {
            A.row(j) = beta.segment(n+j*n,n);
            for (int i=0;i<n;i++) {
                int start = n+n*n +i*n*n + n*j;
                Q[i].row(j) = beta.segment(start,n);
            }
        }
    
        //Compute each row => each data value
        for (int row=0; row<nData; row++) {
            MatrixXd paramJacobian(n,nParam);
            paramJacobian.setZero(n,nParam);
            //The Quad Model
            dx = x[row] - xs;
            for (int i=0;i<n;i++) partialQdx += Q[i]*dx[i];
            xdot = A*dx + 0.5*partialQdx.transpose()*dx;
            //Chi^2 or error
            diff = xdot - y[row];
    
            //--- HARD CODED FOR 2D!! ---//
            // ->Future:  Try complex step derivatives for nD case
            //Fill entries of parameter Jacobian (dxdot/dbeta)
            paramJacobian(0,0) = -A(0,0) -Q[0](0,0)*dx[0] - 0.5*(Q[0](1,0) + Q[1](0,0))*dx[1];
            paramJacobian(0,1) = -A(0,1) -Q[1](1,0)*dx[1] - 0.5*(Q[0](1,0) + Q[1](0,0))*dx[0];
            paramJacobian(1,0) = -A(1,0) -Q[0](0,1)*dx[0] - 0.5*(Q[0](1,1) + Q[1](0,1))*dx[1];
            paramJacobian(1,1) = -A(1,1) -Q[1](1,1)*dx[1] - 0.5*(Q[0](1,1) + Q[1](0,1))*dx[0];
            for (int j=0; j<n; j++) {
                for (int i=0;i<n;i++) paramJacobian(j,j*n+i) = dx[i];
            }
            paramJacobian(0,6) = 0.5*dx[0]*dx[0];
            paramJacobian(1,7) = 0.5*dx[0]*dx[0];
    
            paramJacobian(0,8) = 0.5*dx[0]*dx[1];
            paramJacobian(1,9) = 0.5*dx[0]*dx[1];
            paramJacobian(0,10) = 0.5*dx[0]*dx[1];
            paramJacobian(1,11) = 0.5*dx[0]*dx[1];
    
            paramJacobian(0,12) = 0.5*dx[1]*dx[1];
            paramJacobian(1,13) = 0.5*dx[1]*dx[1];
            //--- End Hard coded section ---//
    
            //Fill out Jacobian row
            fjac.row(row) = 2.0*diff.transpose()*paramJacobian;
    
        }
        //success
        return 0;
    }*/
    
    ///Evaluate based on parameters
    VectorXd eval(const VectorXd& beta, const VectorXd& xVal)
    {
        assert(beta.size() == nParam);
        assert(xVal.size() == n);
        //beta is parameter vector, consisting of (xm,A.row(0),A.row(1),...)
        VectorXd xs(n), xdot(n), dx(n);
        MatrixXd A(n,n), partialQdx(n,n);
        //Q->Tensor or n nxn matrices
        MatrixXd Jx(n,n), Jy(n,n);
        partialQdx.setZero(n,n);
        Jx.setZero();
        Jy.setZero();
        xs = beta.head(n);
        for (int j=0; j<n; j++) {
            A.row(j) =  beta.segment(n+j*n          , n);
            Jx.row(j) = beta.segment(n+n*n +j*n     , n);
            Jy.row(j) = beta.segment(n+n*n +n*n +j*n, n);
        }
        //The Quad Model
        //For each value provided, compute fvec value (=chi^2)
        dx = xVal - xs;
        partialQdx = Jx*dx[0] + Jy*dx[1];
        xdot = A*dx + 0.5*partialQdx.transpose()*dx;
        return xdot;
    }
};


} //End xavier

#endif
