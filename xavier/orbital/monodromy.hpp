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


//Fixed-point classification for high-dim systems
// - Overrides some functionality in fixpoints.hpp
// - This assumes that you have the STM embedded in your RHS file.
//Author:  Wayne Schlei
//Date:  1/14/2013
//
//Mod: 10/30/2013 WRS- Modified classification of stable/unstable
//to use the trace(Monodromy) instead of numerical eigenvalues.
//Mod: 11/5/2016 WRS- Sometimes, eigenvalues don't come out as
//correct pairings and send out "nan" as eigenvectors. Trying fix:
//
#ifndef __MONODROMY_HPP__
#define __MONODROMY_HPP__

#include <iostream>
#include <vector>
#include <complex>

//Map Analysis
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/angle.hpp>
#include <maps/fixpoints.hpp>
#include <maps/mapNd.hpp>
#include <maps/poincare_map.hpp>
#include <maps/section.hpp>

//Eigen API
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;
using namespace xavier;


namespace orbital {

typedef nvis::fixed_matrix<double,6>         mat6;
typedef nvis::fixed_vector<double,6>         vec6;
typedef state_info<double, 6>        return_state;

/// Conversion function nvis::fixed_matrix<> to MatrixXd for Eigen
template<typename MAT>
MatrixXd mat2MatrixXd(const MAT& A)
{
    //Assumes we are transitioning a nvis::fixed_matrix<>
    const int rowDim = A.getNumRows();
    const int colDim = A.getNumCols();
    MatrixXd theMatrix(rowDim,colDim);
    for(int i=0; i<rowDim; i++)
        for(int j=0; j<colDim; j++) {
            theMatrix(i,j) = A[i][j];
        }
        
    return theMatrix;
};

/// MatrixXd to nvis::fixed_matrix<>
template<typename MAT>
MAT MatrixXd2Mat(const MatrixXd& theMatrix)
{
    MAT A(0);
    const int rowDim = A.getNumRows();
    const int colDim = A.getNumCols();
    for(int i=0; i<rowDim; i++)
        for(int j=0; j<colDim; j++) {
            A[i][j] = theMatrix(i,j);
        }
    return A;
}

/// Balancing a matrix as preparation for eigenvalue/vector computation (preserves eigenvalues and returns scales)
/// Apply scale output as eigvec[i] *= scale[i] or with balbak() function
template<class MAT>
VectorXd balanceInPlace(MAT& a)
{
    const double RADIX = std::numeric_limits<double>::radix;
    bool done = false;
    double sqrdx = RADIX*RADIX;
    //Assume STM is square
    const int n = a.cols();
    //Setup scale vector
    VectorXd scale = VectorXd::Ones(n);
    
    while (!done) {
        done = true;
        //Calculate row and column norms
        for (int i=0; i<n; i++) {
            double r=0.0,c=0.0;
            for(int j=0; j<n; j++) {
                if(j!=i) {
                    c += fabs( a(j,i) );
                    r += fabs( a(i,j) );
                }
            }
            //if both are nonzero,
            if (c!=0.0 && r!=0.0) {
                double g=r/RADIX;
                double f=1.0;
                double s=c+r;
                //Find the interger power of the machine radix that comes closest to balancing matrix
                while (c<g) {
                    f *= RADIX;
                    c *= sqrdx;
                }
                g = r*RADIX;
                while (c>g) {
                    f /= RADIX;
                    c /= sqrdx;
                }
                
                //Apply similarity transformation
                if((c+r)/f < 0.95*s) {
                    done = false;
                    g=1.0/f;
                    scale(i) *= f;
                    for(int j=0; j<n; j++) {
                        a(i,j) *= g;
                    }
                    for(int j=0; j<n; j++) {
                        a(j,i) *= f;
                    }
                }
            }
        }
    }
    return scale;
}

/// Apply the back transformation to the eigenvector matrix
template <class MAT>
MAT balbak(const VectorXd& scale, const MAT& evecs)
{
    const int n = evecs.cols();
    MatrixXcd zz = evecs;
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++) {
            zz(i,j) *= scale(i);
        }
    return zz;
}

/// Projecting an Eigen space vector onto a section (using the map)
template<typename MAP>
bool projectEigenVectorToSection(
    const MAP& theMap,
    const typename MAP::lvec_type& pos,       //Fixed point Map Coord
    const typename MAP::state_type& fullEvec, //Full Eigen Vector (normalized)
    typename MAP::lvec_type& evec,      //Section coord (output)
    const double perturbation,                //Amount to perturb from
    bool verbose=false)
{
    typedef typename MAP::section_type  Section;
    typedef typename MAP::lvec_type        lvec_type;
    typedef typename MAP::lmat_type        lmat_type;
    typedef typename MAP::state_type       StateType;
    typedef typename MAP::return_type      ReturnType;
    typedef typename MAP::rhs_type         RHSType;
    typedef typename Section::gvec_type    gvec_type;
    static const int N = RHSType::dimension;
    
    //The Section
    const Section& theSection = theMap.section();
    //The Fixed point
    StateType x0 = theSection.mapToState(pos);
    x0 += perturbation*fullEvec;
    gvec_type y0(0);
    for(int i=0; i<N; i++) {
        y0[i]=x0[i];
    }
    //Poincare map output
    std::vector<ReturnType> rInfo;
    int iter = 1;
    if (theSection.distance(y0) > 0) {
        iter = -1;
    }
    try {
        //Use StateType to ReturnType poincare_map function
        theMap.map(x0,rInfo,iter);
        evec = rInfo[0].x - pos;
        evec /= nvis::norm(evec);
    } catch(...) {
        ///Bad failure
        std::cerr << __FILE__ << ": ProjectEigenVectorToSection failure during mapping!\n";
        evec = lvec_type(0);
        return false;
    }
    return true;
}

/// Function to test if monodromy analysis is trustworthy (namely, the eigenvalues/vectors are usuable)
template<typename MAP>
void monodromyTrustTest(
    const MAP& theMap,
    const MatrixXd& stm,
    std::vector<fixpoint>& fps,
    bool verbose=false)
{
    bool trustResults = true;
    int period = fps[0].K;
    const double nanValue = std::numeric_limits<double>::quiet_NaN();
    vec_type nanVec(nanValue,nanValue);
    //Test determinant == 1
    if (fabs(stm.determinant() - 1.0) > 0.1) {
        trustResults = false;
    }
    
    //Test if we have NANs anywhere
    if (trustResults) {
        for(int i=0; i<period; i++) {
            if (ISNANCALL( fps[i].eval[0]) || ISNANCALL( fps[i].eval[1] ) ||
                    ISNANCALL( fps[i].evec[0][0] ) || ISNANCALL( fps[i].evec[0][1] ) ||
                    ISNANCALL( fps[i].evec[1][0] ) || ISNANCALL( fps[i].evec[1][1] ) ) {
                trustResults = false;
            }
        }
    }
    
    //Test we have reciprocal pairs in fps
    if (trustResults) {
        if (fabs( fps[0].eval[0]*fps[0].eval[1] - 1.0 ) > 0.1) {
            trustResults = false;
        }
    }
    
    //User prompt for debugging
    if(verbose) {
        std::cerr << " Trusted Monodromy Results = " << trustResults << "\n";
        if(!trustResults) {
            std::cerr << " Bad eigenvalues = [" << fps[0].eval[0] << ", " << fps[0].eval[1] << "]\n";
        }
    }
    
    //If we don't trust the results, lets reset values to indicate later
    if (!trustResults) {
        //Reset the minimum eigenvalue (as we can typically expect the largest to be correct)
        double lamMin = 1.0 / fps[0].eval[1];
        for(int i=0; i<period; i++) {
            fps[i].eval[0] = lamMin;
            //Set eigenvectors to nan so we can pickup faults later.
            fps[i].evec[0] = nanVec;
            fps[i].evec[1] = nanVec;
        }
        if(verbose) std::cerr << " Switch to eigenvalues = [" << fps[0].eval[0] << ", " << fps[0].eval[1]
                                  << "] and evecs = " << nanVec << "\n";
    }
}

/// Compute eigenVals/Vecs of Monodromy Matrix with Eigen Package (Return on-Map evecs)
/// Note: this version mixes analytical and numerical results, meaning results may not be the best
template<typename MAP>
bool computeEigen(const MAP& theMap, MatrixXd& numericalSTM,
                  std::vector<return_state>& stateInfo,
                  std::vector<fixpoint>& fps, bool verbose=false)
{
    typedef typename MAP::section_type  Section;
    typedef typename MAP::lvec_type  lvec_type;
    typedef typename MAP::lmat_type  lmat_type;
    typedef typename MAP::state_type StateType;
    typedef typename Section::gvec_type gvec_type;
    //const int N = Section::global_dimension;
    
    //The Section
    const Section& theSection = theMap.section();
    
    //Monodromy matrix is input numerical STM
    int period = (int) stateInfo.size(); //Total period = num elements
    
    double dp = 1.e-9; //Make a param?
    
    //Fixed Point Structure
    fixpoint fp = fps[0];
    std::pair<lvec_type,lmat_type> projection;
    //projection = theSection.project(stateInfo[0].getState());
    //fp.pos = projection.first;
    fp.K = period;
    fp.saddle = false; //Assume it isn't a saddle or it didn't work
    //Common Fixed point parameters
    fp.timePeriod = stateInfo.back().t;
    fp.t = 0.0;
    
    //Solve for Eigen vals/vecs (complex)
    ComplexEigenSolver<MatrixXd> ces(numericalSTM);
    bool didItWork = (ces.info() == Success ? true : false);
    if(verbose) {
        std::cerr << " Eigen Value Decomposition : ";
        if (didItWork) {
            std::cerr << "SUCCESS\n";
            std::cerr << "  Eigen Values :\n " << ces.eigenvalues() << "\n";
            std::cerr << "  Eigen Vectors:\n " << ces.eigenvectors() << "\n";
        } else {
            std::cerr << "FAILED\n";
        }
    }
    if (!didItWork) {
        fp.evec[0] = nvis::vec2(0,0);
        fp.evec[1] = nvis::vec2(0,0);
        fp.eval[0] = fp.eval[1] = 0.0;//0->Invalid, so it didn't work
        return false;
    }
    
    //Pair into recipricals (or conjugates if complex)
    std::map<int,int> valueMap; //Index to index
    double maxVal = 0.0;
    int maxID = 0;
    for (int i=0; i<6; i++) {
        for (int j=0; j<6; j++) {
            if (j!=i) {
                //Test if reciprical pair
                std::complex<double> temp = ces.eigenvalues()[i]*ces.eigenvalues()[j];
                double tempReal = std::real(temp);
                double tempImag = std::imag(temp);
                if ( (std::fabs( std::fabs(tempReal) - 1.0) < 1.e-8) && (std::fabs(tempImag) < 1.e-8) ) {
                    valueMap[i] = j;
                }
                //Conjugate? - Conjugate & reciprocal are the same if complex number is on unit circle
            }
        }
        
        //Find the stable/unstable pair
        // - Pick largest value pair for now (What do we do when there are 2 unstable pairs?)
        double realVal = std::real(ces.eigenvalues()[i]);
        if (abs(realVal) > maxVal) {
            realVal = maxVal;
            maxID = i;
        }
    }
    
    //Trace of the monodromy matrix for Stability information ------------------------------------
    double traceValue = numericalSTM.trace();
    // trace(numericalSTM) = 1+1 + stabilityIndex1 + stabilityIndex2; //SI1 is planar(xy) in CR3BP
    // HARD CODED:  The following applies only to the CR3B Problem
    MatrixXd planarMonodromy = MatrixXd::Zero(4,4);
    planarMonodromy << numericalSTM.block<2,2>(0,0),    numericalSTM.block<2,2>(0,3),
                    numericalSTM.block<2,2>(3,0),    numericalSTM.block<2,2>(3,3);
    double planarTraceValue = planarMonodromy.trace();
    double planarSI = planarTraceValue/2.0 - 1.0;
    double otherSI = (traceValue - planarTraceValue)/2.0;
    bool planarStability = (fabs(planarSI) > 1.0) ? false : true;
    //bool otherStability = (fabs(otherSI) >1.0) ? false : true;
    // Future Work:  For general dynamical system, we might need some rhs_coupling mechanism to
    // detach stability information in uncoupled systems. Maybe a class that will outline
    // which dimensions are coupled or not...
    // -------------------------------------------------------------------------------------------
    
    //Classify as stable or unstable fixed point
    //if ( fabs(std::real(ces.eigenvalues()[maxID])) >= 1.001) {
    if ( !planarStability ) {
        //Saddle - use stable/unstable pair
        fp.eval[0] = std::real( ces.eigenvalues()[valueMap[maxID]] ); //Stable
        fp.eval[1] = std::real( ces.eigenvalues()[maxID] ); //Unstable
        
        //Format Eigen vector
        gvec_type uvec(0.0), svec(0.0);
        for (int i=0; i<6; i++) {
            fp.fullEvec[0][i] = std::real( ces.eigenvectors().col(valueMap[maxID])[i] );
            svec[i] = std::real( ces.eigenvectors().col(valueMap[maxID])[i] );
        }
        for (int i=0; i<6; i++) {
            fp.fullEvec[1][i] = std::real( ces.eigenvectors().col(maxID)[i] );
            uvec[i] = std::real( ces.eigenvectors().col(maxID)[i] );
        }
        //Standardize directions
        if (fp.fullEvec[0][0]<0.0) {
            fp.fullEvec[0] *= -1.0;
        }
        if (fp.fullEvec[1][0]<0.0) {
            fp.fullEvec[1] *= -1.0;
        }
        
        /*//Straight Projection from state space onto section
        projection = theSection.project(svec);
        fp.evec[0] = projection.first; //Stable
        fp.evec[0] = (fp.evec[0][0]<0.0?-1:1)*fp.evec[0]/nvis::norm(fp.evec[0]); //Return normalized (and v[0]>0)
        projection = theSection.project(uvec);
        fp.evec[1] = projection.first; //Unstable
        fp.evec[1] = (fp.evec[1][0]<0.0?-1:1)*fp.evec[1]/nvis::norm(fp.evec[1]);*/
        
        //Project via mapping (takes a step off section and energy but in eigenvector direction)
        //Compute Stable direction with mapping
        bool sok = projectEigenVectorToSection(theMap,fp.pos,fp.fullEvec[0],fp.evec[0],dp);
        //Compute Unstable direction with mapping
        bool uok = projectEigenVectorToSection(theMap,fp.pos,fp.fullEvec[1],fp.evec[1],dp);
        if (!sok || !uok) {
            didItWork = false;
        }
        
        //Indicate as Saddle
        fp.saddle = true;
    } else {
        //Unable to determine if saddle due to numerical precision
        // - Likely a center (or extremely weak saddle)
        fp.evec[0] = nvis::vec2(0,0);
        fp.evec[1] = nvis::vec2(0,0);
        fp.eval[0] = fp.eval[1] = 1.0; //Center Eigen values are 1.0
    }
    fp.si = planarSI;
    fp.otherSI = otherSI;
    //Store the first fixed point in the chain
    fps.clear();
    fps.reserve(period);
    fps.push_back(fp);
    
    //The subsequent returns
    for (int i=0; i<(period-1); i++) {
        fixpoint fp_i(fp); //Copy
        //Set the new point
        projection = theSection.project(stateInfo[i].getState());
        fp_i.pos = projection.first;
        //Time
        fp_i.t = stateInfo[i].t;
        if ( fp.saddle ) {
            //STM transformation of stable/unstable eigenvectors
            MatrixXd stm = mat2MatrixXd<mat6>(stateInfo[i].J);
            VectorXcd sT = stm * ces.eigenvectors().col(valueMap[maxID]);
            VectorXcd uT = stm * ces.eigenvectors().col(maxID);
            gvec_type uvec(0.0), svec(0.0);
            for (int j=0; j<6; j++) {
                svec[j] = std::real( sT[j] );
                fp_i.fullEvec[0][j] = std::real( sT[j] );
                uvec[j] = std::real( uT[j] );
                fp_i.fullEvec[1][j] = std::real( uT[j] );
            }
            //Standardize directions
            if (fp_i.fullEvec[0][0]<0.0) {
                fp_i.fullEvec[0] *= -1.0;
            }
            if (fp_i.fullEvec[1][0]<0.0) {
                fp_i.fullEvec[1] *= -1.0;
            }
            /*//Set the fix point's eigenvectors
            projection = theSection.project(svec);
            fp_i.evec[0] = projection.first;
            fp_i.evec[0] = (fp_i.evec[0][0]<0.0?-1:1)*fp_i.evec[0]/nvis::norm(fp_i.evec[0]);
            projection = theSection.project(uvec);
            fp_i.evec[1] = projection.first;
            fp_i.evec[1] = (fp_i.evec[1][0]<0.0?-1:1)*fp_i.evec[1]/nvis::norm(fp_i.evec[1]);*/
            
            //Compute Stable direction with mapping
            bool sok = projectEigenVectorToSection(theMap,fp_i.pos,fp_i.fullEvec[0],fp_i.evec[0],dp);
            //Compute Unstable direction with mapping
            bool uok = projectEigenVectorToSection(theMap,fp_i.pos,fp_i.fullEvec[1],fp_i.evec[1],dp);
            if (!sok || !uok) {
                didItWork = false;
            }
        }
        //Store
        fps.push_back(fp_i);
    }
    
    return didItWork;
}

/// Compute eigenVals/Vecs of Monodromy Matrix with Eigen Package (Return on-Map evecs)
template<typename MAP>
bool computeEigen(const MAP& theMap, std::vector<return_state>& stateInfo,
                  std::vector<fixpoint>& fps, bool verbose=false)
{
    typedef typename MAP::section_type  Section;
    typedef typename MAP::lvec_type  lvec_type;
    typedef typename MAP::lmat_type  lmat_type;
    typedef typename MAP::state_type StateType;
    typedef typename Section::gvec_type gvec_type;
    //const int N = Section::global_dimension;
    
    //The Section
    const Section& theSection = theMap.section();
    
    //Get Monodromy matrix
    mat6 M = stateInfo.back().J; //Last fp's STM - analytical linear
    int period = (int) stateInfo.size(); //Total period = num elements
    double dp = 1.e-9; //Make a param?
    
    //mat6 -> MatrixXd
    MatrixXd monodromy = mat2MatrixXd<mat6>(M);
    
    //Fixed Point Structure
    fixpoint fp = fps[0];
    std::pair<lvec_type,lmat_type> projection;
    //projection = theSection.project(stateInfo[0].getState());
    //fp.pos = projection.first;
    fp.K = period;
    fp.saddle = false; //Assume it isn't a saddle or it didn't work
    //Common Fixed point parameters
    fp.timePeriod = stateInfo.back().t;
    fp.t = 0.0;
    
    //Solve for Eigen vals/vecs (complex)
    ComplexEigenSolver<MatrixXd> ces(monodromy);
    bool didItWork = (ces.info() == Success ? true : false);
    if(verbose) {
        std::cerr << " Eigen Value Decomposition : ";
        if (didItWork) {
            std::cerr << "SUCCESS\n";
            std::cerr << "  Eigen Values :\n " << ces.eigenvalues() << "\n";
            std::cerr << "  Eigen Vectors:\n " << ces.eigenvectors() << "\n";
        } else {
            std::cerr << "FAILED\n";
        }
    }
    if (!didItWork) {
        fp.evec[0] = nvis::vec2(0,0);
        fp.evec[1] = nvis::vec2(0,0);
        fp.eval[0] = fp.eval[1] = 0.0;//0->Invalid, so it didn't work
        return false;
    }
    
    //Pair into recipricals (or conjugates if complex)
    std::map<int,int> valueMap; //Index to index
    double maxVal = 0.0;
    int maxID = 0;
    for (int i=0; i<6; i++) {
        for (int j=0; j<6; j++) {
            if (j!=i) {
                //Test if reciprical pair
                std::complex<double> temp = ces.eigenvalues()[i]*ces.eigenvalues()[j];
                double tempReal = std::real(temp);
                double tempImag = std::imag(temp);
                if ( (std::fabs( std::fabs(tempReal) - 1.0) < 1.e-8) && (std::fabs(tempImag) < 1.e-8) ) {
                    valueMap[i] = j;
                }
                //Conjugate? - Conjugate & reciprocal are the same if complex number is on unit circle
            }
        }
        
        //Find the stable/unstable pair
        // - Pick largest value pair for now (What do we do when there are 2 unstable pairs?)
        double realVal = std::real(ces.eigenvalues()[i]);
        if (abs(realVal) > maxVal) {
            realVal = maxVal;
            maxID = i;
        }
    }
    
    //Trace of the monodromy matrix for Stability information ------------------------------------
    double traceValue = monodromy.trace();
    // trace(monodromy) = 1+1 + stabilityIndex1 + stabilityIndex2; //SI1 is planar(xy) in CR3BP
    // HARD CODED:  The following applies only to the CR3B Problem
    MatrixXd planarMonodromy = MatrixXd::Zero(4,4);
    planarMonodromy << monodromy.block<2,2>(0,0),    monodromy.block<2,2>(0,3),
                    monodromy.block<2,2>(3,0),    monodromy.block<2,2>(3,3);
    double planarTraceValue = planarMonodromy.trace();
    double planarSI = planarTraceValue/2.0 - 1.0;
    double otherSI = (traceValue - planarTraceValue)/2.0;
    bool planarStability = (fabs(planarSI) > 1.0) ? false : true;
    //bool otherStability = (fabs(otherSI) >1.0) ? false : true;
    // Future Work:  For general dynamical system, we might need some rhs_coupling mechanism to
    // detach stability information in uncoupled systems. Maybe a class that will outline
    // which dimensions are coupled or not...
    // -------------------------------------------------------------------------------------------
    
    //Classify as stable or unstable fixed point
    //if ( fabs(std::real(ces.eigenvalues()[maxID])) >= 1.001) {
    if ( !planarStability ) {
        //Saddle - use stable/unstable pair
        fp.eval[0] = std::real( ces.eigenvalues()[valueMap[maxID]] ); //Stable
        fp.eval[1] = std::real( ces.eigenvalues()[maxID] ); //Unstable
        
        //Format Eigen vector
        gvec_type uvec(0.0), svec(0.0);
        for (int i=0; i<6; i++) {
            fp.fullEvec[0][i] = std::real( ces.eigenvectors().col(valueMap[maxID])[i] );
            svec[i] = std::real( ces.eigenvectors().col(valueMap[maxID])[i] );
        }
        for (int i=0; i<6; i++) {
            fp.fullEvec[1][i] = std::real( ces.eigenvectors().col(maxID)[i] );
            uvec[i] = std::real( ces.eigenvectors().col(maxID)[i] );
        }
        //Standardize directions
        if (fp.fullEvec[0][0]<0.0) {
            fp.fullEvec[0] *= -1.0;
        }
        if (fp.fullEvec[1][0]<0.0) {
            fp.fullEvec[1] *= -1.0;
        }
        
        /*//Project from state space onto section
        projection = theSection.project(svec);
        fp.evec[0] = projection.first; //Stable
        fp.evec[0] = (fp.evec[0][0]<0.0?-1:1)*fp.evec[0]/nvis::norm(fp.evec[0]); //Return normalized (and v[0]>0)
        projection = theSection.project(uvec);
        fp.evec[1] = projection.first; //Unstable
        fp.evec[1] = (fp.evec[1][0]<0.0?-1:1)*fp.evec[1]/nvis::norm(fp.evec[1]);*/
        
        //Project via mapping (takes a step off section and energy but in eigenvector direction)
        //Compute Stable direction with mapping
        bool sok = projectEigenVectorToSection(theMap,fp.pos,fp.fullEvec[0],fp.evec[0],dp);
        //Compute Unstable direction with mapping
        bool uok = projectEigenVectorToSection(theMap,fp.pos,fp.fullEvec[1],fp.evec[1],dp);
        if (!sok || !uok) {
            didItWork = false;
        }
        fp.saddle = true;
    } else {
        //Unable to determine if saddle due to numerical precision
        // - Likely a center (or extremely weak saddle)
        fp.evec[0] = nvis::vec2(0,0);
        fp.evec[1] = nvis::vec2(0,0);
        fp.eval[0] = fp.eval[1] = 1.0; //Center Eigen values are 1.0
    }
    fp.si = planarSI;
    fp.otherSI = otherSI;
    //Store the first fixed point in the chain
    fps.clear();
    fps.reserve(period);
    fps.push_back(fp);
    
    //The subsequent returns
    for (int i=0; i<(period-1); i++) {
        fixpoint fp_i(fp); //Copy
        //Set the new point
        projection = theSection.project(stateInfo[i].getState());
        fp_i.pos = projection.first;
        //Time
        fp_i.t = stateInfo[i].t;
        if ( fp.saddle ) {
            //STM transformation of stable/unstable eigenvectors
            MatrixXd stm = mat2MatrixXd<mat6>(stateInfo[i].J);
            VectorXcd sT = stm * ces.eigenvectors().col(valueMap[maxID]);
            VectorXcd uT = stm * ces.eigenvectors().col(maxID);
            gvec_type uvec(0.0), svec(0.0);
            for (int j=0; j<6; j++) {
                svec[j] = std::real( sT[j] );
                fp_i.fullEvec[0][j] = std::real( sT[j] );
                uvec[j] = std::real( uT[j] );
                fp_i.fullEvec[1][j] = std::real( uT[j] );
            }
            //Standardize directions
            if (fp_i.fullEvec[0][0]<0.0) {
                fp_i.fullEvec[0] *= -1.0;
            }
            if (fp_i.fullEvec[1][0]<0.0) {
                fp_i.fullEvec[1] *= -1.0;
            }
            /*//Set the fix point's eigenvectors
            projection = theSection.project(svec);
            fp_i.evec[0] = projection.first;
            fp_i.evec[0] = (fp_i.evec[0][0]<0.0?-1:1)*fp_i.evec[0]/nvis::norm(fp_i.evec[0]);
            projection = theSection.project(uvec);
            fp_i.evec[1] = projection.first;
            fp_i.evec[1] = (fp_i.evec[1][0]<0.0?-1:1)*fp_i.evec[1]/nvis::norm(fp_i.evec[1]);*/
            
            //Compute Stable direction with mapping
            bool sok = projectEigenVectorToSection(theMap,fp_i.pos,fp_i.fullEvec[0],fp_i.evec[0],dp);
            //Compute Unstable direction with mapping
            bool uok = projectEigenVectorToSection(theMap,fp_i.pos,fp_i.fullEvec[1],fp_i.evec[1],dp);
            if (!sok || !uok) {
                didItWork = false;
            }
        }
        //Store
        fps.push_back(fp_i);
    }
    
    return didItWork;
}



/// Compute eigenVals/Vecs of Monodromy Matrix with Eigen Package (Return on-Map evecs) [specifically for numerical]
template<typename MAP>
bool computeEigen(const MAP& theMap,
                  std::vector<return_state>& stateInfo,
                  std::vector<MatrixXd>& mats,
                  std::vector<fixpoint>& fps, bool verbose=false)
{
    typedef typename MAP::section_type  Section;
    typedef typename MAP::lvec_type  lvec_type;
    typedef typename MAP::lmat_type  lmat_type;
    typedef typename MAP::state_type StateType;
    typedef typename Section::gvec_type gvec_type;
    //const int N = Section::global_dimension;
    
    //The Section
    const Section& theSection = theMap.section();
    
    //Monodromy matrix is the last input numerical STM
    int period = (int) stateInfo.size(); //Total period = num elements
    MatrixXd numericalSTM = mats.back();
    double dp = 1.e-9; //Make a param?
    
    //Fixed Point Structure
    fixpoint fp = fps[0];
    std::pair<lvec_type,lmat_type> projection;
    //projection = theSection.project(stateInfo[0].getState());
    //fp.pos = projection.first;
    fp.K = period;
    fp.saddle = false; //Assume it isn't a saddle or it didn't work
    //Common Fixed point parameters
    fp.timePeriod = stateInfo.back().t;
    fp.t = 0.0;
    
    //BIG NOTE:  Some STMs and monodromy matrixes will experience problems with scaling
    //(particularly in Saturn-Titan, Saturn-Enceladus, and systems with low mu values (<1e-4)).
    //This is usually due to extremely elements within the matrix (>O(10^5)) with other elements
    //that are typically O(10).
    // ->Counter this problem by "balancing" the matrix, which is a process of similarity
    //   transformations that preserve the eigenvalues but scale the eigenvectors.  This
    //   will increase accuracy of the eigenvalue/vector extraction especially in these
    //   unsymmetric matrices where a small deviation in an element can have drastic
    //   changes in the solution.
    // ->Balancing is referenced in Numerical Recipes in Section 11.6.
    
    
    //Solve for Eigen vals/vecs (complex)
    ComplexEigenSolver<MatrixXd> ces(numericalSTM);
    bool didItWork = (ces.info() == Success ? true : false);
    if(verbose) {
        std::cerr << " Eigen Value Decomposition : ";
        if (didItWork) {
            std::cerr << "SUCCESS\n";
            std::cerr << "  Eigen Values :\n " << ces.eigenvalues() << "\n";
            std::cerr << "  Eigen Vectors:\n " << ces.eigenvectors() << "\n";
        } else {
            std::cerr << "FAILED\n";
        }
    }
    if (!didItWork) {
        fp.evec[0] = nvis::vec2(0,0);
        fp.evec[1] = nvis::vec2(0,0);
        fp.eval[0] = fp.eval[1] = 0.0;//0->Invalid, so it didn't work
        return false;
    }
    
    //Pair into recipricals (or conjugates if complex)
    std::map<int,int> valueMap; //Index to index
    double maxVal = 0.0;
    int maxID = 0;
    for (int i=0; i<6; i++) {
        for (int j=0; j<6; j++) {
            if (j!=i) {
                //Test if reciprical pair
                std::complex<double> temp = ces.eigenvalues()[i]*ces.eigenvalues()[j];
                double tempReal = std::real(temp);
                double tempImag = std::imag(temp);
                if ( (std::fabs( std::fabs(tempReal) - 1.0) < 1.e-8) && (std::fabs(tempImag) < 1.e-8) ) {
                    valueMap[i] = j;
                }
                //Conjugate? - Conjugate & reciprocal are the same if complex number is on unit circle
            }
        }
        
        //Find the stable/unstable pair
        // - Pick largest value pair for now (What do we do when there are 2 unstable pairs?)
        double realVal = std::real(ces.eigenvalues()[i]);
        if (abs(realVal) > maxVal) {
            realVal = maxVal;
            maxID = i;
        }
    }
    
    //Trace of the monodromy matrix for Stability information ------------------------------------
    double traceValue = numericalSTM.trace();
    // trace(numericalSTM) = 1+1 + stabilityIndex1 + stabilityIndex2; //SI1 is planar(xy) in CR3BP
    // HARD CODED:  The following applies only to the CR3B Problem
    MatrixXd planarMonodromy = MatrixXd::Zero(4,4);
    planarMonodromy << numericalSTM.block<2,2>(0,0),    numericalSTM.block<2,2>(0,3),
                    numericalSTM.block<2,2>(3,0),    numericalSTM.block<2,2>(3,3);
    double planarTraceValue = planarMonodromy.trace();
    double planarSI = planarTraceValue/2.0 - 1.0;
    double otherSI = (traceValue - planarTraceValue)/2.0;
    bool planarStability = (fabs(planarSI) > 1.0) ? false : true;
    //bool otherStability = (fabs(otherSI) >1.0) ? false : true;
    // Future Work:  For general dynamical system, we might need some rhs_coupling mechanism to
    // detach stability information in uncoupled systems. Maybe a class that will outline
    // which dimensions are coupled or not...
    // -------------------------------------------------------------------------------------------
    
    //Classify as stable or unstable fixed point
    //if ( fabs(std::real(ces.eigenvalues()[maxID])) >= 1.001) {
    if ( !planarStability ) {
        //Saddle - use stable/unstable pair
        fp.eval[0] = std::real( ces.eigenvalues()[valueMap[maxID]] ); //Stable
        fp.eval[1] = std::real( ces.eigenvalues()[maxID] ); //Unstable
        
        //Format Eigen vector
        gvec_type uvec(0.0), svec(0.0);
        for (int i=0; i<6; i++) {
            fp.fullEvec[0][i] = std::real( ces.eigenvectors().col(valueMap[maxID])[i] );
            svec[i] = std::real( ces.eigenvectors().col(valueMap[maxID])[i] );
        }
        for (int i=0; i<6; i++) {
            fp.fullEvec[1][i] = std::real( ces.eigenvectors().col(maxID)[i] );
            uvec[i] = std::real( ces.eigenvectors().col(maxID)[i] );
        }
        //Standardize directions
        if (fp.fullEvec[0][0]<0.0) {
            fp.fullEvec[0] *= -1.0;
        }
        if (fp.fullEvec[1][0]<0.0) {
            fp.fullEvec[1] *= -1.0;
        }
        
        /*//Project from state space onto section
        projection = theSection.project(svec);
        fp.evec[0] = projection.first; //Stable
        fp.evec[0] = (fp.evec[0][0]<0.0?-1:1)*fp.evec[0]/nvis::norm(fp.evec[0]); //Return normalized (and v[0]>0)
        projection = theSection.project(uvec);
        fp.evec[1] = projection.first; //Unstable
        fp.evec[1] = (fp.evec[1][0]<0.0?-1:1)*fp.evec[1]/nvis::norm(fp.evec[1]);*/
        
        //Project via mapping (takes a step off section and energy but in eigenvector direction)
        //Compute Stable direction with mapping
        bool sok = projectEigenVectorToSection(theMap,fp.pos,fp.fullEvec[0],fp.evec[0],dp);
        //Compute Unstable direction with mapping
        bool uok = projectEigenVectorToSection(theMap,fp.pos,fp.fullEvec[1],fp.evec[1],dp);
        if (!sok || !uok) {
            didItWork = false;
        }
        fp.saddle = true;
    } else {
        //Unable to determine if saddle due to numerical precision
        // - Likely a center (or extremely weak saddle)
        fp.evec[0] = nvis::vec2(0,0);
        fp.evec[1] = nvis::vec2(0,0);
        fp.eval[0] = fp.eval[1] = 1.0; //Center Eigen values are 1.0
    }
    fp.si = planarSI;
    fp.otherSI = otherSI;
    //Store the first fixed point in the chain
    fps.clear();
    fps.reserve(period);
    fps.push_back(fp);
    
    //The subsequent returns
    for (int i=0; i<(period-1); i++) {
        fixpoint fp_i(fp); //Copy
        //Set the new point
        projection = theSection.project(stateInfo[i].getState());
        fp_i.pos = projection.first;
        //Time
        fp_i.t = stateInfo[i].t;
        if ( fp.saddle ) {
            //STM transformation of stable/unstable eigenvectors
            //MatrixXd stm = mat2MatrixXd<mat6>(stateInfo[i].J);
            MatrixXd stm = mats[i];
            VectorXcd sT = stm * ces.eigenvectors().col(valueMap[maxID]);
            VectorXcd uT = stm * ces.eigenvectors().col(maxID);
            gvec_type uvec(0.0), svec(0.0);
            for (int j=0; j<6; j++) {
                svec[j] = std::real( sT[j] );
                fp_i.fullEvec[0][j] = std::real( sT[j] );
                uvec[j] = std::real( uT[j] );
                fp_i.fullEvec[1][j] = std::real( uT[j] );
            }
            //Standardize directions
            if (fp_i.fullEvec[0][0]<0.0) {
                fp_i.fullEvec[0] *= -1.0;
            }
            if (fp_i.fullEvec[1][0]<0.0) {
                fp_i.fullEvec[1] *= -1.0;
            }
            
            /*//Set the fix point's eigenvectors
            projection = theSection.project(svec);
            fp_i.evec[0] = projection.first;
            fp_i.evec[0] = (fp_i.evec[0][0]<0.0?-1:1)*fp_i.evec[0]/nvis::norm(fp_i.evec[0]);
            projection = theSection.project(uvec);
            fp_i.evec[1] = projection.first;
            fp_i.evec[1] = (fp_i.evec[1][0]<0.0?-1:1)*fp_i.evec[1]/nvis::norm(fp_i.evec[1]);*/
            
            //Compute Stable direction with mapping
            bool sok = projectEigenVectorToSection(theMap,fp_i.pos,fp_i.fullEvec[0],fp_i.evec[0],dp);
            //Compute Unstable direction with mapping
            bool uok = projectEigenVectorToSection(theMap,fp_i.pos,fp_i.fullEvec[1],fp_i.evec[1],dp);
            if (!sok || !uok) {
                didItWork = false;
            }
        }
        //Store
        fps.push_back(fp_i);
    }
    
    return didItWork;
}

/// Compute eigenVals/Vecs of Monodromy Matrix with Eigen Package (Return on-Map evecs) With BALANCING
template<typename MAP>
bool computeEigenWithBalancing(
    const MAP& theMap,
    std::vector<return_state>& stateInfo,
    std::vector<MatrixXd>& mats,
    std::vector<fixpoint>& fps, bool verbose=false)
{
    typedef typename MAP::section_type  Section;
    typedef typename MAP::lvec_type  lvec_type;
    typedef typename MAP::lmat_type  lmat_type;
    typedef typename MAP::state_type StateType;
    typedef typename Section::gvec_type gvec_type;
    //const int N = Section::global_dimension;
    
    //The Section
    const Section& theSection = theMap.section();
    
    //Monodromy matrix is the last input numerical STM
    int period = (int) stateInfo.size(); //Total period = num elements
    MatrixXd numericalSTM = mats.back();
    double dp = 1.e-9; //Make a param?
    
    //Fixed Point Structure
    fixpoint fp = fps[0];
    std::pair<lvec_type,lmat_type> projection;
    //projection = theSection.project(stateInfo[0].getState());
    //fp.pos = projection.first;
    fp.K = period;
    fp.saddle = false; //Assume it isn't a saddle or it didn't work
    //Common Fixed point parameters
    fp.timePeriod = stateInfo.back().t;
    fp.t = 0.0;
    
    //BIG NOTE:  Some STMs and monodromy matrixes will experience problems with scaling
    //(particularly in Saturn-Titan, Saturn-Enceladus, and systems with low mu values (<1e-4)).
    //This is usually due to extremely elements within the matrix (>O(10^5)) with other elements
    //that are typically O(10).
    // ->Counter this problem by "balancing" the matrix, which is a process of similarity
    //   transformations that preserve the eigenvalues but scale the eigenvectors.  This
    //   will increase accuracy of the eigenvalue/vector extraction especially in these
    //   unsymmetric matrices where a small deviation in an element can have drastic
    //   changes in the solution.
    // ->Balancing is referenced in Numerical Recipes in Section 11.6.
    
    // Apply balancing in place
    MatrixXd nSTM_Bal = numericalSTM;
    VectorXd scale = balanceInPlace(nSTM_Bal);
    
    //Debug output:
    if(verbose) {
        std::cerr << " Balanced STM Matrix:\n" << nSTM_Bal << "\n";
    }
    
    //Solve for Eigen vals/vecs (complex)
    ComplexEigenSolver<MatrixXd> ces(nSTM_Bal);
    bool didItWork = (ces.info() == Success ? true : false);
    //Rescale the eigenvectors to correct size
    MatrixXcd eVecsM = balbak( scale, ces.eigenvectors() );
    if(verbose) {
        std::cerr << " Eigen Value Decomposition (BALANCED matrix): ";
        if (didItWork) {
            std::cerr << "SUCCESS\n";
            std::cerr << "  Eigen Values :\n " << ces.eigenvalues() << "\n";
            std::cerr << "  Eigen Vectors (scaled):\n " << ces.eigenvectors() << "\n";
            std::cerr << "  Eigen Vectors :\n" << eVecsM << "\n";
        } else {
            std::cerr << "FAILED\n";
        }
    }
    if (!didItWork) {
        fp.evec[0] = nvis::vec2(0,0);
        fp.evec[1] = nvis::vec2(0,0);
        fp.eval[0] = fp.eval[1] = 0.0;//0->Invalid, so it didn't work
        return false;
    }
    
    //Pair into recipricals (or conjugates if complex)
    std::map<int,int> valueMap; //Index to index
    double maxVal = 0.0;
    int maxID = 0;
    for (int i=0; i<6; i++) {
        for (int j=0; j<6; j++) {
            if (j!=i) {
                //Test if reciprical pair
                std::complex<double> temp = ces.eigenvalues()[i]*ces.eigenvalues()[j];
                double tempReal = std::real(temp);
                double tempImag = std::imag(temp);
                if ( (std::fabs( std::fabs(tempReal) - 1.0) < 1.e-8) && (std::fabs(tempImag) < 1.e-8) ) {
                    valueMap[i] = j;
                }
                //Conjugate? - Conjugate & reciprocal are the same if complex number is on unit circle
            }
        }
        
        //Find the stable/unstable pair
        // - Pick largest value pair for now (What do we do when there are 2 unstable pairs?)
        double realVal = std::real(ces.eigenvalues()[i]);
        if (abs(realVal) > maxVal) {
            realVal = maxVal;
            maxID = i;
        }
    }
    
    //Trace of the monodromy matrix for Stability information ------------------------------------
    double traceValue = numericalSTM.trace();
    // trace(numericalSTM) = 1+1 + stabilityIndex1 + stabilityIndex2; //SI1 is planar(xy) in CR3BP
    // HARD CODED:  The following applies only to the CR3B Problem
    MatrixXd planarMonodromy = MatrixXd::Zero(4,4);
    planarMonodromy << numericalSTM.block<2,2>(0,0),    numericalSTM.block<2,2>(0,3),
                    numericalSTM.block<2,2>(3,0),    numericalSTM.block<2,2>(3,3);
    double planarTraceValue = planarMonodromy.trace();
    double planarSI = planarTraceValue/2.0 - 1.0;
    double otherSI = (traceValue - planarTraceValue)/2.0;
    bool planarStability = (fabs(planarSI) > 1.0) ? false : true;
    //bool otherStability = (fabs(otherSI) >1.0) ? false : true;
    // Future Work:  For general dynamical system, we might need some rhs_coupling mechanism to
    // detach stability information in uncoupled systems. Maybe a class that will outline
    // which dimensions are coupled or not...
    // -------------------------------------------------------------------------------------------
    
    //Classify as stable or unstable fixed point
    //if ( fabs(std::real(ces.eigenvalues()[maxID])) >= 1.001) {
    if ( !planarStability ) {
        //Saddle - use stable/unstable pair
        fp.eval[0] = std::real( ces.eigenvalues()[valueMap[maxID]] ); //Stable
        fp.eval[1] = std::real( ces.eigenvalues()[maxID] ); //Unstable
        
        //Format Eigen vector
        gvec_type uvec(0.0), svec(0.0);
        for (int i=0; i<6; i++) {
            fp.fullEvec[0][i] = std::real( eVecsM.col(valueMap[maxID])[i] );
            svec[i] = std::real( eVecsM.col(valueMap[maxID])[i] );
        }
        for (int i=0; i<6; i++) {
            fp.fullEvec[1][i] = std::real( eVecsM.col(maxID)[i] );
            uvec[i] = std::real( eVecsM.col(maxID)[i] );
        }
        //Standardize directions
        if (fp.fullEvec[0][0]<0.0) {
            fp.fullEvec[0] *= -1.0;
        }
        if (fp.fullEvec[1][0]<0.0) {
            fp.fullEvec[1] *= -1.0;
        }
        
        /*//Project from state space onto section
        projection = theSection.project(svec);
        fp.evec[0] = projection.first; //Stable
        fp.evec[0] = (fp.evec[0][0]<0.0?-1:1)*fp.evec[0]/nvis::norm(fp.evec[0]); //Return normalized (and v[0]>0)
        projection = theSection.project(uvec);
        fp.evec[1] = projection.first; //Unstable
        fp.evec[1] = (fp.evec[1][0]<0.0?-1:1)*fp.evec[1]/nvis::norm(fp.evec[1]);*/
        
        //Project via mapping (takes a step off section and energy but in eigenvector direction)
        //Compute Stable direction with mapping
        bool sok = projectEigenVectorToSection(theMap,fp.pos,fp.fullEvec[0],fp.evec[0],dp);
        //Compute Unstable direction with mapping
        bool uok = projectEigenVectorToSection(theMap,fp.pos,fp.fullEvec[1],fp.evec[1],dp);
        if (!sok || !uok) {
            didItWork = false;
        }
        fp.saddle = true;
    } else {
        //Unable to determine if saddle due to numerical precision
        // - Likely a center (or extremely weak saddle)
        fp.evec[0] = nvis::vec2(0,0);
        fp.evec[1] = nvis::vec2(0,0);
        fp.eval[0] = fp.eval[1] = 1.0; //Center Eigen values are 1.0
    }
    fp.si = planarSI;
    fp.otherSI = otherSI;
    //Store the first fixed point in the chain
    fps.clear();
    fps.reserve(period);
    fps.push_back(fp);
    
    //The subsequent returns
    for (int i=0; i<(period-1); i++) {
        fixpoint fp_i(fp); //Copy
        //Set the new point
        projection = theSection.project(stateInfo[i].getState());
        fp_i.pos = projection.first;
        //Time
        fp_i.t = stateInfo[i].t;
        if ( fp.saddle ) {
            //STM transformation of stable/unstable eigenvectors
            //MatrixXd stm = mat2MatrixXd<mat6>(stateInfo[i].J);
            MatrixXd stm = mats[i];
            VectorXcd sT = stm * eVecsM.col(valueMap[maxID]);
            VectorXcd uT = stm * eVecsM.col(maxID);
            gvec_type uvec(0.0), svec(0.0);
            for (int j=0; j<6; j++) {
                svec[j] = std::real( sT[j] );
                fp_i.fullEvec[0][j] = std::real( sT[j] );
                uvec[j] = std::real( uT[j] );
                fp_i.fullEvec[1][j] = std::real( uT[j] );
            }
            //Standardize directions
            if (fp_i.fullEvec[0][0]<0.0) {
                fp_i.fullEvec[0] *= -1.0;
            }
            if (fp_i.fullEvec[1][0]<0.0) {
                fp_i.fullEvec[1] *= -1.0;
            }
            
            /*//Set the fix point's eigenvectors
            projection = theSection.project(svec);
            fp_i.evec[0] = projection.first;
            fp_i.evec[0] = (fp_i.evec[0][0]<0.0?-1:1)*fp_i.evec[0]/nvis::norm(fp_i.evec[0]);
            projection = theSection.project(uvec);
            fp_i.evec[1] = projection.first;
            fp_i.evec[1] = (fp_i.evec[1][0]<0.0?-1:1)*fp_i.evec[1]/nvis::norm(fp_i.evec[1]);*/
            
            //Compute Stable direction with mapping
            bool sok = projectEigenVectorToSection(theMap,fp_i.pos,fp_i.fullEvec[0],fp_i.evec[0],dp);
            //Compute Unstable direction with mapping
            bool uok = projectEigenVectorToSection(theMap,fp_i.pos,fp_i.fullEvec[1],fp_i.evec[1],dp);
            if (!sok || !uok) {
                didItWork = false;
            }
        }
        //Store
        fps.push_back(fp_i);
    }
    
    return didItWork;
}

/// Compute eigenVals/Vecs of PLANAR Monodromy Matrix with Eigen Package (Return on-Map evecs) [specifically for numerical]
template<typename MAP>
bool computeEigenPlanar(
    const MAP& theMap,
    std::vector<return_state>& stateInfo,
    std::vector<MatrixXd>& mats,
    std::vector<fixpoint>& fps, bool verbose=false)
{
    typedef typename MAP::section_type  Section;
    typedef typename MAP::lvec_type  lvec_type;
    typedef typename MAP::lmat_type  lmat_type;
    typedef typename MAP::state_type StateType;
    typedef typename Section::gvec_type gvec_type;
    //const int N = Section::global_dimension;
    
    //The Section
    const Section& theSection = theMap.section();
    
    //Monodromy matrix is the last input numerical STM
    int period = (int) stateInfo.size(); //Total period = num elements
    MatrixXd numericalSTM = mats.back();
    double dp = 1.e-9; //Make a param?
    
    //Fixed Point Structure
    fixpoint fp = fps[0];
    std::pair<lvec_type,lmat_type> projection;
    //projection = theSection.project(stateInfo[0].getState());
    //fp.pos = projection.first;
    fp.K = period;
    fp.saddle = false; //Assume it isn't a saddle or it didn't work
    //Common Fixed point parameters
    fp.timePeriod = stateInfo.back().t;
    fp.t = 0.0;
    
    //BIG NOTE:  Some STMs and monodromy matrixes will experience problems with scaling
    //(particularly in Saturn-Titan, Saturn-Enceladus, and systems with low mu values (<1e-4)).
    //This is usually due to extremely elements within the matrix (>O(10^5)) with other elements
    //that are typically O(10).
    // ->Counter this problem by "balancing" the matrix, which is a process of similarity
    //   transformations that preserve the eigenvalues but scale the eigenvectors.  This
    //   will increase accuracy of the eigenvalue/vector extraction especially in these
    //   unsymmetric matrices where a small deviation in an element can have drastic
    //   changes in the solution.
    // ->Balancing is referenced in Numerical Recipes in Section 11.6.
    
    //Trace of the monodromy matrix for Stability information ------------------------------------
    double traceValue = numericalSTM.trace();
    // trace(numericalSTM) = 1+1 + stabilityIndex1 + stabilityIndex2; //SI1 is planar(xy) in CR3BP
    // HARD CODED:  The following applies only to the CR3B Problem
    MatrixXd planarMonodromy = MatrixXd::Zero(4,4);
    planarMonodromy << numericalSTM.block<2,2>(0,0),    numericalSTM.block<2,2>(0,3),
                    numericalSTM.block<2,2>(3,0),    numericalSTM.block<2,2>(3,3);
    double planarTraceValue = planarMonodromy.trace();
    double planarSI = planarTraceValue/2.0 - 1.0;
    double otherSI = (traceValue - planarTraceValue)/2.0;
    bool planarStability = (fabs(planarSI) > 1.0) ? false : true;
    //bool otherStability = (fabs(otherSI) >1.0) ? false : true;
    // Future Work:  For general dynamical system, we might need some rhs_coupling mechanism to
    // detach stability information in uncoupled systems. Maybe a class that will outline
    // which dimensions are coupled or not...
    // -------------------------------------------------------------------------------------------
    
    //Solve for Eigen vals/vecs (complex)
    ComplexEigenSolver<MatrixXd> ces(planarMonodromy);
    bool didItWork = (ces.info() == Success ? true : false);
    if(verbose) {
        std::cerr << " Eigen Value Decomposition (PLANAR) : ";
        if (didItWork) {
            std::cerr << "SUCCESS\n";
            std::cerr << "  Eigen Values :\n " << ces.eigenvalues() << "\n";
            std::cerr << "  Eigen Vectors:\n " << ces.eigenvectors() << "\n";
        } else {
            std::cerr << "FAILED\n";
        }
    }
    if (!didItWork) {
        fp.evec[0] = nvis::vec2(0,0);
        fp.evec[1] = nvis::vec2(0,0);
        fp.eval[0] = fp.eval[1] = 0.0;//0->Invalid, so it didn't work
        return false;
    }
    
    //Pair into recipricals (or conjugates if complex)
    std::map<int,int> valueMap; //Index to index
    double maxVal = 0.0;
    int maxID = 0;
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            if (j!=i) {
                //Test if reciprical pair
                std::complex<double> temp = ces.eigenvalues()[i]*ces.eigenvalues()[j];
                double tempReal = std::real(temp);
                double tempImag = std::imag(temp);
                if ( (std::fabs( std::fabs(tempReal) - 1.0) < 1.e-8) && (std::fabs(tempImag) < 1.e-8) ) {
                    valueMap[i] = j;
                }
                //Conjugate? - Conjugate & reciprocal are the same if complex number is on unit circle
            }
        }
        
        //Find the stable/unstable pair
        // - Pick largest value pair for now (What do we do when there are 2 unstable pairs?)
        double realVal = std::real(ces.eigenvalues()[i]);
        if (abs(realVal) > maxVal) {
            realVal = maxVal;
            maxID = i;
        }
    }
    
    
    //Classify as stable or unstable fixed point
    //if ( fabs(std::real(ces.eigenvalues()[maxID])) >= 1.001) {
    if ( !planarStability ) {
        //Saddle - use stable/unstable pair
        fp.eval[0] = std::real( ces.eigenvalues()[valueMap[maxID]] ); //Stable
        fp.eval[1] = std::real( ces.eigenvalues()[maxID] ); //Unstable
        
        //Format Eigen vector
        gvec_type uvec(0.0), svec(0.0);
        for (int i=0; i<4; i++) {
            int idx = (i>1)? i+1 : i;
            fp.fullEvec[0][idx] = std::real( ces.eigenvectors().col(valueMap[maxID])[i] );
            svec[i] = std::real( ces.eigenvectors().col(valueMap[maxID])[i] );
        }
        for (int i=0; i<4; i++) {
            int idx = (i>1)? i+1 : i;
            fp.fullEvec[1][idx] = std::real( ces.eigenvectors().col(maxID)[i] );
            uvec[i] = std::real( ces.eigenvectors().col(maxID)[i] );
        }
        //Standardize directions
        if (fp.fullEvec[0][0]<0.0) {
            fp.fullEvec[0] *= -1.0;
        }
        if (fp.fullEvec[1][0]<0.0) {
            fp.fullEvec[1] *= -1.0;
        }
        
        /*//Project from state space onto section
        projection = theSection.project(svec);
        fp.evec[0] = projection.first; //Stable
        fp.evec[0] = (fp.evec[0][0]<0.0?-1:1)*fp.evec[0]/nvis::norm(fp.evec[0]); //Return normalized (and v[0]>0)
        projection = theSection.project(uvec);
        fp.evec[1] = projection.first; //Unstable
        fp.evec[1] = (fp.evec[1][0]<0.0?-1:1)*fp.evec[1]/nvis::norm(fp.evec[1]);*/
        
        //Project via mapping (takes a step off section and energy but in eigenvector direction)
        //Compute Stable direction with mapping
        bool sok = projectEigenVectorToSection(theMap,fp.pos,fp.fullEvec[0],fp.evec[0],dp);
        //Compute Unstable direction with mapping
        bool uok = projectEigenVectorToSection(theMap,fp.pos,fp.fullEvec[1],fp.evec[1],dp);
        if (!sok || !uok) {
            didItWork = false;
        }
        fp.saddle = true;
    } else {
        //Unable to determine if saddle due to numerical precision
        // - Likely a center (or extremely weak saddle)
        fp.evec[0] = nvis::vec2(0,0);
        fp.evec[1] = nvis::vec2(0,0);
        fp.eval[0] = fp.eval[1] = 1.0; //Center Eigen values are 1.0
    }
    fp.si = planarSI;
    fp.otherSI = otherSI;
    //Store the first fixed point in the chain
    fps.clear();
    fps.reserve(period);
    fps.push_back(fp);
    
    //The subsequent returns
    for (int i=0; i<(period-1); i++) {
        fixpoint fp_i(fp); //Copy
        //Set the new point
        projection = theSection.project(stateInfo[i].getState());
        fp_i.pos = projection.first;
        //Time
        fp_i.t = stateInfo[i].t;
        if ( fp.saddle ) {
            //STM transformation of stable/unstable eigenvectors
            //MatrixXd stm = mat2MatrixXd<mat6>(stateInfo[i].J);
            MatrixXd stm = mats[i];
            MatrixXd pSTM = MatrixXd::Zero(4,4);
            pSTM << stm.block<2,2>(0,0), stm.block<2,2>(0,3), stm.block<2,2>(3,0), stm.block<2,2>(3,3);
            VectorXcd sT = pSTM * ces.eigenvectors().col(valueMap[maxID]);
            VectorXcd uT = pSTM * ces.eigenvectors().col(maxID);
            gvec_type uvec(0.0), svec(0.0);
            for (int j=0; j<6; j++) {
                int jdx = (j>1)? j+1 : j;
                svec[jdx] = std::real( sT[j] );
                fp_i.fullEvec[0][jdx] = std::real( sT[j] );
                uvec[jdx] = std::real( uT[j] );
                fp_i.fullEvec[1][jdx] = std::real( uT[j] );
            }
            //Standardize directions
            if (fp_i.fullEvec[0][0]<0.0) {
                fp_i.fullEvec[0] *= -1.0;
            }
            if (fp_i.fullEvec[1][0]<0.0) {
                fp_i.fullEvec[1] *= -1.0;
            }
            
            /*//Set the fix point's eigenvectors
            projection = theSection.project(svec);
            fp_i.evec[0] = projection.first;
            fp_i.evec[0] = (fp_i.evec[0][0]<0.0?-1:1)*fp_i.evec[0]/nvis::norm(fp_i.evec[0]);
            projection = theSection.project(uvec);
            fp_i.evec[1] = projection.first;
            fp_i.evec[1] = (fp_i.evec[1][0]<0.0?-1:1)*fp_i.evec[1]/nvis::norm(fp_i.evec[1]);*/
            
            //Compute Stable direction with mapping
            bool sok = projectEigenVectorToSection(theMap,fp_i.pos,fp_i.fullEvec[0],fp_i.evec[0],dp);
            //Compute Unstable direction with mapping
            bool uok = projectEigenVectorToSection(theMap,fp_i.pos,fp_i.fullEvec[1],fp_i.evec[1],dp);
            if (!sok || !uok) {
                didItWork = false;
            }
        }
        //Store
        fps.push_back(fp_i);
    }
    
    return didItWork;
}

/// Classify stability characteristics of a fixed point (Output: chain -> fixed-points for a particular periodic orbit)
template<typename MAP>
bool monodromy_analysis(const MAP& theMap, const fixpoint& fp0, std::vector<fixpoint>& chain, bool linearSTM = true, bool verbose=false)
{
    typedef typename MAP::rhs_type      RHStype;
    typedef typename MAP::xstate_type   xstate_type;
    static const int n = RHStype::dimension;
    const typename MAP::section_type theSection = theMap.section();
    //Compute the monodromy matrix and data per return
    std::vector<return_state> stateDataPerReturn;
    nvis::vec2 x = fp0.pos;
    int period = fp0.K;
    std::vector<double> tCrossings;
    //Map for linearSTM all starting from x for full period (Really bad approximation)
    theMap.map_complete(x,stateDataPerReturn,period);
    
    //Lets check if fp0 passes a sanity check:   abs(fp0.pos - stateDataPerReturn[period-1].x) = 0.0
    nvis::vec2 test = fp0.pos - theSection.project( stateDataPerReturn.back().getState() ).first;
    if (nvis::norm(test) >= 1.e-4) {
        if (verbose) {
            std::cerr << " Fixed Point Sanity Check FAILED!\n";
            std::cerr << " Fixed Point : " << fp0.pos << " with p=" << fp0.K << "\n";
            std::cerr << "   Propagated to : " << stateDataPerReturn.back().x << "\n";
            std::cerr << "   with error : " << test << "\n";
        }
        return false;
    }
    for(int i=0; i<period; i++) {
        tCrossings.push_back( stateDataPerReturn[i].t );
    }
    MatrixXd numericalSTM = MatrixXd::Zero(n,n);
    std::vector<MatrixXd> stmPerReturn;
    //bool tryBalancing = false;
    //bool tryPlanar = false;
    bool revertToLinearSTM = false;
    
    if(linearSTM) {
        //Map:  Resets the STM at each crossing (Better Approx, but will be inconsistent)
        //theMap.map_STMperReturn(x,stateDataPerReturn,period); //Need to multiply STMs to get true value
        //Map: Runs Map accumulating STMs with maximum time spans (Best Linear Approx)
        //double dtmax = 2.0; //Param.linearSTM_dtMax?
        //theMap.map_LinearMonodromy(x,dtmax,tCrossings,stateDataPerReturn);
    } else {
        //Computing the numerical monodromy matrix - Technically need the numerical value for all subpoints too!
        //return_state finalInfo = stateDataPerReturn.back();
        //numericalSTM.resize(n,n);
        //numericalSTM.setZero();
        //Need a matrix for each fixed point except last
        stmPerReturn.assign(period,numericalSTM);
        //Propagation times
        std::vector<double> dt;
        for(int i=0; i<period; i++) {
            dt.push_back( tCrossings[i] - ((i==0)?0.0:tCrossings[i-1]) );
        }
        for(int i=0; i<n; i++) {
            //Perturb each state
            double h = 1e-10; //Perturbation
            //Central difference
            xstate_type yp = theSection.unproject(fp0.pos); //IC
            xstate_type ym = theSection.unproject(fp0.pos); //IC
            yp[i] += h;
            ym[i] -= h;
            return_state pInfo,mInfo;
            //Loop for each time step
            for (int k=0; k<period; k++) {
                try {
                    //contemporaneous variation (i.e., same times)
                    pInfo = theMap.integrate_state(yp,dt[k]);
                    mInfo = theMap.integrate_state(ym,dt[k]);
                } catch(...) {
                    //Unsaveable error - (shouldn't happen, BUT IT DOES!!!)
                    //throw std::runtime_error("Failed to propagate pertubation for Numerical Monodromy Matrix");
                    //In the event of this happening, we should recover linearSTM from propagation as a filler
                    revertToLinearSTM = true;
                    //OR use forward difference OR a modified new step???
                }
                //Fill out corresponding STM row with numerical values
                xstate_type yfp =  pInfo.getState();
                xstate_type yfm =  mInfo.getState();
                //phi  = dxf/dx0 = dxf/h -> column of numerical monodromy
                for(int j=0; j<n; j++) {
                    numericalSTM(j,i) = (yfp[j] - yfm[j])/(2.0*h);
                }
                for(int j=0; j<n; j++) {
                    stmPerReturn[k](j,i) = (yfp[j] - yfm[j])/(2.0*h);
                }
                //Update state (perturbed for next propagation)
                yp = yfp;
                ym = yfm;
            }
        }
        //Debug prompt:  Output the numericalSTM
        if(verbose) {
            std::cerr << "NumericalSTM: \n" << numericalSTM << "\n";
        }
        if(verbose) {
            std::cerr << "det(NumericalSTM) : " << numericalSTM.determinant() << "\n";
        }
        //std::cerr << "stmPerReturn.back(): \n" << stmPerReturn.back() << "\n";
    }
    
    //Unfortunately, map_complete() does not store the first point,
    //   so we save it as the first point in the chain
    chain.push_back(fp0); //This will be rewritten in computeEigen() with the rest of the information.
    
    //Observe the numerical STM behavior if activated
    bool ok = false;
    if (!linearSTM) {
        //Numerical Monodromy : try the determinant first
        if (fabs(numericalSTM.determinant() - 1.0) > 0.1) {
            ok = false; //Just revert back to linear STM
            if(verbose) {
                std::cerr << "NumericalSTM has poor matrix (det != 1). Reverting to linearSTM.\n";
            }
        } else {
            ok = true;
        }
        
        //Numerical STM's version
        if(ok) {
            ok = computeEigen<MAP>(theMap, stateDataPerReturn, stmPerReturn, chain, verbose);
        }
        //std::cerr << " Numerical Monodromy = \n" << stmPerReturn.back() << "\n\n";
        
        //Now, test for nan values in eigenvalues or eigenvectors
        // Note:  some numericalSTM create erroneous values that generate nans
        if(!ok) {
            revertToLinearSTM = true;
        } else {
            for(int i=0; i<period; i++) {
                if (ISNANCALL( chain[i].eval[0]) || ISNANCALL( chain[i].eval[1] ) ||
                        ISNANCALL( chain[i].evec[0][0] ) || ISNANCALL( chain[i].evec[0][1] ) ||
                        ISNANCALL( chain[i].evec[1][0] ) || ISNANCALL( chain[i].evec[1][1] ) ) {
                    //tryBalancing = true;
                    //tryPlanar = true;
                    revertToLinearSTM = true;
                }
            }
        }
        //Prompt user
        //if(verbose) std::cerr << "tryBalancing: " << tryBalancing << "\n";
        //if(verbose) std::cerr << "tryPlanar: " << tryPlanar << "\n";
        if(verbose) {
            std::cerr << "revertToLinearSTM: " << revertToLinearSTM << "\n";
        }
    }
    
    /*//Run a balancing operation to attempt to fix issues
    if(tryBalancing) {
      //Numerical STM's version with balancing applied
      ok = computeEigenWithBalancing<MAP>(theMap, stateDataPerReturn, stmPerReturn, chain, verbose);
      //Now, test for nan values in eigenvalues or eigenvectors
      // Note:  some numericalSTM create erroneous values that generate nans
      if(!ok) {revertToLinearSTM = true;}
      else {
        for(int i=0;i<period;i++) {
        if (ISNANCALL( chain[i].eval[0]) || ISNANCALL( chain[i].eval[1] ) ||
            ISNANCALL( chain[i].evec[0][0] ) || ISNANCALL( chain[i].evec[0][1] ) ||
            ISNANCALL( chain[i].evec[1][0] ) || ISNANCALL( chain[i].evec[1][1] ) ) {
          revertToLinearSTM = true;
        }
        }
      }
      //Prompt user
      if(verbose) std::cerr << "revertToLinearSTM: " << revertToLinearSTM << "\n";
    }*/
    
    /*//Try a version with the planar Monodromy matrix (buggy, but doesn't help)
    if(tryPlanar) {
      //Numerical STM's version for planar Monodromy
      ok = computeEigenPlanar<MAP>(theMap, stateDataPerReturn, stmPerReturn, chain, verbose);
      //Now, test for nan values in eigenvalues or eigenvectors
      // Note:  some numericalSTM create erroneous values that generate nans
      if(!ok) {revertToLinearSTM = true;}
      else {
        for(int i=0;i<period;i++) {
        if (ISNANCALL( chain[i].eval[0]) || ISNANCALL( chain[i].eval[1] ) ||
            ISNANCALL( chain[i].evec[0][0] ) || ISNANCALL( chain[i].evec[0][1] ) ||
            ISNANCALL( chain[i].evec[1][0] ) || ISNANCALL( chain[i].evec[1][1] ) ) {
          revertToLinearSTM = true;
        }
        }
      }
      //Prompt user
      if(verbose) std::cerr << "revertToLinearSTM: " << revertToLinearSTM << "\n";
    }*/
    
    //Switch between linearSTM with equations or numericalSTM
    if (linearSTM || revertToLinearSTM) {
        ok = computeEigen<MAP>(theMap, stateDataPerReturn, chain, verbose);
        //MatrixXd mon = mat2MatrixXd(stateDataPerReturn.back().J);
        //std::cerr << " Analytical (Linear) Monodromy =\n" << mon << "\n\n";
    }
    
    if (!ok) {
        throw std::runtime_error("Error in computing eigenvalues and eigenvectors.");
    }
    
    //Can we trust these results?
    //  -If no, we reset eigenvalues to be reciprocal pairs and set eigenvectors as NANs.
    if (!linearSTM && !revertToLinearSTM) {
        //Numerical STM test:
        monodromyTrustTest<MAP>(theMap, numericalSTM, chain, verbose);
    } else {
        //Linear model STM test:
        MatrixXd monoM = mat2MatrixXd<mat6>(stateDataPerReturn.back().J);
        monodromyTrustTest<MAP>(theMap, monoM, chain, verbose);
    }
    
    //Successfully computed eigenvalues & eigenvectors
    return true;
}

/// A function for checking if a fixed point is unique compared to other members of a list
template<typename InputIterator, typename PLANEMET>
InputIterator epsilon_find(const InputIterator first, const InputIterator last,
                           const nvis::vec2& x, const PLANEMET& _metric, double eps)
{
    // return element corresponding to min distance under epsilon to
    // reference position
    std::map<double, InputIterator> norm_to_it;
    for (InputIterator i=first ; i!=last ; ++i) {
        const nvis::vec2& y = i->pos;
        norm_to_it[_metric.distance(x, y)] = i;
    }
    if (norm_to_it.begin()->first < eps) {
        return norm_to_it.begin()->second;
    }
    
    return last;
}

}// orbital
#endif //__MONODROMY_HPP__
