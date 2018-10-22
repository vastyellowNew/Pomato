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


/// Orbit Convolution
/// Convolution process functions for maps
// Author: Wayne Schlei

#ifndef ORBIT_CONVOLUTION_HPP
#define ORBIT_CONVOLUTION_HPP

#include <iostream>
#include <vector>
#include <complex>
#include <sstream>
#include <boost/random.hpp>

#include <util/wall_timer.hpp>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif
//#include <mpi.h> //Usually only needed for propagation

namespace OrbitConvolution {

/// Converting a color from hsv to rgb assuming each are 3-element vectors with []-operator
template<typename COLOR>
COLOR hsvtorgb( COLOR hsv )
{
    if( hsv[2] == 0 ) {
        return COLOR( 0, 0, 0 );
    } else if( hsv[1] == 0 ) {
        return COLOR( hsv[2], hsv[2], hsv[2] );
    }
    
    const float hf = hsv[0] * 6.0;
    const int    i = (int)floor( hf );
    const float  f = hf - i;
    const float pv = hsv[2] * ( 1 - hsv[1] );
    const float qv = hsv[2] * ( 1 - hsv[1] * f );
    const float tv = hsv[2] * ( 1 - hsv[1] * ( 1 - f ) );
    
    switch( i ) {
        case 0:
            return COLOR( hsv[2], tv, pv );
        case 1:
            return COLOR( qv, hsv[2], pv );
        case 2:
            return COLOR( pv, hsv[2], tv );
        case 3:
            return COLOR( pv, qv, hsv[2] );
        case 4:
            return COLOR( tv, pv, hsv[2] );
        case 5:
            return COLOR( hsv[2], pv, qv );
    }
}

/** Map Data to NRRD Conversion function
 *  - Converts return data from the map to NRRD formatted data (and saves file)
 *  - NRRD formatted data is an array of ints indicating what index in the image each return is
 *     unless the return is non-existent which is signified by a -1
 *  - DATASET is usually a std::map<>
 *  - RETURN_DATA is typically std::vector<return_type>  [ref. maps/mapNd.hpp at return_map_info()]
 */
template<class DATASET, class RETURN_DATA>
void convertMapDataToNRRD(
    const nvis::ivec2& res, const int niter, const nvis::bbox2& imageBounds,
    DATASET& mapData, const char* nDataFilename)
{
    typedef typename DATASET::iterator      MapIterator;
    typedef typename RETURN_DATA::const_iterator    DataIterator;
    const unsigned int resx = res[0];
    const unsigned int resy = res[1];
    std::vector<int> result(resx*resy*niter,-1);
    
    
    #pragma omp parallel for schedule(dynamic,1)
    for(int k=0; k<resx*resy; ++k) {
        const int ix = k % resx;
        const int iy = k / resx;
        //Find initial condition on image
        nvis::vec2 s  = imageBounds.interpolate( nvis::vec2( (ix+0.5)/resx, (iy+0.5)/resy ) );
        //Current Data Location
        MapIterator mapIT = mapData.find( k );
        //If not found, the point is invalid (shouldn't happen)
        if (mapIT == mapData.end()) {
            std::cerr << "Can't find data for point " << k << " =  Image position " << s << "\n";
            continue;
        }
        
        const RETURN_DATA& dataVector = mapIT->second;
        std::vector<int>::iterator mi = result.begin() + niter*k;
        
        //Loop through returns of kth initial state (typically stored in vector)
        for (DataIterator dataIT = dataVector.cbegin(); dataIT!=dataVector.cend(); ++dataIT) {
            //Find the normalized coordinate of the current iterate's map position
            nvis::vec2 p = imageBounds.normalizeCoordinate( dataIT->x );
            
            //Add to data only if within bounds
            if( p[0] < 0 || p[0] > 1 || p[1] < 0 || p[1] > 1) {
                //mi++;//Don't increment pointer until we have valid data
                continue;
            }
            
            //Find the positions
            int x = floor( p[0]*resx + 0.5 );
            int y = floor( p[1]*resy + 0.5 );
            //if on edge, trace back
            if( x == resx ) {
                --x;
            }
            if( y == resy ) {
                --y;
            }
            if( x < 0 || x >= resx || y < 0 || y >= resy ) {
                printf( "error out of bounds\n" );
                std::cerr << " Values (x,y) = (" << x << "," << y
                          << ") are not in range with max (" << resx << "," << resy << ")\n";
                exit( EXIT_FAILURE );
            }
            //Store
            *(mi) = x + y*resx;
            //Increment
            mi++;
        }
    }
    
    //Write to nrrd object
    size_t size[3] = {(unsigned int)std::abs(niter),resx,resy};
    Nrrd* nout = nrrdNew();
    if ( nrrdWrap_nva(nout, &result.front(), nrrdTypeInt, 3, size) ||
            nrrdSave(nDataFilename, nout, NULL) ) {
        std::cout << "ERROR while exporting map data:" << biffGetDone(NRRD) << "\n";
        nrrdNix( nout );
    }
    
}

/** NRRD - Convolve Map Base function
 *  - Runs the convolution process for a specified number of passes.
 *  - Image is assumed to be the same size as your sampling domain
 *  NRRD input object
 */
void nrrdConvolvePoincareMap(Nrrd* nin, const char* output,
                             unsigned int& maxi, const unsigned int npasses,
                             const unsigned int seed=0, const bool sharpenEachPass=false)
{
    //Assume nin is set before this function is called
    const unsigned int niter = nin->axis[0].size;
    const unsigned int resx  = nin->axis[1].size;
    const unsigned int resy  = nin->axis[2].size;
    
    if( maxi > niter ) {
        maxi = niter;
    }
    
    boost::mt19937 rng;
    if (seed > 0) {
        rng.seed( seed );    //specified seeding
    }
    boost::uniform_01<> dist;
    fprintf( stderr, "%d %d %d\n", resx, resy, niter );
    
    const int* iterdata = (const int*)nin->data;
    
    // initialize input color noise image
    std::vector<nvis::fvec3> image( resx*resy );
    
    for( unsigned int i=0; i<resx*resy; ++i ) {
        image[i][0] = dist(rng);//drand48();
        image[i][1] = dist(rng);
        image[i][2] = dist(rng);
        // image[i] = hsvtorgb( nvis::fvec3( dist(rng), 1.0, 1.0 ) );
    }
    
    // perform convolution passes
    for( unsigned int pass=0; pass<npasses; ++pass ) {
        fprintf( stderr, "convolution pass %d\n", pass );
        
        std::vector<nvis::fvec3> conv( resx*resy, 0.0 );
        
        #pragma omp parallel for schedule(dynamic,1)
        for( unsigned int k=0; k<resx*resy; ++k ) {
            const int* iterptr = iterdata + niter*k;
            
            nvis::fvec3 accumv = image[k];
            int   accumn = 0;
            
            for( unsigned int i=0; i<maxi; ++i, ++iterptr ) {
                if( *iterptr < 0 ) {
                    break;
                }
                
                accumv += image[*iterptr];
                accumn += 1;
            }
            
            if( accumn < 1 ) {
                conv[k] = (pass == npasses-1) ? 1.0 : 0.5;
            } else {
                conv[k] = accumv / (accumn + 1);
            }
        }
        
        //Set the image as the new convolution image
        image.swap( conv );
        
        //High-pass filtering
        // -> Teem: unu resample -i tmp.nrrd -k gauss:2,2 -s = x1 x1 | unu 2op x- -0.5 | unu 2op + tmp.nrrd - -o output.nrrd
        
    }
    
    //Store resulting image as a NRRD
    Nrrd* nout = nrrdNew();
    size_t size[3] = { 3, resx, resy };
    if( nrrdWrap_nva( nout, &image.front(), nrrdTypeFloat, 3, size ) ||
            nrrdSave( output, nout, NULL ) ) {
        std::cerr << "ERROR while exporting IMAGE file: " << biffGetDone(NRRD) << '\n';
        nrrdNix( nout );
    }
}

/** Single orbit convolution pass
 *  - Executing a single orbit convolution pass given an input image and
 *     Poincare map return data within a NRRD formated object
 */
bool singleOrbitConvolutionPass(
    const unsigned int pass, const unsigned int numPasses, unsigned int& maxi,
    Nrrd* nin, std::vector<nvis::fvec3>& image, int shortFallPenalty=0)
{
    //Assume nin is set before this function is called
    const unsigned int niter = nin->axis[0].size;
    const unsigned int resx  = nin->axis[1].size;
    const unsigned int resy  = nin->axis[2].size;
    if( maxi > niter ) {
        maxi = niter;
    }
    //Data pointer
    const int* iterdata = (const int*)nin->data;
    
    //Initialized Convolved image
    std::vector<nvis::fvec3> conv( resx*resy, 0.0 );
    
    //Perform convolution pass
    #pragma omp parallel for schedule(dynamic,1)
    for( unsigned int k=0; k<resx*resy; ++k ) {
        const int* iterptr = iterdata + niter*k;
        
        nvis::fvec3 accumv = image[k];
        int   accumn = 0;
        
        for( unsigned int i=0; i<maxi; ++i, ++iterptr ) {
            if( *iterptr < 0 ) {
                break;
            }
            
            accumv += image[*iterptr];
            accumn += 1;
        }
        
        if( accumn < 1 )
            //No available data
        {
            conv[k] = (pass == numPasses-1) ? 1.0 : 0.5;
        } else if ( accumn < maxi ) {
            //Insufficient number of returns
            int numOffset = maxi - accumn - 1;
            float offRange = numOffset/(maxi);
            float graySkew = 0.8 - offRange*0.3;
            //Options
            switch (shortFallPenalty) {
                case 0 :
                    conv[k] = (pass == numPasses-1) ? 1.0 : 0.5;
                    break;
                case 1 :
                    //Penalty
                    conv[k] = (accumv + numOffset*0.5)/(maxi+1);
                    break;
                case 2 :
                    //Percentage gray penalty
                    conv[k] = (accumv + numOffset*graySkew)/(maxi+1);
                    //Gray skew indicates low num as whitish, higher number as darker
                    break;
                default :
                    //Only use the current number of returns
                    conv[k] = accumv / (accumn + 1); //Like original
                    break;
            }
        } else
            //Correct number of returns
        {
            conv[k] = accumv / (accumn + 1);
        }
    }
    
    //Set the image as the new convolution image
    image.swap( conv );
    
    //High-pass filtering after a convolution pass
    // -> Teem: unu resample -i tmp.nrrd -k gauss:2,2 -s = x1 x1
    //           | unu 2op x- -0.5 | unu 2op + tmp.nrrd - -o output.nrrd
    
    return true;
}


} // end convolution

#endif
