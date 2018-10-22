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

/*
  * OrbitAvgStandardMap.cu
  *  - Performs orbit averaging texture method for Standard Map
  *
  * Author:  Wayne Schlei
  * Date: 4/24/2012
  * Mod:  10/13/2014 - Update to CUDA 6.5
  *
  * Note:  Assumes a 2D grid. To Use:  call cudaOrbitAvgStandardMap();
*/

#include <stdio.h>
#include <stdlib.h>

// Include statements
//#include <cuda.h>
#include <cuda_runtime.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <driver_types.h>
#include <driver_functions.h>
//Obsolete CUDA Utilities - deprecated in CUDA 6.5
//#include <cutil_inline.h>
//#include <cutil_math.h>
#include <helper_cuda.h>
#include <helper_functions.h>
#include <helper_math.h>
#include <helper_double_math.h>

// Constants
#define PI 3.1415926535897932

// GPU Textures
texture<float, 2,cudaReadModeElementType>  texWhiteNoise;
texture<float4, 2, cudaReadModeElementType> texRGBANoise;

// GPU constant memory
//Bounding Box(bbox.x = min_x, bbox.y = min_y, bbox.z = max_x, bbox.w = max_y)
__constant__ float min_x;
__constant__ float min_y;
__constant__ float max_x;
__constant__ float max_y;
__constant__ float k; //Chaos parameter
__constant__ float stretch; //A color stetching factor

//Debug helper message
void checkCUDAError(const char *msg) {
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err) {
    fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
    exit(EXIT_FAILURE);
  }
}

//Setting Device constants
extern "C"
void cudaSetStandardMapConstants(float4 bounds, float chaosParam, float stretchFactor) {
	float minx,miny,maxx,maxy;
	minx = bounds.x; miny = bounds.y;
	maxx = bounds.z; maxy = bounds.w;
	//Copy constants from __host__ to __device__
	checkCudaErrors ( cudaMemcpyToSymbol (min_x,&minx,sizeof(float)) );
	checkCudaErrors ( cudaMemcpyToSymbol (min_y,&miny,sizeof(float)) );
	checkCudaErrors ( cudaMemcpyToSymbol (max_x,&maxx,sizeof(float)) );
	checkCudaErrors ( cudaMemcpyToSymbol (max_y,&maxy,sizeof(float)) );
	checkCudaErrors ( cudaMemcpyToSymbol (k,&chaosParam,sizeof(float)) );
	checkCudaErrors ( cudaMemcpyToSymbol (stretch,&stretchFactor,sizeof(float) ));
}

//Identify the position for a given pixel -> maybe not necessary
__device__ float2 getPosition( int2 pixel, const uint2 GridDims) {
	//Evaluate the spacing
	float2 h = make_float2( 0.0, 0.0);
	h.x = (max_x - min_x) / (float) (GridDims.x - 1);
	h.y = (max_y - min_y) / (float) (GridDims.y - 1);

	//Set Position as center of specified cell
	float2 pos = make_float2(
		min_x + (float) pixel.x * h.x + h.x/2.0,
		min_y + (float) pixel.y * h.y + h.y/2.0
	);

	return (pos);
}

//Identify the new pixel location
__device__ int2 getPixel(float2 pos, const uint2 GridDims) {
	//Evaluate the span of the grid
	float2 span = make_float2(max_x-min_x,max_y-min_y);

	//Move position inside of grid (modulo)
	//-> Maybe texture does it for me with wrap mode


	int2 pixel = make_int2(0,0);

	pixel.x = (int) ((pos.x/span.x)*(float)GridDims.x);
	pixel.y = (int) ((pos.y/span.y)*(float)GridDims.y);
	return (pixel);
}

//Running the Map in Forward Time
__device__ float2 advectStandardMapFwd( float2 pos ) {
	float2 newPos = make_float2( 0.0 , 0.0 );
	//EOMS -> But just Difference Equation
	newPos.y = pos.y + k / ( 2.0*PI ) * sin ( 2.0*PI * pos.x);
	newPos.x = pos.x + newPos.y;
	return (newPos);
}

//Modulous distance -  for Repeating maps like the Standard Map
__device__ float2 modulo(float2 pos) {
	//Transform
	float2 val = pos-make_float2(min_x,min_y);
	float2 span = make_float2(max_x-min_x,max_y-min_y);
	float2 r = make_float2(fmod(val.x,span.x),fmod(val.y,span.y)); //Remainder of val/span
	//Check x
	if (r.x==0.0) {
	 	val.x = pos.x - min_x;
	} else {
		val.x = (val.x>=0 ? r.x : span.x + r.x);
	}
	//Check y
	if (r.y==0.0) {
		val.y = pos.y - min_y;
	} else {
		val.y = (val.y>=0 ? r.y : span.y + r.y);
	}
	//Add min back on - Transform back
	val.x += min_x;
	val.y += min_y;
	return val;
}


__global__ void testkernel(float *imageOut, const int numIters, const uint2 GridDims) {
	unsigned int x = threadIdx.x + __umul24(blockIdx.x, blockDim.x);
	unsigned int y = threadIdx.y + __umul24(blockIdx.y, blockDim.y);

	if (x >= GridDims.x || y >= GridDims.y) return;
	int index = x+y*GridDims.x;
	imageOut[index] = ((float) index) * 0.001;
}

//Oribt-Avg Kernel: Gray channel
__global__ void orbit_average_stdMap_kernel(float *imageOutDevice, const int numIters, const uint2 GridDims) {
	unsigned int x = threadIdx.x + __umul24(blockIdx.x, blockDim.x);
	unsigned int y = threadIdx.y + __umul24(blockIdx.y, blockDim.y);

	if (x >= GridDims.x || y >= GridDims.y) return;

	//Get start position (normalized)
	int2 pixel = make_int2(x,y);
	float2 pos = getPosition(pixel, GridDims);

	//Transform Coords - Normalized
	float u = x / (float) GridDims.x;
	float v = y / (float) GridDims.y;


	//Loop for all iterations
	float acc = tex2D(texWhiteNoise,u,v);
	for (int ii=0; ii<numIters; ii++) {
		//Advect Standard Map
		float2 newPos = advectStandardMapFwd( pos );
		//Run Modulous math on position
		newPos = modulo(newPos);
		//Evaluate new pixel location in grid
		pixel = getPixel(newPos,GridDims);
		//Convert (int2) pixel to normalized coords
		u = pixel.x / (float) GridDims.x;
		v = pixel.y / (float) GridDims.y;
		//Get pixel color from white noise
		acc += tex2D(texWhiteNoise,u,v);
		pos = newPos;
	}

	//Evaluate the resulting color
	acc /= (float) (numIters+1);
	acc = (acc-0.5) * stretch + 0.5; //Apply stretch factor

	//Store resulting color
	int index = x+y*GridDims.x;
	imageOutDevice[index] = acc;

	//__syncthreads();
}

//Oribt-Avg Kernel: RGBA 4-channel
__global__ void orbit_average_stdMap_color_kernel(
	float *imageOutDevice, const int numIters, const uint2 GridDims) {
	unsigned int x = threadIdx.x + __umul24(blockIdx.x, blockDim.x);
	unsigned int y = threadIdx.y + __umul24(blockIdx.y, blockDim.y);

	if (x >= GridDims.x || y >= GridDims.y) return;

	//Get start position (normalized)
	int2 pixel = make_int2(x,y);
	float2 pos = getPosition(pixel, GridDims);

	//Transform Coords - Normalized
	float u = x / (float) GridDims.x;
	float v = y / (float) GridDims.y;


	//Loop for all iterations
	float4 acc = tex2D(texRGBANoise,u,v);
	for (int ii=0; ii<numIters; ii++) {
		//Advect Standard Map
		float2 newPos = advectStandardMapFwd( pos );
		//Run Modulous math on position
		newPos = modulo(newPos);
		//Evaluate new pixel location in grid
		pixel = getPixel(newPos,GridDims);
		//Convert (int2) pixel to normalized coords
		u = pixel.x / (float) GridDims.x;
		v = pixel.y / (float) GridDims.y;
		//Get pixel color from white noise
		acc += tex2D(texRGBANoise,u,v);
		pos = newPos;
	}

	//Evaluate the resulting color
	acc /= (float) (numIters+1);
	float4 half = make_float4(0.5f,0.5f,0.5f,0.5f);
	acc = (acc-half) * stretch + half; //Apply stretch factor

	//Store resulting color
	int index = x+y*GridDims.x;
	imageOutDevice[4*index+0] = acc.x;
	imageOutDevice[4*index+1] = acc.y;
	imageOutDevice[4*index+2] = acc.z;
	imageOutDevice[4*index+3] = acc.w;

	//__syncthreads();
}

//Format Descriptor for textures:
//Single value float:
//cudaChannelFormatDesc channelDescFloat =
//		cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
//float4:
//cudaChannelFormatDesc channelDescFloat4 =
//		cudaCreateChannelDesc(32,32,32,32,cudaChannelFormatKindFloat);

//Texture-based Map Visualization for Standard Map:
//-------------------------------------------------------------------------------------
//Intensity only (black/white):
extern "C"
void cudaOrbitAvgStandardMap( float *imageOut, float *noiseData, const uint2 GridDims,const int numIters)
{
	//Block and thread dimensions
	dim3  threads( 16, 16, 1);
	dim3  blocks( (GridDims.x % threads.x !=0) ? (GridDims.x / threads.x + 1) : (GridDims.x / threads.x) ,
		      (GridDims.y % threads.y !=0) ? (GridDims.y / threads.y + 1) : (GridDims.y / threads.y) );

	//Allocate space for output data
	float *imageOutDevice;
	int imageSize = GridDims.x*GridDims.y;
	checkCudaErrors (  cudaMalloc ( (void**)&imageOutDevice, imageSize * sizeof ( float ) ) );

	//Channel Descriptor
	cudaChannelFormatDesc channelDescFloat =
		cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat); //Gray-scale only


	//Set Texture parameters ->Note:  These statements only work with cudaArrays or pitch linear memory (not cudaMalloc)
	texWhiteNoise.addressMode[0] = cudaAddressModeWrap; //Wraps texture so x=1.1 (out of grid) becomes x=0.1
	texWhiteNoise.addressMode[1] = cudaAddressModeWrap;
	//texWhiteNoise.filterMode     = cudaFilterModeLinear; //Linear interpolation between cells.
	texWhiteNoise.filterMode     = cudaFilterModePoint; //Like Nearest-Neighbor
	texWhiteNoise.normalized     = true; //Normalized Coords (u,v on [0-1])


	//Host noise array
	//In this case, noiseData is already appropriately allocated
	//Allocate Noise array on Device
    cudaArray *noiseDataDevice = 0;
    checkCudaErrors(cudaMallocArray (&noiseDataDevice, &channelDescFloat, GridDims.x, GridDims.y));

	//Copy noiseData to device array
	checkCudaErrors (
	  cudaMemcpyToArray(noiseDataDevice, 0, 0, noiseData, sizeof(float)*imageSize, cudaMemcpyHostToDevice)
	);
	//Bind Array to Texture
	checkCudaErrors ( cudaBindTextureToArray(texWhiteNoise, noiseDataDevice, channelDescFloat) );



	//Run Map and Orbit Averaging Procedure -> Returns a float4 per pixel for coloring
	orbit_average_stdMap_kernel<<<blocks,threads>>>(imageOutDevice,numIters,GridDims);
	checkCUDAError("run orbit_avg");
	//testkernel<<<blocks,threads>>>(imageOutDevice,numIters,GridDims);
	//checkCUDAError("run testkernel");

	//Copy Data from Device to Host
	checkCudaErrors ( cudaMemcpy (imageOut, imageOutDevice, imageSize*sizeof(float), cudaMemcpyDeviceToHost));

	//Unbind
	cudaUnbindTexture( texWhiteNoise );
	//Free Device Memory
	checkCudaErrors ( cudaFreeArray(noiseDataDevice) );
	checkCudaErrors ( cudaFree((void*)imageOutDevice) );
}

//Color version -> 4D vector [r g b alpha]
extern "C"
void cudaColorOrbitAvgStandardMap( float *imageOut, float *noiseData, const uint2 GridDims,const int numIters)
{
	//Block and thread dimensions
	dim3  threads( 16, 16, 1);
	dim3  blocks( (GridDims.x % threads.x !=0) ? (GridDims.x / threads.x + 1) : (GridDims.x / threads.x) ,
		      (GridDims.y % threads.y !=0) ? (GridDims.y / threads.y + 1) : (GridDims.y / threads.y) );

	//Allocate space for output data
	float *imageOutDevice;
	int imageSize = GridDims.x*GridDims.y;
	checkCudaErrors (  cudaMalloc ( (void**)&imageOutDevice, 4 * imageSize * sizeof ( float ) ) );

	//Channel Descriptor
	cudaChannelFormatDesc channelDescFloat4 =
		cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindFloat); //4D color+transparency

	//Set Texture parameters
	//Note:  These statements only work with cudaArrays or pitch linear memory (not cudaMalloc)
	texRGBANoise.addressMode[0] = cudaAddressModeWrap; //Wraps texture so x=1.1 (out of grid) becomes x=0.1
	texRGBANoise.addressMode[1] = cudaAddressModeWrap;
	//texRGBANoise.filterMode     = cudaFilterModeLinear; //Linear interpolation between cells.
	texRGBANoise.filterMode     = cudaFilterModePoint; //Like Nearest-Neighbor
	texRGBANoise.normalized     = true; //Normalized Coords (u,v on [0-1])


	//Host noise array
	//In this case, noiseData is already appropriately allocated
	//Allocate Noise array on Device
    	cudaArray *noiseDataDevice = 0;
    	checkCudaErrors(cudaMallocArray (&noiseDataDevice, &channelDescFloat4, GridDims.x, GridDims.y));

	//Copy noiseData to device array
	checkCudaErrors (
	  cudaMemcpyToArray(noiseDataDevice, 0, 0, noiseData, 4*imageSize*sizeof(float), cudaMemcpyHostToDevice)
	);
	//Bind Array to Texture
	checkCudaErrors ( cudaBindTextureToArray(texRGBANoise, noiseDataDevice, channelDescFloat4) );



	//Run Map and Orbit Averaging Procedure -> Returns a float4 per pixel for coloring
	orbit_average_stdMap_color_kernel<<<blocks,threads>>>(imageOutDevice,numIters,GridDims);
	checkCUDAError("run orbit_avg");
	//testkernel<<<blocks,threads>>>(imageOutDevice,numIters,GridDims);
	//checkCUDAError("run testkernel");

	//Copy Data from Device to Host
	checkCudaErrors ( cudaMemcpy (imageOut, imageOutDevice, 4*imageSize*sizeof(float), cudaMemcpyDeviceToHost));

	//Unbind
	cudaUnbindTexture( texRGBANoise );
	//Free Device Memory
	checkCudaErrors ( cudaFreeArray(noiseDataDevice) );
	checkCudaErrors ( cudaFree((void*)imageOutDevice) );
}
