#ifndef _DEVICE_FOURIER_TRANSFORM_H_
#define _DEVICE_FOURIER_TRANSFORM_H_

#include<iostream>
#include<string>
#include<fstream>
#include <sstream>
#include <algorithm>
#include <mpi.h>

#include<getopt.h>



using namespace std;

// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order.
__device__ void get_monomial_bosons_cuda(unsigned long Index, int Sites, int Bosons, int *Monomial);

// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order and also return the normal.
__device__ void get_monomial_bosons_cuda(unsigned long Index, int Sites, int Bosons, int *Monomial, float *Normal);

// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order.
__device__ void get_monomial_fermions_cuda(unsigned long Index, int Sites, int Fermions, int *Monomial) ;

__device__ int permute_cuda(int *Monomial, int Len);

__device__ float configuration_normal_cuda(unsigned long Config, int NumberSites, unsigned long Mask);

__global__ void calculate_FT_elements(int VectorLength, float *VectorElements, unsigned long *BasisElements, int FTBasisDim, float *FTVectorElements, unsigned long *FTBasisElements, int NumberSites, int Filling, int CurrentMomentum, int SpinMomentum, int *fmonomial, int *fmonomialTmp, int *omonomial, float *phaseArrayR);

void fourier_transform_vector_cuda(int VectorLength, float *VectorElements, unsigned long *BasisElements, int FTBasisDim, float *FTVectorElements, unsigned long *FTBasisElements, int NumberSites, int Filling, int CurrentMomentum, int SpinMomentum, bool single, int CudaDevice, string OutputPrefix, int Offset);

//void fourier_transform_vector_cuda2(int VectorLength, float *VectorElements, unsigned long *BasisElements, int FTBasisDim, float *FTVectorElements, unsigned long *FTBasisElements, int NumberSites, int Filling, int CurrentMomentum, int SpinMomentum);
  
#endif