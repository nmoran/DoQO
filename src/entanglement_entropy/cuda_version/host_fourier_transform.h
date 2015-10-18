#ifndef _HOST_FOURIER_TRANSFORM_H_
#define _HOST_FOURIER_TRANSFORM_H_

#include<iostream>
#include<string>
#include<fstream>
#include <sstream>
#include <algorithm>
#include <mpi.h>

#include<getopt.h>



using namespace std;

// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order.
void get_monomial_bosons(unsigned long Index, int Sites, int Bosons, int *Monomial);

// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order and also return the normal.
void get_monomial_bosons(unsigned long Index, int Sites, int Bosons, int *Monomial, float *Normal);


// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order.
void get_monomial_fermions(unsigned long Index, int Sites, int Fermions, int *Monomial) ;

int permute(int *Monomial, int Len);


float configuration_normal(unsigned long Config, int NumberSites, unsigned long Mask);


void fourier_transform_vector(int VectorLength, float *VectorElements, unsigned long *BasisElements, int FTBasisDim, float *FTVectorElements, unsigned long *FTBasisElements, int NumberSites, int Filling, int CurrentMomentum, int SpinMomentum);
  
#endif