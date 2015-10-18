/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 *
 * Header file for solver functions  
 * 
 * Version: $Id$ 
 * 
*/



#ifndef _SOLVER_FUNCS_H_
#define _SOLVER_FUNCS_H_

#include "doqo.h"
#include "utils.h"

int createSolveDisplayProb(Mat *A, EPS *eps, int eigPairs ,PetscReal tol ,int maxits);
int createEigenSolver(struct parameter_struct *parameters,Mat *A, EPS *eps, PetscInt eigenPairs, PetscReal tol, int maxits);
PetscInt solveEigenProblem(EPS *eps,struct parameter_struct *parameters,struct solver_results *results);
int getEigenValues(EPS *eps, PetscReal *eigenvalues, PetscInt nconv);
int getEigenVector(EPS *eps, Vec eigenvector_real, Vec eigenvector_imag, PetscInt n);
int cleanupEigenSolver(struct parameter_struct *parameters, EPS *eps, Mat *A);
int degeneracy_analysis(struct solver_results *results, PetscReal tol);
int getSolverResults(EPS *eps, struct parameter_struct *parameters, struct solver_results *results); 
PetscInt MyEPSComputeRelativeError(EPS *eps,PetscInt i,PetscReal *error);

#endif
