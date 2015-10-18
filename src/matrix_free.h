/*
 * Copyright (C) 2008    Niall Moran
 * 
 * Header file for matrix free method code.
 * 
 * Version: $Id$ 
 * 
*/

#ifndef _MATRIX_FREE_H_
#define _MATRIX_FREE_H_

#include "doqo.h"
#include "utils.h"


PetscErrorCode MatSpinSys_Mult( Mat A, Vec x, Vec y );
int prepare_matrix_free(struct shell_matrix_context *shell_context);
int cleanup_matrix_free(struct shell_matrix_context *shell_context);
int getrow_matrix_free(struct shell_matrix_context *shell_context, PetscInt Istart, PetscInt Iend , PetscInt Istart_other, PetscInt Iend_other, PetscScalar *px, PetscScalar *buf) ;



#endif
