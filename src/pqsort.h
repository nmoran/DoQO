/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 *
 * Header file for parallel quick sort function. 
 * 
 * Version: $Id$ 
 * 
*/


#ifndef _PQSORT_H_
#define _PQSORT_H_

#include "doqo.h"

int parallel_qsort_uPetscInt(struct parameter_struct *parameters, MPI_Comm *subcomm, uPetscInt *array, PetscInt llength);
#endif
