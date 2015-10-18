/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 *
 * Header file for spectral flow functions. 
 * 
 * Version: $Id$ 
 * 
*/


#ifndef _SPECTRAL_FLOW_H_
#define _SPECTRAL_FLOW_H_

#include "doqo.h"
#include "hamiltonian.h"
#include "utils.h"

int setup_spectral_flow_info(struct parameter_struct *parameters,struct matrix_meta_data *matrix_data) ;
uPetscInt translate_using_masks_fermions_spectral(struct parameter_struct *parameters, uPetscInt basis, PetscInt *amount, PetscScalar *phase);
PetscReal  momentum_group_normal_spectral(struct parameter_struct *parameters, uPetscInt index);
#endif
