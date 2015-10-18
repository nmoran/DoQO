/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 *
 * Header file for utility functions for dealing with momentum sectors.
 * 
 * Version: $Id$ 
 * 
*/



#ifndef _MOMENTUM_H_
#define _MOMENTUM_H_

#include "doqo.h"
#include "hamiltonian.h"
#include "utils.h"
#include "conservation.h"

int parse_momentum_file(struct parameter_struct *parameters);
void print_momentum_information(struct parameter_struct *parameters);
uPetscInt translate_using_masks_spins(struct parameter_struct *parameters, uPetscInt basis, PetscInt *amount,PetscScalar *phase);
uPetscInt translate_using_masks_fermions(struct parameter_struct *parameters, uPetscInt basis, PetscInt *amount,PetscScalar *phase);
int calculate_translation_masks(struct parameter_struct *parameters);
PetscScalar get_translation_phase(struct parameter_struct *parameters,PetscInt *amount, PetscInt *sector);
PetscScalar* calculate_translation_phases(struct parameter_struct *parameters);
PetscTruth is_representative(struct parameter_struct *parameters, uPetscInt basis, PetscReal *normal, PetscTruth phase_check, PetscInt *sector);
uPetscInt get_representative(struct parameter_struct *parameters, uPetscInt basis, PetscTruth phase_check, PetscScalar *phase, PetscInt *sector);
PetscInt MatrixDistributeMomentumSector(Mat *A,struct parameter_struct *parameters, struct matrix_meta_data *matrix_data, PetscInt *mom_sector);
PetscInt MatrixDistributeMomentumConservationSector(Mat *A,struct parameter_struct *parameters, struct matrix_meta_data *matrix_data, PetscInt *mom_sector, PetscInt cons_sector);
void test_momentum_tools(struct parameter_struct *parameters);
PetscTruth dummy_representative_check(struct parameter_struct *parameters, uPetscInt index); 
uPetscInt  dummy_group_representative(struct parameter_struct *parameters, uPetscInt basis, PetscTruth phase_check, PetscScalar *phase, PetscTruth *group_exists); 
PetscReal  dummy_group_normal(struct parameter_struct *parameters, uPetscInt index); 
PetscTruth momentum_representative_check(struct parameter_struct *parameters, uPetscInt index);
uPetscInt  momentum_group_representative(struct parameter_struct *parameters, uPetscInt basis, PetscTruth phase_check, PetscScalar *phase, PetscTruth *group_exists) ;
PetscReal  momentum_group_normal(struct parameter_struct *parameters, uPetscInt index);
#endif

