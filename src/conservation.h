/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 *
 * Header file for symmetry functions. These functions implement the symmetries. 
 * 
 * Version: $Id$ 
 * 
*/


#ifndef _CONSERVATION_H_
#define _CONSERVATION_H_

#include "doqo.h"
#include "hamiltonian.h"
#include "utils.h"
#include "symmetry.h"

uPetscInt get_size_of_conservation_sector(struct parameter_struct *parameters, uPetscInt sector);
uPetscInt get_mask_filling(struct parameter_struct *parameters, uPetscInt sector, int mask_number);
int search_for_conservation_sectors(struct matrix_meta_data *md, struct parameter_struct *parameters);
uPetscInt get_full_conservation_basis_index(struct parameter_struct *parameters, uPetscInt reduced_index);
uPetscInt get_reduced_conservation_basis_index(struct parameter_struct *parameters, uPetscInt full_index, PetscTruth *exists);
uPetscInt get_full_basis_index_passthru(struct parameter_struct *parameters, uPetscInt reduced_index);
uPetscInt get_reduced_basis_index_passthru(struct parameter_struct *parameters, uPetscInt full_index,PetscTruth *exists);
#endif
