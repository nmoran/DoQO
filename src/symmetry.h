/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 *
 * Header file for symmetry functions. These functions implement the symmetries. 
 * 
 * Version: $Id$ 
 * 
*/


#ifndef _SYMMETRY_H_
#define _SYMMETRY_H_

#define SWAP(a,b,c) { c = a; a = b; b = c; }

#include "doqo.h"
#include "hamiltonian.h"
#include "utils.h"

int search_for_parity_sectors(struct matrix_meta_data *md, struct parameter_struct *parameters) ;
int read_parity_masks(char *filename, uPetscInt **masks, int *number);
int get_mask_parity(int sector, int mask_number);
uPetscInt num_permutations_at_parity(struct parameter_struct *parameters, uPetscInt mask,int parity);
int free_parity_info(struct parameter_struct *parameters);
uPetscInt NchooseC(int n , int c );
uPetscInt NchooseC_safe(int n , int c );
int factorial(int from, int to);
uPetscInt get_full_parity_basis_index(struct parameter_struct *parameters, uPetscInt reduced_index);
uPetscInt get_reduced_parity_basis_index(struct parameter_struct *parameters, uPetscInt full_index, PetscTruth *exists);
uPetscInt get_index_with_rank_parity(struct parameter_struct *parameters, uPetscInt rank, int mask, int parity);
uPetscInt get_index_n_choose_c(uPetscInt rank, int n, int c);
uPetscInt get_rank_n_choose_c(uPetscInt index, int n, int c);
uPetscInt get_rank_with_index_parity(uPetscInt index, int parity, int count );
uPetscInt get_size_of_parity_sector(struct parameter_struct *parameters, uPetscInt sector);
#endif
