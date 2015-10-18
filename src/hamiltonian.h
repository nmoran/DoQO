/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala.
 *
 * Header file for hamiltonian functions. 
 * 
 * Version: $Id$ 
 * 
*/



#ifndef _HAMILTONIAN_H_
#define _HAMILTONIAN_H_

#include <string.h>
#include "doqo.h"
#include "utils.h"
#include "spectral_flow.h"
#include "bst.h"
#include "pqsort.h"


//these macros will be used in the fermionic getrow and prepare_term functions.
#define ANNIHILATION_OPERATOR 0
#define CREATION_OPERATOR 1
#define SPIN_ONE_Z_OPERATOR 2
#define SPIN_ONE_RAISING 3
#define SPIN_ONE_LOWERING 4

int setup_operators(struct parameter_struct *parameters);
int setup_operator(struct parameter_struct *parameters , struct matrix_meta_data *matrix_data );
int expand_format(char **in, char **out, int columns, int rows);
int prepare_matrix_data(struct parameter_struct *parameters,struct matrix_meta_data *matrix_data);
int getrow_spin_half(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters, int nz );
int getrowlocal_spin_half(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters,int nz,PetscInt Istart, PetscInt Iend );
int getrow_fermionic(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters, int nz );
int getrowlocal_fermionic(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters,int nz,PetscInt Istart, PetscInt Iend );
int getrow_spin_one(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters, int nz);
int getrowlocal_spin_one(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters, int nz,PetscInt Istart, PetscInt Iend );
int allocate_row_space(struct matrix_row_values *row, struct matrix_meta_data *matrix_data);
int deallocate_row_space(struct matrix_row_values *row);
int deallocate_matrix_data(struct matrix_meta_data *matrix_data);
int allocate_matrix_data(struct matrix_meta_data *matrix_data);
//int print_matrix_info(struct matrix_meta_data *md );
//int get_parameter_name_idx(struct matrix_meta_data *md , char *str);
int get_term(char *str, int *spin, char *type, int *pos);
int prepare_term_data(struct parameter_struct *parameters, struct matrix_meta_data *matrix_data, int term );
PetscInt build_local_basis_lists(struct parameter_struct *parameters, uPetscInt **local_reps_array, uPetscInt **ending_idx_for_proc);
PetscInt build_local_basis_lists_using_disk(struct parameter_struct *parameters, uPetscInt **local_reps_array, uPetscInt **ending_idx_for_proc);
PetscInt MatrixDistribute(Mat *A,struct parameter_struct *parameters, struct matrix_meta_data *matrix_data);
PetscInt exchange_required_basis_indices(struct parameter_struct *parameters,struct matrix_meta_data *matrix_data,uPetscInt *local_reps_array,uPetscInt *ending_idx_for_proc,PetscInt **d_nnz, PetscInt **o_nnz,uPetscInt **reps, uPetscInt **indices,uPetscInt *num_unique_cols);
PetscInt set_matrix_elements(Mat *A,struct parameter_struct *parameters,struct matrix_meta_data *matrix_data, uPetscInt *local_reps_array,uPetscInt *ending_idx_for_proc,uPetscInt *reps, uPetscInt *indices,uPetscInt num_unique_cols);
PetscInt create_matrix_communicator(struct parameter_struct *parameters);
#endif
