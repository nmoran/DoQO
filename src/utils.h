/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 *
 * Header file for utility functions  
 * 
 * Version: $Id$ 
 * 
*/



#ifndef _UTILS_H_
#define _UTILS_H_

#include "doqo.h"
#include "hamiltonian.h"
#include "conservation.h"
#include "momentum.h"
#include "tasks.h"
#include "rotation.h"

/*#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif*/

int 	initialise(int argc, char **argv , struct parameter_struct *parameters );
int 	parse_main_input_file(struct parameter_struct *parameters);
int 	read_input_files(struct parameter_struct *parameters);
void 	write_doqo_heading(struct parameter_struct *parameters);
int 	assemble_filename(struct parameter_struct *parameters);
int 	write_symmetry_information(struct parameter_struct *parameters, struct symmetry_information *symmetry_info);
PetscInt get_rank_chunk_size(PetscInt n, PetscInt rank, PetscInt size); 
PetscInt get_max_chunk_size(PetscInt n, PetscInt size); 
PetscInt get_rank_chunk_start(PetscInt n, PetscInt rank, PetscInt size);
char* 	dec2bin(uPetscInt a, int len, char *num);
uPetscInt bin2dec(char *bnum);
//void getUsedMem(long long *mem);
PetscInt find_in_sorted_array(uPetscInt *array,PetscInt len, uPetscInt item); //find index of item in sorted array
int 	insert_into_sorted_list(list<uPetscInt> *sorted_list, uPetscInt item);
void 	output_memory_usage();
int trim(const string str, string *trimmed) ;
PetscInt get_range_with_item(PetscInt n, PetscInt item, PetscInt size);
int 	parse_adjacency_info_file(struct parameter_struct *parameters, char *filename);
int 	deallocate_adjacency_info(struct parameter_struct *parameters);
PetscTruth check_adjacency(struct parameter_struct *parameters, uPetscInt index);
PetscTruth dummy_basis_check(struct parameter_struct *parameters, uPetscInt index);
int 	save_matrix(struct parameter_struct *parameters,Mat *A, const char *filename);
int 	tokenize_string(string str, string delim, vector<string> *parts);
PetscScalar str2scalar(string str); 
int find_valid_configs(struct parameter_struct *parameters, ofstream *basis_file, PetscInt *local_reps);
int find_valid_configs_with_this_subset(struct parameter_struct *parameters, int set_size, int set_filling, uPetscInt set, ofstream *basis_file, PetscInt *local_reps );
int parse_prod_wf_overlap_file(struct parameter_struct *parameters, char *prod_wf_file);
int deallocate_prod_wf_info(struct parameter_struct *parameters);
int create_prod_wf_state(struct parameter_struct *parameters, Vec *x);
PetscTruth susy_staggered_subspace_check(struct parameter_struct *parameters, uPetscInt index);
uPetscInt str_to_mask(const char *buf);
#endif

