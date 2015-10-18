/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 * 
 * Header file for main solver code. 
 * 
 * Version: $Id$ 
 * 
*/

/*! \file */ 

#ifndef _DOQO_H_
#define _DOQO_H_

#include "slepceps.h"
#include "tinyxml.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <algorithm>

//define this constant to turn on debugging
//#define DEBUG
//#define SUPERDEBUG
#define _BOOST_BST_
#define DOQO_VERSION "1.0.0"
using namespace std;


#define 	CODE_REV "$Rev$"
#define 	MAX_STRING 200 /*!< The maximum length of a character array. Used at various places for buffers and filenames. */
#define 	MAX_PARM 20    /*!< The maximum length of the name of a parameter. */
#define 	MAX_TERM 200  /*!< The maximum length in characters of a line defining a term in the input file. */

#if defined(PETSC_USE_64BIT_INDICES)
typedef 	unsigned long long uPetscInt;
#define 	PETSC_MPI_INT MPI_LONG_LONG_INT
#define 	PETSC_UNSIGNED_MPI_INT MPI_UNSIGNED_LONG_LONG
#else
typedef 	unsigned int uPetscInt;
#define 	PETSC_MPI_INT MPI_INT
#define 	PETSC_UNSIGNED_MPI_INT MPI_UNSIGNED
#endif

#if defined(PETSC_USE_SINGLE)
#define 	PETSC_MPI_REAL MPI_REAL 
#else
#define 	PETSC_MPI_REAL MPI_DOUBLE
#endif

//minimum macro used in momentum.C
#define		my_min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define		two_pow(a)  ( ((uPetscInt)1) << a)
#define 	log2_uPetscInt(num) ((uPetscInt)log2((double)num))


//#if (PETSC_VERSION_MAJOR > 3) && (PETSC_VERSION_MINOR > 1)
typedef PetscBool PetscTruth;
//#endif

//define some enumerations
enum ModelTypeEnum {SPIN_HALF, FERMIONIC, SPIN_ONE};

/*! \brief Structure that stores an array of the edges connecting the sites of the lattice or vertices of the graph. */
struct adjacency_information {
	int 		number_masks;
	uPetscInt 	*masks;
};

struct wf_conf {
	PetscScalar coefficient;
	uPetscInt	mask;
};

struct wf_group {
	uPetscInt	mask;
	PetscInt	number_confs;
	struct wf_conf *confs;
};

struct wf_subterm {
	PetscScalar coefficient;
	PetscInt	number_groups;
	struct wf_group *groups;
};

struct wf_term {
	PetscScalar coefficient;
	PetscInt	number_subterms;
	struct wf_subterm *subterms;
};

/*! \brief Structure that stores an array of the edges connecting the sites of the lattice or vertices of the graph. */
struct prod_wf_information {
	PetscInt	number_terms; 
	struct wf_term *terms;
};

/*! \brief This structure stores the information necessary to implement the conservation of parity and is also used for conservation of filling as the same data fields can be used. */
struct parity_information {
	int 		number_parity_sectors;
	int 		number_masks;
	uPetscInt 	*masks ; //the mask for each parity sector.
	uPetscInt 	*permutations; //number of permutations for each mask.
	int 		*counts; // number of 1's in each mask
	int 		*relevant_sectors, num_relevant_sectors; //array of sectors to use
};

/*! \brief This structure stores the fields necessary to exploit translational invariance. */
struct momentum_information { 
	PetscInt 	lattice_vector[2][2]; //!< Stores the matrix vector along each direction in terms of the lattice spacings.
	PetscReal 	normalised_vector[2][2]; //stores the vector along each direction normalised to one.
	PetscInt 	num_sectors[2];  //the number of sectors in each direction.
	PetscInt 	lattice_dimensions[2]; //stores the size of the lattice in each direction.
	int 		*relevant_sectors, num_relevant_sectors; //array of sectors to use
	PetscScalar *phases; //array which will store the phase corresponding to each translation in each sector. Will save recalculating each time.
	uPetscInt 	***translation_masks;
	PetscInt 	current_sector[2];
};

/*! \brief This structure stores the fields necessary to exploit rotational invariance. */
struct rotation_information { 
	PetscInt 	num_sectors;  //the number of sectors in each direction.
	int 		*relevant_sectors, num_relevant_sectors; //array of sectors to use
	PetscScalar *phases; //array which will store the phase corresponding to each rotation in each sector. Will save recalculating each time.
	uPetscInt 	**mapping_masks;
	PetscInt 	current_sector,real_sector;
	PetscInt    sector_modulo, sector_divisor;
};

/*! \brief Structure to hold data on how to generate the matrix for a given operator. 
	Instances are referenced from the instance of the 
	parameter_struct. There is one instance for the principle operator and one for each of the additional operators. */
struct matrix_meta_data {
	int 		number_terms; /*!< The number of terms in this operator.*/
	char 		**terms;     /*!< Terms in character format. */ 
	//array of flip masks. One per part
	uPetscInt 	*flips; /*!< The flip mask for retrieving new column index. For over 32 particles 64 bit indices needs to be enables. */ 
	//array of y masks. One per y term per part
	uPetscInt 	**ys; /*!< Masks storking position of sigma y in each term. For over 32 particles 64 bit indices needs to be enables. */ 
	uPetscInt 	**zs; /*!< Masks storking position of sigma z in each term. For over 32 particles 64 bit indices needs to be enables. */ 
	//array of y counts per part
	int 		*y_count; 
	//array of z counts per part
	int 		*z_count;
	//array to keep count of X's per part 	
	int 		*x_count;
	PetscScalar *term_multipliers;
	int 		nz ; /*!< Number of non zero entries per row. */  
	int 		number_parameters;  /*!< The number of parameters. */
	int 		*parameter_indices; /*!< Stores the index of the parameter in the list of parameter names stored in the task_parameters struct. */
	int 		*number_term_params; /*!< The number of parameters associated with each term. */
	char	 	filename[MAX_STRING]; /*!< The input filename for the operator.*/
	int 		**term_parameters;   /*!< This is a 2D array with an array for each term. */
	PetscTruth  is_diagonal;
};

/*! \brief Structure to store parameter names and values for each task. */ 
struct task_parameters {
	char 		filename[MAX_STRING]; /*!< Filename of task input file.*/
	int 		number_tasks, start_task;  /*!< Number of tasks to run. */
	int 		number_parameters;
	PetscScalar	**task_parameter_values;
	vector<string> parameter_names;  	 /*!< Vector storing the names of all the parameters. */
	
};

/*! \brief This structure is used for storing the information needed to perform spectral flow calculations. This includes the name of the alpha parameter and conjugate, number
of points to use between 0 and 1 for alpha and which values to use. */
struct spectral_flow_info { 
	PetscInt 	number_points[2], current_point[2], real_current_point[2], number_relevant_points[2], *relevant_points[2];
	char 		alpha_param[2][MAX_STRING], alpha_conj_param[2][MAX_STRING];
	PetscReal	alpha[2];
	PetscScalar e_alpha[2], e_alpha_c[2], e_alpha_f_l[2], eigenvalues[2],eigenvalues_exponent[2];
};

/*! \brief This structure contains pointers to the basis information. In certain cases when it is needed later it is kept. */
struct basis_information { 
	uPetscInt *local_reps_array,*ending_idx_for_proc ;	
};

/*! \brief This structure contains pointers to the basis information. In certain cases when it is needed later it is kept. */
struct expectation_values_information {
	PetscInt 	number_spaces;
	PetscInt 	*number_vectors_per_space;
	PetscInt 	number_vectors;
	PetscInt 	momentum_sector; //not implemented yet
	PetscInt 	rotation_sector;
	PetscInt 	filling_sector;
	PetscInt 	parity_sector; // not implemented yet
	string**	vectors;
	string 		overlap_vec;
};


/*! \brief Structure to store input parameters, settings. There is only on instance of this structure per process. */
struct parameter_struct {
  int 			number_particles;
  PetscInt  	eigPairs , 
  				use_sectors, 
  				max_its , 
  				conservation_filling,
  				parity_sector, 
  				conservation_sector, 
  				momentum_sector, 
  				current_basis_size, 
  				real_parity_sector, 
  				real_conservation_sector, 
  				real_momentum_sector,
  				Istart, 
  				Iend, 
  				verbosity, 
  				current_task,
  				number_ops,
  				numberBasisStates,
  				number_rotation_ops,
  				number_rotation_sectors,
  				rotation_number_relevant_sectors,
  				rotation_current_sector,
  				*rotation_relevant_sectors;
  PetscMPIInt 	size, 
  				rank,
  				restricted_size;		
  PetscReal 	tol, deg_tol, phase_tol; 
  PetscTruth	specific_point,  
				use_shell, 
				use_sz_symmetries, 
				prepare, 
				distribute, 
				solver, 
				benchmark, 
				save_states,
				save_states_ascii,
				save_states_matlab,
				save_states_real,
				multiply_test, 
				use_parity_sectors, 
				use_conservation_sectors, 
				use_momentum_sectors,
				nn_exclusion,
				nn_recursive_algorithm,
				save_matrix,
				save_matrix_ascii,
				save_matrix_matlab,
				use_hermitian,
				use_disk,
				use_bst,
				calculate_correlations,
				in_communicator,
				use_spectral_flow,
				prod_wf_overlap,
				other_wf_overlap,
				use_rotation_invariance,
				staggered_subspace,
				save_basis,
				calculate_expectation_values;
  struct 		matrix_meta_data hamiltonian,
  				*ops_data;
  struct 		task_parameters tasks;
  struct 		momentum_information momentum_info;
  struct 		rotation_information *rotation_info; 
  struct 		adjacency_information adjacency_info;
  struct 		parity_information parity_info ;
  struct 		spectral_flow_info 	spectral_info;
  struct 		basis_information basis_info;
  struct 		prod_wf_information prod_wf_info;
  struct		expectation_values_information exp_vals_info;
  char 			op_dat_filename[MAX_STRING], 
				op_meta_filename[MAX_STRING], 
				output_prefix[MAX_STRING], 
				parity_masks[MAX_STRING], 
				method[MAX_STRING], 
				conservation_masks[MAX_STRING], 
				momentum_info_file[MAX_STRING],
				rotation_info_file[MAX_STRING],
				input_file[MAX_STRING],
				current_output_additional_prefix[MAX_STRING],
				prod_wf_file[MAX_STRING],
				other_wf_file[MAX_STRING],
				expectation_values_file[MAX_STRING]; 
  ModelTypeEnum model_type;
  MPI_Comm   	mat_solver_comm;
  
  //These are the function references. Functions they reference depend on the input parameters.
  int 			(*getrow)(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters, int nz );
  int 			(*getrowlocal)(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters,int nz,PetscInt Istart, PetscInt Iend );
  uPetscInt 	(*get_full_basis_index)(struct parameter_struct *parameters, uPetscInt reduced_index);
  uPetscInt 	(*get_reduced_basis_index)(struct parameter_struct *parameters, uPetscInt full_index,PetscTruth *exists);
  PetscTruth 	(*simple_basis_check)(struct parameter_struct *parameters, uPetscInt index); //this function pointer can be used for cases where there is a simple check to see if a basis element is in or out.
  PetscTruth 	(*representative_check)(struct parameter_struct *parameters, uPetscInt index); //this function pointer can be used for cases where only one basis out of a group of connected basis elements is used. TODO: lowest representative must also satisfy other checks. 
  uPetscInt  	(*group_representative)(struct parameter_struct *parameters, uPetscInt basis, PetscTruth phase_check, PetscScalar *phase,PetscTruth *group_exists); //TODO: lowest representative must also satisfy other checks. 
  PetscReal  	(*group_normal)(struct parameter_struct *parameters, uPetscInt index); //return the normal to use which is 1/sqrt(N) where n is the number of times the same basis appears in the group.
  uPetscInt 	(*translate_using_masks)(struct parameter_struct *parameters, uPetscInt basis, PetscInt *amount, PetscScalar *phase);
  uPetscInt 	(*rotate_using_masks)(struct parameter_struct *parameters, uPetscInt basis, int *amount, int rot_op_idx, PetscScalar *phase);
};


/*! \brief This is the structure that acts as the context structure when using matrix free methods. 
	this context is to contain enough information about the operator to perform the matrix vector multiplication.
	Still at the experimental stage and not scaling very well. */
struct shell_matrix_context {
  struct 		parameter_struct *parameters;
  struct 		matrix_meta_data *matrix_data;
  PetscScalar 	*buf1, *buf2 ; /*!< Buffers for sending and receiving partial results.*/
  int 			number_relevant_terms;
  int 			final_number_terms;
  //array of flip masks. One per term
  uPetscInt 	*flips; 
  int 			*sub_term_numbers;
  int 			**subterms;  /*!< Array of arrays of indexes of terms.*/
};

/*! \brief This structure is used for returning the calculated non zero values for a row of the matrix. Values are stored along with row and column index. */
struct matrix_row_values {
  uPetscInt 	row;
  PetscInt 		*cols;
#if defined(PETSC_USE_COMPLEX)
  PetscScalar 	*values;
#else	
  PetscReal 	*real_vals;
  PetscReal 	*imag_vals;
#endif	
  PetscInt val_count; 
};

/*! \brief Structure for storing the results from a single diagonalisation before they are written to disk. */
struct solver_results {
  PetscReal  	*eigenvalues, *errors, *overlaps, *other_overlaps;
  PetscScalar   ***correlations; //pointer to an array of matrices of correlation values. 
  PetscScalar 	*trial_correlations; //point to correlations for the trial wavefunction 
  PetscInt   	eigenvalues_converged, *degeneracies, degenerate_bands, iterations_taken;
  PetscReal 	time_taken;
};

#endif
