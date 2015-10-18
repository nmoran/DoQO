/*
 *  p_entanglement_entropy.h
 *  
 *
 *  Created by Niall Moran on 16/11/2010.
 *  Copyright 2010 All rights reserved.
 *
 */

#include "slepceps.h"
#include "slepcsvd.h"
#include "tinyxml.h"

#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>

#if defined(PETSC_USE_64BIT_INDICES)
typedef 	unsigned long long uPetscInt;
#define 	PETSC_MPI_INT MPI_LONG_LONG_INT
#define 	PETSC_UNSIGNED_MPI_INT MPI_UNSIGNED_LONG_LONG
#else
typedef 	unsigned int uPetscInt;
#define 	PETSC_MPI_INT MPI_INT
#define 	PETSC_UNSIGNED_MPI_INT MPI_UNSIGNED
#endif


#define MAX_STRING 	200
#define TOLERANCE  	1e-30
#define NSV_LIMIT	200
#define ERROR_MARGIN	0.05
#define INPUT_VEC_TOLERANCE 1e-20

#define 	log2_uPetscInt(num) ((uPetscInt)log2((double)num))

using namespace std;

typedef PetscBool PetscTruth;

class my_uPetscInt {
public: 
	uPetscInt num;
	static uPetscInt mask;
	
	bool operator==(const my_uPetscInt &other) const {
    	return (num & this->mask) == (other.num & this->mask);  
	}
	
	bool operator<(const my_uPetscInt &other) const {
    	return (num & this->mask) < (other.num & this->mask);  
	}
	
	bool operator>(const my_uPetscInt &other) const {
    	return (num & this->mask) > (other.num & this->mask);  
	}
};

//class to store the vector value along with the associated vector configuration.
class my_vec_entry { 
public: 
	uPetscInt mask;
	uPetscInt conf;
	PetscScalar val;
	
	bool operator==(const my_vec_entry &other) const {
    	return (conf & mask) == (other.conf & mask);  
	}
	
	bool operator<(const my_vec_entry &other) const {
    	return (conf & mask) < (other.conf & mask);  
	}
	
	bool operator>(const my_vec_entry &other) const {
    	return (conf & mask) > (other.conf & mask);  
	}
};


//class to store the vector value along with the associated vector configuration.
class my_vec_element { 
public: 
	uPetscInt basis_element;
	uPetscInt sorting_token;
	PetscScalar val;
	
	bool operator==(const my_vec_element &other) const {
    	return sorting_token == (other.sorting_token);  
	}
	
	bool operator<(const my_vec_element &other) const {
    	return sorting_token < (other.sorting_token);  
	}
	
	bool operator>(const my_vec_element &other) const {
    	return sorting_token > (other.sorting_token);  
	}
};

struct adjacency_information {
	int 		number_masks;
	uPetscInt 	*masks;
};

//info necessary for generating all configs in a list
struct conf_list_info {
	struct 		adjacency_information adj_info;
	PetscInt 	num_sites, filling, max_filling;
};


/*! \brief This structure stores the fields necessary to exploit rotational invariance. */
struct rotation_information { 
	PetscInt 	num_sectors;  //the number of sectors in each direction.
	int 		*relevant_sectors, num_relevant_sectors; //array of sectors to use
	PetscScalar *phases; //array which will store the phase corresponding to each rotation in each sector. Will save recalculating each time.
	uPetscInt 	**mapping_masks;
	PetscInt 	current_sector,real_sector;
	char		filename[MAX_STRING];
	
};

struct parameters {
	PetscMPIInt rank, size;
	uPetscInt *AB_confs,*A_confs,*B_confs;
	PetscInt sites, filling, A_filling, B_filling, l_A_confs, g_A_confs, l_B_confs, g_B_confs, l_AB_confs, g_AB_confs, A_momentum, B_momentum;
	PetscInt AB_sites, A_sites, B_sites, number_particles; 
	uPetscInt A_mask, B_mask, AB_mask; 
	PetscInt total_momentum, A_total_momentum, B_total_momentum;
	struct adjacency_information adj_info;
	Vec  psi; 
	Mat  W;
	PetscTruth nn_exclusion,ascii, rot_symmetry;
	bool spin_one;
	class my_vec_entry *my_vec_entries;	
	struct rotation_information rotation_info;
	//this function pointer can be used for cases where only one basis out of a group of connected basis elements is used. TODO: lowest representative must also satisfy other checks. 
	PetscTruth 	(*representative_check)(struct parameters *params, uPetscInt index); 
	PetscInt step;
};





