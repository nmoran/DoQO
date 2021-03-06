/*
 *  p_entanglement_entropy.h
 *  
 *
 *  Created by Niall Moran on 16/11/2010.
 *  Copyright 2010 All rights reserved.
 *
 */


#include "tinyxml.h"

#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <algorithm>

#define MAX_STRING 	200
#define TOLERANCE  	1e-30
#define NSV_LIMIT	200
#define ERROR_MARGIN	0.05
#define INPUT_VEC_TOLERANCE 1e-20

using namespace std;

typedef unsigned long int uPetscInt;
typedef long long int PetscInt;
typedef complex<double> PetscScalar;


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





