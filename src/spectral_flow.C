/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 *
 * Functions that implement spectral flow features. 
 * 
 * Version: $Id$ 
 * 
*/

#include "spectral_flow.h"

/*! \brief Function that sets up variables used for spectral flow calculation.*/
int setup_spectral_flow_info(struct parameter_struct *parameters,struct matrix_meta_data *matrix_data) {
#if defined(PETSC_USE_COMPLEX)
	struct spectral_flow_info *spectral_info;
	spectral_info = &(parameters->spectral_info);
	
	spectral_info->real_current_point[0] = spectral_info->current_point[0];
	spectral_info->real_current_point[1] = spectral_info->current_point[1];
	
	//work out which point we are using.
	if ( spectral_info->number_relevant_points[0] < spectral_info->number_points[0] ) {
		spectral_info->real_current_point[0] = spectral_info->relevant_points[0][spectral_info->current_point[0]];
	}
	if ( spectral_info->number_relevant_points[1] < spectral_info->number_points[1] ) {
		spectral_info->real_current_point[1] = spectral_info->relevant_points[1][spectral_info->current_point[1]];
	}
	
	//calculatte the value of alpha and the exponentiated versions.
	spectral_info->alpha[0] = (PetscReal)spectral_info->real_current_point[0] * 1.0/(PetscReal)(spectral_info->number_points[0]-1);
	spectral_info->alpha[1] = (PetscReal)spectral_info->real_current_point[1] * 1.0/(PetscReal)( (spectral_info->number_points[1] == 1) ? 1 : spectral_info->number_points[1]-1);
	spectral_info->e_alpha[0] = PetscExpScalar(2.0*PETSC_i*PETSC_PI*spectral_info->alpha[0]);
	spectral_info->e_alpha[1] = PetscExpScalar(2.0*PETSC_i*PETSC_PI*spectral_info->alpha[1]);
	spectral_info->e_alpha_c[0] = PetscConj(spectral_info->e_alpha[0]);
	spectral_info->e_alpha_c[1] = PetscConj(spectral_info->e_alpha[1]);
	spectral_info->e_alpha_f_l[0] = PetscExpScalar(2.0*PETSC_i*PETSC_PI*spectral_info->alpha[0]*(PetscReal)parameters->real_conservation_sector/(PetscReal)parameters->number_particles);
	spectral_info->e_alpha_f_l[1] = PetscExpScalar(2.0*PETSC_i*PETSC_PI*spectral_info->alpha[1]*(PetscReal)parameters->real_conservation_sector/(PetscReal)parameters->number_particles);
	if (parameters->verbosity >= 10 ) {
		stringstream ss; 
		ss << "Alpha[0]: " << spectral_info->alpha[0] << endl;
		ss << "Alpha[1]: " << spectral_info->alpha[1] << endl;
		ss << "e^(2pi alpha[0] i): " << spectral_info->e_alpha[0] << endl;
		ss << "e^(2pi alpha[1] i): " << spectral_info->e_alpha[1] << endl;
		ss << "e^(-2pi alpha[0] i): " << spectral_info->e_alpha_c[0] << endl;
		ss << "e^(-2pi alpha[1] i): " << spectral_info->e_alpha_c[1] << endl;
		ss << "e^(2pi alpha[0] i F/L): " << spectral_info->e_alpha_f_l[0] << endl;
		ss << "e^(2pi alpha[1] i F/L): " << spectral_info->e_alpha_f_l[1] << endl;
		PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
	}

	//find the index of the parameters that are on the edge in the list of parameters.
	int alpha_param_idx[2], alpha_conj_param_idx[2];
	alpha_param_idx[0] = -1;alpha_param_idx[1] = -1;alpha_conj_param_idx[0] = -1;alpha_conj_param_idx[1] = -1;
	for ( int i = 0 ; i < parameters->tasks.number_parameters ; i++ ) {
		if ( strcmp(spectral_info->alpha_param[0],parameters->tasks.parameter_names[i].c_str()) == 0 ) {
			alpha_param_idx[0] = i;
		} else if ( strcmp(spectral_info->alpha_conj_param[0],parameters->tasks.parameter_names[i].c_str()) == 0 ) {
			alpha_conj_param_idx[0] = i;
		}	
		if ( strcmp(spectral_info->alpha_param[1],parameters->tasks.parameter_names[i].c_str()) == 0 ) {
			alpha_param_idx[1] = i;
		} else if ( strcmp(spectral_info->alpha_conj_param[1],parameters->tasks.parameter_names[i].c_str()) == 0 ) {
			alpha_conj_param_idx[1] = i;
		}	
	}
	
	//now find which terms have this parameter and work out the multiplier again using the e_alpha and e_alpha_c as values for these.
	for ( int i = 0 ; i < matrix_data->number_terms ; i++ ) {
		for ( int j = 0 ; j < matrix_data->number_term_params[i] ; j++ ) {
			if ( matrix_data->term_parameters[i][j] == alpha_param_idx[0] ) { //check if alpha_param is one of the coefficients of the ith term. 
				matrix_data->term_multipliers[i] = spectral_info->e_alpha[0];	
				for ( int k = 0 ; k < matrix_data->number_term_params[i] ; k++ ) {
					if (matrix_data->term_parameters[i][j] != alpha_param_idx[0] ) { 
						matrix_data->term_multipliers[i] *=  parameters->tasks.task_parameter_values[parameters->current_task][matrix_data->term_parameters[i][k]] ;
					}
				}
			}
			if ( matrix_data->term_parameters[i][j] == alpha_conj_param_idx[0] ) { //check if alpha_conj_param is one of the coefficients of the ith term. 
				matrix_data->term_multipliers[i] = spectral_info->e_alpha_c[0];	
				for ( int k = 0 ; k < matrix_data->number_term_params[i] ; k++ ) {
					if (matrix_data->term_parameters[i][j] != alpha_conj_param_idx[0] ) { 
						matrix_data->term_multipliers[i] *=  parameters->tasks.task_parameter_values[parameters->current_task][matrix_data->term_parameters[i][k]] ;
					}
				}
			}
			if ( matrix_data->term_parameters[i][j] == alpha_param_idx[1] ) { //check if alpha_param is one of the coefficients of the ith term. 
				matrix_data->term_multipliers[i] = spectral_info->e_alpha[1];	
				for ( int k = 0 ; k < matrix_data->number_term_params[i] ; k++ ) {
					if (matrix_data->term_parameters[i][j] != alpha_param_idx[1] ) { 
						matrix_data->term_multipliers[i] *=  parameters->tasks.task_parameter_values[parameters->current_task][matrix_data->term_parameters[i][k]] ;
					}
				}
			}
			if ( matrix_data->term_parameters[i][j] == alpha_conj_param_idx[1] ) { //check if alpha_conj_param is one of the coefficients of the ith term. 
				matrix_data->term_multipliers[i] = spectral_info->e_alpha_c[1];	
				for ( int k = 0 ; k < matrix_data->number_term_params[i] ; k++ ) {
					if (matrix_data->term_parameters[i][j] != alpha_conj_param_idx[1] ) { 
						matrix_data->term_multipliers[i] *=  parameters->tasks.task_parameter_values[parameters->current_task][matrix_data->term_parameters[i][k]] ;
					}
				}
			}
		}
	}	
	
	//work out the eigenvalues of the translation operators in each direction for this eigenspace. Note that we are only looking at the spectal flow in the 0 direction for now. 
	spectral_info->eigenvalues_exponent[0] = 2.0*PETSC_PI*PETSC_i*(PetscScalar)(  (PetscReal)parameters->momentum_info.current_sector[0]/(PetscReal)parameters->momentum_info.num_sectors[0] + (spectral_info->alpha[0]*(PetscScalar)parameters->real_conservation_sector)/(PetscReal)parameters->momentum_info.num_sectors[0]);
	spectral_info->eigenvalues_exponent[1] = 2.0*PETSC_PI*PETSC_i*(PetscScalar)(  (PetscReal)parameters->momentum_info.current_sector[1]/(PetscReal)parameters->momentum_info.num_sectors[1] + (spectral_info->alpha[1]*(PetscScalar)parameters->real_conservation_sector)/(PetscReal)parameters->momentum_info.num_sectors[1]);
	spectral_info->eigenvalues[0] = PetscExpScalar(spectral_info->eigenvalues_exponent[0]);
	spectral_info->eigenvalues[1] = PetscExpScalar(spectral_info->eigenvalues_exponent[1]);
	if (parameters->verbosity >= 10 ) {
		stringstream ss; 
		ss << "Eigenvalue[0]: " << spectral_info->eigenvalues[0] << endl;
		ss << "Eigenvalue[1]: " << spectral_info->eigenvalues[1] << endl;
		ss << "Eigenvalue exponent[0] pi: " << spectral_info->eigenvalues_exponent[0]/PETSC_PI << endl;
		ss << "Eigenvalue exponent[1] pi: " << spectral_info->eigenvalues_exponent[1]/PETSC_PI << endl;

		PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
	}

	
#endif

	return 0;
}	


#undef __FUNCT__
#define __FUNCT__ "translate_using_masks_fermions_spectral"
/*! \brief Function which uses translation masks to perform the translation faster. Function is meant for fermiions when using spectral flow feature 
and keeps track of the sign and phase in 'phase' argument. */
uPetscInt translate_using_masks_fermions_spectral(struct parameter_struct *parameters, uPetscInt basis, PetscInt *amount, PetscScalar *phase){
	uPetscInt new_basis;
	struct momentum_information *mom_info;
	
	mom_info = &parameters->momentum_info;
	
	*phase = 1.0; //initialize the phase to one
	list<int>  creation_ops;  //list to store creation operators 
	int from;
	//PetscPrintf(PETSC_COMM_WORLD,"Translate by %d, %d\n", translate_vector[0], translate_vector[1] ) ;
	new_basis = 0; 
	for ( int spin = 0 ; spin < parameters->number_particles ; spin++ ){ //loop over all the particles. 
		if ( (uPetscInt)(basis & mom_info->translation_masks[amount[0]][amount[1]][spin]) > (uPetscInt)0 ) { //if the basis has a fermion in the position that the mask has it in then move the fermion to the positon of spin.
			new_basis += ((uPetscInt)1) << spin; //this adds a fermion in the position indicated by 'spin'.
			from = (int)log2((double)mom_info->translation_masks[amount[0]][amount[1]][spin]); //this is the position it was moved from. 
			list<int>::iterator it  = creation_ops.begin(); //going to populate a list of the creation operators. Sorted in desending order.
			if ( creation_ops.begin() == creation_ops.end() || from > *(it)) { //its empty or belongs at the beginning we just insert.
				creation_ops.insert(it,from);
			} else { //otherwise have to pass it through to its proper place
				it  = creation_ops.begin(); //start at the beginning. 
				while ( from < *(it) && it != creation_ops.end() ) { //keep going forward till we find the suitable place or till we come to the end.
					(*phase) *= -1.0; //flip the sign due to the anticommutation relation.
					it++;  //move onto the op. 
				}
				creation_ops.insert(it,from); //insert before the position of the iterator which should now be at the correct location.
			}
			//note this only works for rectangular lattice geomentries at the moment. 
			//if a fermion has been moved over the boundary then we multiply by the phase e^(2pi i alpha). 
			if ( from % mom_info->lattice_dimensions[0] > spin % mom_info->lattice_dimensions[0]) {
				(*phase) *= parameters->spectral_info.e_alpha[0];
			}
			if ( from / mom_info->lattice_dimensions[0] > spin / mom_info->lattice_dimensions[0]) {
				(*phase) *= parameters->spectral_info.e_alpha[1];
			}
		}
	}
	return new_basis;
} 



#undef __FUNCT__
#define __FUNCT__ "momentum_group_normal_spectral"
/*! \brief Function returns the normal for a given configuration. The normal is the square root of the number of times the same configuration is encoutered via translations. 
	The configuration 0001 would have normal sqrt(1) and the configuration 0101 would have normal sqrt(2) under translation by one site left or right. */
PetscReal  momentum_group_normal_spectral(struct parameter_struct *parameters, uPetscInt index) {
	PetscScalar normal = 0.0;
	struct momentum_information *mom_info;
	uPetscInt new_basis;
	PetscInt amount[2];
	PetscScalar sign;
	
	mom_info = &parameters->momentum_info;
	for ( amount[0] = 0 ; amount[0] < mom_info->num_sectors[0] ; amount[0]++ ) {
		for ( amount[1] = 0 ; amount[1] < mom_info->num_sectors[1] ; amount[1]++ ) {
			new_basis = (*parameters->translate_using_masks)(parameters,index,amount,&sign);		
			if	( new_basis == index ) {
				//normal += 1.0;
				normal += sign*parameters->momentum_info.phases[amount[1]*parameters->momentum_info.lattice_dimensions[0]+amount[0]];
			}
		}	
	}	
	normal = PetscSqrtScalar(normal);
	return PetscRealPart(normal);
}





