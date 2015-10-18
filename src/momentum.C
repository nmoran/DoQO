/* 
 * Copyright (C) 2010    Niall Moran, Graham Kells, Jiri Vala
 * 
 * 
*/

#include "momentum.h"

using namespace std;


#undef __FUNCT__
#define __FUNCT__ "parse_momentum_file"
/*! \brief Function that reads the information about the translation invariance in the lattice.
 	To provide information about the nature of the lattice and the translations possible an information file is provided. This function parses that file and populates the momentum information
 	structure with this information.
 */
int parse_momentum_file(struct parameter_struct *parameters) {
	FILE *fp;
	struct momentum_information *mom_info;
	
	mom_info = &parameters->momentum_info;
	
	if (parameters->use_momentum_sectors == PETSC_TRUE ) {
		fp = fopen(parameters->momentum_info_file,"r");
		if ( fp == NULL ) return -1;
#if defined(PETSC_USE_64BIT_INDICES)
		if ( fscanf(fp,"LATTICE_VECTOR1 = %lld,%lld\n",&mom_info->lattice_vector[0][0], &mom_info->lattice_vector[0][1]) != 2 ) return -1;
		if ( fscanf(fp,"LATTICE_VECTOR2 = %lld,%lld\n",&mom_info->lattice_vector[1][0], &mom_info->lattice_vector[1][1]) != 2 ) return -1;
		if ( fscanf(fp,"NORM_VECTOR1 = %lf,%lf\n",&mom_info->normalised_vector[0][0], &mom_info->normalised_vector[0][1])!= 2 ) return -1;
		if ( fscanf(fp,"NORM_VECTOR2 = %lf,%lf\n",&mom_info->normalised_vector[1][0], &mom_info->normalised_vector[1][1]) != 2 ) return -1;
		if ( fscanf(fp,"LATTICE_DIMENSIONS = %lld,%lld",&mom_info->lattice_dimensions[0], &mom_info->lattice_dimensions[1]) != 2 ) return -1;
#else
		if ( fscanf(fp,"LATTICE_VECTOR1 = %d,%d\n",&mom_info->lattice_vector[0][0], &mom_info->lattice_vector[0][1]) != 2 ) return -1;
		if ( fscanf(fp,"LATTICE_VECTOR2 = %d,%d\n",&mom_info->lattice_vector[1][0], &mom_info->lattice_vector[1][1]) != 2 ) return -1;
		if ( fscanf(fp,"NORM_VECTOR1 = %lf,%lf\n",&mom_info->normalised_vector[0][0], &mom_info->normalised_vector[0][1])!= 2 ) return -1;
		if ( fscanf(fp,"NORM_VECTOR2 = %lf,%lf\n",&mom_info->normalised_vector[1][0], &mom_info->normalised_vector[1][1]) != 2 ) return -1;
		if ( fscanf(fp,"LATTICE_DIMENSIONS = %d,%d",&mom_info->lattice_dimensions[0], &mom_info->lattice_dimensions[1]) != 2 ) return -1;
#endif
		fclose(fp);
	}
	
	//now that we have the translation vectors and the lattice dimensions we can work out the number of sectors we have in each direction.
	mom_info->num_sectors[0] = my_min((mom_info->lattice_vector[0][0] == 0 ? mom_info->lattice_dimensions[0]*mom_info->lattice_dimensions[1] : mom_info->lattice_dimensions[0]/mom_info->lattice_vector[0][0]), 
	(mom_info->lattice_vector[0][1] == 0 ? mom_info->lattice_dimensions[0]*mom_info->lattice_dimensions[1]  : mom_info->lattice_dimensions[1]/mom_info->lattice_vector[0][1]));
	mom_info->num_sectors[1] = my_min((mom_info->lattice_vector[1][0] == 0 ? mom_info->lattice_dimensions[0]*mom_info->lattice_dimensions[1]  : mom_info->lattice_dimensions[0]/mom_info->lattice_vector[1][0]), 
	(mom_info->lattice_vector[1][1] == 0 ? mom_info->lattice_dimensions[0]*mom_info->lattice_dimensions[1]  : mom_info->lattice_dimensions[1]/mom_info->lattice_vector[1][1]));

	
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "print_momentum_information"
/*! \brief Function that simply prints the momentum information read from the input file.*/
void print_momentum_information(struct parameter_struct *parameters) {
	if ( parameters->verbosity > 2 ) { 
		PetscPrintf(PETSC_COMM_WORLD,"Momentum Information.\n");
#if defined(PETSC_USE_64BIT_INDICES)
		PetscPrintf(PETSC_COMM_WORLD,"Vector 1: %lld, %lld\n",parameters->momentum_info.lattice_vector[0][0],parameters->momentum_info.lattice_vector[0][1]); 
		PetscPrintf(PETSC_COMM_WORLD,"Vector 2: %lld, %lld\n",parameters->momentum_info.lattice_vector[1][0],parameters->momentum_info.lattice_vector[1][1]);
		PetscPrintf(PETSC_COMM_WORLD,"Norm Vector 1: %lf, %lf\n",parameters->momentum_info.normalised_vector[0][0],parameters->momentum_info.normalised_vector[0][1]);
		PetscPrintf(PETSC_COMM_WORLD,"Norm Vector 2: %lf, %lf\n",parameters->momentum_info.normalised_vector[1][0],parameters->momentum_info.normalised_vector[1][1]);
		PetscPrintf(PETSC_COMM_WORLD,"Lattice dimensions: %lld, %lld\n",parameters->momentum_info.lattice_dimensions[0],parameters->momentum_info.lattice_dimensions[1]);
		PetscPrintf(PETSC_COMM_WORLD,"Sector sizes: %lld, %lld\n",parameters->momentum_info.num_sectors[0],parameters->momentum_info.num_sectors[1]);
#else		
		PetscPrintf(PETSC_COMM_WORLD,"Vector 1: %d, %d\n",parameters->momentum_info.lattice_vector[0][0],parameters->momentum_info.lattice_vector[0][1]); 
		PetscPrintf(PETSC_COMM_WORLD,"Vector 2: %d, %d\n",parameters->momentum_info.lattice_vector[1][0],parameters->momentum_info.lattice_vector[1][1]);
		PetscPrintf(PETSC_COMM_WORLD,"Norm Vector 1: %lf, %lf\n",parameters->momentum_info.normalised_vector[0][0],parameters->momentum_info.normalised_vector[0][1]);
		PetscPrintf(PETSC_COMM_WORLD,"Norm Vector 2: %lf, %lf\n",parameters->momentum_info.normalised_vector[1][0],parameters->momentum_info.normalised_vector[1][1]);
		PetscPrintf(PETSC_COMM_WORLD,"Lattice dimensions: %d, %d\n",parameters->momentum_info.lattice_dimensions[0],parameters->momentum_info.lattice_dimensions[1]);
		PetscPrintf(PETSC_COMM_WORLD,"Sector sizes: %d, %d\n",parameters->momentum_info.num_sectors[0],parameters->momentum_info.num_sectors[1]);
#endif
	}
}

#undef __FUNCT__
#define __FUNCT__ "translate"
/*! \brief Function which translates a given configuration by the specified about according to the translation vectors being used.*/
uPetscInt translate(struct parameter_struct *parameters, uPetscInt basis, PetscInt *amount){
	uPetscInt new_basis;
	uPetscInt translate_vector[2];
	int x, y , new_x, new_y;
	struct momentum_information *mom_info;
	
	mom_info = &parameters->momentum_info;
	
	translate_vector[0] = amount[0] * mom_info->lattice_vector[0][0] ;
	translate_vector[0] += amount[1] * mom_info->lattice_vector[1][0] ;
	translate_vector[1] = amount[0] * mom_info->lattice_vector[0][1] ;
	translate_vector[1] += amount[1] * mom_info->lattice_vector[1][1] ;
	
	//PetscPrintf(PETSC_COMM_WORLD,"Translate by %d, %d\n", translate_vector[0], translate_vector[1] ) ;
	new_basis = 0;
	for ( x = 0 ; x < mom_info->lattice_dimensions[0] ; x++ ){
		for ( y = 0 ; y < mom_info->lattice_dimensions[1] ; y++ ){
			if ( (basis & (uPetscInt)pow(2.0,(double)(y*mom_info->lattice_dimensions[0]+x))) > (uPetscInt)0 ) {
				new_x = (x + translate_vector[0]) % mom_info->lattice_dimensions[0];
				new_y = (y + translate_vector[1]) % mom_info->lattice_dimensions[1];
				new_basis += (uPetscInt)pow(2.0,(double)(new_y*mom_info->lattice_dimensions[0] + new_x));
			} 
		}
	}
	return new_basis;
} 

#undef __FUNCT__
#define __FUNCT__ "translate_using_masks_spins"
/*! \brief Function which uses translation masks to perform the translation faster. Function is meant for spin systems. */
uPetscInt translate_using_masks_spins(struct parameter_struct *parameters, uPetscInt basis, PetscInt *amount, PetscScalar *phase){
	uPetscInt new_basis;
	struct momentum_information *mom_info;
	
	mom_info = &parameters->momentum_info;
	
	*phase = 1.0;
	//PetscPrintf(PETSC_COMM_WORLD,"Translate by %d, %d\n", translate_vector[0], translate_vector[1] ) ;
	new_basis = 0;
	for ( int spin = 0 ; spin < parameters->number_particles ; spin++ ){
		if ( (uPetscInt)(basis & mom_info->translation_masks[amount[0]][amount[1]][spin]) > (uPetscInt)0 ) {
			new_basis += ((uPetscInt)1) << spin;
		}
	}
	return new_basis;
} 

#undef __FUNCT__
#define __FUNCT__ "translate_using_masks_fermions"
/*! \brief Function which uses translation masks to perform the translation faster. Function is meant for fermiions and keeps track of the sign in 'sign'. */
uPetscInt translate_using_masks_fermions(struct parameter_struct *parameters, uPetscInt basis, PetscInt *amount, PetscScalar *phase){
	uPetscInt new_basis;
	struct momentum_information *mom_info;
	
	mom_info = &parameters->momentum_info;
	
	*phase = 1.0;
	list<int>  creation_ops;
	int from;
	//PetscPrintf(PETSC_COMM_WORLD,"Translate by %d, %d\n", translate_vector[0], translate_vector[1] ) ;
	new_basis = 0;
	for ( int spin = 0 ; spin < parameters->number_particles ; spin++ ){
		if ( (uPetscInt)(basis & mom_info->translation_masks[amount[0]][amount[1]][spin]) > (uPetscInt)0 ) {
			new_basis += ((uPetscInt)1) << spin;
			from = (int)log2((double)mom_info->translation_masks[amount[0]][amount[1]][spin]);
			list<int>::iterator it  = creation_ops.begin(); //going to create a list of the creation operators. Sorted in desending order.
			if ( creation_ops.begin() == creation_ops.end() || from > *(it)) { //its empty or belongs at the beginning we just insert.
				creation_ops.insert(it,from);
			} else { //otherwise have to pass it through to its proper place
				it  = creation_ops.begin(); //start at the beginning. 
				while ( from < *(it) && it != creation_ops.end() ) { //keep going forward till we find the suitable place or till we come to the end.
					(*phase) *= -1.0; //flip the sign due to the anticommutation relation.
					it++;  //move onto the next op. 
				}
				creation_ops.insert(it,from); //insert before the position of the iterator which should now be at the correct location.
			}
		}
	}
	return new_basis;
} 

#undef __FUNCT__
#define __FUNCT__ "calculate_translation_masks"
/*! \brief  This function calculates the translation masks for each translation to make the translation operation more efficient. */
int calculate_translation_masks(struct parameter_struct *parameters){
	//first allocate space for the translations
	PetscErrorCode ierr;
	int amount[2];
	struct momentum_information *mom_info;
	mom_info = &parameters->momentum_info;
	
	ierr = PetscMalloc(parameters->momentum_info.num_sectors[0]*sizeof(uPetscInt**),&(parameters->momentum_info.translation_masks));CHKERRQ(ierr);
	for ( int i = 0 ; i < parameters->momentum_info.num_sectors[0] ; i++ ) {
		ierr = PetscMalloc(parameters->momentum_info.num_sectors[1]*sizeof(uPetscInt*),&(parameters->momentum_info.translation_masks[i]));CHKERRQ(ierr);
		for ( int j = 0 ; j < parameters->momentum_info.num_sectors[1] ; j++ ) {
			ierr = PetscMalloc(parameters->number_particles*sizeof(uPetscInt),&(parameters->momentum_info.translation_masks[i][j]));CHKERRQ(ierr);
		}
	}
	
	uPetscInt translate_vector[2];
	for ( amount[0] = 0 ; amount[0] < parameters->momentum_info.num_sectors[0]; amount[0]++){
		for ( amount[1] = 0 ; amount[1] < parameters->momentum_info.num_sectors[1]; amount[1]++){
			//we do a reverse translation to get where the final value of each comes from.
			translate_vector[0] =  amount[0] * mom_info->lattice_vector[0][0] ;
			translate_vector[0] += amount[1] * mom_info->lattice_vector[1][0] ;
			translate_vector[1] =  amount[0] * mom_info->lattice_vector[0][1] ;
			translate_vector[1] += amount[1] * mom_info->lattice_vector[1][1] ;
			
			//PetscPrintf(PETSC_COMM_WORLD,"Translate by %d, %d\n", translate_vector[0], translate_vector[1] ) ;
			for ( int spin = 0 ; spin < parameters->number_particles ; spin++ ) {
				int x = spin % mom_info->lattice_dimensions[0];
				int y = spin / mom_info->lattice_dimensions[0];
				int new_x = ((x + mom_info->lattice_dimensions[0]) - translate_vector[0]) % mom_info->lattice_dimensions[0];
				int new_y = ((y + mom_info->lattice_dimensions[1]) - translate_vector[1]) % mom_info->lattice_dimensions[1];
				mom_info->translation_masks[amount[0]][amount[1]][spin] = ((uPetscInt)( 1 ) << (new_y*mom_info->lattice_dimensions[0] + new_x));
			}			
		}	
	}		

	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "get_translation_phase"
/*! \brief Gets the phase associated with the specified translation.
 	Given the amount to translate along each of the translation vectors and the sector we are currently in this function returns the phases
	associated with the translation for a given momentum sector. This is usually a complex number and so a complex build of PETSc should be used. */ 
PetscScalar get_translation_phase(struct parameter_struct *parameters,PetscInt *amount, PetscInt *sector){
	PetscScalar phase;

	if ( parameters->use_spectral_flow == PETSC_FALSE) {
		PetscReal k[2];
		PetscReal theta;
		PetscReal r[2];
		struct momentum_information *mom_info;
		
		mom_info = &parameters->momentum_info;
		
		k[0] = (PetscReal)sector[0]*2.0*PETSC_PI/(PetscReal)mom_info->num_sectors[0];
		k[1] = (PetscReal)sector[1]*2.0*PETSC_PI/(PetscReal)mom_info->num_sectors[1];
		
		r[0] = (PetscReal)amount[0] * mom_info->normalised_vector[0][0] ;
		r[0] += (PetscReal)amount[1] * mom_info->normalised_vector[1][0] ;
		r[1] = (PetscReal)amount[0] * mom_info->normalised_vector[0][1] ;
		r[1] += (PetscReal)amount[1] * mom_info->normalised_vector[1][1] ;
		
		theta = r[0] * k[0] + r[1] * k[1];
		
	
		phase = cos(theta);
	#if defined(PETSC_USE_COMPLEX)
		phase += sin(theta)*PETSC_i;
	#endif
	} else {
		phase = PetscExpScalar(-parameters->spectral_info.eigenvalues_exponent[0]*(PetscScalar)amount[0]);
		phase *= PetscExpScalar(-parameters->spectral_info.eigenvalues_exponent[1]*(PetscScalar)amount[1]);
	}
//	PetscPrintf(PETSC_COMM_WORLD,"Phase calc: sector %d, %d, amount %d %d, momentum: %lf , %lf ,R: %lf, %lf,  theta %lf, phase %lf + i %lf \n",sector[0],sector[1],amount[0],amount[1],k[0],k[1],r[0],r[1],theta,PetscRealPart(phase),PetscImaginaryPart(phase));
	
	return phase;
}

#undef __FUNCT__
#define __FUNCT__ "calculate_translation_phases"
/*! \brief Calculate and store all possible translation phases so that they can just be used later without calculation. */
PetscScalar* calculate_translation_phases(struct parameter_struct *parameters){
	PetscScalar *phases;
	int number_translations = parameters->momentum_info.num_sectors[0]*parameters->momentum_info.num_sectors[1];
	PetscInt amount[2];
	
	//allocate space to store phases.
	PetscMalloc(number_translations*sizeof(PetscScalar),&phases);//allocate space to store phases
	
	//go through translations and set the phase for each. 
	for ( int i = 0 ; i < number_translations ; i++ ) {
		amount[0] = i / parameters->momentum_info.num_sectors[1];
		amount[1] = i % parameters->momentum_info.num_sectors[1];
		phases[i] = get_translation_phase(parameters,amount,parameters->momentum_info.current_sector);
	}
	
	return phases;
}

#undef __FUNCT__
#define __FUNCT__ "is_representative"
/*! \brief This function loops over all the translations and checks to see if the given basis element is the lowest of the connected ones. If it isnt false is returned. If 
	it is the normal is calculated and returned. If phase_check is true it also checks to make sure the phase does not dissappear.*/
PetscTruth is_representative(struct parameter_struct *parameters, uPetscInt basis, PetscReal *normal, PetscTruth phase_check, PetscInt *sector){
	struct momentum_information *mom_info;
	uPetscInt new_basis;
	PetscInt amount[2];
	PetscScalar phase_total;
	PetscScalar phase_tmp;
	PetscScalar sign;
	
	phase_total = 0.0;
	*normal = 0.0;
	mom_info = &parameters->momentum_info;
	for ( amount[0] = 0 ; amount[0] < mom_info->num_sectors[0] ; amount[0]++ ) {
		for ( amount[1] = 0 ; amount[1] < mom_info->num_sectors[1] ; amount[1]++ ) {
			//new_basis = translate(parameters,basis,amount);
			new_basis = (*parameters->translate_using_masks)(parameters,basis,amount,&sign);
//			PetscPrintf(PETSC_COMM_WORLD,"%d: New basis %d, sector %d, %d\n", basis,new_basis,sector[0], sector[1] );
			if (new_basis < basis ) { 
				return PETSC_FALSE;
			} else if ( new_basis == basis ) {
				(*normal) += 1.0;
				if ( phase_check == PETSC_TRUE) {
					//phase_tmp = get_translation_phase(parameters,amount,sector);
					//phase_total += get_translation_phase(parameters,amount,sector);
					phase_total += sign * parameters->momentum_info.phases[amount[0]*parameters->momentum_info.num_sectors[1]+ amount[1]];
				}
//				if ( amount[0] != 0 || amount[1] != 0 ) PetscPrintf(PETSC_COMM_WORLD,"Match: basis %d and new basis : %d and amount is %d, %d\n",basis,new_basis,amount[0],amount[1]);
			}
		}	
	}
	
	//PetscPrintf(PETSC_COMM_WORLD,"Phase %lf + %lf with abs %.8e\n",PetscRealPart(phase_total),PetscImaginaryPart(phase_total),PetscAbsScalar(phase_total));
	
	if ( phase_check == PETSC_TRUE) {
		if ( PetscAbsScalar(phase_total) < parameters->phase_tol ) {
			return PETSC_FALSE;
		}
	}
	
	*normal = sqrt(*normal);
	
	return PETSC_TRUE;
}

#undef __FUNCT__
#define __FUNCT__ "get_representative"
/*! \brief This function loops over all the translations and finds the representative of the current basis element. If phase check is true it also returns the phase 
	of the operation to bring the basis to its representative. */
uPetscInt get_representative(struct parameter_struct *parameters, uPetscInt basis, PetscTruth phase_check, PetscScalar *phase, PetscInt *sector){
	struct momentum_information *mom_info;
	uPetscInt new_basis, lowest_basis;
	PetscInt amount[2];
	PetscScalar phase_total = 0 ;
	PetscScalar sign;
	
	lowest_basis = basis;
	*phase = 1.0;
	mom_info = &parameters->momentum_info;
	for ( amount[0] = 0 ; amount[0] < mom_info->num_sectors[0] ; amount[0]++ ) {
		for ( amount[1] = 0 ; amount[1] < mom_info->num_sectors[1] ; amount[1]++ ) {
			//new_basis = translate(parameters,basis,amount);
			new_basis = (*parameters->translate_using_masks)(parameters,basis,amount,&sign);
//			PetscPrintf(PETSC_COMM_WORLD,"%d: New basis %d, sector %d, %d\n", basis,new_basis,sector[0], sector[1] );
			if ( new_basis == lowest_basis ) { //must be equal
				//if ( phase_check == PETSC_TRUE ) phase_total += get_translation_phase(parameters,amount,sector);
				if ( phase_check == PETSC_TRUE ) phase_total += sign* parameters->momentum_info.phases[amount[0]*parameters->momentum_info.num_sectors[1]+ amount[1]];
			} //end if 
			if (new_basis > lowest_basis ) { 
				continue;
			} else if ( new_basis < lowest_basis ) {
				lowest_basis = new_basis;
				if ( phase_check == PETSC_TRUE ) {
					//*phase = get_translation_phase(parameters,amount,sector);
					*phase = sign * parameters->momentum_info.phases[amount[0]*parameters->momentum_info.num_sectors[1]+ amount[1]];
				} //end if
			} 
		}	//end for amount[1]
	} //end for amount[0]
	
	if ( phase_check == PETSC_TRUE && lowest_basis < basis ) {  //in this case we should return the conjugate of the phase instead.
		*phase = PetscConj(*phase);
	}
	
	if ( phase_check == PETSC_TRUE && PetscAbsScalar(phase_total) < parameters->phase_tol ) return -1;	

	return lowest_basis;
}

#undef __FUNCT__
#define __FUNCT__ "test_momentum_tools"
/*! \brief Function used during debugging to exercise the other momentum functions. */
void test_momentum_tools(struct parameter_struct *parameters){
	struct momentum_information *mom_info;
	uPetscInt  basis;
	PetscReal normal;
	PetscInt sector[2];
	//PetscScalar phase;
	
	mom_info = &parameters->momentum_info;
	
	sector[0] = 0;
	sector[1] = 0;
	
	for ( sector[0] = 0 ; sector[0] < mom_info->num_sectors[0]; sector[0] ++ ) { 
		for ( sector[1] = 0 ; sector[1] < mom_info->num_sectors[1]; sector[1] ++ ) { 
			PetscPrintf(PETSC_COMM_WORLD,"Sector : %d %d\n-----------------------\n",sector[0],sector[1]);
			for ( basis = 0 ; basis < (uPetscInt)parameters->numberBasisStates ; basis++ ) {
//			for ( basis = 0 ; basis < 5 ; basis++ ) {
				if ( (is_representative(parameters,basis,&normal,PETSC_TRUE,sector)) == PETSC_TRUE ) {
					PetscPrintf(PETSC_COMM_WORLD,"%d : %lf\n", basis,normal);
				}
			}
		/*	for ( basis = 0 ; basis < parameters->hamiltonian.numberBasisStates ; basis++ ) {
//			for ( basis = 0 ; basis < 5 ; basis++ ) {
					PetscPrintf(PETSC_COMM_WORLD,"%d represented by %d\n", basis,get_representative(parameters,basis,PETSC_TRUE,&phase,sector));
			}*/
		}
	}
}


#undef __FUNCT__
#define __FUNCT__ "momentum_group_representative"
/*! \brief Function returns the configuration with the lowest index which is related via the allows translations to the given configuration. */
uPetscInt  momentum_group_representative(struct parameter_struct *parameters, uPetscInt basis, PetscTruth phase_check, PetscScalar *phase, PetscTruth *group_exists) {
	struct momentum_information *mom_info;
	uPetscInt new_basis, lowest_basis;
	PetscInt amount[2];
	PetscScalar phase_total = 0 ;
	PetscTruth exists;
	PetscScalar sign;
	
	*group_exists = PETSC_TRUE;
	lowest_basis = basis;
	*phase = 1.0;
	mom_info = &parameters->momentum_info;
	for ( amount[0] = 0 ; amount[0] < mom_info->num_sectors[0] ; amount[0]++ ) {
		for ( amount[1] = 0 ; amount[1] < mom_info->num_sectors[1] ; amount[1]++ ) {
			//new_basis = translate(parameters,basis,amount);
			new_basis = (*parameters->translate_using_masks)(parameters,basis,amount,&sign);
//			PetscPrintf(PETSC_COMM_WORLD,"%d: New basis %d, sector %d, %d\n", basis,new_basis,sector[0], sector[1] );
			(*parameters->get_reduced_basis_index)(parameters,new_basis,&exists);
			if ( exists == PETSC_TRUE ) {
				if ( new_basis == lowest_basis ) { //must be equal
					//if ( phase_check == PETSC_TRUE ) phase_total += get_translation_phase(parameters,amount,sector);
					phase_total += parameters->momentum_info.phases[amount[0]*parameters->momentum_info.num_sectors[1]+ amount[1]] * sign;
				} //end if 
				if (new_basis > lowest_basis ) { 
					continue;
				} else if ( new_basis < lowest_basis ) {
					lowest_basis = new_basis;
					phase_total = parameters->momentum_info.phases[amount[0]*parameters->momentum_info.num_sectors[1]+ amount[1]] * sign;
				} 
			}
		}	//end for amount[1]
	} //end for amount[0]
	
	if ( lowest_basis < basis ) {
		phase_total = PetscConj(phase_total);
	}
	
	if ( PetscAbsScalar(phase_total) > parameters->phase_tol ) {
		*phase = phase_total / sqrt(PetscRealPart(phase_total * PetscConj(phase_total)));
	} else phase_total = 0.0 ;
	
	if ( phase_check == PETSC_TRUE && PetscAbsScalar(phase_total) < parameters->phase_tol ) *group_exists = PETSC_FALSE;	

	return lowest_basis;
}

#undef __FUNCT__
#define __FUNCT__ "momentum_group_normal"
/*! \brief Function returns the normal for a given configuration. The normal is the square root of the number of times the same configuration is encoutered via translations. 
	The configuration 0001 would have normal sqrt(1) and the configuration 0101 would have normal sqrt(2) under translation by one site left or right. */
PetscReal  momentum_group_normal(struct parameter_struct *parameters, uPetscInt index) {
	PetscReal normal = 0.0;
	struct momentum_information *mom_info;
	uPetscInt new_basis;
	PetscInt amount[2];
	PetscScalar sign;
	
	mom_info = &parameters->momentum_info;
	for ( amount[0] = 0 ; amount[0] < mom_info->num_sectors[0] ; amount[0]++ ) {
		for ( amount[1] = 0 ; amount[1] < mom_info->num_sectors[1] ; amount[1]++ ) {
			new_basis = (*parameters->translate_using_masks)(parameters,index,amount,&sign);		
			if	( new_basis == index ) {
					normal += 1.0;
			}
		}	
	}	
	normal = sqrt(normal);
	return normal;
}

#undef __FUNCT__
#define __FUNCT__ "momentum_representative_check"
/*! \brief Check if this particular index is in the current momentum sector. */
PetscTruth momentum_representative_check(struct parameter_struct *parameters, uPetscInt index){
	struct momentum_information *mom_info;
	uPetscInt new_basis;
	PetscInt amount[2];
	PetscScalar phase_total;
	PetscScalar phase_tmp;
	PetscTruth exists;
	PetscScalar sign;
	
	phase_total = 0.0;
	mom_info = &parameters->momentum_info;
	for ( amount[0] = 0 ; amount[0] < mom_info->num_sectors[0] ; amount[0]++ ) {
		for ( amount[1] = 0 ; amount[1] < mom_info->num_sectors[1] ; amount[1]++ ) {
			//new_basis = translate(parameters,basis,amount);
			new_basis = (*parameters->translate_using_masks)(parameters,index,amount,&sign);
//			PetscPrintf(PETSC_COMM_WORLD,"%d: New basis %d, sector %d, %d\n", basis,new_basis,sector[0], sector[1] );
			(*parameters->get_reduced_basis_index)(parameters,new_basis,&exists);
			if ( exists == PETSC_TRUE ) {
				if (new_basis < index ) { 
					return PETSC_FALSE;
				} else if ( new_basis == index ) {
					//(*normal) += 1.0;
					//if ( phase_check == PETSC_TRUE) {
						//phase_tmp = get_translation_phase(parameters,amount,sector);
						//phase_total += get_translation_phase(parameters,amount,sector);
						phase_total += sign * parameters->momentum_info.phases[amount[0]*parameters->momentum_info.num_sectors[1]+amount[1]];
					//}
	//				if ( amount[0] != 0 || amount[1] != 0 ) PetscPrintf(PETSC_COMM_WORLD,"Match: basis %d and new basis : %d and amount is %d, %d\n",basis,new_basis,amount[0],amount[1]);
				}
			}
		}	
	}
	
	//PetscPrintf(PETSC_COMM_WORLD,"Phase %lf + %lf with abs %.8e\n",PetscRealPart(phase_total),PetscImaginaryPart(phase_total),PetscAbsScalar(phase_total));
	
	//if ( phase_check == PETSC_TRUE) {
	if ( PetscAbsScalar(phase_total) < parameters->phase_tol ) {
		return PETSC_FALSE;
	}
	//}
	
	//*normal = sqrt(*normal);
	
	return PETSC_TRUE;
}
