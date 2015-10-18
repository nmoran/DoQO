/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 *
 * Source file for rotation functions  
 * 
 * Version: $Id$ 
 * 
*/

#include "rotation.h"

#undef __FUNCT__
#define __FUNCT__ "parse_rotation_info"
int parse_rotation_info(struct parameter_struct *parameters,TiXmlHandle *docHandle){
	PetscErrorCode ierr;
	//section to read in product wavefunction over lap details
	TiXmlElement* element = docHandle->FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SYMMETRIES").FirstChild( "ROTATION").Element();
	if ( element)
	{	
		parameters->use_rotation_invariance = PETSC_TRUE;
		parameters->representative_check = &rotation_representative_check; 
		parameters->group_representative = &rotation_group_representative;
		parameters->group_normal = &rotation_group_normal; 
		parameters->rotation_number_relevant_sectors = -1;

		char buf2[MAX_STRING];
		if ( element->Attribute("file") == NULL ) {
			PetscPrintf(PETSC_COMM_WORLD,"A file must be supplied to describe the rotation when using rotational invariance.\n");
			return 1;
		}
		strcpy(parameters->rotation_info_file,element->Attribute("file"));
		
		element = docHandle->FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SYMMETRIES").FirstChild("ROTATION").FirstChild("RELEVANT_SECTORS").Element();
		if ( element ) {
			if ( element->Attribute("number") == NULL ) {
				PetscPrintf(PETSC_COMM_WORLD,"The number of relevant sectors must be specified.\n");
				return 1;
			}
			parameters->rotation_number_relevant_sectors = atoi(element->Attribute("number"));
			ierr = PetscMalloc(parameters->rotation_number_relevant_sectors*sizeof(PetscInt),&parameters->rotation_relevant_sectors);CHKERRQ(ierr);
			for ( int i = 0 ; i < parameters->rotation_number_relevant_sectors ; i++ ) {
				element = docHandle->FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SYMMETRIES").FirstChild("ROTATION").FirstChild("RELEVANT_SECTORS").Child("SECTOR",i).Element();
				if ( element ) {
					parameters->rotation_relevant_sectors[i] = atoi(element->GetText());
				} else {
					PetscPrintf(PETSC_COMM_WORLD,"Number of relevant rotation sectors incorrect.\n");
					return 1;
				}
			}
			if ( parameters->verbosity >= 3) {
				stringstream ss ;
				ss << "\tNumber of rotation sectors to use: " << parameters->rotation_number_relevant_sectors << " (" ; 
				PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
				for ( int i = 0 ; i < parameters->rotation_number_relevant_sectors ; i++ ) {
					PetscPrintf(PETSC_COMM_WORLD,"%d ", (int)parameters->rotation_relevant_sectors[i]);
				}
				PetscPrintf(PETSC_COMM_WORLD,")\n");
			}
		}	
	}	
	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "read_rotation_information"
/*! \brief This function reads the input file containing the details of the transformation being used.*/
int read_rotation_information(struct parameter_struct *parameters){
	string filename(parameters->rotation_info_file);
	PetscErrorCode ierr;
	
	//open rotation file as xml document.
	TiXmlDocument doc(filename);
	bool loadOkay = doc.LoadFile();
	if ( !loadOkay ) 
	{
		PetscPrintf(PETSC_COMM_WORLD,"Could not load rotation info file '%s'. Error='%s'. Exiting.\n", parameters->rotation_info_file, doc.ErrorDesc() );
		return 1;
	}
	TiXmlHandle docHandle(&doc);
	
	//find ROTATION_OP tag.
	
	parameters->number_rotation_sectors = 1;
	
	TiXmlElement* parent_element = docHandle.FirstChild("ROTATION_OPS").Element(); 	
	if ( parent_element) 
	  {
	    if (parent_element->Attribute("number") == NULL ) 
	      { 
		PetscPrintf(PETSC_COMM_WORLD,"You must specify the number of rotation operators.\n");
		return 1;
	      }
	    parameters->number_rotation_ops = atoi(parent_element->Attribute("number"));
	    
	    ierr = PetscMalloc(parameters->number_rotation_ops *sizeof(struct rotation_information), &parameters->rotation_info); CHKERRQ(ierr);
	    
	    for ( int rot_idx = 0 ; rot_idx < parameters->number_rotation_ops ; rot_idx ++ ) 
	      {
			TiXmlElement* element = docHandle.FirstChild("ROTATION_OPS").Child("ROTATION_OP", rot_idx).Element(); //tags surrounding the mapping. Has attribute "number" which contains the 
														  //number of times the rotation op can be applied before getting back to the start.
			TiXmlElement* element2;
			if ( element) 
			{	
				if (element->Attribute("number") == NULL ) { 
					PetscPrintf(PETSC_COMM_WORLD,"You must specify a number attribute to the rotation operation.\n");
					return 1;
				}

				parameters->rotation_info[rot_idx].num_sectors = atoi(element->Attribute("number"));
				parameters->rotation_info[rot_idx].num_relevant_sectors = parameters->rotation_info[rot_idx].num_sectors;
				parameters->number_rotation_sectors *= parameters->rotation_info[rot_idx].num_sectors;
				
				ierr = PetscMalloc( parameters->rotation_info[rot_idx].num_sectors * sizeof(PetscScalar),&parameters->rotation_info[rot_idx].phases);CHKERRQ(ierr);
				
				//allocate space for the mapping masks.
				ierr = PetscMalloc( parameters->rotation_info[rot_idx].num_sectors * sizeof(uPetscInt*),&parameters->rotation_info[rot_idx].mapping_masks);CHKERRQ(ierr);
				for ( int i = 0 ; i < parameters->rotation_info[rot_idx].num_sectors ; i++ ) {
					ierr = PetscMalloc(parameters->number_particles * sizeof(uPetscInt), &parameters->rotation_info[rot_idx].mapping_masks[i]);CHKERRQ(ierr);
				}
				
				//applying 0 rotations will give the same stage so add masks that do not move any particles.
				for ( int i = 0 ; i < parameters->number_particles ; i++ ) {
					parameters->rotation_info[rot_idx].mapping_masks[0][i] = ((uPetscInt)1) << i;
				}
				
				//read in the mappings to position 1 which is a single rotation with the op.
				for ( int i = 0 ; i < parameters->number_particles ; i++ ) {
					element2 = docHandle.FirstChild("ROTATION_OPS").Child("ROTATION_OP", rot_idx).Child("MAPPING",i).Element();
					if ( element2 ) {
						parameters->rotation_info[rot_idx].mapping_masks[1][atoi(element2->Attribute("to"))-1] = ((uPetscInt)1) << (atoi(element2->Attribute("from"))-1);
					} else {
						PetscPrintf(PETSC_COMM_WORLD,"Mapping %d does not exist.\n",i);
						return 1;
					}
				}
				
				//Now iterate this rotation to fill in the masks for the other possible mutiples of the rotation.
				for ( int i = 2 ; i < parameters->rotation_info[rot_idx].num_sectors ; i++ ) {
					for ( int j = 0 ; j < parameters->number_particles ; j++ ) {
						parameters->rotation_info[rot_idx].mapping_masks[i][j] = parameters->rotation_info[rot_idx].mapping_masks[1][log2_uPetscInt(parameters->rotation_info[rot_idx].mapping_masks[i-1][j])];
					}
				}
												
				//print out the mapping for each number of rotations that can be applied. For debugging purposes. 
				/*for (int i = 0 ; i < parameters->rotation_info.num_sectors ; i++ ) {
					PetscPrintf(PETSC_COMM_WORLD,"Applying %d rotations:\n",i);
					for ( int j = 0 ; j < parameters->number_particles ; j++ ) {
						PetscPrintf(PETSC_COMM_WORLD,"Site %d moved to site %d.\n",j,(int)log2_uPetscInt(parameters->rotation_info.mapping_masks[i][j]));
					}
				}*/
				
				parameters->rotation_info[rot_idx].sector_modulo = 1;
				parameters->rotation_info[rot_idx].sector_divisor = 1;
				for ( int j = 0 ; j <= rot_idx ; j++ )
				{
					parameters->rotation_info[rot_idx].sector_modulo *= parameters->rotation_info[j].num_sectors;
					if ( j <  rot_idx) parameters->rotation_info[rot_idx].sector_divisor *= parameters->rotation_info[j].num_sectors;
				}
				
				
			} else {
				PetscPrintf(PETSC_COMM_WORLD,"No ROTATION_OP %d tag in rotation info file.\n",rot_idx);
				return 1;
			}
	      }
	  } 
	else 
	  {
	    PetscPrintf(PETSC_COMM_WORLD,"No ROTATION_OPS tag in rotation info file.\n");
	    return 1;
	  }
	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "apply_rotation"
uPetscInt apply_rotation(struct parameter_struct *parameters, uPetscInt basis, int *number, int rot_op_idx)
{
	if ( parameters->model_type != SPIN_ONE ) {
		if ( *number > 0 ) {
			uPetscInt new_basis;
			struct rotation_information *rot_info;
			rot_info = &parameters->rotation_info[rot_op_idx];	
			new_basis = 0;
			for ( int spin = 0 ; spin < parameters->number_particles ; spin++ ){
				if ( (uPetscInt)(basis & rot_info->mapping_masks[*number][spin]) > (uPetscInt)0 ) {
					new_basis += ((uPetscInt)1) << spin;
				}
			}
			return new_basis;
		} else {
			return basis;
		}
	} else {
		if ( *number > 0 ) {
			uPetscInt new_basis;
			struct rotation_information *rot_info;
			rot_info = &parameters->rotation_info[rot_op_idx];	
			new_basis = 0;
			for ( int spin = 0 ; spin < parameters->number_particles ; spin++ ){
				if ( (uPetscInt)(basis & rot_info->mapping_masks[*number][spin]) > (uPetscInt)0 ) {
					new_basis += ((uPetscInt)1) << spin;
				}
				if ( (uPetscInt)(basis & (rot_info->mapping_masks[*number][spin] << (parameters->number_particles))) > (uPetscInt)0 ) {
					new_basis += ((uPetscInt)1) << (spin + parameters->number_particles);
				}
			}
			return new_basis;
		} else {
			return basis;
		}
	  
	}
} 

#undef __FUNCT__
#define __FUNCT__ "apply_rotation"
uPetscInt apply_rotation(struct parameter_struct *parameters, uPetscInt basis, int *number, int rot_op_idx, PetscScalar *phase)
{
	if ( parameters->model_type != SPIN_ONE ) {	  
		*phase = 1.0;
		if ( number[0] > 0 ) {
			list<int>  creation_ops;
			int from;
			uPetscInt new_basis;
			struct rotation_information *rot_info;
			rot_info = &parameters->rotation_info[rot_op_idx];	
			new_basis = 0;
			for ( int spin = 0 ; spin < parameters->number_particles ; spin++ ){
				if ( (uPetscInt)(basis & rot_info->mapping_masks[number[0]][spin]) > (uPetscInt)0 ) {
					new_basis += ((uPetscInt)1) << spin;
					from = (int)log2((double)rot_info->mapping_masks[number[0]][spin]);
					list<int>::iterator it  = creation_ops.begin(); //going to create a list of the creation operators. Sorted in desending order.
					if ( creation_ops.begin() == creation_ops.end() || from > *(it)) { //its empty or belongs at the beginning we just insert.
						creation_ops.insert(it,from);
					} else { //otherwise have to pass it through to its proper place					
						it  = creation_ops.begin(); //start at the beginning. 
						while ( from < *(it) && it != creation_ops.end() ) { //keep going forward till we find the suitable place or till we come to the end.
							if ( parameters->model_type == FERMIONIC ) (*phase) *= -1.0; //flip the sign due to the anticommutation relation.
							it++;  //move onto the next op. 
						}					 
						creation_ops.insert(it,from); //insert before the position of the iterator which should now be at the correct location.
					}				
				}
			}
			return new_basis;
		} else {
			return basis;
		}
	} else {	  
		*phase = 1.0;
		if ( number[0] > 0 ) {
			list<int>  creation_ops;
			int from;
			uPetscInt new_basis;
			struct rotation_information *rot_info;
			rot_info = &parameters->rotation_info[rot_op_idx];	
			new_basis = 0;
			for ( int spin = 0 ; spin < parameters->number_particles ; spin++ ){
				if ( (uPetscInt)(basis & rot_info->mapping_masks[number[0]][spin]) > (uPetscInt)0 ) {
					new_basis += ((uPetscInt)1) << spin;
					from = (int)log2((double)rot_info->mapping_masks[number[0]][spin]);
					list<int>::iterator it  = creation_ops.begin(); //going to create a list of the creation operators. Sorted in desending order.
					if ( creation_ops.begin() == creation_ops.end() || from > *(it)) { //its empty or belongs at the beginning we just insert.
						creation_ops.insert(it,from);
					} else { //otherwise have to pass it through to its proper place					
						it  = creation_ops.begin(); //start at the beginning. 
						while ( from < *(it) && it != creation_ops.end() ) { //keep going forward till we find the suitable place or till we come to the end.
							if ( parameters->model_type == FERMIONIC ) (*phase) *= -1.0; //flip the sign due to the anticommutation relation.
							it++;  //move onto the next op. 
						}					 
						creation_ops.insert(it,from); //insert before the position of the iterator which should now be at the correct location.
					}				
				}
				if ( (uPetscInt)(basis & (rot_info->mapping_masks[*number][spin] << (parameters->number_particles))) > (uPetscInt)0 ) {
					new_basis += ((uPetscInt)1) << (spin + parameters->number_particles);
				}
			}
			return new_basis;
		} else {
			return basis;
		}
	}
}


#undef __FUNCT__
#define __FUNCT__ "get_rotation_phase"
/*! \brief Gets the phase associated with the specified rotation.
 	Given the amount to rotate and the sector we are currently in this function returns the phases
	associated with the rotation for a given rotation sector. This is usually a complex number and so a complex build of PETSc should be used. */ 
PetscScalar get_rotation_phase(struct parameter_struct *parameters,PetscInt sector, int rot_op_idx)
{
	PetscScalar phase;
	PetscScalar k;
	struct rotation_information *rot_info;
	rot_info = &parameters->rotation_info[rot_op_idx];
    // This feature requires a complex build so only adding directive so that DoQO will compile but this feature should not be used without complex support.
    #if defined(PETSC_USE_COMPLEX)
	k = ((PetscScalar)sector)*2.0*PETSC_PI*((PetscScalar)rot_info->real_sector)/(PetscScalar)rot_info->num_sectors*PETSC_i;
    #endif
	phase = PetscExpScalar(k);
	return phase;
}


#undef __FUNCT__
#define __FUNCT__ "calculate_rotation_phases"
/*! \brief Calculate and store all possible rotation phases so that they can just be used later without calculation. */
void calculate_rotation_phases(struct parameter_struct *parameters, int rot_op_idx)
{
	int number_rotations = parameters->rotation_info[rot_op_idx].num_sectors;
	
	//allocate space to store phases.
	//PetscMalloc(number_rotations*sizeof(PetscScalar),&phases);//allocate space to store phases
	
	//go through translations and set the phase for each. 
	for ( PetscInt i = 0 ; i < number_rotations ; i++ ) {
		parameters->rotation_info[rot_op_idx].phases[i] = get_rotation_phase(parameters,i, rot_op_idx);
	}
}

#undef __FUNCT__
#define __FUNCT__ "rotation_group_representative"
/*! \brief Function returns the configuration with the lowest index which is related via the allowed rotation to the given configuration. */
uPetscInt  rotation_group_representative(struct parameter_struct *parameters, uPetscInt basis, PetscTruth phase_check, PetscScalar *phase, PetscTruth *group_exists) 
{
	struct rotation_information *rot_info;
	uPetscInt new_basis, lowest_basis;
	int amount;
	PetscScalar phase_total = 0 , phase_tmp;
	PetscTruth exists;
	PetscScalar sign;
	
	*group_exists = PETSC_TRUE;
	lowest_basis = basis;
	*phase = 1.0;
	
	PetscInt total_sectors = 1;
	for ( int rot_op_idx = 0 ; rot_op_idx < parameters->number_rotation_ops ; rot_op_idx++ ) 
	{
		total_sectors *= parameters->rotation_info[rot_op_idx].num_sectors;		
	}
	
	for ( PetscInt sector = 0 ; sector < total_sectors ; sector++ ) 
	{	
		new_basis = basis;
		phase_tmp = 1.0;
		for ( int rot_op_idx = 0 ; rot_op_idx < parameters->number_rotation_ops ; rot_op_idx++ )
		{
			rot_info = &parameters->rotation_info[rot_op_idx];
			amount = (sector % rot_info->sector_modulo) / rot_info->sector_divisor;
			//new_basis = translate(parameters,basis,amount);
			new_basis = (*parameters->rotate_using_masks)(parameters,new_basis,&amount,rot_op_idx, &sign);
			phase_tmp *= sign * rot_info->phases[amount];
		}
		//PetscPrintf(PETSC_COMM_WORLD,"%d: New basis %d, sector %d, %d\n", basis,new_basis,sector[0], sector[1] );
		(*parameters->get_reduced_basis_index)(parameters,new_basis,&exists);
		if ( exists == PETSC_TRUE ) {
			if ( new_basis == lowest_basis ) { //must be equal
				phase_total += phase_tmp;
			} //end if 
			if (new_basis > lowest_basis ) { 
				continue;
			} else if ( new_basis < lowest_basis ) {
				lowest_basis = new_basis;
				phase_total = phase_tmp;
			} 
		}
	}	

	
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
#define __FUNCT__ "rotation_group_normal"
/*! \brief Function returns the normal for a given configuration. The normal is the square root of the number of times the same configuration is encoutered via rotations. 
	The configuration 0001 would have normal sqrt(1) and the configuration 0101 would have normal sqrt(2) under rotation by one site left or right. */
PetscReal  rotation_group_normal(struct parameter_struct *parameters, uPetscInt index) 
{
	PetscReal normal = 0.0;
	struct rotation_information *rot_info;
	uPetscInt new_basis;
	int amount;
	PetscScalar sign;
	
	PetscInt total_sectors = 1;
	for ( int rot_op_idx = 0 ; rot_op_idx < parameters->number_rotation_ops ; rot_op_idx++ ) 
	{
		total_sectors *= parameters->rotation_info[rot_op_idx].num_sectors;		
	}
	
	for ( PetscInt sector = 0 ; sector < total_sectors ; sector++ ) 
	{	
		new_basis = index;
		for (int rot_op_idx = 0 ; rot_op_idx < parameters->number_rotation_ops ; rot_op_idx++ )
		{
			rot_info = &parameters->rotation_info[rot_op_idx];
			amount = (sector % rot_info->sector_modulo) / rot_info->sector_divisor;
			new_basis = (*parameters->rotate_using_masks)(parameters,new_basis,&amount,rot_op_idx,&sign);
		}
		if	( new_basis == index ) {
				normal += 1.0;
		}
	}
	
	normal = sqrt(normal);
	return normal;
}

#undef __FUNCT__
#define __FUNCT__ "rotation_representative_check"
/*! \brief Check if this particular index is in the current rotation sector. */
PetscTruth rotation_representative_check(struct parameter_struct *parameters, uPetscInt index)
{
	struct rotation_information *rot_info;
	uPetscInt new_basis;
	int amount;
	PetscScalar phase_total;
	PetscScalar phase_tmp;
	PetscTruth exists;
	PetscScalar sign;
	
	PetscInt total_sectors = 1;
	for ( int rot_op_idx = 0 ; rot_op_idx < parameters->number_rotation_ops ; rot_op_idx++ ) 
	{
		total_sectors *= parameters->rotation_info[rot_op_idx].num_sectors;		
	}
	
	phase_total = 0.0;
	
	for ( PetscInt sector = 0 ; sector < total_sectors ; sector++ ) 
	{	
		new_basis = index; 
		phase_tmp = 1.0;
		for (int  rot_op_idx = 0 ; rot_op_idx < parameters->number_rotation_ops ; rot_op_idx++ )
		{
			rot_info = &parameters->rotation_info[rot_op_idx];
			amount = (sector % rot_info->sector_modulo) / rot_info->sector_divisor;
			new_basis = (*parameters->rotate_using_masks)(parameters,new_basis,&amount,rot_op_idx,&sign);
			phase_tmp *= rot_info->phases[amount] * sign;
		}
		(*parameters->get_reduced_basis_index)(parameters,new_basis,&exists);
		if ( exists == PETSC_TRUE ) {
			if (new_basis < index ) { 
				return PETSC_FALSE;
			} else if ( new_basis == index ) {
				phase_total += phase_tmp;
			}
		}
	}
	
	if ( PetscAbsScalar(phase_total) < parameters->phase_tol ) {
		return PETSC_FALSE;
	}
	
	return PETSC_TRUE;
}


void set_rotation_sector(struct parameter_struct *parameters, PetscInt sector)
{
	PetscInt modulo = 1, divisor = 1;
	for ( int rot_op_idx = 0 ; rot_op_idx < parameters->number_rotation_ops; rot_op_idx ++ ) 
	{
		struct rotation_information *rot_info;
		rot_info = &parameters->rotation_info[rot_op_idx];
		modulo *= rot_info->num_sectors;
		rot_info->real_sector = (sector % modulo) / divisor;
		divisor *= rot_info->num_sectors;
	}
}

void deallocate_rotation_memory(struct parameter_struct *parameters)
{
	for ( int rot_op_idx = 0 ; rot_op_idx < parameters->number_rotation_ops ; rot_op_idx++ )
	{
		for ( int i = 0 ; i < parameters->rotation_info[rot_op_idx].num_sectors ; i++ ) 
		{
			PetscFree(parameters->rotation_info[rot_op_idx].mapping_masks[i]);
		}
		PetscFree(parameters->rotation_info[rot_op_idx].mapping_masks);
		PetscFree(parameters->rotation_info[rot_op_idx].phases);
	}
	if ( parameters->rotation_number_relevant_sectors > 0 ) 
		PetscFree(parameters->rotation_relevant_sectors);
	
	PetscFree(parameters->rotation_info);
}
