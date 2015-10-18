/* 
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala.
 *
 * Utility functions for implementing conservation of particle number or sz. 
 * 
 * Version: $Id$
 * 
 * 
*/


#include "conservation.h"



/*! \brief For a given sector and subset of sites described by the mask with index 'mask_number' the filling on those sites is returned.
 	For two masks 0 would correspond to both empty, 1: first containing one and second empty and so on.*/
uPetscInt get_mask_filling(struct parameter_struct *parameters, uPetscInt sector, int mask_number){
	struct parity_information *cons_info = &parameters->parity_info;
	uPetscInt filling;
	
	for ( int i = cons_info->number_masks-1 ; i >= mask_number ; i-- ) {
		if ( i > 0 ) { 
			uPetscInt sectors_at_mask = 1;
			for ( int j = 0 ; j < i ; j++ ) {
				sectors_at_mask *= (cons_info->counts[j] + 1);
			}
			filling = sector / sectors_at_mask ;  
			sector = sector % sectors_at_mask;
		} else {
			filling = sector;
		}
	}
	return filling; 
}

/*! \brief This function returns the number of valid configuration in the given filling sector.*/
uPetscInt get_size_of_conservation_sector(struct parameter_struct *parameters, uPetscInt sector){
	uPetscInt sector_size = 1;
	struct parity_information *cons_info = &parameters->parity_info; //reusing parity structure as same information required and 
	//do not envisage needing both parity and conservation at the same time.
	if ( parameters->model_type == SPIN_ONE ) {
		for ( int i = 0 ; i < cons_info->number_masks ; i++ ) {	 
			int filling = get_mask_filling(parameters,sector,i);
			int sites = cons_info->counts[i];
			int itr = 0;
			cons_info->permutations[i] = 0;
			while ( (filling - (itr*2)) >= 0 && (sites - itr) >= 0 ) {
			      cons_info->permutations[i] += NchooseC_safe(sites, itr) * NchooseC_safe(sites - itr, filling - (itr*2));
			      itr++;			      
			}			
			sector_size *= cons_info->permutations[i];
			if ( i > 0 ) {
				cons_info->permutations[i] *= cons_info->permutations[i-1];
			} 
		}
	} else {
		for ( int i = 0 ; i < cons_info->number_masks ; i++ ) {	 
			cons_info->permutations[i] = NchooseC(cons_info->counts[i],get_mask_filling(parameters,sector,i)); 
			sector_size *= cons_info->permutations[i];
			if ( i > 0 ) {
				cons_info->permutations[i] *= cons_info->permutations[i-1];
			} 
		}
	}
	
	return sector_size; 
}

/*! \brief This Function works out the amount of different fillings possible.*/
int search_for_conservation_sectors(struct matrix_meta_data *md, struct parameter_struct *parameters) {
	int i,j ;
	char buf[MAX_STRING];
	struct parity_information *cons_info;
	int count;
	PetscErrorCode ierr;
	
	//read_parity_masks(parameters->parity_masks, &(md->parity_info.masks), &(md->parity_info.number_masks));
	
	cons_info = &(parameters->parity_info);
	ierr = PetscMalloc(sizeof(uPetscInt) * cons_info->number_masks,&cons_info->permutations);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(int) * cons_info->number_masks,&cons_info->counts);CHKERRQ(ierr);
	
	cons_info->number_parity_sectors = 1;
	//count the number of 1s in the mask and store in counts array
	for ( i = 0 ; i < cons_info->number_masks ; i++ ) {
		dec2bin(cons_info->masks[i],parameters->number_particles, buf);
		count = 0;
		for ( j = 0 ; j < parameters->number_particles ; j++ ){
			if ( buf[j] == '1' ) count++;
		}
		cons_info->counts[i] = count;
		if ( parameters->model_type != SPIN_ONE ) {
			cons_info->number_parity_sectors *= (count+1);
		} else if ( parameters->model_type == SPIN_ONE ) {
			cons_info->number_parity_sectors *= ((2*count)+1);
		}
	}
	
	
	
	//PetscFPrintf(PETSC_COMM_WORLD,stderr,"%d parity sectors and masks %d.\n", cons_info->number_parity_sectors, cons_info->number_masks);
	/*for ( i = 0 ; i < cons_info->number_masks ; i++ ) {
#if defined(PETSC_USE_64BIT_INDICES)
		PetscFPrintf(PETSC_COMM_WORLD,stderr,"Mask %d is given by %llu.\n", i, cons_info->masks[i]);
#else
		PetscFPrintf(PETSC_COMM_WORLD,stderr,"Mask %d is given by %u.\n", i, cons_info->masks[i]);
#endif
	}*/
	return 0;
}

/*! \brief  Given an index in the reduced basis set and information the conservation sector this function will return the corresponding basis element from the full basis set. */
uPetscInt get_full_conservation_basis_index(struct parameter_struct *parameters, uPetscInt reduced_index) {
	uPetscInt full_index = 0;
	PetscInt sector = parameters->real_conservation_sector;	
	uPetscInt tmp;
	
	if ( parameters->model_type != SPIN_ONE )
	{
		for ( int i = 0 ; i < parameters->parity_info.number_masks ; i++ ) {
			uPetscInt mask_rank = reduced_index; 
			if ( i != (parameters->parity_info.number_masks - 1 ) ) {
				mask_rank %= parameters->parity_info.permutations[i+1];  
			}
			if ( i > 0 ) {
				mask_rank /= parameters->parity_info.permutations[i-1];
			}
			
			//PetscPrintf(PETSC_COMM_WORLD,"Getting n choose c for %d choose %d for rank %u.\n", parameters->parity_info.counts[i],(int)get_mask_filling(parameters,(uPetscInt)sector,i),mask_rank);
			tmp = get_index_n_choose_c(mask_rank,parameters->parity_info.counts[i],(int)get_mask_filling(parameters,(uPetscInt)sector,i));
			int local_count = 0 ;
			for ( int j = 0 ; j < parameters->number_particles ; j++ ) {
	//			if ( ( (uPetscInt)pow(2.0,(double)j) & parameters->parity_info.masks[i]) > 0 ) {
				if ( ( (((uPetscInt)1) << j) & parameters->parity_info.masks[i]) > 0 ) {
	//				if  ( (tmp & (uPetscInt)pow(2.0,(double)local_count)) > 0 ) {
					if  ( (tmp & (((uPetscInt)1) << local_count)) > 0 ) {
	//					full_index += (uPetscInt)pow(2.0,(double)j);
						full_index += ((uPetscInt)1) << j;
					}
					local_count += 1 ;
				}
			}
		}
	} else if ( parameters->model_type == SPIN_ONE ) 
	{
		for ( int i = 0 ; i < parameters->parity_info.number_masks ; i++ ) {
			uPetscInt mask_rank = reduced_index; 
			if ( i != (parameters->parity_info.number_masks - 1 ) ) {
				mask_rank %= parameters->parity_info.permutations[i+1];  
			}
			if ( i > 0 ) {
				mask_rank /= parameters->parity_info.permutations[i-1];
			}
			
			int filling = get_mask_filling(parameters,sector,i);
			int sites = parameters->parity_info.counts[i];
			int itr = 0;			
			uPetscInt accumulation ; 
			int pairs = 0;
			while ( (filling - (itr*2)) >= 0 && (sites - itr) >= 0  ) {
			      accumulation = NchooseC(sites, itr) * NchooseC(sites - itr, filling - (itr*2));
			      if ( accumulation <= mask_rank ) {
				      pairs++;
				      mask_rank -= accumulation; 
			      } else break;
			      itr++;			      
			}		
			
			/*if ( itr > 0 ) {
				 pairs = itr -1;
			} else pairs = 0;*/
			

			uPetscInt pair_rank = mask_rank / NchooseC(sites - pairs, filling - (pairs*2));			
			uPetscInt pair_pos = get_index_n_choose_c(pair_rank,sites,pairs);		
			tmp = pair_pos << parameters->number_particles;
			
			uPetscInt single_rank = mask_rank % NchooseC(sites - pairs, filling - (pairs*2));												
			uPetscInt single_pos = get_index_n_choose_c(single_rank,sites-pairs,filling - (pairs*2));
			
			int local_count = 0 ;
			for ( int j = 0 ; j < sites ; j++ ) {	
				if ( ( (((uPetscInt)1) << j) & pair_pos) == (uPetscInt)0 ) {
					if  ( (single_pos & (((uPetscInt)1) << local_count)) == (uPetscInt)0 ) {	
						tmp += ((uPetscInt)1) << j;
					}
					local_count += 1 ;
				}								  
			}			
			
			local_count = 0 ;
			for ( int j = 0 ; j < parameters->number_particles ; j++ ) {	
				if ( ( (((uPetscInt)1) << j) & parameters->parity_info.masks[i]) > 0 ) {	
					if  ( (tmp & (((uPetscInt)1) << local_count)) > 0 ) {	
						full_index += ((uPetscInt)1) << j;
					}
					local_count += 1 ;
				}
			}
			
			local_count = 0 ;
			for ( int j = parameters->number_particles ; j < (parameters->number_particles << 1) ; j++ ) {	
				if ( ( (((uPetscInt)1) << (j-parameters->number_particles)) & parameters->parity_info.masks[i]) > 0 ) {	
					if  ( (tmp & (((uPetscInt)1) << (local_count+parameters->number_particles))) > 0 ) {	
						full_index += ((uPetscInt)1) << j;
					}
					local_count += 1 ;
				}
			}
		}	  
	}
	return full_index;
}

/*! \brief  This function works out the index of the configuration and says whether it is valid or not via the exists flag. */
uPetscInt get_reduced_conservation_basis_index(struct parameter_struct *parameters, uPetscInt full_index, PetscTruth *exists) { //TODO enforce conservation is conserved and same with parity. Check full index and return -1 if it is not part of the current sector. Done
	uPetscInt reduced_index = 0 ;
	PetscInt sector = parameters->real_conservation_sector;
	*exists = PETSC_TRUE;
	
	if ( parameters->model_type != SPIN_ONE ) 
	{	
		for ( int i = 0 ; i < parameters->parity_info.number_masks ; i++ ) {
			uPetscInt mask_index = 0, mask_rank ;
			int local_count = 0 ;
			int filling = 0;
			for ( int j = 0 ; j < parameters->number_particles ; j++ ) {
				if ( ( (uPetscInt)pow(2.0,(double)j) & parameters->parity_info.masks[i]) > 0 ) {
					if ( (full_index & (uPetscInt)pow(2.0,(double)j)) > 0 ) {
						mask_index += (uPetscInt)pow(2.0,(double)local_count);
						filling++;
					}
					local_count++;
				}
			}
			if ( filling != (int)get_mask_filling(parameters,(uPetscInt)sector,i) ) {
				*exists = PETSC_FALSE;
				return -1; //returning an unsigned int so this doesnt make sense.
			}
			mask_rank = get_rank_n_choose_c(mask_index,parameters->parity_info.counts[i],(int)get_mask_filling(parameters,(uPetscInt)sector,i));
			if ( i > 0 ) {
				reduced_index += mask_rank * parameters->parity_info.permutations[i-1];
			} else {
				reduced_index += mask_rank; 
			}
		}
	} else if ( parameters->model_type == SPIN_ONE ) //for now only works with a single conservation mask
	{
		if ( parameters->parity_info.number_masks > 1 ) PetscFPrintf(PETSC_COMM_WORLD,stderr,"Spin one systems and more than one filling mask is not currently supported.\n");
		for ( int i = 0 ; i < parameters->parity_info.number_masks ; i++ ) {
			uPetscInt mask_index = 0, mask_rank = 0;
			int local_count = 0 ;
			int filling = 0;
			int pair_filling = 0;
			int single_filling = 0;
			for ( int j = 0 ; j < (parameters->number_particles << 1) ; j++ ) {
				//if ( ( (uPetscInt)pow(2.0,(double)j) & parameters->parity_info.masks[i]) > 0 ) {
					if ( (full_index & (uPetscInt)pow(2.0,(double)j)) > 0 ) {
						mask_index += (uPetscInt)pow(2.0,(double)local_count);
						if ( j >= parameters->number_particles ) { 
							filling+=2;
							pair_filling++;
						}
					} else {
						if ( j < parameters->number_particles ) {
							if ( (((uPetscInt)1 << (j+parameters->number_particles)) & full_index) == (uPetscInt)0 ) {
								filling++;
								single_filling++;
							}									
						} 
					}
					local_count++;
				//}
			}
			if ( filling != (int)get_mask_filling(parameters,(uPetscInt)sector,i) ) {
				*exists = PETSC_FALSE;
				return -1; //returning an unsigned int so this doesnt make sense.
			}
			
			//take into account configurations with less pairs than current config.
			for ( int j = 0 ; j < pair_filling ; j++ ) {
				mask_rank += NchooseC(parameters->parity_info.counts[i],j)*NchooseC(parameters->parity_info.counts[i]-j,filling-2*j);
			}
			uPetscInt pair_pos = mask_index >> parameters->number_particles;
			mask_rank += get_rank_n_choose_c(pair_pos, parameters->parity_info.counts[i],pair_filling) * NchooseC(parameters->parity_info.counts[i]-pair_filling, filling-2*pair_filling);
			
			uPetscInt single_pos = (uPetscInt)0;
			local_count = 0;
			for ( int j = 0 ; j < parameters->parity_info.counts[i] ; j++ ) {
				if ( (pair_pos & ( (uPetscInt)1 << j )) == (uPetscInt)0 ) { // if there is no pair in this position then its a valid position for single excitation.
					if ( (mask_index & ( (uPetscInt)1 << j )) == (uPetscInt)0 ) { //if sz is -1 here then set bit in string
						  single_pos += (uPetscInt)1 << local_count;						  
					}
					local_count++;
				}
			}
			
			mask_rank += get_rank_n_choose_c(single_pos,parameters->parity_info.counts[i] - pair_filling, filling - 2*pair_filling);
			
			if ( i > 0 ) {
				reduced_index += mask_rank * parameters->parity_info.permutations[i-1];
			} else {
				reduced_index += mask_rank; 
			}
		}  
	  
	}
	
	return reduced_index;
}

/*! \brief  This function is a dummy function for when no conservation of filling has been enabled. */
uPetscInt get_full_basis_index_passthru(struct parameter_struct *parameters, uPetscInt reduced_index){
	if ( parameters->model_type != SPIN_ONE ) {
		return reduced_index; //dummy function that just passthru the value as reduced basis is the same as full basis.
	} else  {
		uPetscInt full_index = 0;
		int spin_idx = parameters->number_particles;
		while ( spin_idx > 0 ) {
			if ( reduced_index >= 2*((uPetscInt)pow(3.0, (double)spin_idx-1)) ) {
				  reduced_index -= 2*((uPetscInt)pow(3.0, (double)spin_idx-1));
				  full_index += ((uPetscInt)1) << (spin_idx - 1 + parameters->number_particles);
			} else if ( reduced_index >= ((uPetscInt)pow(3.0, (double)spin_idx-1)) ) {
				  reduced_index -= ((uPetscInt)pow(3.0, (double)spin_idx-1));
			} else {
				  full_index += ((uPetscInt)1) << (spin_idx-1);
			}
			spin_idx--;
		}		
		return full_index;
	}
}

/*! \brief  This function is a dummy function for when no conservation of filling has been enabled. */
uPetscInt get_reduced_basis_index_passthru(struct parameter_struct *parameters, uPetscInt full_index,PetscTruth *exists){
	if ( parameters->model_type != SPIN_ONE ) {
		*exists = PETSC_TRUE;
		return full_index; //dummy function that just passthru the value as reduced basis is the same as full basis.
	} else  {
		uPetscInt reduced_index = 0;
		int spin_idx = parameters->number_particles;
		while ( spin_idx > 0 ) {
			if ( (full_index & (((uPetscInt)1) << (spin_idx-1+parameters->number_particles))) != ((uPetscInt)0) ) {
				  reduced_index += 2*((uPetscInt)pow(3.0, (double)spin_idx-1));
			} else if ( (full_index & (((uPetscInt)1) << (spin_idx-1))) == ((uPetscInt)0) ) {
				  reduced_index += ((uPetscInt)pow(3.0, (double)spin_idx-1));
			} 
			spin_idx--;
		}	
		*exists = PETSC_TRUE;
		return reduced_index;
	}
}

