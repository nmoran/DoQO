/* 
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 *
 * Functions for analysing and working out symmetries. 
 * 
 * Version: $Id$
 * 
 * 
*/


#include "symmetry.h"



#undef __FUNCT__
#define __FUNCT__ "read_parity_masks"
/*! \brief Read parity mask from a file named filename. These are one per line binary strings.*/
int read_parity_masks(char *filename, uPetscInt **masks, int *number) {
	FILE *fp;
	int c;
	char buf[MAX_STRING];
	int  bufidx, idx;
	PetscErrorCode ierr;
	
	fp = fopen(filename, "r" ) ;
	
	if ( fp == NULL ) { 
		PetscPrintf(PETSC_COMM_SELF, "Error opening file %s.\n", filename);
		return -1;
	}
	
	*number = 0;
	//open file and first count the lines
	int char_count = 0;
	while( (c = fgetc(fp)) != EOF ) {
		if ( c == '\n' && char_count > 0 ) {
			(*number) = (*number) + 1;
			char_count = 0;
		} else if ( c != ' ')  {
			char_count++; 
		}
	}
	if (char_count > 0 ) (*number) = (*number) + 1;
	fclose(fp);
	//PetscFPrintf(PETSC_COMM_WORLD,stderr,"%d masks.\n", *number);
	
	//allocate space for the masks
	//*masks = (uPetscInt *) malloc(sizeof(uPetscInt) * (*number));
	ierr = PetscMalloc(sizeof(uPetscInt) * (*number),masks);CHKERRQ(ierr);
	
	//this time populate the array
	fp = fopen(filename, "r" ) ;
	
	if ( fp == NULL ) { 
		PetscPrintf(PETSC_COMM_SELF, "Error opening file %s.\n", filename);
		return -1;
	}
	
	idx = 0;
	bufidx = 0;
	char_count = 0;
	while( (c = fgetc(fp)) != EOF ) {
		if ( c == '0' || c == '1' ) {
			buf[bufidx++] = c;
			char_count++;
		} 
		if ( c == '\n' && char_count > 0 ) {
			buf[bufidx] = '\0';

/*#if defined(PETSC_USE_64BIT_INDICES)
			PetscFPrintf(PETSC_COMM_WORLD,stderr,"Another mask %s (%llu) set by index %d.\n" , buf, bin2dec(buf),idx);
#else
			PetscFPrintf(PETSC_COMM_WORLD,stderr,"Another mask %s (%u) set by index %d.\n" , buf, bin2dec(buf),idx);
#endif*/
			(*masks)[idx] = bin2dec(buf);
			idx++;
			bufidx = 0;
			char_count = 0 ;
		}
	}
	if ( char_count > 0 ) {
		buf[bufidx] = '\0';
		(*masks)[idx] = bin2dec(buf);
		idx++;
		bufidx = 0;
		char_count = 0 ;
	}
	fclose(fp);
	
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "search_for_parity_sectors"
/*! \brief Process the parity information read from file and work out number of sites in each subset of sites over which parity applies.*/
int search_for_parity_sectors(struct matrix_meta_data *md, struct parameter_struct *parameters) {
	int i,j ;
	char buf[MAX_STRING];
	struct parity_information *parity_info;
	int count;
	PetscErrorCode ierr;
	
	//read_parity_masks(parameters->parity_masks, &(md->parity_info.masks), &(md->parity_info.number_masks));
	
	parity_info = &(parameters->parity_info);
	ierr = PetscMalloc(sizeof(uPetscInt) * parameters->parity_info.number_masks,&parity_info->permutations);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(int) * parity_info->number_masks,&parity_info->counts);CHKERRQ(ierr);
	
	//count the number of 1s in the mask and store in counts array
	for ( i = 0 ; i < parity_info->number_masks ; i++ ) {
		dec2bin(parity_info->masks[i],parameters->number_particles, buf);
		count = 0;
		for ( j = 0 ; j < parameters->number_particles ; j++ ){
			if ( buf[j] == '1' ) count++;
		}
		parity_info->counts[i] = count; 
	}
	
	
	//PetscFPrintf(PETSC_COMM_WORLD,stderr,"%d parity sectors and masks %d.\n", parity_info->number_parity_sectors, parity_info->number_masks);
	for ( i = 0 ; i < parity_info->number_masks ; i++ ) {
/*#if defined(PETSC_USE_64BIT_INDICES)
		PetscFPrintf(PETSC_COMM_WORLD,stderr,"Mask %d is given by %llu.\n", i, parity_info->masks[i]);
#else
		PetscFPrintf(PETSC_COMM_WORLD,stderr,"Mask %d is given by %u.\n", i, parity_info->masks[i]);
#endif*/
	}
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "get_size_of_parity_sector"
/*! \brief Work out the size of a given parity sector.*/
uPetscInt get_size_of_parity_sector(struct parameter_struct *parameters, uPetscInt sector){
	uPetscInt sector_size = 1;
	struct parity_information *parity_info;
	int parity;
	
	parity_info = &parameters->parity_info;
	// first task is to find the number of basis states in the sector. This is gotten from counting the permutations
	// in each parity sector and multiplying
	for ( int i = 0 ; i < parity_info->number_masks ; i++ ) {
		parity = get_mask_parity((int)sector,i);
		//PetscPrintf(PETSC_COMM_WORLD,"Permutations of mask %d with parity %d : %d\n",i, parity,num_permutations_at_parity(parameters,parity_info->masks[i], parity));
		parity_info->permutations[i] = sector_size; 
		sector_size *= num_permutations_at_parity(parameters,parity_info->masks[i], parity);
	}

	return sector_size;
}


#undef __FUNCT__
#define __FUNCT__ "get_mask_parity"
/*! \brief Returns the parity that each mask will have in each sector. For two masks 0 would correspond to even even, 1: even odd, 2: odd, even, 3: odd odd. */
int get_mask_parity(int sector, int mask_number){
	uPetscInt tmp_mask;
	
	tmp_mask = (uPetscInt)pow((double)2,(double)mask_number);
	if ( (((uPetscInt)sector) & tmp_mask)  > 0 ) {
		return 1;
	} else return 0; 
}

#undef __FUNCT__
#define __FUNCT__ "num_permutations_at_parity"
/*! \brief Returns the number of permutations of this mask that can be achieved with the given parity. */
uPetscInt num_permutations_at_parity(struct parameter_struct *parameters, uPetscInt mask, int parity){
	char buf[MAX_STRING];
	int count ,i;
	uPetscInt permutations;
	
	//count the number of 1s in the mask
	dec2bin(mask,parameters->number_particles, buf);
	count = 0;
	for ( i = 0 ; i < parameters->number_particles ; i++ ){
		if ( buf[i] == '1' ) count++;
	}
	//PetscPrintf(PETSC_COMM_WORLD,"Number of ones %d.\n", count);
	
	permutations = 0;
	for ( i = parity ; i <= count ; i+= 2 ) { //conserve parity in i 
		permutations += NchooseC(count,i);
		//PetscPrintf(PETSC_COMM_WORLD,"%d choose %d is : %d\n", count, i , NchooseC(count,i));
	}
	
	return permutations;
}

#undef __FUNCT__
#define __FUNCT__ "NchooseC"
/*! \brief Calculate the number of ways of choosing c elements from n. */
uPetscInt NchooseC(int n , int c ) {
	if ( c == n || c == 0 ) return 1;
	if ( c == 1 ) return n;
	return NchooseC(n,c-1) * (n-(c-1))/c;
}

#undef __FUNCT__
#define __FUNCT__ "NchooseC"
/*! \brief Calculate the number of ways of choosing c elements from n. */
uPetscInt NchooseC_safe(int n , int c ) {
	if ( c < 0 ) return 0;
	if ( n < 0 ) return 0;
	if ( c == n || c == 0 ) return 1;
	if ( c == 1 ) return n;
	return NchooseC(n,c-1) * (n-(c-1))/c;
}

#undef __FUNCT__
#define __FUNCT__ "free_parity_info"
/*! \brief Free memory used by parity information structure. */
int free_parity_info(struct parameter_struct *parameters){
	PetscErrorCode ierr;
	
	//free((md->parity_info.masks));
	//free((md->parity_info.permutations));
	//free(md->parity_info.counts);
	ierr = PetscFree(parameters->parity_info.masks);CHKERRQ(ierr);
	ierr = PetscFree(parameters->parity_info.permutations);CHKERRQ(ierr);
	ierr = PetscFree(parameters->parity_info.counts);CHKERRQ(ierr);
	
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "factorial"
/*! \brief Function which calculates factorials and partial factorials. */
int factorial(int from, int to){
	uPetscInt ans;
	int i ;
	
	ans = 1;
	for ( i = from ; i <= to ; i++ ) {
		ans *= i ;  
	}
	return ans;
}

#undef __FUNCT__
#define __FUNCT__ "get_full_parity_basis_index"
/*! \brief This function given a basis index in the reduced basis will return the full basis index corresponding to it in the given sector.*/
uPetscInt get_full_parity_basis_index(struct parameter_struct *parameters, uPetscInt reduced_index) {
	struct parity_information *parity_info;
	uPetscInt full_index;
	uPetscInt *rank_in_sector;
	//uPetscInt rank_in_sector[50];
	int i ;
	PetscErrorCode ierr;
	PetscInt sector = parameters->real_parity_sector;
	

	parity_info = &parameters->parity_info ;

	//rank_in_sector = (uPetscInt *) malloc(sizeof(uPetscInt) * parity_info->number_masks); 
	ierr = PetscMalloc(sizeof(uPetscInt) * parity_info->number_masks,&rank_in_sector);CHKERRQ(ierr);
	
	

	full_index = 0;
	for( i = parity_info->number_masks - 1 ; i >= 0 ; i-- ) {
		rank_in_sector[i] = reduced_index / parity_info->permutations[i] ; 
		reduced_index = reduced_index % parity_info->permutations[i];
//		PetscPrintf(PETSC_COMM_WORLD,"Rank in %d is %d and full is %d.\n",i, rank_in_sector[i],get_index_with_rank_parity(parameters, rank_in_sector[i],i,get_mask_parity(sector,i)));
		full_index += get_index_with_rank_parity(parameters, rank_in_sector[i],i,get_mask_parity(sector,i)); 
	}
	
	//free(rank_in_sector);
	ierr = PetscFree(rank_in_sector);CHKERRQ(ierr);
	
	return full_index;
}


#undef __FUNCT__
#define __FUNCT__ "get_index_with_rank_parity"
/*! \brief Given the rank, the sector mask and the parity. Want to return the full basis index.*/
uPetscInt get_index_with_rank_parity(struct parameter_struct *parameters, uPetscInt rank, int mask, int parity) {
	struct parity_information *parity_info;
	char buf[MAX_STRING];
	int i ,count;
	uPetscInt *values, index, counter, permutations, raw_index;
	PetscErrorCode ierr;
	
	parity_info = &parameters->parity_info ; 
	//count the number of 1s in the mask
	dec2bin(parity_info->masks[mask],parameters->number_particles, buf);
	count = 0;
	for ( i = 0 ; i < parameters->number_particles ; i++ ){
		if ( buf[i] == '1' ) count++;
	}
	//PetscPrintf(PETSC_COMM_WORLD,"Number of ones %d.\n", count);
	//values = (uPetscInt *) malloc(sizeof(uPetscInt) * count); //assign space for the values
	ierr = PetscMalloc(sizeof(uPetscInt) * count,&values);CHKERRQ(ierr);
	//loop through and assign the values
	count = 0;
	for ( i = 0 ; i < parameters->number_particles ; i++ ){
		if ( buf[parameters->number_particles - 1 - i] == '1' ) {
			values[count] = (uPetscInt)pow((double)2,(double)i);
			count++;
		}
	}
	
	/*for ( i = 0 ; i < count ; i++ ) { 
		PetscPrintf(PETSC_COMM_WORLD,"Mask %s(%d): value %d is %d.\n", buf,parity_info->masks[mask],i, values[i]);
	}*/
	
	counter = 0;
	index = 0;
	for ( i = parity ; i <= count ; i+= 2 ) { //conserve parity in i 
		permutations = NchooseC(count,i);
		if ( rank < (counter + permutations) ) { //if will have i number of 1's
			//PetscPrintf(PETSC_COMM_WORLD,"Rank %d with parity %d contained within %d choose %d as rank %d.\n", rank, parity, count, i, rank - counter);
			//PetscPrintf(PETSC_COMM_WORLD,"Index rank %d from %d choose %d: %d.\n", rank-counter, count, i, get_index_n_choose_c(rank-counter, count, i ) );
			raw_index = get_index_n_choose_c( rank - counter, count, i );
			break;
		} else {
			counter += permutations; 
		}
	}
	
	
	dec2bin(raw_index,count, buf);
	for ( i = 0 ; i < count ; i++ ) {
		if ( buf[count-1-i] == '1' ) {
			index += values[i];
		}
	}
	//PetscPrintf(PETSC_COMM_WORLD,"Raw index for rank %d with parity %d : %d and full %d.\n",rank, parity , raw_index,index);
	
	//free(values);
	ierr = PetscFree(values);CHKERRQ(ierr);
	
	return index;
}

#undef __FUNCT__
#define __FUNCT__ "get_index_n_choose_c"
/*! \brief Iterative function which returns the index with given rank with c 1's out of n.*/
uPetscInt get_index_n_choose_c(uPetscInt rank, int n, int c){
	uPetscInt index, accummulation, current, tmp_add;
	int i ;
	
	//PetscPrintf(PETSC_COMM_WORLD,"- rank %d from %d choose %d.\n", rank, n, c);
	
	index = 0;
	accummulation = 0;
	
	if ( c == 1 ) { 
		//PetscPrintf(PETSC_COMM_WORLD,"--- adding %d to index %d (raised %d).\n", (uPetscInt)pow(2,rank), index,rank );
		index += (uPetscInt)pow((double)2,(double)rank) ;
		
		return index;
	} else if ( c == 0 ) {
		return 0;
	}
	
	for ( i = 0; i <= n-c ; i++ ) {
		current = NchooseC(n-1-i,c-1) ;
		if ( (accummulation + current) > rank ) { // if this rank has a 1 in the (i+1)th position
			//PetscPrintf(PETSC_COMM_WORLD,"-- adding %d to index %d.\n", (uPetscInt)pow(2,i), index);
			index += (uPetscInt)pow((double)2,(double)i);
			tmp_add = (get_index_n_choose_c(rank - accummulation, n-1-i, c - 1)) << (i+1);
			//PetscPrintf(PETSC_COMM_WORLD,"---- adding %d to index %d (raised 2^%d).\n", tmp_add, index, i+1 );
			index += tmp_add;
			return index;
		}
		accummulation += current;
	}
	
	return index;
}

#undef __FUNCT__
#define __FUNCT__ "get_reduced_parity_basis_index"
/*! \brief This function when given an index in the full basis will return the corresponding reduced basis index.*/
uPetscInt get_reduced_parity_basis_index(struct parameter_struct *parameters, uPetscInt full_index, PetscTruth *exists) { // TODO: get exists variable checking if parity is correct.
	struct parity_information *parity_info;
	uPetscInt reduced_index, rank, raw_index;
	char buf[MAX_STRING], buf2[MAX_STRING];
	int i ,j, power, ns;
	PetscInt sector = parameters->real_parity_sector;
	
	*exists = PETSC_TRUE;  
	parity_info = &parameters->parity_info ;
	reduced_index = 0; 
	ns = parameters->number_particles;
	
	//find the bit pattern in each mask
	for ( i = 0 ; i < parity_info->number_masks ; i++ ) {
		raw_index = 0;
		rank = full_index & parity_info->masks[i] ; 
		dec2bin(rank,parameters->number_particles, buf); // convert to a character string
		dec2bin(parity_info->masks[i],parameters->number_particles, buf2); // convert mask to a character string
		power = 0 ;
		int count = 0;
		for ( j = ns - 1 ; j >= 0 ; j-- ) {
			if ( buf2[j] == '1' ) {
				if ( buf[j] == '1' ) {
				count++;
					raw_index += (uPetscInt)pow((double)2,(double)power);
				}
				power++;
			}
		}
		if ( (count % 2 ) != get_mask_parity(sector,i) ) {
			*exists = PETSC_FALSE;
		}
		//PetscPrintf(PETSC_COMM_WORLD,"From full_index %d to raw index %d using mask %s and string %s.\n",full_index, raw_index, buf2,buf);
		reduced_index += parity_info->permutations[i] * get_rank_with_index_parity(raw_index, get_mask_parity(sector,i), parity_info->counts[i] );
	}

	return reduced_index;
}

#undef __FUNCT__
#define __FUNCT__ "get_rank_with_index_parity"
/*! \brief Return the rank of an index with given size in bits and given parity.*/
uPetscInt get_rank_with_index_parity(uPetscInt index, int parity, int count ){
	int ones, i;
	char buf[MAX_STRING];
	uPetscInt rank;
	
	ones = 0; 
	dec2bin(index, count, buf); //convert to a character string to count the ones
	for ( i = 0 ; i < count ; i++ ) { 
		if ( buf[i] == '1' ) { 
			ones++;
		}	
	}
	
	//depending on the number of ones we can already add an offset to the rank. If no ones then rank is just 0. 
	rank = 0;
	
	//this will work up to the rank maintaining the parity and adding the offset each time.
	for ( i = parity ; i < ones ; i += 2 ) {
		rank += NchooseC(count, i); //add the number of choices of i elements from a set of size count.
	}
	
	//PetscPrintf(PETSC_COMM_WORLD,"Buf %s, Index %d, count %d, parity %d starting %d, ones %d.\n",buf, index, count, parity, rank, ones);
	//now that all of the combinations with number of less that we have have been taken care of we want to find the rank
	//within the n choose ones.
	if ( ones > 0 ) { //does not work for no ones at all rank will just be 0.
		rank += get_rank_n_choose_c(index, count, ones); 
	}
	
	return rank;
}


#undef __FUNCT__
#define __FUNCT__ "get_rank_n_choose_c"
/*! \brief Iterative function which returns the rank with given the index with c 1's out of n.*/
uPetscInt get_rank_n_choose_c(uPetscInt index, int n, int c){
	uPetscInt rank; 
	int i ;
	
	rank = 0 ;
	for ( i = 0 ; i < n ; i++ ) {
		if ( (index &  (uPetscInt)pow((double)2,(double)i) ) == 0 ) { //if there is a 1 in the ith position from the right
			//missed out on a whole (n-1-i) choose (c-1) so add this
			if ( c > 1) {
				rank += NchooseC(n-1-i,c-1);
	//			PetscPrintf(PETSC_COMM_WORLD,"No one again so add %d choose %d: %d (total %d).\n" , n-1-i, c-1, NchooseC(n-1-i,c-1), rank);
			}
		} else { //otherwise if it was one then we iterate
			if ( c == 1 ) {
				rank += i ; 
				break;
			} else {
				rank += get_rank_n_choose_c( (index >> (i+1)), n - i - 1, c - 1 );	
			}
			break;
		}
	}
	//PetscPrintf(PETSC_COMM_WORLD,"- index %d from %d choose %d. return %d\n", index	, n, c,rank);
	return rank; 
}
