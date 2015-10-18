#include "p_entanglement_entropy_fourier_cut.h"


PetscInt min(PetscInt A,PetscInt B){ 
	if ( A < B ){
		return A;
	} else {
		return B;
	}
}

//utility function to convert decimal to binary 
char *dec2bin(uPetscInt a, int len, char *num) {
	//char *num ;
	int i;
	uPetscInt tmp;
	
	//allocate space for each character
	//num = (char *)malloc((len + 1 )* sizeof(char));
	
	//Start from the highest values 
	for ( i = len-1 ; i >= 0 ; i-- ) {
		tmp = (uPetscInt)pow((double)2,(double)i);
		if ( a >= tmp ) {
			num[len - 1 -i] = '1';
			a -= tmp;
		} else {
			num[len - 1 -i] = '0';
		}
	}
	num[len] = '\0';
	
	return num;
}


int parse_adjacency_info_file(string filename, struct adjacency_information *adj_info,PetscInt *sites){
	//string filename(adjacency_file);
	PetscErrorCode ierr;
	
	*sites = 0 ;
	
	TiXmlDocument doc(filename);
	
	bool loadOkay = doc.LoadFile();
	
	if ( !loadOkay )
	{
		PetscPrintf(PETSC_COMM_WORLD,"Could not load adjacency file '%s'. Error='%s'. Exiting.\n", filename.c_str(), doc.ErrorDesc());
		exit(1);
	}
	
	TiXmlHandle docHandle( &doc );
	
	TiXmlElement* element = docHandle.FirstChild("EDGES").Element();
	TiXmlElement* element2;
	if ( element) {	
		adj_info->number_masks = atoi(element->Attribute("number"));
		if ( adj_info->number_masks <= 0 ){
			PetscPrintf(PETSC_COMM_WORLD,"Invalid number of edges %d.\n",adj_info->number_masks);
			exit(1);
		}
		ierr = PetscMalloc( adj_info->number_masks * sizeof(uPetscInt), &(adj_info->masks));CHKERRQ(ierr);
		for ( int i = 0 ; i < adj_info->number_masks ; i++ ) {
			adj_info->masks[i] = 0;
			element2 = docHandle.FirstChild("EDGES").Child("EDGE",i).Element();
			if ( element2 ) {
				adj_info->masks[i] += (uPetscInt)pow(2.0,double(atoi(element2->Attribute("from"))-1));
				adj_info->masks[i] += (uPetscInt)pow(2.0,double(atoi(element2->Attribute("to"))-1));
				if ( atoi(element2->Attribute("from")) > *sites) *sites = atoi(element2->Attribute("from"));
				if ( atoi(element2->Attribute("to")) > *sites) *sites = atoi(element2->Attribute("to"));
			} else {
				PetscPrintf(PETSC_COMM_WORLD,"Edge %d does not exist.\n",i);
				exit(1);
			}
			
		}
		/*for ( int i = 0 ; i < adj_info->number_masks ; i++ ) {
			char buf[MAX_STRING];
			dec2bin( adj_info->masks[i], *sites, buf );
			PetscPrintf(PETSC_COMM_WORLD,"Mask %d: %s \n",i,buf);			
		}*/
	}
	
	return 0;
}

int parse_a_sites_file(string filename, struct parameters *params){
	
	TiXmlDocument doc(filename);
	
	bool loadOkay = doc.LoadFile();
	
	if ( !loadOkay )
	{
		PetscPrintf(PETSC_COMM_WORLD,"Could not read A sites file '%s'. Error='%s'. Exiting.\n", filename.c_str(), doc.ErrorDesc() );
		exit(1);
	}
	
	TiXmlHandle docHandle( &doc );
	
	TiXmlElement* element = docHandle.FirstChild("SITES").Element();
	TiXmlElement* element2;
	if ( element) {	
		int num_sites = atoi(element->Attribute("number"));
		if ( num_sites <= 0 ){
			PetscPrintf(PETSC_COMM_WORLD,"Invalid number of sites %d.\n",num_sites);
			exit(1);
		}
		params->A_sites = num_sites;
		params->B_sites = params->AB_sites - num_sites;
		params->A_mask = 0; 
		for ( int i = 0 ; i < num_sites ; i++ ) {
			element2 = docHandle.FirstChild("SITES").Child("SITE",i).Element();
			if ( element2 ) {
				params->A_mask += (uPetscInt)pow(2.0,double(atoi(element2->GetText())-1));
			} else {
				PetscPrintf(PETSC_COMM_WORLD,"Site %d does not exist.\n",i);
				exit(1);
			}
		}
		
		//set B mask by first filling with ones and then taking XOR with A mask.
		params->B_mask = 0;
		for ( int i = 0 ; i < params->AB_sites ; i++ ) {
			params->B_mask += (uPetscInt)pow(2.0,double(i));
		}
		params->B_mask ^= params->A_mask; 
		params->AB_mask = params->A_mask | params->B_mask;
	}
	
	return 0;
}

//this function checks the adjacency information for the given index and returns false if it contains adjacent particles.
PetscTruth check_adjacency(struct parameters *params, uPetscInt index, struct conf_list_info *c_info) {
	for ( int i = 0 ; i < c_info->adj_info.number_masks; i++ ) {
		if ( (index & c_info->adj_info.masks[i] ) == c_info->adj_info.masks[i] ) {
			return PETSC_FALSE;
		}
	}
	return PETSC_TRUE;
}

PetscInt get_rank_chunk_size(PetscInt n, PetscInt rank, PetscInt size) {
	PetscInt  max_chunk, min_chunk, cut_off; 
	
	min_chunk = (PetscInt)floor((double)n / (double)size); 
	max_chunk = (PetscInt)ceil((double)n / (double)size);
	cut_off = n - (size * min_chunk); 
	
	if ( rank < cut_off ) {
		return max_chunk;
	} else {
		return min_chunk;
	}
}

PetscInt get_rank_chunk_start(PetscInt n, PetscInt rank, PetscInt size){
	PetscInt  max_chunk, min_chunk, cut_off; 
	
	min_chunk = (PetscInt)floor((double)n / (double)size); 
	max_chunk = (PetscInt)ceil((double)n / (double)size);
	cut_off = n - (size * min_chunk); 

	if ( rank < cut_off ) {
		return rank * max_chunk;
	} else {
		return (cut_off * max_chunk) + ((rank - cut_off) * min_chunk);
	}
}

/*! \brief This function is a dummy function which will always return true. When the standard basis is being used the relvant function reference points to this function.*/
PetscTruth dummy_basis_check(struct parameters *params, uPetscInt index){
	return PETSC_TRUE;
}

uPetscInt apply_rotation(struct parameters *params, uPetscInt basis,int *number,PetscScalar *phase){
	*phase = 1.0;
	if ( number[0] > 0 ) {
		list<int>  creation_ops;
		int from;
		uPetscInt new_basis;
		struct rotation_information *rot_info;
		rot_info = &params->rotation_info;	
		new_basis = 0;
		for ( int spin = 0 ; spin < params->AB_sites ; spin++ ){
			if ( (uPetscInt)(basis & rot_info->mapping_masks[number[0]][spin]) > (uPetscInt)0 ) {
				new_basis += ((uPetscInt)1) << spin;
				from = (int)log2((double)rot_info->mapping_masks[number[0]][spin]);
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
	} else {
		return basis;
	}
}

/*! \brief Gets the phase associated with the specified rotation.
 	Given the amount to rotate and the sector we are currently in this function returns the phases
	associated with the rotation for a given rotation sector. This is usually a complex number and so a complex build of PETSc should be used. */ 
PetscScalar get_rotation_phase(struct parameters *params,PetscInt sector){
	PetscScalar phase;
	PetscScalar k;
	struct rotation_information *rot_info;
	rot_info = &params->rotation_info;
	k = ((PetscScalar)sector)*2.0*PETSC_PI*((PetscScalar)rot_info->real_sector)/(PetscScalar)rot_info->num_sectors*PETSC_i;
	phase = PetscExpScalar(k);
	return phase;
}

/*! \brief Calculate and store all possible rotation phases so that they can just be used later without calculation. */
int calculate_rotation_phases(struct parameters *params){
	int number_rotations = params->rotation_info.num_sectors;
	
	//allocate space to store phases.
	PetscMalloc(number_rotations*sizeof(PetscScalar),&params->rotation_info.phases);//allocate space to store phases
	
	//go through translations and set the phase for each. 
	for ( PetscInt i = 0 ; i < number_rotations ; i++ ) {
		params->rotation_info.phases[i] = get_rotation_phase(params,i);
		PetscPrintf(PETSC_COMM_WORLD,"Phase %d: %lf + i%lf\n",i,PetscRealPart(params->rotation_info.phases[i]),PetscImaginaryPart(params->rotation_info.phases[i]));
	}
	return 0;
}


/*! \brief Check if this particular index is in the current rotation sector. */
PetscTruth rotation_representative_check(struct parameters *params, uPetscInt index){
	struct rotation_information *rot_info;
	uPetscInt new_basis;
	int amount;
	PetscScalar phase_total;
	//PetscScalar phase_tmp;
	//PetscTruth exists;
	PetscScalar sign;
	//char buf[MAX_STRING];
	//char buf2[MAX_STRING];
	
	phase_total = 0.0;
	rot_info = &params->rotation_info;
	for ( amount = 0 ; amount < rot_info->num_sectors ; amount++ ) {
		new_basis = apply_rotation(params,index,&amount,&sign);
		//if ( exists == PETSC_TRUE ) {
			if (new_basis < index ) { 
				/*dec2bin(index, params->AB_sites, buf );
				dec2bin(new_basis, params->AB_sites, buf );
				PetscPrintf(PETSC_COMM_SELF,"Not using %s as %s is less.\n",buf, buf2);*/
				return PETSC_FALSE;
			} else if ( new_basis == index ) {
				phase_total += sign * rot_info->phases[amount];
			}
		//}
	}
	
	if ( PetscAbsScalar(phase_total) < 1.0e-9 ) {
		/*dec2bin(index,params->AB_sites,buf);
		PetscPrintf(PETSC_COMM_SELF,"Phase is 0 for %s.\n",buf);*/
		return PETSC_FALSE;
	}
	
	return PETSC_TRUE;
}

int find_valid_configs_with_this_subset_list(struct parameters *params, struct conf_list_info *c_info, int set_size, int set_filling, uPetscInt set , list<uPetscInt> *ll) {
	//first check the stop condition if set size is already the size of the full graph or there is not enough space left to reach required filling.
	if ( set_size == c_info->num_sites ) {
		if ( c_info->filling == -1 || set_filling == c_info->filling ) { //if there is a filling limit and we are at it just add the configuration.
			if ( (*params->representative_check)(params,set) ) { ll->push_back(set); } 
			return 0;
		} else {
			return 0;
		}
	} else if ( c_info->filling != -1 ) {
		if ((c_info->filling - set_filling) > (c_info->num_sites - set_size) )  { //if filling limit and will not be able to make it just return.
			return 0 ;
		} else if ( set_filling == c_info->filling ) {
			if ( (*params->representative_check)(params,set) ) { ll->push_back(set); } 
			return 0;
		}
	} else if ( c_info->filling == -1 ) { //we check that everything is alright with respect to max filling.
		if ( set_filling == c_info->max_filling ) {
			if ( (*params->representative_check)(params,set) ) { ll->push_back(set); } 
			return 0;
		} else if ( set_filling > c_info->max_filling ) { //if its gone over the max then just return. 
			return 0;
		}
	}
	
	//have two branches one with 0 at site set_size and one with 1 at set_size
	//now add a 1 to position set_size, check if valid and increment filling.
	uPetscInt new_set = set + (((uPetscInt)1) << set_size) ; 
	if ( check_adjacency(params,new_set,c_info) ) {
		find_valid_configs_with_this_subset_list(params, c_info,set_size+1 , set_filling+1, new_set, ll);
	}
	
	
	//do not need to check the addition of a 0 as will already be full.
	find_valid_configs_with_this_subset_list(params,c_info,set_size+1 , set_filling, set, ll);
	return 0;
}

int find_valid_configs_with_this_subset(struct parameters *params,struct conf_list_info *c_info,int set_size, int set_filling, uPetscInt set , ofstream *basis_file, PetscInt *local_reps ) {
	uPetscInt set2;
	//first check the stop condition if set size is already the size of the full graph or there is not enough space left to reach required filling.
	if ( set_size == c_info->num_sites ) {
		if ( c_info->filling == -1 || set_filling == c_info->filling ) {
			if ( (*params->representative_check)(params,set) ) {   
				set2 = set;
				basis_file->write((char*)&set2,sizeof(uPetscInt));	
				//cout << "Found " << set << endl;
				(*local_reps)++;
			} /*else {
					char buf[MAX_STRING];
					dec2bin(set,params->AB_sites,buf);
					PetscPrintf(PETSC_COMM_SELF,"Rejected as not representative: %s.\n",buf);
			}*/
			return 0;
		} else {
			return 0;
		}
	} else if ( c_info->filling != -1 ) {
		if ((c_info->filling - set_filling) > (c_info->num_sites - set_size)) {
			return 0 ;
		} else if ( set_filling == c_info->filling ) {
			if ( (*params->representative_check)(params,set) ) {   
				set2 = set;
				basis_file->write((char*)&set2,sizeof(uPetscInt));	
				//cout << "Found " << set << endl;
				(*local_reps)++;
			} 	/*else {
					char buf[MAX_STRING];
					dec2bin(set,params->AB_sites,buf);
					PetscPrintf(PETSC_COMM_SELF,"Rejected as not representative: %s.\n",buf);
			}*/
			return 0;
		}
	} else if ( c_info->filling == -1 ) {
		if ( set_filling == c_info->max_filling ) {
			//ll->push_back(set);
			if ( (*params->representative_check)(params,set) ) {   
				set2 = set;
				basis_file->write((char*)&set2,sizeof(uPetscInt));	
				//cout << "Found " << set << endl;
				(*local_reps)++;
			} /*else {
				char buf[MAX_STRING];
				dec2bin(set,params->AB_sites,buf);
				PetscPrintf(PETSC_COMM_SELF,"Rejected as not representative: %s.\n",buf);
			}*/
			return 0;
		} else if ( set_filling > c_info->max_filling ) { //if its gone over the max then just return. 
			return 0;
		}
	}
	
	//have two branches one with 0 at site set_size and one with 1 at set_size
	//now add a 1 to position set_size, check if valid and increment filling.
	uPetscInt new_set = set + (((uPetscInt)1) << set_size) ; 
	if ( check_adjacency(params,new_set,c_info) ) {
		find_valid_configs_with_this_subset(params,c_info,set_size+1 , set_filling+1, new_set,basis_file,local_reps);
	}
	
	//do not need to check the addition of a 0 as will already be full.
	find_valid_configs_with_this_subset(params, c_info, set_size+1 , set_filling, set,basis_file, local_reps);
	return 0;
}

int find_valid_configs_list(struct parameters *params, struct conf_list_info *c_info, list<uPetscInt> *ll){
	// will work across the tree breadth first until we arrive at a depth which has as many valid configurations on sub lattices as 
	// processes and then will distribute divide the valid configurations among the processes.
	vector<uPetscInt> current_configs;
	//vector<uPetscInt> next_configs;
	vector<int> current_fillings;
	//vector<int> next_fillings;
	int set_size = 0;
	
	//current_configs = new vector<uPetscInt>(); 
	//current_fillings = new vector<int>();	
	current_configs.push_back((uPetscInt)0);
	current_fillings.push_back(0);
	
	find_valid_configs_with_this_subset_list(params,c_info, set_size,(int)(current_fillings)[0],(uPetscInt)(current_configs)[0],ll);	
	int count = 0;
	for ( list<uPetscInt>::iterator iter = ll->begin();  iter != ll->end() ; iter++ ) {
		count++;
	}
	//cout << "Rank : " << params->rank << " number : " << count << endl;
	
	return 0;
}

// this function is overloaded and this version writes the configurations to disk and works for parallel implemenetations.
int find_valid_configs(struct parameters *params, struct conf_list_info *c_info, ofstream *basis_file, PetscInt *local_reps){
	// will work across the tree breadth first until we arrive at a depth which has as many valid configurations on sub lattices as 
	// processes and then will distribute divide the valid configurations among the processes.
	vector<uPetscInt> *current_configs;
	vector<uPetscInt> *next_configs;
	vector<int> *current_fillings;
	vector<int> *next_fillings;
	int set_size = 0;
	
	current_configs = new vector<uPetscInt>(); 
	current_fillings = new vector<int>();	
	current_configs->push_back((uPetscInt)0);
	current_fillings->push_back(0);
	
	//count number of current configs.
	int count = 0;
	for ( vector<uPetscInt>::iterator iter = current_configs->begin() ; iter != current_configs->end(); iter++ ) {
		count++;
	}
	
	//while there are less valid configs than processes 
	while ( count < params->size && set_size < c_info->num_sites ) {
		next_configs = new vector<uPetscInt>(); 
		next_fillings = new vector<int>();
		int i = 0 ;
		int new_count = 0; 
		while ( i < count ) {		
			if ( (c_info->filling == -1 && ((*current_fillings)[i] <= c_info->max_filling ) )  || ((c_info->filling - (*current_fillings)[i]) <= (c_info->num_sites - set_size)) ) {
				next_configs->push_back((*current_configs)[i]);
				next_fillings->push_back((*current_fillings)[i]);
				new_count ++;
			}
				
			if ( ( c_info->filling == -1 && ((*current_fillings)[i] < c_info->max_filling ) ) || ((*current_fillings)[i] < c_info->filling) ) {
				uPetscInt tmp;
				tmp = (*current_configs)[i] + (((uPetscInt)1) << set_size) ;
				if ( check_adjacency(params,tmp,c_info)) {
					next_configs->push_back(tmp);
					next_fillings->push_back((*current_fillings)[i]+1);
					new_count++;
				}	
			}
			i++;
		}	
		set_size++; 
		count = new_count;
		delete current_configs;
		delete current_fillings;
		current_configs = next_configs;
		current_fillings = next_fillings;
	}
	
	*local_reps = 0 ;
	if ( count >= params->size && set_size <= c_info->num_sites ) {
		//cout << "Count " << count << " rank " << params->rank << " size " << get_rank_chunk_size(count,params->rank,params->size) << " and start " <<  get_rank_chunk_start(count,params->rank,params->size)  << endl;
		for ( PetscInt idx = 0; idx <  get_rank_chunk_size(count, params->rank,params->size) ; idx++ ) {
			//cout << "Rank " << params->rank << " getting started with " << (*current_configs)[get_rank_chunk_start(count, params->rank,params->size)+idx] << endl;
			find_valid_configs_with_this_subset(params,c_info,set_size,(*current_fillings)[get_rank_chunk_start(count, params->rank,params->size)+idx],  (*current_configs)[get_rank_chunk_start(count, params->rank,params->size)+idx],basis_file,local_reps);	
		}		
	} else if ( count < params->size && set_size <= c_info->num_sites ) {
		if ( params->rank < count ) {
			find_valid_configs_with_this_subset(params,c_info,set_size,(*current_fillings)[params->rank],  (*current_configs)[params->rank],basis_file,local_reps);	
		}
	}
	
	delete current_configs;
	delete current_fillings;
	//cout << "Rank : " << params->rank << " number : " << *local_reps << endl;
	
	return 0;
}

PetscInt find_in_sorted_array(uPetscInt *array,PetscInt len, uPetscInt item) {
	if ( len == 0 || item < array[0] || item > array[len-1]) {
		return -1;
	}
	
	PetscInt index, l,r;
	l = 0;
	r = len-1;
	index = len/2;
	if ( array[l] == item ) return l;
	if ( array[r] == item ) return r;
	
	while ( r - l > 0 ) {
		if ( array[l+(r-l)/2] == item ) {
			return l+(r-l)/2;
		}
		if ( array[l+(r-l)/2] > item ) {
			r = l+(r-l)/2;
		//} else if ( array[l+(r-l)/2] < item ) {
		} else { 
			//l = array[l+(r-l)/2];
			l = l+(r-l)/2 + 1;
		}
	}
	return -1;
}

PetscInt find_in_sorted_array(PetscInt *array,PetscInt len, PetscInt item) {
	if ( len == 0 || item < array[0] || item > array[len-1]) {
		return -1;
	}
	
	PetscInt index, l,r;
	l = 0;
	r = len-1;
	index = len/2;
	if ( array[l] == item ) return l;
	if ( array[r] == item ) return r;
	
	while ( r - l > 0 ) {
		if ( array[l+(r-l)/2] == item ) {
			return l+(r-l)/2;
		}
		if ( array[l+(r-l)/2] > item ) {
			r = l+(r-l)/2;
		//} else if ( array[l+(r-l)/2] < item ) {
		} else { 
			//l = array[l+(r-l)/2];
			l = l+(r-l)/2 + 1;
		}
	}
	return -1;
}

/*! \brief Function to perform parallel quicksort to sort an array of unsigned PetscInts.*/
int parallel_qsort_uPetscInt(MPI_Comm *subcomm,uPetscInt *array, PetscInt llength) {
	//get size of sub communicator as well as the rank of the current process within it.
	PetscMPIInt rank, size ;
	PetscErrorCode ierr;
	ierr = MPI_Comm_rank(*subcomm,&rank);CHKERRQ(ierr);
	ierr = MPI_Comm_size(*subcomm,&size);CHKERRQ(ierr);
	
	//PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, size %d and length %d.\n",rank, size, llength);	
	
	//first we sort the local portions of the array
	sort(array,array+llength);
	
	//If the communicator is only made up of a single process then we are done.
	if ( size > 1 ) {  
		PetscInt *lengths; //array to store the lengths of the other chunks.
		ierr = PetscMalloc(size*sizeof(PetscInt),&lengths);CHKERRQ(ierr);
		MPI_Allgather(&llength,1,PETSC_MPI_INT,lengths,1,PETSC_MPI_INT,*subcomm); //gather lengths of all chunks on each process.
		PetscInt total_length = 0;
		for ( int i = 0 ; i < size ; i++ ) total_length += lengths[i]; //find the total length of the subset.
	
		//select the pivot as the central item in the sorted list.
		PetscReal pivot = 0 ;
		if ( llength > 0 ) {
			pivot = (PetscReal)array[llength/2] ;
		}  
		
		//sum the pivots among all the processes in the communicator.
		PetscReal sum = 0.0; 
		MPI_Allreduce(&pivot, &sum, 1,MPIU_REAL, MPI_SUM, *subcomm); 
		//divide by the number of processes to get the mean.
		//pivot = sum/total_length;
		PetscReal real_pivot = (PetscReal)sum / (PetscReal)size;
		
		PetscInt counts[2] = {0,0}, *g_below, *g_above; //array containing counts of array elements below and above the pivot.
		ierr = PetscMalloc(size*sizeof(PetscInt),&g_below);CHKERRQ(ierr);  //allocat space to store the numbers above and below the pivot on each process.
		ierr = PetscMalloc(size*sizeof(PetscInt),&g_above);CHKERRQ(ierr);
		{
			PetscInt i = 0;
			while ( i < llength && ((PetscReal)array[i]) < real_pivot ) {
				if ( ((PetscReal)array[i]) < real_pivot ) counts[0]++;
				i++;
			}
			counts[1] = llength - i;
		}
		MPI_Allgather(&counts[0],1,PETSC_MPI_INT,g_below,1,PETSC_MPI_INT,*subcomm); 
		MPI_Allgather(&counts[1],1,PETSC_MPI_INT,g_above,1,PETSC_MPI_INT,*subcomm); 
		
		//find the total number of elements less than the pivot and the total amount more than the pivot.
		PetscInt total_below = 0, total_above = 0 ;
		for ( int i = 0 ; i < size ; i++ ) {
			total_below += g_below[i];
			total_above += g_above[i];
		}
		
		//find the rank of the process which contains the pivot.
		int pivot_process = 0;
		PetscInt pivot_pos_within_chunk = 0;
		PetscInt l = 0;
		for ( int i = 0 ; i < size ; i++ ) {
			if ( total_below <= (l + lengths[i]) ){
				pivot_process = i;
				pivot_pos_within_chunk = total_below - l; 
				break;
			}
			l += lengths[i];
		}
		
		//check if there are any elements to exchange and if not we can finish here.
		/*PetscInt number_to_exchange_below = 0 ;
		for ( int i = 0 ; i < pivot_process ; i++ ) {
			number_to_exchange_below += g_above[i];
		}
		PetscInt number_to_exchange_above = 0 ;
		for ( int i = size-1 ; i > pivot_process ; i-- ) {
			number_to_exchange_above += g_below[i];
		}
		
		if ( number_to_exchange_above == 0 && number_to_exchange_below == 0  ) {
			ierr = PetscFree(lengths);CHKERRQ(ierr);
			ierr = PetscFree(g_below);CHKERRQ(ierr);
			ierr = PetscFree(g_above);CHKERRQ(ierr);
			return 0;
		}*/
//		PetscPrintf(PETSC_COMM_SELF,"Pivot process %d.\n",pivot_process);
		
		if ( rank == pivot_process ) { //if this is the pivot process.
			//find the total length of the subset
			l = 0 ;
			for ( int i = 0 ; i < pivot_process ; i++ ) {
				l += lengths[i];
			}
			PetscInt end_of_lower = total_below - l; //the local index of the end of the lower partition.
			//find the total length of the sites below the pivot
			l = 0;
			for ( int i = size-1 ; i > pivot_process ; i-- ) {
				l += lengths[i];
			}
			//PetscInt start_of_upper = llength - (total_above - l); //the local index of the start of the upper index.
			//if ( start_of_upper != end_of_lower) cout << "Ends dont meet : " <<  end_of_lower << ", " << start_of_upper << endl;
			
			if ( end_of_lower > counts[0] ) { // if this process has too many elements that are greater than or equal to the pivot. Need to exchange these elements for ones that are below the pivot.
				PetscInt to_exchange = end_of_lower - counts[0]; //the number of items there are to exchange.
				PetscInt to_exchange_first = 0;
				for ( int i = 0 ; i < pivot_process ; i++ ) to_exchange_first += g_above[i]; //count up number each process below has to share first.
				uPetscInt *recv_ptr = &(array[counts[0]]);  
				for ( int proc = pivot_process+1 ; proc < size ; proc++ ) { //loop through processes until we find the one to exchange with.
					if ( (to_exchange_first - g_below[proc]) >= 0 ) {
						to_exchange_first -= g_below[proc] ; //subtract anything to be exchanged between other processes.
					} else {
						uPetscInt *send_buf;
						int send_size; 
						if ( to_exchange > 0 ) {
							if ( g_below[proc] - to_exchange_first > to_exchange ) { //if there is more to give away than is needed. 
								send_size = (int)(to_exchange);					
							} else {
								send_size = (int)(g_below[proc] - to_exchange_first);
							}
							if ( to_exchange_first > 0 ) { 
								//recv_ptr += (int)to_exchange_first;
								to_exchange_first = 0;
							}
							ierr = PetscMalloc(send_size * sizeof(uPetscInt),&send_buf);CHKERRQ(ierr);
							memcpy(send_buf,recv_ptr,send_size*sizeof(uPetscInt)); 
							MPI_Request send_req; MPI_Status send_stat;
							MPI_Isend(send_buf,send_size,PETSC_UNSIGNED_MPI_INT,proc,1,*subcomm,&send_req);//send using non blocking communication.
	//						PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, sending and receiving %d from %d (0).\n",rank, send_size, proc);
							MPI_Recv(recv_ptr,send_size,PETSC_UNSIGNED_MPI_INT,proc,1,*subcomm,&send_stat);//wait for receive. 
							MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
							PetscFree(send_buf);
							to_exchange -= send_size;
							recv_ptr += send_size;
						} //end if no_exchange > 0 
					} //ebd if 
				}  //end for over processes
			} else if ( end_of_lower < counts[0] ) { //otherwise if this process has too many elements that are lower than the pivot. Need to exchange these elements for ones that are above the pivot.
				PetscInt to_exchange = counts[0] - end_of_lower; //the number of items there are to exchange.
				uPetscInt *recv_ptr = &(array[end_of_lower]);
				for ( int proc = 0 ; proc < pivot_process ; proc++ ) { //loop through processes before this one.
					uPetscInt *send_buf;
					int send_size; 
					if ( to_exchange > 0  && g_above[proc] > 0) {
						if ( g_above[proc] > to_exchange ) { //if there is more to give away than is needed. 
							send_size = (int)to_exchange;					
						} else {
							send_size = (int)g_above[proc];
						}
						ierr = PetscMalloc(send_size * sizeof(uPetscInt),&send_buf);CHKERRQ(ierr);
						memcpy(send_buf,recv_ptr,send_size*sizeof(uPetscInt)); 
						MPI_Request send_req; MPI_Status send_stat;
						MPI_Isend(send_buf,send_size,PETSC_UNSIGNED_MPI_INT,proc,1,*subcomm,&send_req);//send using non blocking communication.
//						PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, sending and receiving %d from %d (1).\n",rank, send_size, proc);
						MPI_Recv(recv_ptr,send_size,PETSC_UNSIGNED_MPI_INT,proc,1,*subcomm,&send_stat);//wait for receive. 
						MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
						PetscFree(send_buf);
						to_exchange -= send_size;
						recv_ptr += send_size;
					}
				}
			}
		} else if ( rank < pivot_process ) { // means that everything above must be exchanged for elements below.
			PetscInt to_exchange = counts[1];
			uPetscInt *recv_ptr = &(array[counts[0]]);
			PetscInt to_exchange_first = 0;
			for ( int i = 0 ; i < rank ; i++ ) to_exchange_first += g_above[i];
			
			//first check the pivot process as it is a bit anomalous.
			if ( g_below[pivot_process] - pivot_pos_within_chunk > 0 ) {
				if ( (to_exchange_first - (g_below[pivot_process] - pivot_pos_within_chunk) ) >= 0 ) { 
					to_exchange_first -= (g_below[pivot_process] - pivot_pos_within_chunk); //subtract anything to be exchanged.
				} else {
					uPetscInt *send_buf;
					int send_size ;	
					if ( to_exchange > 0 ) {
						if ( ((g_below[pivot_process] - pivot_pos_within_chunk) - to_exchange_first) > to_exchange ) { //if there is more to give away than is needed. 
							send_size = (int)(to_exchange);	
						} else {
							send_size = (int)(g_below[pivot_process]-pivot_pos_within_chunk-to_exchange_first);
							to_exchange_first = 0;
						}
						ierr = PetscMalloc(send_size * sizeof(uPetscInt),&send_buf);CHKERRQ(ierr);
						memcpy(send_buf,recv_ptr,send_size*sizeof(uPetscInt)); 
						MPI_Request send_req; MPI_Status send_stat;
						MPI_Isend(send_buf,send_size,PETSC_UNSIGNED_MPI_INT,pivot_process,1,*subcomm,&send_req);//send using non blocking communication.
	//					PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, sending and receiving %d from %d (2).\n",rank, send_size, pivot_process);
						MPI_Recv(recv_ptr,send_size,PETSC_UNSIGNED_MPI_INT,pivot_process,1,*subcomm,&send_stat);//wait for receive. 
						MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
						PetscFree(send_buf);
						to_exchange -= send_size;
						recv_ptr += send_size ;
					} // end if to_exchange > 0 
				} // end if else there is enough to cover what needs to be exchanged first
			} //end of check of pivot process
			
			//now go throug the rest of the processes. 
			for ( int proc = pivot_process+1; proc < size ; proc ++ ) { 
				uPetscInt *send_buf;
				int send_size;	
				if ( to_exchange > 0 ) {
					if ( to_exchange_first - g_below[proc] >= 0 ) {
						to_exchange_first -= g_below[proc];
					} else {
						if ( g_below[proc] - to_exchange_first >= to_exchange ) { //if there is more to give away than is needed. 
							send_size = (int)(to_exchange);					
						} else {
							send_size = (int)(g_below[proc] - to_exchange_first);
						}
					
						if ( to_exchange_first > 0 ) { 
							//recv_ptr += (int)to_exchange_first;
							to_exchange_first = 0;
						}
						ierr = PetscMalloc(send_size * sizeof(uPetscInt),&send_buf);CHKERRQ(ierr);
						memcpy(send_buf,recv_ptr,send_size*sizeof(uPetscInt)); 
						MPI_Request send_req; MPI_Status send_stat;
						MPI_Isend(send_buf,send_size,PETSC_UNSIGNED_MPI_INT,proc,1,*subcomm,&send_req);//send using non blocking communication.
		//				PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, sending and receiving %d from %d (3).\n",rank, send_size, proc);				
						MPI_Recv(recv_ptr,send_size,PETSC_UNSIGNED_MPI_INT,proc,1,*subcomm,&send_stat);//wait for receive. 
						MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
						PetscFree(send_buf);
						to_exchange -= send_size;
						recv_ptr += send_size;
					} // end if 
				} //end if exchange > 0 
			} //end for 
		} else if ( rank > pivot_process ) {
			PetscInt to_exchange = counts[0];
			uPetscInt *recv_ptr = array;
			PetscInt to_exchange_first = 0;
			if (  g_below[pivot_process]-pivot_pos_within_chunk > 0 ) to_exchange_first += g_below[pivot_process] - pivot_pos_within_chunk ;
			for ( int i = pivot_process+1 ; i < rank ; i++ ) to_exchange_first += g_below[i];
			
			
			//now go throug the rest of the processes. 
			for ( int proc = 0; proc < pivot_process ; proc ++ ) { 
				uPetscInt *send_buf;
				int send_size;	
				if ( to_exchange > 0 ) {
					if ( to_exchange_first - g_above[proc] >= 0 ) {
						to_exchange_first -= g_above[proc];
					} else {
						if ( g_above[proc] - to_exchange_first > to_exchange ) { //if there is more to give away than is needed. 
							send_size = (int)(to_exchange);					
						} else {
							send_size = (int)(g_above[proc] - to_exchange_first);
						}
						
						if ( to_exchange_first > 0 ) { 
							//recv_ptr += (int)to_exchange_first;
							to_exchange_first = 0;
						}
						ierr = PetscMalloc(send_size * sizeof(uPetscInt),&send_buf);CHKERRQ(ierr);
						memcpy(send_buf,recv_ptr,send_size*sizeof(uPetscInt)); 
						MPI_Request send_req; MPI_Status send_stat;
						MPI_Isend(send_buf,send_size,PETSC_UNSIGNED_MPI_INT,proc,1,*subcomm,&send_req);//send using non blocking communication.
	//					PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, sending and receiving %d from %d (4).\n",rank, send_size, proc);
						MPI_Recv(recv_ptr,send_size,PETSC_UNSIGNED_MPI_INT,proc,1,*subcomm,&send_stat);//wait for receive. 
						MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
						PetscFree(send_buf);
						to_exchange -= send_size;
						recv_ptr += send_size;
					}
				}
			}
			
			
			PetscInt pivot_exchange = pivot_pos_within_chunk - g_below[pivot_process] ;
			//first check the pivot process as it is a bit anomalous.
			if ( pivot_exchange > 0 ) {
				if ( (to_exchange_first - pivot_exchange) >= 0 ) {
					to_exchange_first -= pivot_exchange; //subtract anything to be exchanged.
				} else {
					uPetscInt *send_buf;
					int send_size ;	
					if ( to_exchange > 0 ) {
						if ( (pivot_exchange - to_exchange_first) > to_exchange ) { //if there is more to give away than is needed. 
							send_size = (int)(to_exchange);	
						} else {
							send_size = (int)(pivot_exchange-to_exchange_first);
						}
						
						if ( to_exchange_first > 0 ) { 
							//recv_ptr += (int)to_exchange_first;
							to_exchange_first = 0;
						}
						ierr = PetscMalloc(send_size * sizeof(uPetscInt),&send_buf);CHKERRQ(ierr);
						memcpy(send_buf,recv_ptr,send_size*sizeof(uPetscInt)); 
						MPI_Request send_req; MPI_Status send_stat;
						MPI_Isend(send_buf,send_size,PETSC_UNSIGNED_MPI_INT,pivot_process,1,*subcomm,&send_req);//send using non blocking communication.
		//				PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, sending and receiving %d from %d (5).\n",rank, send_size, pivot_process);
						MPI_Recv(recv_ptr,send_size,PETSC_UNSIGNED_MPI_INT,pivot_process,1,*subcomm,&send_stat);//wait for receive. 
						MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
						PetscFree(send_buf);
						to_exchange -= send_size;
						recv_ptr += send_size ;
					}
				}
			}
		} //conditional on whether rank is equal , less or more than the pivot process rank and partitioning performed accordingly.
		
		
		//create new communicators
		MPI_Group processes_being_used;
		MPI_Group global_group;
		MPI_Comm  comm_world = *subcomm;
		MPI_Comm left_subcomm, right_subcomm;
		int number_relevant_ranks_left, number_relevant_ranks_right ;
		int *relevant_ranks;
		
		if ( pivot_pos_within_chunk == 0 ) {
			number_relevant_ranks_left = pivot_process;
		}else {
			number_relevant_ranks_left = pivot_process+1;
		}
		
		
		if ( number_relevant_ranks_left > 0 ) {
			PetscMalloc(sizeof(int)*number_relevant_ranks_left,&relevant_ranks);
			for ( int proc = 0 ; proc < pivot_process ; proc++ ) {
				relevant_ranks[proc] = proc;
			}
			if ( pivot_pos_within_chunk != 0 ) {
				relevant_ranks[pivot_process] = pivot_process;
			}
			
			MPI_Comm_group(comm_world,&global_group);
			MPI_Group_incl(global_group,number_relevant_ranks_left,relevant_ranks,&processes_being_used);
			MPI_Comm_create(comm_world,processes_being_used,&left_subcomm);
			MPI_Group_free(&processes_being_used);
			MPI_Group_free(&global_group);
			PetscFree(relevant_ranks);
		}
	
		
		if ( pivot_pos_within_chunk == lengths[pivot_process] ) {	
			number_relevant_ranks_right = size - pivot_process - 1;
		} else {
			number_relevant_ranks_right = size - pivot_process ;
		}
		if ( number_relevant_ranks_right > 0 ) {
			PetscMalloc(sizeof(int)*number_relevant_ranks_right,&relevant_ranks);
			int rank_idx = 0 ;
			if ( pivot_pos_within_chunk != lengths[pivot_process] ) relevant_ranks[rank_idx++] = pivot_process;
			for ( int proc = pivot_process + 1 ; proc < size ; proc++ ) {
				relevant_ranks[rank_idx++] = proc;
			}
			MPI_Comm_group(comm_world,&global_group);
			MPI_Group_incl(global_group,number_relevant_ranks_right,relevant_ranks,&processes_being_used);
			MPI_Comm_create(comm_world,processes_being_used,&right_subcomm);
			MPI_Group_free(&processes_being_used);
			MPI_Group_free(&global_group);
			PetscFree(relevant_ranks);
		}
		
		
		//partitioning done so must create new communicators and call functions recursively.
		if ( rank == pivot_process) {
			if ( number_relevant_ranks_left ) {
				if ( pivot_pos_within_chunk != 0 ) {
// 					cout << "Pivot process going to left partition of size " << number_relevant_ranks_left << " and local array size " << pivot_pos_within_chunk << ": "; 
// 					for ( int i = 0 ; i < pivot_pos_within_chunk ; i++ ) cout << array[i] << " " ;
// 					cout << endl;
					parallel_qsort_uPetscInt(&left_subcomm,array,pivot_pos_within_chunk);
				}
			}
			if ( number_relevant_ranks_right ) { 
				if ( pivot_pos_within_chunk != lengths[pivot_process] ) {
				// 	cout << "Pivot process going to right partition of size " << number_relevant_ranks_right << " and local array size " << llength-pivot_pos_within_chunk << ": ";
// 					for ( int i = 0 ; i < llength - pivot_pos_within_chunk ; i++ ) cout << array[pivot_pos_within_chunk + i] << " " ;
// 					cout << endl;
					parallel_qsort_uPetscInt(&right_subcomm,&(array[pivot_pos_within_chunk]),llength-pivot_pos_within_chunk);
				}
			}
		} else if ( rank < pivot_process) {
			if ( number_relevant_ranks_left ) {
				//cout << "Process left of pivot partition of size " << number_relevant_ranks_left << " and local array size " << llength << endl;
				parallel_qsort_uPetscInt(&left_subcomm,array,llength);
			}
		} else if ( rank > pivot_process ) {
			if ( number_relevant_ranks_right ) { 
				//cout << "Process right of pivot partition of size " << number_relevant_ranks_right << " and local array size " << llength << ": " ;
				//for ( int i = 0 ; i < llength  ; i++ ) cout << array[i] << " " ;
				//cout << endl;
				parallel_qsort_uPetscInt(&right_subcomm,array,llength);
			}
		}
		
		if ( left_subcomm != MPI_COMM_NULL ) {
			MPI_Barrier(left_subcomm);
			MPI_Comm_free(&left_subcomm);
		}
		if ( right_subcomm != MPI_COMM_NULL ) {
			MPI_Barrier(right_subcomm);
			MPI_Comm_free(&right_subcomm);
		}
		
		ierr = PetscFree(lengths);CHKERRQ(ierr);
		ierr = PetscFree(g_below);CHKERRQ(ierr);
		ierr = PetscFree(g_above);CHKERRQ(ierr);
	} //end if there is more than one process.	
	return 0;
}

/*! \brief Function to perform parallel quicksort to sort an array of unsigned PetscInts.*/
template <class T>
int parallel_qsort(MPI_Comm *subcomm,T *array, PetscInt llength) {
	//get size of sub communicator as well as the rank of the current process within it.
	PetscMPIInt rank, size ;
	PetscErrorCode ierr;
	ierr = MPI_Comm_rank(*subcomm,&rank);CHKERRQ(ierr);
	ierr = MPI_Comm_size(*subcomm,&size);CHKERRQ(ierr);
	
	//PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, size %d and length %d.\n",rank, size, llength);	
	
	//first we sort the local portions of the array
	sort(array,array+llength);
	
	//If the communicator is only made up of a single process then we are done.
	if ( size > 1 ) {  
		PetscInt *lengths; //array to store the lengths of the other chunks.
		ierr = PetscMalloc(size*sizeof(PetscInt),&lengths);CHKERRQ(ierr);
		MPI_Allgather(&llength,1,PETSC_MPI_INT,lengths,1,PETSC_MPI_INT,*subcomm); //gather lengths of all chunks on each process.
		PetscInt total_length = 0;
		for ( int i = 0 ; i < size ; i++ ) total_length += lengths[i]; //find the total length of the subset.
		
		//select the pivot as the central item in the sorted list.
		PetscReal pivot = 0 ;
		if ( llength > 0 ) 
		  {
		    pivot = (PetscReal)array[llength/2].sorting_token ;
		  }  
		
		//sum the pivots among all the processes in the communicator.
		PetscReal sum = 0.0; 
		MPI_Allreduce(&pivot, &sum, 1,MPIU_REAL, MPI_SUM, *subcomm); 
		//divide by the number of processes to get the mean.
		//pivot = sum/total_length;
		PetscReal real_pivot = (PetscReal)sum / (PetscReal)size;
		
		PetscInt counts[2] = {0,0}, *g_below, *g_above; //array containing counts of array elements below and above the pivot.
		ierr = PetscMalloc(size*sizeof(PetscInt),&g_below);CHKERRQ(ierr);  //allocat space to store the numbers above and below the pivot on each process.
		ierr = PetscMalloc(size*sizeof(PetscInt),&g_above);CHKERRQ(ierr);
		{
			PetscInt i = 0;
			while ( i < llength && ((PetscReal)array[i].sorting_token) < real_pivot ) {
				if ( ((PetscReal)array[i].sorting_token) < real_pivot ) counts[0]++;
				i++;
			}
			counts[1] = llength - i;
		}
		MPI_Allgather(&counts[0],1,PETSC_MPI_INT,g_below,1,PETSC_MPI_INT,*subcomm); 
		MPI_Allgather(&counts[1],1,PETSC_MPI_INT,g_above,1,PETSC_MPI_INT,*subcomm); 
		
		//find the total number of elements less than the pivot and the total amount more than the pivot.
		PetscInt total_below = 0, total_above = 0 ;
		for ( int i = 0 ; i < size ; i++ ) {
			total_below += g_below[i];
			total_above += g_above[i];
		}
		
		//find the rank of the process which contains the pivot.
		int pivot_process = 0;
		PetscInt pivot_pos_within_chunk = 0;
		PetscInt l = 0;
		for ( int i = 0 ; i < size ; i++ ) {
			if ( total_below <= (l + lengths[i]) ){
				pivot_process = i;
				pivot_pos_within_chunk = total_below - l; 
				break;
			}
			l += lengths[i];
		}
		
		//check if there are any elements to exchange and if not we can finish here.
		/*PetscInt number_to_exchange_below = 0 ;
		 for ( int i = 0 ; i < pivot_process ; i++ ) {
		 number_to_exchange_below += g_above[i];
		 }
		 PetscInt number_to_exchange_above = 0 ;
		 for ( int i = size-1 ; i > pivot_process ; i-- ) {
		 number_to_exchange_above += g_below[i];
		 }
		 
		 if ( number_to_exchange_above == 0 && number_to_exchange_below == 0  ) {
		 ierr = PetscFree(lengths);CHKERRQ(ierr);
		 ierr = PetscFree(g_below);CHKERRQ(ierr);
		 ierr = PetscFree(g_above);CHKERRQ(ierr);
		 return 0;
		 }*/
		//		PetscPrintf(PETSC_COMM_SELF,"Pivot process %d.\n",pivot_process);
		
		if ( rank == pivot_process ) { //if this is the pivot process.
			//find the total length of the subset
			l = 0 ;
			for ( int i = 0 ; i < pivot_process ; i++ ) {
				l += lengths[i];
			}
			PetscInt end_of_lower = total_below - l; //the local index of the end of the lower partition.
			//find the total length of the sites below the pivot
			l = 0;
			for ( int i = size-1 ; i > pivot_process ; i-- ) {
				l += lengths[i];
			}
			PetscInt start_of_upper = llength - (total_above - l); //the local index of the start of the upper index.
			//if ( start_of_upper != end_of_lower) cout << "Ends dont meet : " <<  end_of_lower << ", " << start_of_upper << endl;
			
			if ( end_of_lower > counts[0] ) { // if this process has too many elements that are greater than or equal to the pivot. Need to exchange these elements for ones that are below the pivot.
				PetscInt to_exchange = end_of_lower - counts[0]; //the number of items there are to exchange.
				PetscInt to_exchange_first = 0;
				for ( int i = 0 ; i < pivot_process ; i++ ) to_exchange_first += g_above[i]; //count up number each process below has to share first.
				T *recv_ptr = &(array[counts[0]]);  
				for ( int proc = pivot_process+1 ; proc < size ; proc++ ) { //loop through processes until we find the one to exchange with.
					if ( (to_exchange_first - g_below[proc]) >= 0 ) {
						to_exchange_first -= g_below[proc] ; //subtract anything to be exchanged between other processes.
					} else {
						T *send_buf;
						int send_size; 
						if ( to_exchange > 0 ) {
							if ( g_below[proc] - to_exchange_first > to_exchange ) { //if there is more to give away than is needed. 
								send_size = (int)(to_exchange);					
							} else {
								send_size = (int)(g_below[proc] - to_exchange_first);
							}
							if ( to_exchange_first > 0 ) { 
								//recv_ptr += (int)to_exchange_first;
								to_exchange_first = 0;
							}
							ierr = PetscMalloc(send_size * sizeof(T),&send_buf);CHKERRQ(ierr);
							memcpy(send_buf,recv_ptr,send_size*sizeof(T)); 
							MPI_Request send_req; MPI_Status send_stat;
							MPI_Isend((char *)send_buf,send_size*sizeof(T),MPI_BYTE,proc,1,*subcomm,&send_req);//send using non blocking communication.
							//						PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, sending and receiving %d from %d (0).\n",rank, send_size, proc);
							MPI_Recv((char *)recv_ptr,send_size*sizeof(T),MPI_BYTE,proc,1,*subcomm,&send_stat);//wait for receive. 
							MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
							PetscFree(send_buf);
							to_exchange -= send_size;
							recv_ptr += send_size;
						} //end if no_exchange > 0 
					} //ebd if 
				}  //end for over processes
			} else if ( end_of_lower < counts[0] ) { //otherwise if this process has too many elements that are lower than the pivot. Need to exchange these elements for ones that are above the pivot.
				PetscInt to_exchange = counts[0] - end_of_lower; //the number of items there are to exchange.
				uPetscInt *recv_ptr = &(array[end_of_lower].sorting_token);
				for ( int proc = 0 ; proc < pivot_process ; proc++ ) { //loop through processes before this one.
					T *send_buf;
					int send_size; 
					if ( to_exchange > 0  && g_above[proc] > 0) {
						if ( g_above[proc] > to_exchange ) { //if there is more to give away than is needed. 
							send_size = (int)to_exchange;					
						} else {
							send_size = (int)g_above[proc];
						}
						ierr = PetscMalloc(send_size * sizeof(T),&send_buf);CHKERRQ(ierr);
						memcpy(send_buf,recv_ptr,send_size*sizeof(T)); 
						MPI_Request send_req; MPI_Status send_stat;
						MPI_Isend((char*)send_buf,send_size*sizeof(T),MPI_BYTE,proc,1,*subcomm,&send_req);//send using non blocking communication.
						//						PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, sending and receiving %d from %d (1).\n",rank, send_size, proc);
						MPI_Recv((char*)recv_ptr,send_size*sizeof(T),MPI_BYTE,proc,1,*subcomm,&send_stat);//wait for receive. 
						MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
						PetscFree(send_buf);
						to_exchange -= send_size;
						recv_ptr += send_size;
					}
				}
			}
		} else if ( rank < pivot_process ) { // means that everything above must be exchanged for elements below.
			PetscInt to_exchange = counts[1];
			T *recv_ptr = &(array[counts[0]]);
			PetscInt to_exchange_first = 0;
			for ( int i = 0 ; i < rank ; i++ ) to_exchange_first += g_above[i];
			
			//first check the pivot process as it is a bit anomalous.
			if ( g_below[pivot_process] - pivot_pos_within_chunk > 0 ) {
				if ( (to_exchange_first - (g_below[pivot_process] - pivot_pos_within_chunk) ) >= 0 ) { 
					to_exchange_first -= (g_below[pivot_process] - pivot_pos_within_chunk); //subtract anything to be exchanged.
				} else {
					T *send_buf;
					int send_size ;	
					if ( to_exchange > 0 ) {
						if ( ((g_below[pivot_process] - pivot_pos_within_chunk) - to_exchange_first) > to_exchange ) { //if there is more to give away than is needed. 
							send_size = (int)(to_exchange);	
						} else {
							send_size = (int)(g_below[pivot_process]-pivot_pos_within_chunk-to_exchange_first);
							to_exchange_first = 0;
						}
						ierr = PetscMalloc(send_size * sizeof(T),&send_buf);CHKERRQ(ierr);
						memcpy(send_buf,recv_ptr,send_size*sizeof(T)); 
						MPI_Request send_req; MPI_Status send_stat;
						MPI_Isend((char*)send_buf,send_size*sizeof(T),MPI_BYTE,pivot_process,1,*subcomm,&send_req);//send using non blocking communication.
						//					PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, sending and receiving %d from %d (2).\n",rank, send_size, pivot_process);
						MPI_Recv((char*)recv_ptr,send_size*sizeof(T),MPI_BYTE,pivot_process,1,*subcomm,&send_stat);//wait for receive. 
						MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
						PetscFree(send_buf);
						to_exchange -= send_size;
						recv_ptr += send_size ;
					} // end if to_exchange > 0 
				} // end if else there is enough to cover what needs to be exchanged first
			} //end of check of pivot process
			
			//now go throug the rest of the processes. 
			for ( int proc = pivot_process+1; proc < size ; proc ++ ) { 
				T *send_buf;
				int send_size;	
				if ( to_exchange > 0 ) {
					if ( to_exchange_first - g_below[proc] >= 0 ) {
						to_exchange_first -= g_below[proc];
					} else {
						if ( g_below[proc] - to_exchange_first >= to_exchange ) { //if there is more to give away than is needed. 
							send_size = (int)(to_exchange);					
						} else {
							send_size = (int)(g_below[proc] - to_exchange_first);
						}
						
						if ( to_exchange_first > 0 ) { 
							//recv_ptr += (int)to_exchange_first;
							to_exchange_first = 0;
						}
						ierr = PetscMalloc(send_size * sizeof(T),&send_buf);CHKERRQ(ierr);
						memcpy(send_buf,recv_ptr,send_size*sizeof(T)); 
						MPI_Request send_req; MPI_Status send_stat;
						MPI_Isend((char*)send_buf,send_size*sizeof(T),MPI_BYTE,proc,1,*subcomm,&send_req);//send using non blocking communication.
						//				PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, sending and receiving %d from %d (3).\n",rank, send_size, proc);				
						MPI_Recv((char*)recv_ptr,send_size*sizeof(T),MPI_BYTE,proc,1,*subcomm,&send_stat);//wait for receive. 
						MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
						PetscFree(send_buf);
						to_exchange -= send_size;
						recv_ptr += send_size;
					} // end if 
				} //end if exchange > 0 
			} //end for 
		} else if ( rank > pivot_process ) {
			PetscInt to_exchange = counts[0];
			T *recv_ptr = array;
			PetscInt to_exchange_first = 0;
			if (  g_below[pivot_process]-pivot_pos_within_chunk > 0 ) to_exchange_first += g_below[pivot_process] - pivot_pos_within_chunk ;
			for ( int i = pivot_process+1 ; i < rank ; i++ ) to_exchange_first += g_below[i];
			
			
			//now go throug the rest of the processes. 
			for ( int proc = 0; proc < pivot_process ; proc ++ ) { 
				T *send_buf;
				int send_size;	
				if ( to_exchange > 0 ) {
					if ( to_exchange_first - g_above[proc] >= 0 ) {
						to_exchange_first -= g_above[proc];
					} else {
						if ( g_above[proc] - to_exchange_first > to_exchange ) { //if there is more to give away than is needed. 
							send_size = (int)(to_exchange);					
						} else {
							send_size = (int)(g_above[proc] - to_exchange_first);
						}
						
						if ( to_exchange_first > 0 ) { 
							//recv_ptr += (int)to_exchange_first;
							to_exchange_first = 0;
						}
						ierr = PetscMalloc(send_size * sizeof(T),&send_buf);CHKERRQ(ierr);
						memcpy(send_buf,recv_ptr,send_size*sizeof(T)); 
						MPI_Request send_req; MPI_Status send_stat;
						MPI_Isend((char*)send_buf,send_size*sizeof(T),MPI_BYTE,proc,1,*subcomm,&send_req);//send using non blocking communication.
						//					PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, sending and receiving %d from %d (4).\n",rank, send_size, proc);
						MPI_Recv((char*)recv_ptr,send_size*sizeof(T),MPI_BYTE,proc,1,*subcomm,&send_stat);//wait for receive. 
						MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
						PetscFree(send_buf);
						to_exchange -= send_size;
						recv_ptr += send_size;
					}
				}
			}
			
			
			PetscInt pivot_exchange = pivot_pos_within_chunk - g_below[pivot_process] ;
			//first check the pivot process as it is a bit anomalous.
			if ( pivot_exchange > 0 ) {
				if ( (to_exchange_first - pivot_exchange) >= 0 ) {
					to_exchange_first -= pivot_exchange; //subtract anything to be exchanged.
				} else {
					T *send_buf;
					int send_size ;	
					if ( to_exchange > 0 ) {
						if ( (pivot_exchange - to_exchange_first) > to_exchange ) { //if there is more to give away than is needed. 
							send_size = (int)(to_exchange);	
						} else {
							send_size = (int)(pivot_exchange-to_exchange_first);
						}
						
						if ( to_exchange_first > 0 ) { 
							//recv_ptr += (int)to_exchange_first;
							to_exchange_first = 0;
						}
						ierr = PetscMalloc(send_size * sizeof(T),&send_buf);CHKERRQ(ierr);
						memcpy(send_buf,recv_ptr,send_size*sizeof(T)); 
						MPI_Request send_req; MPI_Status send_stat;
						MPI_Isend((char*)send_buf,send_size*sizeof(T),MPI_BYTE,pivot_process,1,*subcomm,&send_req);//send using non blocking communication.
						//				PetscFPrintf(PETSC_COMM_SELF,stderr,"Rank %d, sending and receiving %d from %d (5).\n",rank, send_size, pivot_process);
						MPI_Recv((char*)recv_ptr,send_size*sizeof(T),MPI_BYTE,pivot_process,1,*subcomm,&send_stat);//wait for receive. 
						MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
						PetscFree(send_buf);
						to_exchange -= send_size;
						recv_ptr += send_size ;
					}
				}
			}
		} //conditional on whether rank is equal , less or more than the pivot process rank and partitioning performed accordingly.
		
		
		//create new communicators
		MPI_Group processes_being_used;
		MPI_Group global_group;
		MPI_Comm  comm_world = *subcomm;
		MPI_Comm left_subcomm, right_subcomm;
		int number_relevant_ranks_left, number_relevant_ranks_right ;
		int *relevant_ranks;
		
		if ( pivot_pos_within_chunk == 0 ) {
			number_relevant_ranks_left = pivot_process;
		}else {
			number_relevant_ranks_left = pivot_process+1;
		}
		
		
		if ( number_relevant_ranks_left > 0 ) {
			PetscMalloc(sizeof(int)*number_relevant_ranks_left,&relevant_ranks);
			for ( int proc = 0 ; proc < pivot_process ; proc++ ) {
				relevant_ranks[proc] = proc;
			}
			if ( pivot_pos_within_chunk != 0 ) {
				relevant_ranks[pivot_process] = pivot_process;
			}
			
			MPI_Comm_group(comm_world,&global_group);
			MPI_Group_incl(global_group,number_relevant_ranks_left,relevant_ranks,&processes_being_used);
			MPI_Comm_create(comm_world,processes_being_used,&left_subcomm);
			MPI_Group_free(&processes_being_used);
			MPI_Group_free(&global_group);
			PetscFree(relevant_ranks);
		}
		
		
		if ( pivot_pos_within_chunk == lengths[pivot_process] ) {	
			number_relevant_ranks_right = size - pivot_process - 1;
		} else {
			number_relevant_ranks_right = size - pivot_process ;
		}
		if ( number_relevant_ranks_right > 0 ) {
			PetscMalloc(sizeof(int)*number_relevant_ranks_right,&relevant_ranks);
			int rank_idx = 0 ;
			if ( pivot_pos_within_chunk != lengths[pivot_process] ) relevant_ranks[rank_idx++] = pivot_process;
			for ( int proc = pivot_process + 1 ; proc < size ; proc++ ) {
				relevant_ranks[rank_idx++] = proc;
			}
			MPI_Comm_group(comm_world,&global_group);
			MPI_Group_incl(global_group,number_relevant_ranks_right,relevant_ranks,&processes_being_used);
			MPI_Comm_create(comm_world,processes_being_used,&right_subcomm);
			MPI_Group_free(&processes_being_used);
			MPI_Group_free(&global_group);
			PetscFree(relevant_ranks);
		}
		
		
		//partitioning done so must create new communicators and call functions recursively.
		if ( rank == pivot_process) {
			if ( number_relevant_ranks_left ) {
				if ( pivot_pos_within_chunk != 0 ) {
					// 					cout << "Pivot process going to left partition of size " << number_relevant_ranks_left << " and local array size " << pivot_pos_within_chunk << ": "; 
					// 					for ( int i = 0 ; i < pivot_pos_within_chunk ; i++ ) cout << array[i] << " " ;
					// 					cout << endl;
					parallel_qsort(&left_subcomm,array,pivot_pos_within_chunk);
				}
			}
			if ( number_relevant_ranks_right ) { 
				if ( pivot_pos_within_chunk != lengths[pivot_process] ) {
					// 	cout << "Pivot process going to right partition of size " << number_relevant_ranks_right << " and local array size " << llength-pivot_pos_within_chunk << ": ";
					// 					for ( int i = 0 ; i < llength - pivot_pos_within_chunk ; i++ ) cout << array[pivot_pos_within_chunk + i] << " " ;
					// 					cout << endl;
					parallel_qsort(&right_subcomm,&(array[pivot_pos_within_chunk]),llength-pivot_pos_within_chunk);
				}
			}
		} else if ( rank < pivot_process) {
			if ( number_relevant_ranks_left ) {
				//cout << "Process left of pivot partition of size " << number_relevant_ranks_left << " and local array size " << llength << endl;
				parallel_qsort(&left_subcomm,array,llength);
			}
		} else if ( rank > pivot_process ) {
			if ( number_relevant_ranks_right ) { 
				//cout << "Process right of pivot partition of size " << number_relevant_ranks_right << " and local array size " << llength << ": " ;
				//for ( int i = 0 ; i < llength  ; i++ ) cout << array[i] << " " ;
				//cout << endl;
				parallel_qsort(&right_subcomm,array,llength);
			}
		}
		
		if ( left_subcomm != MPI_COMM_NULL ) {
			MPI_Barrier(left_subcomm);
			MPI_Comm_free(&left_subcomm);
		}
		if ( right_subcomm != MPI_COMM_NULL ) {
			MPI_Barrier(right_subcomm);
			MPI_Comm_free(&right_subcomm);
		}
		
		ierr = PetscFree(lengths);CHKERRQ(ierr);
		ierr = PetscFree(g_below);CHKERRQ(ierr);
		ierr = PetscFree(g_above);CHKERRQ(ierr);
	} //end if there is more than one process.	
	return 0;
}


uPetscInt translate(uPetscInt conf, PetscInt *site_indices, int count) { 
	uPetscInt new_conf = 0; 
	
	for ( int i = 0 ; i < count ; i++ ) {
		if ( (conf & (((uPetscInt)1) << site_indices[i])) > (uPetscInt)0 ) {
			new_conf += (((uPetscInt)1) << i);
		}
	}	
	return new_conf;
}

uPetscInt translate_back(uPetscInt conf, PetscInt *site_indices, int count) { 
	uPetscInt new_conf = 0; 
	
	for ( int i = 0 ; i < count ; i++ ) {
		if ( (conf & (((uPetscInt)1) << i)) > (uPetscInt)0 ) {
			new_conf += (((uPetscInt)1) << site_indices[i]);
		}
	}	
	return new_conf;
}

PetscInt build_config_list(struct parameters *params, uPetscInt **array, PetscTruth parallel, PetscInt *l_cnt, PetscInt *g_cnt, uPetscInt mask, PetscInt filling, PetscInt max_filling) {
	struct conf_list_info c_info; //to store config info for this particular list	
	c_info.filling = filling;	//filling to use or -1 if filling not defined on the subsection.
	c_info.max_filling = max_filling; //if filling not set then use max filling to keep the size down. 
	PetscErrorCode ierr;
//	c_info.num_sites = __builtin_popcount(mask); //this builtin function does a fast count on the number of 1 bits in the mask.
	c_info.num_sites = 0;
	for ( int i = 0 ; i < params->AB_sites ; i++ ) {
		if ( (mask & ((uPetscInt)1) << i ) > (uPetscInt) 0 ) {
			c_info.num_sites++;
		}
	}
	
	int count = 0;
	PetscInt *site_indices;
	ierr = PetscMalloc(c_info.num_sites*sizeof(PetscInt),&site_indices);CHKERRQ(ierr);
	int pos = 0;
	while ( count < c_info.num_sites && pos < params->AB_sites) {
		if ( (mask & (((uPetscInt)1) << pos) ) > (uPetscInt)0 )  {
			site_indices[count++] = pos;
		}
		pos++;
	}
	
	if ( params->nn_exclusion ) {
		//count the number of masks that apply here. 
		list<uPetscInt> relevant_masks;
		c_info.adj_info.number_masks = 0;
		for ( int i = 0 ; i < params->adj_info.number_masks ; i++ ) {
			if ( (((params->adj_info.masks[i] & mask) & ((params->adj_info.masks[i] & mask)-1)) > (uPetscInt)0) ) {	
				/*char buf[MAX_STRING];
				dec2bin(params->adj_info.masks[i],params->AB_sites, buf);
				PetscPrintf(PETSC_COMM_WORLD,"Before: %s\n",buf);*/
				relevant_masks.push_back(translate(params->adj_info.masks[i],site_indices,count)); c_info.adj_info.number_masks++;
				/*dec2bin(translate(params->adj_info.masks[i],site_indices,count),params->AB_sites, buf);
				PetscPrintf(PETSC_COMM_WORLD,"After: %s\n",buf);*/
			} /*else {
				char buf[MAX_STRING];
				dec2bin( params->adj_info.masks[i], params->AB_sites, buf );
				PetscPrintf(PETSC_COMM_WORLD,"Unsuitable mask: %s\n",buf);
			}*/
		}
		struct adjacency_information *adj_info = &c_info.adj_info ;
		ierr = PetscMalloc(adj_info->number_masks * sizeof(uPetscInt), &(adj_info->masks) ); CHKERRQ(ierr);
		for ( int i = 0 ; i < c_info.adj_info.number_masks ; i++ ) {
			c_info.adj_info.masks[i] = relevant_masks.front() ; relevant_masks.pop_front();
			/*char buf[MAX_STRING];
			dec2bin(c_info.adj_info.masks[i], params->AB_sites, buf );
			PetscPrintf(PETSC_COMM_WORLD,"Mask: %s %d\n",buf, c_info.adj_info.number_masks);*/
		}
	
		if ( parallel ) {
			string basis_file_prefix("ee_basis_list_portion_");
			stringstream ss;
			ss << basis_file_prefix << params->rank << ".bas" ; 
			string basis_filename = ss.str();
			ofstream basis_file( basis_filename.c_str(),ios::out | ios::binary);
			find_valid_configs(params, &c_info, &basis_file, l_cnt);
			//PetscPrintf(PETSC_COMM_SELF,"Rank %d found %d configs.\n",params->rank,*l_cnt);
			basis_file.close();		
			MPI_Allreduce(l_cnt, g_cnt, 1, PETSC_MPI_INT, MPI_SUM, PETSC_COMM_WORLD); //sum the number local numbers. 
			
			if (*g_cnt > 0 ) {
				PetscInt *all_counts;
				//PetscPrintf(PETSC_COMM_SELF,"Rank %d total configs found: %d.\n",params->rank,*g_cnt);
				ierr = PetscMalloc((params->size+1)*sizeof(PetscInt), &all_counts);CHKERRQ(ierr);
				MPI_Allgather(l_cnt,1,PETSC_MPI_INT,&all_counts[1],1,PETSC_MPI_INT,PETSC_COMM_WORLD); //gather the number of reps found on each process.
				ierr = PetscMalloc(get_rank_chunk_size(*g_cnt,params->rank,params->size)*sizeof(uPetscInt),array);CHKERRQ(ierr);
				all_counts[0] = 0; 
				for ( int i = 1 ; i <= params->size ; i++ ) {
					all_counts[i] += all_counts[i-1];
				}
				PetscInt offset = 0;
				//need to now find which process and offset to start at.			
				PetscInt to_fill = get_rank_chunk_size(*g_cnt,params->rank,params->size); //keep account of how many array places are left to fill. 
				PetscInt proc = 0; //start at the first process.
				while (to_fill > 0 ) {
					while ( (get_rank_chunk_start(*g_cnt,params->rank,params->size)+offset) >= all_counts[proc+1] ) { //while the starting index is not within this range move on.
						proc++; 
					}
					//open the file corresponding to processes proc
					ss.str("");
					ss << basis_file_prefix << proc << ".bas" ;
					basis_filename = ss.str();
					ifstream basis_file_in ( basis_filename.c_str(),ios::in | ios::binary );
					//seek if we need to
					PetscInt available = all_counts[proc+1] - all_counts[proc];
					if ( (get_rank_chunk_start(*g_cnt,params->rank,params->size)+offset) > all_counts[proc] ) {
						basis_file_in.seekg((get_rank_chunk_start(*g_cnt,params->rank,params->size)+offset-all_counts[proc]) * sizeof(uPetscInt),ios_base::beg);
						available -= (get_rank_chunk_start(*g_cnt,params->rank,params->size)+offset-all_counts[proc]);
					}
					//now read in as many elements as we want
					if ( to_fill < available ) {
						basis_file_in.read((char *) &(*array)[offset], sizeof(uPetscInt) * to_fill );
						to_fill -= to_fill;
						offset += to_fill;
					} else {
						basis_file_in.read((char *) &(*array)[offset], sizeof(uPetscInt) * available );
						offset += available;
						to_fill -= available;
					}
					basis_file_in.close();
				}
				
				*l_cnt = get_rank_chunk_size(*g_cnt,params->rank,params->size);
				
				MPI_Group processes_being_used;
				MPI_Group global_group;
				MPI_Comm  comm_world = PETSC_COMM_WORLD;
				MPI_Comm new_subcomm;
				int number_relevant_ranks;
				int *relevant_ranks;
				
				if ( *g_cnt >= params->size ) {
					number_relevant_ranks = params->size;
				} else {
					number_relevant_ranks = *g_cnt;
				}
				
				PetscMalloc(sizeof(int)*number_relevant_ranks,&relevant_ranks);
				for ( int proc = 0 ; proc < number_relevant_ranks ; proc++ ) {
					relevant_ranks[proc] = proc;
				}
				MPI_Comm_group(comm_world,&global_group);
				MPI_Group_incl(global_group,number_relevant_ranks,relevant_ranks,&processes_being_used);
				MPI_Comm_create(comm_world,processes_being_used,&new_subcomm);
				MPI_Group_free(&processes_being_used);
				MPI_Group_free(&global_group);
				PetscFree(relevant_ranks);
		
				if ( *l_cnt > 0 ) {
					parallel_qsort_uPetscInt(&new_subcomm,*array,*l_cnt);
				}
				
				if ( new_subcomm != MPI_COMM_NULL ) {
					MPI_Barrier(new_subcomm);
					MPI_Comm_free(&new_subcomm);
				}
				PetscFree(all_counts);
			}
			
		} else {
			list<uPetscInt> ll;
			find_valid_configs_list(params, &c_info,&ll);
			*l_cnt = ll.size();
			*g_cnt = *l_cnt;
			ierr = PetscMalloc(*g_cnt * sizeof(uPetscInt),array);CHKERRQ(ierr);
			for ( PetscInt i = 0 ; i < *g_cnt; i++ ) {
				(*array)[i] = ll.front() ; ll.pop_front();
			}
			sort(*array, *array + *l_cnt );
		}
		//going to translate configurations back according to the mask. 
		for ( PetscInt i = 0 ; i < *l_cnt; i++ ) {
			(*array)[i] = translate_back((*array)[i], site_indices, count);
		}
	
	} else {
		//come back to this later when needed.
	}

	PetscFree(site_indices);
	
	return 1; 
}

//find the start of the range by using divide and conquer approach to find a single element and then descend to beginning of the range.
PetscInt find_start_of_range(class my_vec_entry *entries, PetscInt count, uPetscInt lower_limit, uPetscInt upper_limit, uPetscInt mask, PetscInt *size) {
	PetscInt bottom = 0 , top = count;
	PetscInt pos =  count / 2; 
	PetscInt start;
	
	//keep dividing by 2 until we are in the vincinity.
	while ( ((entries[pos].conf & mask) <= lower_limit 
			 			|| (entries[pos].conf & mask) > upper_limit) && (top - bottom) > 1 ) { //if its not in the range
		if ((entries[pos].conf & mask) <= lower_limit ) {
			bottom = pos;
			pos = (top + bottom)/2; 
		} else {
			top = pos;
			pos = (top + bottom)/2; 
		}
	}
	
	if (   ((lower_limit == 0 && (entries[pos].conf & mask) >= lower_limit) || (entries[pos].conf & mask) > lower_limit ) 
		&& (entries[pos].conf & mask) <= upper_limit ) {
		start = pos;
		while (start > 0 &&  ( ( (entries[start-1].conf & mask) > lower_limit) ||  (lower_limit == 0 && ( (entries[start-1].conf & mask) >= lower_limit)) )  ) {
			start--;
		}
		*size = pos - start + 1 ;
		while ( pos < (count-1) && (entries[pos+1].conf & mask) <= upper_limit ) {
			pos++;
			(*size) += 1;
		}
		return start; 
	} else if (((lower_limit == 0 && (entries[bottom].conf & mask) >= lower_limit) || (entries[bottom].conf & mask) > lower_limit ) 
		&& (entries[bottom].conf & mask) <= upper_limit ) {
		start = bottom;
		while (start > 0 &&  ( ( (entries[start-1].conf & mask) > lower_limit) ||  (lower_limit == 0 && ( (entries[start-1].conf & mask) >= lower_limit)) )  ) {
			start--;
		}
		*size = bottom - start + 1 ;
		while ( bottom < (count-1) && (entries[bottom+1].conf & mask) <= upper_limit ) {
			bottom++;
			(*size) += 1;
		}
		return start; 
	} else {
		*size = 0;
	}
	return 0;
}


/*! \brief This function reads the input file containing the details of the transformation being used.*/
int read_rotation_information(struct parameters *params){
	string filename(params->rotation_info.filename);
	PetscErrorCode ierr;
	
	//open rotation file as xml document.
	TiXmlDocument doc(filename);
	bool loadOkay = doc.LoadFile();
	if ( !loadOkay ) 
	{
		PetscPrintf(PETSC_COMM_WORLD,"Could not load rotation info file '%s'. Error='%s'. Exiting.\n", params->rotation_info.filename, doc.ErrorDesc() );
		return 1;
	}
	TiXmlHandle docHandle(&doc);
	
	//find ROTATION_OP tag.
	TiXmlElement* element = docHandle.FirstChild("ROTATION_OP").Element(); //tags surrounding the mapping. Has attribute "number" which contains the 																   //number of times the rotation op can be applied before getting back to the start.												
	TiXmlElement* element2;
	if ( element) 
	{	
		if (element->Attribute("number") == NULL ) { 
			PetscPrintf(PETSC_COMM_WORLD,"You must specify a number attribute to the rotation operation.\n");
			return 1;
		}
		params->rotation_info.num_sectors = atoi(element->Attribute("number"));
		params->rotation_info.num_relevant_sectors = params->rotation_info.num_sectors;
		
		//allocate space for the mapping masks.
		ierr = PetscMalloc( params->rotation_info.num_sectors * sizeof(uPetscInt*),&params->rotation_info.mapping_masks);CHKERRQ(ierr);
		for ( int i = 0 ; i < params->rotation_info.num_sectors ; i++ ) {
			ierr = PetscMalloc(params->number_particles * sizeof(uPetscInt), &params->rotation_info.mapping_masks[i]);CHKERRQ(ierr);
		}
		
		//applying 0 rotations will give the same stage so add masks that do not move any particles.
		for ( int i = 0 ; i < params->number_particles ; i++ ) {
			params->rotation_info.mapping_masks[0][i] = ((uPetscInt)1) << i;
		}
		
		//read in the mappings to position 1 which is a single rotation with the op.
		for ( int i = 0 ; i < params->number_particles ; i++ ) {
			element2 = docHandle.FirstChild("ROTATION_OP").Child("MAPPING",i).Element();
			if ( element2 ) {
				params->rotation_info.mapping_masks[1][atoi(element2->Attribute("to"))-1] = ((uPetscInt)1) << (atoi(element2->Attribute("from"))-1);
			} else {
				PetscPrintf(PETSC_COMM_WORLD,"Mapping %d does not exist.\n",i);
				return 1;
			}
		}
		
		//Now iterate this rotation to fill in the masks for the other possible mutiples of the rotation.
		for ( int i = 2 ; i < params->rotation_info.num_sectors ; i++ ) {
			for ( int j = 0 ; j < params->number_particles ; j++ ) {
				params->rotation_info.mapping_masks[i][j] = params->rotation_info.mapping_masks[1][log2_uPetscInt(params->rotation_info.mapping_masks[i-1][j])];
			}
		}
		
		//print out the mapping for each number of rotations that can be applied. For debugging purposes. 
		/*for (int i = 0 ; i < params->rotation_info.num_sectors ; i++ ) {
		 PetscPrintf(PETSC_COMM_WORLD,"Applying %d rotations:\n",i);
		 for ( int j = 0 ; j < params->number_particles ; j++ ) {
		 PetscPrintf(PETSC_COMM_WORLD,"Site %d moved to site %d.\n",j,(int)log2_uPetscInt(params->rotation_info.mapping_masks[i][j]));
		 }
		 }*/
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"No ROTATION_OP tag in rotation info file.\n");
		return 1;
	}
	return 0;
}

int add_rotationally_related(struct parameters *params,uPetscInt conf,PetscScalar val,class my_vec_entry **curr_entry) {
	struct rotation_information *rot_info;
	rot_info = &params->rotation_info;
	class my_vec_entry *l_vec_ptr = *curr_entry;
	PetscScalar sign;
	int num = 0; 
	PetscTruth found;
	uPetscInt new_basis;
	
	for ( int amount = 0 ; amount < rot_info->num_sectors ; amount++ ){ //we go through the rotation sectors
		new_basis = apply_rotation(params,conf,&amount,&sign); //apply the rotation.
		found = PETSC_FALSE;
		for ( int i = 0 ; i < num ; i++ ) { //go through the previous ones added to make sure this configuration is not already present.
			if ( new_basis == l_vec_ptr[i].conf ) { //if its found we add the phase to the existing value
				l_vec_ptr[i].val += sign * rot_info->phases[amount]; //the val field is used to hold the phase on the first time around.
				found = PETSC_TRUE; 	  //indicate that entry already added. 
			}
		}
		if ( found == PETSC_FALSE ) { //if its not found we add a new entry.
			l_vec_ptr[num].conf = new_basis;	//add new basis
			l_vec_ptr[num].mask = (uPetscInt)0;	//set mask to 0 for now.
			l_vec_ptr[num++].val = sign * rot_info->phases[amount];		//val used to hold phase so initialise to initial phase.
		}
	}
	
	//we now go over this and normalise
	PetscReal norm = 1.0/sqrt((PetscReal)num); //global normalisation is the square root of the number of individual terms.
	//for ( int amount = 0 ; amount < rot_info->num_sectors ; amount++ ){ 
	//	new_basis = apply_rotation(params,conf,&amount,&sign);
		for ( int i = 0 ; i < num ; i++ ) {
	//		if ( new_basis == l_vec_ptr[i].conf && l_vec_ptr[i].mask == (uPetscInt)0 ) {
				if (PetscAbsScalar(l_vec_ptr[i].val) > 1.0e-5 ) {
					l_vec_ptr[i].val /= PetscAbsScalar(l_vec_ptr[i].val); //normalise the phase
					l_vec_ptr[i].val *= val*norm ;
				}
				l_vec_ptr[i].mask = params->B_mask;
	//		}
		}
	//}
	*curr_entry = &l_vec_ptr[num];
	return num; //return the number of entries added
}


uPetscInt NchooseC(int n , int c ) {
	if ( c == n || c == 0 ) return 1;
	if ( c == 1 ) return n;
	if ( c > n ) return 0;
	return NchooseC(n,c-1) * (n-(c-1))/c;
}


/*! \brief  Given an index in the reduced basis set and information the conservation sector this function will return the corresponding basis element from the full basis set. */
uPetscInt get_full_basis_index(struct parameters *parameters, PetscInt sites, PetscInt filling, uPetscInt reduced_index) {
	uPetscInt full_index = 0;	
	
	PetscInt i = 1;
	uPetscInt l = reduced_index;
	while (i <= sites){ 
	  if ( l >= NchooseC(sites-i,filling) ) {
	      full_index += 1ul << (sites - i);
	      l -= NchooseC(sites-i,filling);
	      filling--;
	  }
	  i++;
	}
	return full_index;
}


// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order.
void get_monomial_bosons(struct parameters *parameters, uPetscInt index, PetscInt sites, PetscInt bosons, int *monomial) {
	int mom = 1;
	int idx = 0;
	int pos = 0;
	while ( pos < (sites + bosons) ) {
	  if ( (index & (1ul << pos) ) > 0ul ) {
	     monomial[idx++] = mom;
	  } else {
	     mom += 1;
	  }	  
	  pos++;
	}  
}

// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order and also return the normal.
void get_monomial_bosons(struct parameters *parameters, uPetscInt index, PetscInt sites, PetscInt bosons, int *monomial, PetscReal *normal) {
	int mom = 1;
	int idx = 0;
	int pos = 0;
	*normal = 1.0;
	double num = 0;
	while ( pos < (sites + bosons) ) {
	  if ( (index & (1ul << pos) ) > 0ul ) {
	     monomial[idx++] = mom;
	     num += 1.0;
	     *normal *= num;
	  } else {
	     mom += 1;
	     //if ( num > 1.0 ) *normal *= num; 
	     num = 0.0; 
	  }	  
	  pos++;
	}
	//if ( num > 1.0 ) *normal *= num; 
	*normal = sqrt(*normal);
}


// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order.
void get_monomial_fermions(struct parameters *parameters, uPetscInt index, PetscInt sites, PetscInt fermions, int *monomial) {
	int mom = 1;
	int idx = 0;
	int pos = 0;
	while ( pos < sites ) {
	  if ( (index & (1ul << pos) ) > 0ul ) {
	     monomial[idx++] = mom;
	     mom += 1; 
	  } else {
	     mom += 1;
	  }	  
	  pos++;
	}  
}

// calculated the position that each exciton above the ferrmionic spin one state. Does so in ascending order.
void get_monomial_spin_one(struct parameters *parameters, uPetscInt index, PetscInt sites, PetscInt fermions, int *monomial) {
	int mom = 1;
	int idx = 0;
	int pos = 0;
	while ( pos < sites ) 
	  {
	    if ( (index & (1ul << pos) ) == 0ul )  //if it is not in sz = -1
	      {
		if ( (index & (1ul << (pos + sites))) > 0ul ) //if it is in sz=1 then there are two so add one here and one after to monomial list
		  {
		    monomial[idx++] = mom;		    
		  }
	       monomial[idx++] = mom;
	       mom += 1; 
	     }
	   else 
	     {
	       mom += 1;
	     }	  
	   pos++;
	 }  
}

// calculated the position that each exciton above the ferrmionic spin one state. Does so in ascending order.
void get_monomial_spin_one(struct parameters *parameters, uPetscInt index, PetscInt sites, PetscInt fermions, int *monomial, PetscReal *normal) {
	int mom = 1;
	int idx = 0;
	int pos = 0;
	*normal = 1.0;
	double num = 0;
	while ( pos < sites ) 
	  {
	    if ( (index & (1ul << pos) ) == 0ul )  //if it is not in sz = -1
	      {
		if ( (index & (1ul << (pos + sites))) > 0ul ) //if it is in sz=1 then there are two so add one here and one after to monomial list
		  {
		    monomial[idx++] = mom;		
		    num += 1.0;
		    *normal *= num;
		  }
	       monomial[idx++] = mom;
	       mom += 1; 
	       num += 1.0;
	       *normal *= num;
	       num = 0.0;
	     }
	   else 
	     {
	       mom += 1;
	       num = 0.0;
	     }	  
	   pos++;
	 }  
	*normal = sqrt(*normal);
}

PetscInt get_number_permutations(struct parameters *parameters, uPetscInt index, PetscInt sites, PetscInt bosons) 
{
	int pos = 0;
	PetscInt permutations = 1;
	PetscInt num = 0;
	PetscInt choices = bosons;
	while ( pos < (sites + bosons) ) {
	  if ( (index & (1ul << pos) ) > 0ul ) {
	     permutations *= choices-- ;
	     num += 1;
	     permutations /= num;
	  } else {
	     num = 0; 
	  }	  
	  pos++;
	}
	return permutations;
}


/*! \brief  This function works out the index of the configuration and says whether it is valid or not via the exists flag. */
/*uPetscInt get_reduced_conservation_basis_index(struct parameters *parameters, uPetscInt full_index, PetscTruth *exists) { //TODO enforce conservation is conserved and same with parity. Check full index and return -1 if it is not part of the current sector. Done
	uPetscInt reduced_index = 0 ;
	PetscInt sector = parameters->real_conservation_sector;
	*exists = PETSC_TRUE;
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
	return reduced_index;
}*/

/*
 This function calculates the dimension of the boson basis set with given total momentum.
 
 sites : number of sites.
 filling : number of bosons.
 M : total momentum.
*/
PetscInt boson_basis_dimension(int sites, int filling, int M)
{  
  if ( sites > 0 && sites == M && filling == 1 )
    {
      return 1;
    }  
    
  if ( M < filling || M > sites*filling || sites == 0 || filling == 0 ) 
    {
      return 0; 
    } 
  else 
    {
    
      PetscInt dim = 0;  
      //use recursive procedure where can add boson and stay in position, add boson and move on or not add boson and move on.
      if ( sites <= M ) 
	{
	  dim += boson_basis_dimension(sites, filling-1, M-sites);	  
	}
      dim += boson_basis_dimension(sites-1, filling, M);
      return dim; 
    }
}

/*
 This function calculates the dimension of the boson basis set with given total momentum.
 
 sites : number of sites.
 filling : number of bosons.
 M : total momentum.
 states: array of states.
 maxm : the maximum m for a given state.
*/
PetscInt boson_basis_generate(int sites, int filling, int M, int pos, uPetscInt *states, uPetscInt current_state)
{  
  if ( sites > 0 && sites == M && filling == 1 )
    {      
      current_state += 1ul << (sites + filling - 2);      
      states[pos++] = current_state;      
      return pos;
    }  
    
  if ( M < filling || M > sites*filling || sites == 0 || filling == 0 ) 
    {
      return pos; 
    } 
  else 
    {          
      //use recursive procedure where can add boson and stay in position, add boson and move on or not add boson and move on.
      pos = boson_basis_generate(sites-1, filling, M, pos, states, current_state);
      if ( sites <= M ) 
	{
	  current_state += 1ul << (sites + filling - 2);	  
	  pos = boson_basis_generate(sites, filling-1, M-sites, pos, states, current_state);
	}
      
      return pos; 
    }
}

/*
 This function calculates all the different permutations of a given monomial.
 
 monomial : the original monomial to permute.
 filling : number of elements in the monomial.
 number :  the number of permutations.
 returns: 2 dim array of permutted monomials.
*/
int ** get_monomial_permutations(int filling, PetscInt *number)
{
   int **monomial_permutations;
   int *tmp_permutation;
   tmp_permutation = new int[filling];
   *number = 1;
   for ( int i = 1 ; i <= filling ; i++) 
    {
      *number *= i;
      tmp_permutation[i-1] = i-1;
    }
   
   monomial_permutations = new int*[*number];
   for ( PetscInt i = 0 ; i < *number;  i++ ) 
    {    
      monomial_permutations[i] = new int[filling];
    }
   
   int temp;   
   for ( PetscInt permutation = 0 ; permutation < *number ; permutation++ ) 
    {
      // Find largest index j with a[j] < a[j+1]
      int j = filling - 2;
      while (tmp_permutation[j] > tmp_permutation[j+1]) {
	j--;
      }

      // Find index k such that a[k] is smallest integer
      // greater than a[j] to the right of a[j]
      int k = filling - 1;
      while (tmp_permutation[j] > tmp_permutation[k]) {
	k--;
      }

      // Interchange a[j] and a[k]
      temp = tmp_permutation[k];
      tmp_permutation[k] = tmp_permutation[j];
      tmp_permutation[j] = temp;

      // Put tail end of permutation after jth position in increasing order

      int r = filling - 1;
      int s = j + 1;

      while (r > s) {
	temp = tmp_permutation[s];
	tmp_permutation[s] = tmp_permutation[r];
	tmp_permutation[r] = temp;
	r--;
	s++;
      }
      for ( int i = 0 ; i < filling ; i++) 
	{
	  monomial_permutations[permutation][i] = tmp_permutation[i];
	}
    }	
		   
   delete[] tmp_permutation;
   return monomial_permutations;
}


/*
* Convert from monomial representation to the unsigned integer representation.
*
* monomial: array of monomial coefficients.
* sites: the number of sites.
* filling: number of particles.
* return: the integer representation of the monomial.
*/
uPetscInt convert_from_monomial(int *monomial, int sites, int filling)
{
  uPetscInt basis_element = 0ul;
  
  for ( int i = 0 ; i < filling ; i++ ) 
    {
      basis_element += 1ul << (monomial[i] - 1 + i);
    }
  return basis_element;
}

int trim(const string str, string *trimmed) {
	size_t start = str.find_first_not_of(" \t\n\r");
	if(start == string::npos) return 1;
	trimmed->assign(str,start,str.find_last_not_of(" \t\n\r") - start + 1);
	return 0; 
}


uPetscInt* read_spin_basis(string basis_file, PetscInt* spin_basis_dim)
{
  uPetscInt *spin_basis;
  size_t pos;
  string line, tmp, tmp2;
  ifstream myfile (basis_file.c_str());
  if (myfile.is_open())
  {
    getline (myfile,line);
    pos = line.find_first_of("=");
    if ( pos == string::npos ) return NULL;
    tmp = line.substr(pos+1);
    trim(tmp, &tmp2);
    *spin_basis_dim = atoi(tmp2.c_str());
    
    PetscMalloc(*spin_basis_dim * sizeof(uPetscInt), &spin_basis);
    for ( PetscInt idx = 0 ; idx < *spin_basis_dim ; idx++ )
    {
      getline (myfile,line);
      spin_basis[idx] = (uPetscInt)atoi(line.c_str());
    }    
    myfile.close();
  }
  return spin_basis;
}

// This function calculates the normal for a given configuration. This is the inverse square root of the number of configurations that are equal by translation. 
// 
// config: the configuration.
// num_sites: the number of sites.
// mask: bit mask over relevant sites. 
// return: the norm. 
PetscScalar configuration_norm(uPetscInt config, int num_sites, uPetscInt mask, int step=1)
{
    PetscScalar normal = 1.0;    
    
    for ( int i = 1 ; i < num_sites ; i+=step ) 
    {
      if ( config == ( ((config << i) | (config >> (num_sites - i))) & mask ) ) 
	{
	  normal += 1.0;
	}
    }
    return 1.0/sqrt(normal);             
}

// This function calculates the normal for a given spin one configuration. This is the inverse square root of the number of configurations that are equal by translation.
// It is necessary to translate both the parts with bits from 0 to (num_sites-1) and from num_sites to (2*num_sites - 1)
// 
// config: the configuration.
// num_sites: the number of sites.
// mask: bit mask over relevant sites. 
// return: the norm. 
PetscScalar configuration_norm_spin_one(uPetscInt config, int num_sites, uPetscInt mask)
{
    PetscScalar normal = 1.0;    
    uPetscInt mask2 = mask << num_sites;
    for ( int i = 1 ; i < num_sites ; i++ ) 
    {
      if ( config == ( 
			( 
			( (((config & mask) << i) | ((config & mask) >> (num_sites - i))) & mask ) 
			| 
			( (((config & mask2) << i) | ((config & mask2) >> (num_sites - i))) & mask2 )
			)
		     )
	  ) 
	{
	  normal += 1.0;
	}
    }
    return 1.0/sqrt(normal);             
}


inline int permute(int *monomial, int len)
{
    int k,l,tmp;
    
    k = len - 2;
    while (k > 0 && monomial[k] >= monomial[k+1] ) k--;
    
    if ( k == 0 && (monomial[0] >= monomial[1]) ) return 0;
    
    l = len - 1;
    while ( monomial[l] <= monomial[k] ) l--;
    
    tmp = monomial[k];
    monomial[k] = monomial[l];
    monomial[l] = tmp;
    
    l = len - 1;
    k = k + 1;
    
    while ( k < l ) 
    {
      tmp = monomial[k];
      monomial[k] = monomial[l];
      monomial[l] = tmp;
      l--;
      k++;
    }
    return 1;       
}



/*
 This function calculated the fourier transform of vector A and returns it in vector B.
 
 A : local input vector.
 B : pointer to vector that will be global output vector.
 M : total momentum.
*/
int fourier_transform_vector_debug(struct parameters *parameters, Vec* A, Vec *B, int M, uPetscInt *spin_basis, PetscInt spin_basis_dim, PetscInt spin_momentum, uPetscInt *basis, PetscInt basis_dim, bool normalise, PetscInt stop_at) {
	PetscInt sites = parameters->AB_sites;
	PetscInt filling = parameters->filling;
	PetscInt orig_vec_length = spin_basis_dim;
	PetscInt vec_length = basis_dim;
	PetscErrorCode ierr;
	int *fmonomial, *omonomial;	
	ierr = PetscMalloc(filling*sizeof(int),&fmonomial);CHKERRQ(ierr);
	ierr = PetscMalloc(filling*sizeof(int),&omonomial);CHKERRQ(ierr);
	//create global array for fourier transform
	ierr = VecCreateMPI(MPI_COMM_WORLD, get_rank_chunk_size(vec_length, parameters->rank,parameters->size),vec_length,B);CHKERRQ(ierr);

	PetscScalar *original_array, *fourier_array;
	VecGetArray(*A, &original_array);
	VecGetArray(*B, &fourier_array);
	
	//precalculate the phases to save time
	PetscScalar **phase_array;		
	ierr = PetscMalloc((sites+1)*sizeof(PetscScalar*),&phase_array);CHKERRQ(ierr);
	for( int i = 1  ; i <= sites ; i++ ) {
	    ierr = PetscMalloc((sites+1)*sizeof(PetscScalar),&phase_array[i]);CHKERRQ(ierr);	
	}
	
	phase_array[1][0] = 1.0;
	for( int i = 1  ; i <= sites ; i++ ) {	    	    
	    for ( int j = 1; j <= sites ; j++ ) {
		PetscReal angle = ( PETSC_PI * 2.0 * ((PetscReal)(i*j))) / ((PetscReal)sites);
		phase_array[i][j] = cos(angle) + PETSC_i * sin(angle);
		//PetscPrintf(PETSC_COMM_WORLD,"%d, %d: %lf + i%lf\n",i,j,PetscRealPart(phase_array[i][j]), PetscImaginaryPart(phase_array[i][j]));
	    }
	}
		
	PetscScalar tmp_val, tmp_val2;
	PetscScalar norm = sqrt(1.0/pow((PetscReal)sites,(PetscReal)filling));
	PetscScalar momentum_phase = 1.0;
	if ( spin_momentum != -1 ) 
	{
	    for ( int i = 1 ; i < parameters->AB_sites ; i++ )
	    {
		momentum_phase += phase_array[1][(i*(M+spin_momentum))%parameters->AB_sites];
	    }	  
	}
	norm *= sqrt(momentum_phase);
	  		
	int **monomial_permutations;
	PetscInt number;
	if ( filling > 1 ) 
	  {
	    monomial_permutations = get_monomial_permutations(filling, &number);
	  }
	else if ( filling == 1)
	  {
	    number = 1;
	    monomial_permutations = new int*[1];
	    monomial_permutations[0] = new int[1];
	    monomial_permutations[0][0] = 0;
	  }	    
	
	PetscInt fstart = get_rank_chunk_start(vec_length, parameters->rank, parameters->size);
	int phase_idx;  	
	PetscInt configurations_coupled;
	
	PetscReal start_time, end_time; 
	start_time = MPI_Wtime();
	
	PetscInt fend = get_rank_chunk_size(vec_length, parameters->rank, parameters->size);
	if ( stop_at > 0 && stop_at < fend ) fend = stop_at;
	for ( PetscInt fidx = 0 ; fidx < fend ; fidx++ ){	    	      
	    configurations_coupled = 0 ;
	    tmp_val = 0.0; 
	    PetscReal state_normal = 1.0;
	    if ( normalise ) 
	      {
		get_monomial_bosons(parameters, basis[fstart+fidx], sites, filling, fmonomial, &state_normal);
	      }
	    else 
	      {
		get_monomial_bosons(parameters, basis[fstart+fidx], sites, filling, fmonomial);	    	    
	      }	      	    
	      for ( int i = 0 ; i < filling ; i++ ) 
	      {
		PetscPrintf(PETSC_COMM_WORLD,"%d, ", fmonomial[i]);
	      }
	      PetscPrintf(PETSC_COMM_WORLD,"\n");
	    for ( PetscInt oidx = 0 ; oidx < orig_vec_length ; oidx++ )
	      {
		tmp_val2 = 0.0;
		//get_monomial_fermions(parameters, get_full_basis_index(parameters, sites, filling,oidx), sites, filling, omonomial);
		get_monomial_fermions(parameters, spin_basis[oidx], sites, filling, omonomial);		
		PetscPrintf(PETSC_COMM_WORLD,"New fermion config: ");
		for ( int i = 0 ; i < filling ; i++ ) 
		{
		  PetscPrintf(PETSC_COMM_WORLD,"%d, ", omonomial[i]);
		}
		PetscPrintf(PETSC_COMM_WORLD,"\n");
		for ( PetscInt pidx = 0 ; pidx < number ; pidx++ )
		  {		    
		    phase_idx = 0;
		    for ( int i = 0 ; i < filling ; i++ ) 
		      {
			phase_idx += omonomial[i] * fmonomial[monomial_permutations[pidx][i]];
		      }
		    tmp_val2 += phase_array[1][phase_idx % sites];
		    PetscPrintf(PETSC_COMM_WORLD,"P %ld: %ld (mod %ld)\n", pidx, phase_idx, phase_idx%parameters->AB_sites);
		  }
		  PetscPrintf(PETSC_COMM_WORLD,"%lf + i%lf\n", PetscRealPart(tmp_val2), PetscImaginaryPart(tmp_val2));
		//tmp_val += tmp_val2 * original_array[oidx] * configuration_norm(fmonomial, filling, parameters->AB_sites) / state_normal;
		//tmp_val += tmp_val2 * original_array[oidx] / (state_normal * configuration_norm(fmonomial, filling, parameters->AB_sites));
		
		tmp_val += tmp_val2 * original_array[oidx] / state_normal;
		//if (PetscAbsScalar(tmp_val2 * original_array[oidx] / (state_normal * configuration_norm(omonomial, filling, parameters->AB_sites))) > 1e-10 ) configurations_coupled ++;
	      }			     	      
	    fourier_array[fidx] = tmp_val * norm;	    
	    /*for ( int i = 0 ; i < filling ; i++ ) 
	      {
		PetscPrintf(PETSC_COMM_WORLD," %d", fmonomial[i]);
	      }
	    PetscPrintf(PETSC_COMM_WORLD,": %d\n", configurations_coupled);*/
	}		
	
	end_time = MPI_Wtime();
	PetscPrintf(PETSC_COMM_WORLD,"Fourier transform took: %lf seconds.\n", end_time - start_time);
	
	for ( PetscInt pidx = 0 ; pidx < number ; pidx++ )  
	    delete[] monomial_permutations[pidx];
	delete[]monomial_permutations;
	VecRestoreArray(*A, &original_array);
	VecRestoreArray(*B, &fourier_array);	
	
	PetscFree(fmonomial);
	PetscFree(omonomial);
	
	for( int i = 1  ; i <= sites; i++ ) {
	    ierr = PetscFree(phase_array[i]);CHKERRQ(ierr);	
	}
	ierr = PetscFree(phase_array);CHKERRQ(ierr);
	
	return 0;
}

/*
 This function calculated the fourier transform of vector A and returns it in vector B.
 
 A : local input vector.
 B : pointer to vector that will be global output vector.
 M : total momentum.
*/
int fourier_transform_vector_opt(struct parameters *parameters, Vec* A, Vec *B, int M, uPetscInt *spin_basis, PetscInt spin_basis_dim, PetscInt spin_momentum, uPetscInt *basis, PetscInt basis_dim, bool normalise, PetscInt stop_at, PetscInt *end_indices) {
	PetscInt sites = parameters->AB_sites;
	PetscInt filling = parameters->filling;
	PetscInt orig_vec_length = spin_basis_dim;
	PetscInt vec_length = basis_dim;
	PetscErrorCode ierr;
	PetscInt starting_idx = 0;	
	int *fmonomial, *omonomial;
	ierr = PetscMalloc(filling*sizeof(int),&fmonomial);CHKERRQ(ierr);
	ierr = PetscMalloc(filling*sizeof(int),&omonomial);CHKERRQ(ierr);
	//create global array for fourier transform
	
	// if we are using multiple processors and not on the first one then changing the starting index
	if ( parameters->rank != 0 )
	{
		starting_idx = end_indices[parameters->rank-1];
	}
	
	//create the final vector
	ierr = VecCreateMPI(MPI_COMM_WORLD, end_indices[parameters->rank] - starting_idx,vec_length,B);CHKERRQ(ierr);	

	//get direct access to output and input vectors.
	PetscScalar *original_array, *fourier_array;
	VecGetArray(*A, &original_array);
	VecGetArray(*B, &fourier_array);
	
	//precalculate the phases
	PetscScalar **phase_array;		
	ierr = PetscMalloc((sites+1)*sizeof(PetscScalar*),&phase_array);CHKERRQ(ierr);
	for( int i = 1  ; i <= sites ; i++ ) {
	    ierr = PetscMalloc((sites+1)*sizeof(PetscScalar),&phase_array[i]);CHKERRQ(ierr);	
	}		
	phase_array[1][0] = 1.0;
	for( int i = 1  ; i <= sites ; i++ ) {	    	    
	    for ( int j = 1; j <= sites ; j++ ) {
		PetscReal angle = ( PETSC_PI * 2.0 * ((PetscReal)(i*j))) / ((PetscReal)sites);
		phase_array[i][j] = cos(angle) + PETSC_i * sin(angle);
		//PetscPrintf(PETSC_COMM_WORLD,"%d, %d: %lf + i%lf\n",i,j,PetscRealPart(phase_array[i][j]), PetscImaginaryPart(phase_array[i][j]));
	    }
	}
		
	PetscScalar tmp_val, tmp_val2;
	PetscScalar norm = sqrt(1.0/pow((PetscReal)sites,(PetscReal)filling)); //normalisation
	
	PetscScalar momentum_phase = 1.0; // defautls to one if not using translation invariance
	//PetscReal *config_norms;
	if ( spin_momentum != -1 )  //if using translation invariance
	{
	    for ( int i = 1 ; i < parameters->AB_sites ; i++ )
	    {
		momentum_phase += phase_array[1][(i*(M+spin_momentum))%parameters->AB_sites];
	    }		    
	//     ierr = PetscMalloc(spin_basis_dim*sizeof(PetscReal), &config_norms); CHKERRQ(ierr);	    
	}
	norm *= sqrt(momentum_phase);	
				  		
	/*int **monomial_permutations;
	PetscInt number;
	if ( filling > 1 ) 
	  {
	    monomial_permutations = get_monomial_permutations(filling, &number);
	  }
	else if ( filling == 1)
	  {
	    number = 1;
	    monomial_permutations = new int*[1];
	    monomial_permutations[0] = new int[1];
	    monomial_permutations[0][0] = 0;
	  }	    */
	
	
	  
	//PetscInt fstart = get_rank_chunk_start(vec_length, parameters->rank, parameters->size);
	PetscInt fstart = starting_idx;
	int phase_idx;  	
	PetscInt configurations_coupled;
	
	PetscReal start_time, end_time; 
	start_time = MPI_Wtime();
	int k,l,tmp;
	
	//PetscInt fend = get_rank_chunk_size(vec_length, parameters->rank, parameters->size);
	PetscInt fend = end_indices[parameters->rank];
	if ( stop_at > 0 && stop_at < (fend - fstart) ) fend = fstart + stop_at;
	for ( PetscInt fidx = fstart ; fidx < fend ; fidx++ ){	    	      
	    configurations_coupled = 0 ;
	    tmp_val = 0.0; 
	    PetscReal state_normal = 1.0;
	    if ( normalise ) 
	      {
		get_monomial_bosons(parameters, basis[fidx], sites, filling, fmonomial, &state_normal);
	      }
	    else 
	      {
		get_monomial_bosons(parameters, basis[fidx], sites, filling, fmonomial);	    	    
	      }	      	    	      
	    for ( PetscInt oidx = 0 ; oidx < orig_vec_length ; oidx++ )
	      {
		sort(fmonomial, fmonomial+filling); //rest fmonomial to original permutation
		tmp_val2 = 0.0;						
		get_monomial_fermions(parameters, spin_basis[oidx], sites, filling, omonomial);				
		phase_idx = 0;
		for ( int i = 0 ; i < parameters->filling ; i++ ) 
		  {
		    phase_idx += omonomial[i] * fmonomial[i];
		  }
		int permute_flg = 1;
		while ( permute_flg )
		  {
		    tmp_val2 += phase_array[1][phase_idx % sites];
		    k = parameters->filling - 2;
		    while (k > 0 && fmonomial[k] >= fmonomial[k+1] ) k--;
		    
		    if ( k == 0 && (fmonomial[0] >= fmonomial[1]) ) 
		      {
			permute_flg = 0;
		      }
		    else
		      {
			l = parameters->filling - 1;
			while ( fmonomial[l] <= fmonomial[k] ) l--;

			phase_idx += (fmonomial[l] - fmonomial[k])*omonomial[k] + (fmonomial[k] - fmonomial[l])*omonomial[l];
			tmp = fmonomial[k];
			fmonomial[k] = fmonomial[l];
			fmonomial[l] = tmp;
			
			l = parameters->filling - 1;
			k = k + 1;
			
			while ( k < l ) 
			{
			  phase_idx += (fmonomial[l] - fmonomial[k])*omonomial[k] + (fmonomial[k] - fmonomial[l])*omonomial[l];
			  tmp = fmonomial[k];
			  fmonomial[k] = fmonomial[l];
			  fmonomial[l] = tmp;
			  l--;
			  k++;
			}     
			permute_flg = 1;
		      }
		  }
		//commenting this out to put the permute function directly here. 
		/*do 
		{		  		  
		    phase_idx = 0;
		    for ( int i = 0 ; i < filling ; i++ ) 
		      {
			phase_idx += omonomial[i] * fmonomial[i];			
		      }		    
		    tmp_val2 += phase_array[1][phase_idx % sites];
		} while ( permute(fmonomial, filling) );*/
				
		/*for ( PetscInt pidx = 0 ; pidx < number ; pidx++ )
		  {		    
		    phase_idx = 0;
		    for ( int i = 0 ; i < filling ; i++ ) 
		      {
			phase_idx += omonomial[i] * fmonomial[monomial_permutations[pidx][i]];
		      }
		    tmp_val2 += phase_array[1][phase_idx % sites];
		  }*/
		
		/*if ( spin_momentum != -1 ) 
		  tmp_val += tmp_val2 * original_array[oidx] * state_normal;
		else tmp_val += tmp_val2 * original_array[oidx] * state_normal;*/
		tmp_val += tmp_val2 * original_array[oidx] * state_normal;
	      }			     	      
	    fourier_array[fidx - fstart] = tmp_val * norm;	    	    
	}	
	
	end_time = MPI_Wtime();
	PetscPrintf(PETSC_COMM_SELF,"Rank: %d, Fourier transform took: %lf seconds.\n", parameters->rank, end_time - start_time);
	
	/*for ( PetscInt pidx = 0 ; pidx < number ; pidx++ )  
	    delete[] monomial_permutations[pidx];
	delete[]monomial_permutations;*/
	VecRestoreArray(*A, &original_array);
	VecRestoreArray(*B, &fourier_array);	
	
	PetscFree(fmonomial);
	PetscFree(omonomial);
	
	for( int i = 1  ; i <= sites; i++ ) {
	    ierr = PetscFree(phase_array[i]);CHKERRQ(ierr);	
	}
	ierr = PetscFree(phase_array);CHKERRQ(ierr);
	
	return 0;
}




/*
 This function calculates the fourier transform of spin one vector A and returns it in vector B.
 
 A : local input vector.
 B : pointer to vector that will be global output vector.
 M : total momentum.
*/
int fourier_transform_vector_spin_one_opt(struct parameters *parameters, Vec* A, Vec *B, int M, uPetscInt *spin_basis, PetscInt spin_basis_dim, PetscInt spin_momentum, uPetscInt *basis, PetscInt basis_dim, bool normalise, PetscInt stop_at, PetscInt *end_indices) {
	PetscInt sites = parameters->AB_sites;
	PetscInt filling = parameters->filling;
	PetscInt orig_vec_length = spin_basis_dim;
	PetscInt vec_length = basis_dim;
	PetscErrorCode ierr;
	PetscInt starting_idx = 0;	
	int *fmonomial, *omonomial;
	ierr = PetscMalloc(filling*sizeof(int),&fmonomial);CHKERRQ(ierr);
	ierr = PetscMalloc(filling*sizeof(int),&omonomial);CHKERRQ(ierr);
	//create global array for fourier transform
	
	if ( parameters->rank != 0 )
	{
		starting_idx = end_indices[parameters->rank-1];
	}
	
	ierr = VecCreateMPI(MPI_COMM_WORLD, end_indices[parameters->rank] - starting_idx,vec_length,B);CHKERRQ(ierr);	

	PetscScalar *original_array, *fourier_array;
	VecGetArray(*A, &original_array);
	VecGetArray(*B, &fourier_array);
	
	//precalculate the phases to save time
	PetscScalar **phase_array;		
	ierr = PetscMalloc((sites+1)*sizeof(PetscScalar*),&phase_array);CHKERRQ(ierr);
	for( int i = 1  ; i <= sites ; i++ ) {
	    ierr = PetscMalloc((sites+1)*sizeof(PetscScalar),&phase_array[i]);CHKERRQ(ierr);	
	}
	
	phase_array[1][0] = 1.0;
	for( int i = 1  ; i <= sites ; i++ ) {	    	    
	    for ( int j = 1; j <= sites ; j++ ) {
		PetscReal angle = ( PETSC_PI * 2.0 * ((PetscReal)(i*j))) / ((PetscReal)sites);
		phase_array[i][j] = cos(angle) + PETSC_i * sin(angle);
		//PetscPrintf(PETSC_COMM_WORLD,"%d, %d: %lf + i%lf\n",i,j,PetscRealPart(phase_array[i][j]), PetscImaginaryPart(phase_array[i][j]));
	    }
	}
		
	PetscScalar tmp_val, tmp_val2;
	PetscScalar norm = sqrt(1.0/pow((PetscReal)sites,(PetscReal)filling));
	PetscScalar momentum_phase = 1.0;
	if ( spin_momentum != -1 ) 
	{
	    for ( int i = 1 ; i < parameters->AB_sites ; i++ )
	    {
		momentum_phase += phase_array[1][(i*(M+spin_momentum))%parameters->AB_sites];
	    }	  
	}
	norm *= sqrt(momentum_phase);
	  		
	/*int **monomial_permutations;
	PetscInt number;
	if ( filling > 1 ) 
	  {
	    monomial_permutations = get_monomial_permutations(filling, &number);
	  }
	else if ( filling == 1)
	  {
	    number = 1;
	    monomial_permutations = new int*[1];
	    monomial_permutations[0] = new int[1];
	    monomial_permutations[0][0] = 0;
	  }	    */
	
	//PetscInt fstart = get_rank_chunk_start(vec_length, parameters->rank, parameters->size);
	PetscInt fstart = starting_idx;
	int phase_idx;  	
	PetscInt configurations_coupled;
	
	PetscReal start_time, end_time; 
	start_time = MPI_Wtime();
	
	//PetscInt fend = get_rank_chunk_size(vec_length, parameters->rank, parameters->size);
	PetscInt fend = end_indices[parameters->rank];
	if ( stop_at > 0 && stop_at < (fend - fstart) ) fend = fstart + stop_at;
	for ( PetscInt fidx = fstart ; fidx < fend ; fidx++ ){	    	      
	    configurations_coupled = 0 ;
	    tmp_val = 0.0; 
	    PetscReal state_normal = 1.0;
	    if ( normalise ) 
	      {
		get_monomial_bosons(parameters, basis[fidx], sites, filling, fmonomial, &state_normal);
	      }
	    else 
	      {
		get_monomial_bosons(parameters, basis[fidx], sites, filling, fmonomial);	    	    
	      }	      	    	      
	    for ( PetscInt oidx = 0 ; oidx < orig_vec_length ; oidx++ )
	      {
		sort(fmonomial, fmonomial+filling); //restore fmonomial to original permutation
		tmp_val2 = 0.0;	
		PetscReal spin_one_norm;
		get_monomial_spin_one(parameters, spin_basis[oidx], sites, filling, omonomial, &spin_one_norm);		
		do 
		{		  		  
		    phase_idx = 0;
		    for ( int i = 0 ; i < filling ; i++ ) 
		      {
			phase_idx += omonomial[i] * fmonomial[i];			
		      }		    
		    tmp_val2 += phase_array[1][phase_idx % sites];
		} while ( permute(fmonomial, filling) );
				
		/*for ( PetscInt pidx = 0 ; pidx < number ; pidx++ )
		  {		    
		    phase_idx = 0;
		    for ( int i = 0 ; i < filling ; i++ ) 
		      {
			phase_idx += omonomial[i] * fmonomial[monomial_permutations[pidx][i]];
		      }
		    tmp_val2 += phase_array[1][phase_idx % sites];
		  }*/
		
		tmp_val += tmp_val2 * original_array[oidx] * state_normal / spin_one_norm;
	      }			     	      
	    fourier_array[fidx - fstart] = tmp_val * norm ;
	}	
	
	end_time = MPI_Wtime();
	PetscPrintf(PETSC_COMM_SELF,"Rank: %d, Fourier transform took: %lf seconds.\n", parameters->rank, end_time - start_time);
	
	/*for ( PetscInt pidx = 0 ; pidx < number ; pidx++ )  
	    delete[] monomial_permutations[pidx];
	delete[]monomial_permutations;*/
	VecRestoreArray(*A, &original_array);
	VecRestoreArray(*B, &fourier_array);	
	
	PetscFree(fmonomial);
	PetscFree(omonomial);
	
	for( int i = 1  ; i <= sites; i++ ) {
	    ierr = PetscFree(phase_array[i]);CHKERRQ(ierr);	
	}
	ierr = PetscFree(phase_array);CHKERRQ(ierr);
	
	return 0;
}


/*
 This function calculated the fourier transform of vector A and returns it in vector B.
 
 A : fourier transform matrix.
 M : total momentum
*/
int fourier_transform_matrix_opt(struct parameters *parameters, Mat *A, int M, uPetscInt *spin_basis, PetscInt spin_basis_dim, PetscInt spin_momentum, uPetscInt *basis, PetscInt basis_dim, bool normalise, PetscInt stop_at, PetscInt *end_indices) {
	PetscInt sites = parameters->AB_sites;
	PetscInt filling = parameters->filling;
	PetscInt orig_vec_length = spin_basis_dim;
	//PetscInt vec_length = basis_dim;
	PetscErrorCode ierr;
	PetscInt starting_idx = 0;
	//uPetscInt mask = (1ul << parameters->AB_sites) - 1ul; 
	int *fmonomial, *omonomial;
	ierr = PetscMalloc(filling*sizeof(int),&fmonomial);CHKERRQ(ierr);
	ierr = PetscMalloc(filling*sizeof(int),&omonomial);CHKERRQ(ierr);
	//create global array for fourier transform
	
	// if we are using multiple processors and not on the first one then changing the starting index
	if ( parameters->rank != 0 )
	{
		starting_idx = end_indices[parameters->rank-1];
	}
	
	//create the final vector
	//ierr = MatCreateMPIDense(MPI_COMM_WORLD, end_indices[parameters->rank] - starting_idx, PETSC_DECIDE, basis_dim, spin_basis_dim, PETSC_NULL, A);CHKERRQ(ierr);	
	ierr = MatCreateDense(MPI_COMM_WORLD, end_indices[parameters->rank] - starting_idx, PETSC_DECIDE, basis_dim, spin_basis_dim, PETSC_NULL, A);CHKERRQ(ierr);	
	ierr = MatZeroEntries(*A);CHKERRQ(ierr);
		

	//precalculate the phases
	PetscScalar **phase_array;		
	ierr = PetscMalloc((sites+1)*sizeof(PetscScalar*),&phase_array);CHKERRQ(ierr);
	for( int i = 1  ; i <= sites ; i++ ) {
	    ierr = PetscMalloc((sites+1)*sizeof(PetscScalar),&phase_array[i]);CHKERRQ(ierr);	
	}		
	phase_array[1][0] = 1.0;
	for( int i = 1  ; i <= sites ; i++ ) {	    	    
	    for ( int j = 1; j <= sites ; j++ ) {
		PetscReal angle = ( PETSC_PI * 2.0 * ((PetscReal)(i*j))) / ((PetscReal)sites);
		phase_array[i][j] = cos(angle) + PETSC_i * sin(angle);
		//PetscPrintf(PETSC_COMM_WORLD,"%d, %d: %lf + i%lf\n",i,j,PetscRealPart(phase_array[i][j]), PetscImaginaryPart(phase_array[i][j]));
	    }
	}
		
	PetscScalar tmp_val, tmp_val2;
	PetscScalar norm = sqrt(1.0/pow((PetscReal)sites,(PetscReal)filling)); //normalisation
	
	PetscScalar momentum_phase = 1.0; // defautls to one if not using translation invariance

	if ( spin_momentum != -1 )  //if using translation invariance
	{
	    for ( int i = 1 ; i < parameters->AB_sites ; i++ )
	    {
		momentum_phase += phase_array[1][(i*(M+spin_momentum))%parameters->AB_sites];
	    }		    
	}
	norm *= sqrt(momentum_phase);	
				  		

	//PetscInt fstart = get_rank_chunk_start(vec_length, parameters->rank, parameters->size);
	PetscInt fstart = starting_idx;
	int phase_idx;  	
	PetscInt configurations_coupled;
	
	PetscReal start_time, end_time; 
	start_time = MPI_Wtime();
	int k,l,tmp;
	
	//PetscScalar *mat_vals;	
	//MatGetArray(*A,&mat_vals);
	
	//PetscInt fend = get_rank_chunk_size(vec_length, parameters->rank, parameters->size);
	PetscInt fend = end_indices[parameters->rank];
	if ( stop_at > 0 && stop_at < (fend - fstart) ) fend = fstart + stop_at;
	for ( PetscInt fidx = fstart ; fidx < fend ; fidx++ ){	    	      
	    //PetscInt mat_array_offset = (fidx - fstart) * orig_vec_length ;
	    //MatGetRow(*A, fidx, NULL, NULL, (const PetscScalar**)&mat_vals);	  
	    tmp_val = 0.0; 
	    PetscReal state_normal = 1.0;
	    if ( normalise ) 
	      {
		get_monomial_bosons(parameters, basis[fidx], sites, filling, fmonomial, &state_normal);
	      }
	    else 
	      {
		get_monomial_bosons(parameters, basis[fidx], sites, filling, fmonomial);	    	    
	      }	      	    	      
	    for ( PetscInt oidx = 0 ; oidx < orig_vec_length ; oidx++ )
	      {
		sort(fmonomial, fmonomial+filling); //rest fmonomial to original permutation
		tmp_val2 = 0.0;						
		get_monomial_fermions(parameters, spin_basis[oidx], sites, filling, omonomial);				
		phase_idx = 0;
		for ( int i = 0 ; i < parameters->filling ; i++ ) 
		  {
		    phase_idx += omonomial[i] * fmonomial[i];
		  }
		int permute_flg = 1;
		while ( permute_flg )
		  {
		    tmp_val2 += phase_array[1][phase_idx % sites];
		    k = parameters->filling - 2;
		    while (k > 0 && fmonomial[k] >= fmonomial[k+1] ) k--;
		    
		    if ( k == 0 && (fmonomial[0] >= fmonomial[1]) ) 
		      {
			permute_flg = 0;
		      }
		    else
		      {
			l = parameters->filling - 1;
			while ( fmonomial[l] <= fmonomial[k] ) l--;

			phase_idx += (fmonomial[l] - fmonomial[k])*omonomial[k] + (fmonomial[k] - fmonomial[l])*omonomial[l];
			tmp = fmonomial[k];
			fmonomial[k] = fmonomial[l];
			fmonomial[l] = tmp;
			
			l = parameters->filling - 1;
			k = k + 1;
			
			while ( k < l ) 
			{
			  phase_idx += (fmonomial[l] - fmonomial[k])*omonomial[k] + (fmonomial[k] - fmonomial[l])*omonomial[l];
			  tmp = fmonomial[k];
			  fmonomial[k] = fmonomial[l];
			  fmonomial[l] = tmp;
			  l--;
			  k++;
			}     
			permute_flg = 1;
		      }
		  }				
		MatSetValue(*A, fidx, oidx, tmp_val2 * state_normal * norm, ADD_VALUES); 		  
		//mat_vals[oidx] += tmp_val2 * state_normal * norm;
	      }		
	    //MatRestoreRow(*A, fidx, NULL, NULL, (const PetscScalar**)&mat_vals);
	}	
		
	//MatRestoreArray(*A, &mat_vals);
	
	MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);
		
	end_time = MPI_Wtime();
	PetscPrintf(PETSC_COMM_SELF,"Rank: %d, Matrix construction took: %lf seconds.\n", parameters->rank, end_time - start_time);
	
	PetscFree(fmonomial);
	PetscFree(omonomial);
	
	for( int i = 1  ; i <= sites; i++ ) {
	    ierr = PetscFree(phase_array[i]);CHKERRQ(ierr);	
	}
	ierr = PetscFree(phase_array);CHKERRQ(ierr);
	
	return 0;
}

/*
 With n numbers (1..n) to choose from and choosing k populate choice array with idx^th selection of k numbers.
 
 n: number of numbers
 k: number to choose
 idx: index of selection
 choice: array where final selection goes
*/
void get_choosen_configuration(int n, int k, int idx, int *choice, int **binomialCoefs) 
{  
  int i=1, l=idx; int c =k;
  int count = 0;
  while ( i <= n && count < k )
    {
      if ( l >= binomialCoefs[n-i][c] )	
	{
	  choice[count++] = n-i;
	  l -= binomialCoefs[n-i][c];
	  c -= 1;
	}
      i++;
    }
}

inline void get_combinations_from_array(int *array, int n, int *choice, int *others, int c, int l, int **binomialCoefs)
{
    int i = 0;
    int count = 0;
    int other_count = 0;
    while ( i < n )
      {
        if ( l >= binomialCoefs[n-i-1][c] )
          {
            choice[count++] = array[n-i-1];
            l -= binomialCoefs[n-i-1][c];
            c -= 1;
          }
        else
          {
            others[other_count++] = array[n-i-1];
          }
        i++;
      }
}

inline int get_label_for_combination(int *array, int n, int c, int **binomialCoefs)
{
    int label = 0;
    
    for ( int i = 0 ; i < c ; i++ )
      {
        label += binomialCoefs[array[i]][ i+1];
      }
    return label;
}


/* Given array of objects we need to find which combination with in the whole each partition belongs to.

choice: the array of indices
n : total number of choices
c : number in this selection
ln : number in left partition
rn : number in right partition
lidx: pointer to left index
ridx: pointer to right index
idx: idx for picking the ln out of the c elements.
*/
void get_combination_indices(int *choice, int n, int c, int ln, int rn, int *lidx, int *ridx, int idx, int **binomialCoefs)
{
  //first find the idx^th combination of ln elements out of the c elements
  int *leftside; PetscMalloc(ln*sizeof(int), &leftside);
  int *rightside; PetscMalloc(rn*sizeof(int), &rightside);
  get_combinations_from_array(choice, c, leftside, rightside, ln, idx, binomialCoefs);
  
  //have combinations on left and right sides. Now find indices of these out of all particles. 
  *lidx = get_label_for_combination(leftside, n, ln, binomialCoefs);
  *ridx = get_label_for_combination(rightside, n, rn, binomialCoefs);
  PetscFree(leftside);
  PetscFree(rightside);
}

/* Get the numbeer of configurations resulting from a combination of two groups.

configurations: array of arrays of configurations from previous level.
numberConfigurationsPerCollection: array of configuration counts.
leftIndex: index of left group
rightIndex: index of right group (could be -1 to indicate none)
*/

int CombinedConfigurationsCount(uPetscInt **configurations, int *numberConfigurationsPerCollection, int leftIndex, int rightIndex)
{
    if ( rightIndex == -1 )
      {
	return numberConfigurationsPerCollection[leftIndex];
      }
    int ConfigurationCount=0;    
    for ( int i = 0 ; i < numberConfigurationsPerCollection[leftIndex]; i++ )
      {
	for ( int j = 0 ; j < numberConfigurationsPerCollection[rightIndex]; j++ )
	  {
	    if ( (configurations[leftIndex][i] & configurations[rightIndex][j]) == 0ul ) ConfigurationCount++;	    
	  }	
      }          
    return ConfigurationCount;
}

/* Populate the arrays of combined configurations and also the links to their constituent parts.

configurations: array of arrays of configurations from pervious level.
numberConfigurationsPerCollection: array of configuration counts.
leftIndex: index of left group
rightIndex: index of right group (could be -1 to indicate none)
combinedConfigurations: array of combimed configurations to populate
combinedConfigurationConstituentIndices: 2D array of constituents that make up each combimed configuration
*/

int CombineConfigurations(uPetscInt **configurations, int *numberConfigurationsPerCollection, int leftIndex, int rightIndex, uPetscInt *combinedConfigurations,int **combinedConfigurationConstituentIndices)
{
    if ( rightIndex == -1 )
      {
	for ( int i = 0 ; i < numberConfigurationsPerCollection[leftIndex]; i++ ) 
	  {
	    combinedConfigurations[i] = configurations[leftIndex][i];
	    combinedConfigurationConstituentIndices[i][0] = i;
	    combinedConfigurationConstituentIndices[i][1] = -1;
	  }
      }
    int ConfigurationCount=0;    
    for ( int i = 0 ; i < numberConfigurationsPerCollection[leftIndex]; i++ )
      {
	for ( int j = 0 ; j < numberConfigurationsPerCollection[rightIndex]; j++ )
	  {
	    if ( (configurations[leftIndex][i] & configurations[rightIndex][j]) == 0ul ) 
	      {
		combinedConfigurations[ConfigurationCount] = configurations[leftIndex][i] | configurations[rightIndex][j];
		combinedConfigurationConstituentIndices[ConfigurationCount][0] = i;
		combinedConfigurationConstituentIndices[ConfigurationCount][1] = j;
		ConfigurationCount++;
	      }
	  }	
      }          
    return ConfigurationCount;
}


/* Get the numbeer of configurations resulting from a decomposition into two groups.

configurations: the configurations to decompose
numberConfigurations: the number of configurations to decompose.
particles: number of particles in each configuration.
numberConfigurationsLeft: number on left.
numberConfigurationsRight: number on right.
configurationsConstituentIndices: array labelling which configurations each configuration decomposes to.
lowerLevelConfigurations: array of arrys of configurations on the next level. Arrays need to be allocated.
leftCollectionIndex: index on the lower level of the left branch of configurations.
*/

int DecomposeConfiguration(uPetscInt *configurations, int numberConfigurations, int particles, int sites, int *numberConfigurationsLeft, int *numberConfigurationsRight, int **configurationsConstituentIndices, uPetscInt **lowerLevelConfigurations, int leftCollectionIndex)
{
    PetscErrorCode ierr;
    if ( particles == 1 )
      {
	*numberConfigurationsLeft = numberConfigurations;	
	ierr=PetscMalloc(numberConfigurations*sizeof(uPetscInt), &lowerLevelConfigurations[leftCollectionIndex]);CHKERRQ(ierr);
	for ( int i = 0 ; i < numberConfigurations ; i++ ) 
	  {
	    lowerLevelConfigurations[leftCollectionIndex][i] = configurations[i];
	    configurationsConstituentIndices[i][0] = i;
	    configurationsConstituentIndices[i][1] = -1;
	  }
	return 0;
      }
      
    int LeftConfigurationCount=0;
    int RightConfigurationCount=0;    
    map<uPetscInt,int> LeftConfigurations, RightConfigurations;
    map<uPetscInt,int>::iterator lit, rit;
    int LeftParticles, RightParticles;
    uPetscInt LeftConf, RightConf;
    LeftParticles = (particles+1)/2;
    RightParticles = particles - LeftParticles;
    
    for ( int idx = 0 ; idx < numberConfigurations; idx++ ) 
      {
	int count = 0;
	RightConf = 0;
	int pos = 0;
	while ((count < RightParticles) && (pos < sites) )
	  {
	    if ( ((1ul << (sites-pos-1)) & configurations[idx]) > 0ul )
	      {
		RightConf += (1ul << (sites-pos-1));
		count++;
	      }
	    pos++;
	  }
	LeftConf = configurations[idx] ^ RightConf;
	rit = RightConfigurations.find(RightConf);
	lit = LeftConfigurations.find(LeftConf);
	if ( rit == RightConfigurations.end() ) 
	  {	    
	    RightConfigurations.insert(pair<uPetscInt, int>(RightConf, RightConfigurationCount));
	    RightConfigurationCount++;
	  }
	if ( lit == LeftConfigurations.end() ) 
	  {
	    LeftConfigurations.insert(pair<uPetscInt, int>(LeftConf, LeftConfigurationCount));
	    LeftConfigurationCount++;
	  }
      }    
    *numberConfigurationsRight = RightConfigurationCount;        
    *numberConfigurationsLeft = LeftConfigurationCount;
    ierr = PetscMalloc(*numberConfigurationsRight*sizeof(uPetscInt), &lowerLevelConfigurations[leftCollectionIndex+1]);CHKERRQ(ierr);
    ierr = PetscMalloc(*numberConfigurationsLeft*sizeof(uPetscInt), &lowerLevelConfigurations[leftCollectionIndex]);CHKERRQ(ierr);    
    for ( map<uPetscInt, int>::iterator it = RightConfigurations.begin(); it != RightConfigurations.end(); it++ )
      {
	lowerLevelConfigurations[leftCollectionIndex+1][(*it).second] = (*it).first;		
      }
    for ( map<uPetscInt, int>::iterator it = LeftConfigurations.begin(); it != LeftConfigurations.end(); it++ )
      {
	lowerLevelConfigurations[leftCollectionIndex][(*it).second] = (*it).first;	
      }    
      
     for ( int idx = 0 ; idx < numberConfigurations; idx++ ) 
      {  
	int count = 0;
	RightConf = 0;
	int pos = 0;
	while (count < RightParticles )
	  {
	    if ( ((1ul << (sites-pos-1)) & configurations[idx]) > 0ul )
	      {
		RightConf += 1ul << (sites-pos-1);
		count++;
	      }
	    pos++;
	  }
	LeftConf = configurations[idx] ^ RightConf;
	rit = RightConfigurations.find(RightConf);
	lit = LeftConfigurations.find(LeftConf);
	if ( rit != RightConfigurations.end() )
	  {
	    configurationsConstituentIndices[idx][1] = (*rit).second;
	  }
	if ( lit != LeftConfigurations.end() )
	  {
	    configurationsConstituentIndices[idx][0] = (*lit).second;
	  }
      }
	
    return 0;
}


/*
 This function calculated the fourier transform matrix in a more efficient manner with reuse. 
 
 parameters: structure containing basic system information.
 A: the matrix.
 M: total momentum_end.
 spin_basis: pointer to spin basis array.
 spin_basis_dim: dimension of spin_basis.
 spin_momentum: momentum sector that state is in.
 basis: bosonic basis corresponding to fourier transformed vector.
 dim: dimension of bosonic basis
 normalise: flag to say whether normalisation should be applied or not.
 stop_at: the index to stop at
 end_indices: array of ending indices to use on each process to allow for even load balancing
*/
int fourier_transform_fast_method(struct parameters *parameters, Mat *FTMat, Vec* A, Vec *B, int M, uPetscInt *spin_basis, PetscInt spin_basis_dim, PetscInt spin_momentum, uPetscInt *basis, PetscInt basis_dim, bool normalise, PetscInt stop_at, PetscInt *end_indices, bool use_matrix) {
	PetscInt sites = parameters->AB_sites;
	PetscInt filling = parameters->filling;
	PetscInt orig_vec_length = spin_basis_dim;
	//PetscInt vec_length = basis_dim;
	PetscErrorCode ierr;
	PetscInt starting_idx = 0;	
	//uPetscInt mask = (1ul << parameters->AB_sites) - 1ul; 
	//create global array for fourier transform
			
	// if we are using multiple processors and not on the first one then changing the starting index
	if ( parameters->rank != 0 )
	{
		starting_idx = end_indices[parameters->rank-1];
	}
	
	PetscScalar *original_array, *fourier_array;
	//create the final vector
	if ( use_matrix == true )
	  {  
	    //ierr = MatCreateMPIDense(MPI_COMM_WORLD, end_indices[parameters->rank] - starting_idx, PETSC_DECIDE, basis_dim, spin_basis_dim, PETSC_NULL, FTMat);CHKERRQ(ierr);	
	    ierr = MatCreateDense(MPI_COMM_WORLD, end_indices[parameters->rank] - starting_idx, PETSC_DECIDE, basis_dim, spin_basis_dim, PETSC_NULL, FTMat);CHKERRQ(ierr);	
	    ierr = MatZeroEntries(*FTMat);CHKERRQ(ierr);
	  }
	else
	  {
	    ierr = VecCreateMPI(MPI_COMM_WORLD, end_indices[parameters->rank] - starting_idx, basis_dim, B);CHKERRQ(ierr);
	    VecGetArray(*A, &original_array);
	    VecGetArray(*B, &fourier_array);
	  }
	  
		
	//precalculate the phases
	PetscScalar **phase_array;		
	ierr = PetscMalloc((sites+1)*sizeof(PetscScalar*),&phase_array);CHKERRQ(ierr);
	for( int i = 1  ; i <= sites ; i++ ) {
	    ierr = PetscMalloc((sites+1)*sizeof(PetscScalar),&phase_array[i]);CHKERRQ(ierr);	
	}		
	phase_array[1][0] = 1.0;
	for( int i = 1  ; i <= sites ; i++ ) {	    	    
	    for ( int j = 1; j <= sites ; j++ ) {
		PetscReal angle = ( PETSC_PI * 2.0 * ((PetscReal)(i*j))) / ((PetscReal)sites);
		phase_array[i][j] = cos(angle) + PETSC_i * sin(angle);
		//PetscPrintf(PETSC_COMM_WORLD,"%d, %d: %lf + i%lf\n",i,j,PetscRealPart(phase_array[i][j]), PetscImaginaryPart(phase_array[i][j]));
	    }
	}
		
	PetscScalar tmp_val, tmp_val2;
	PetscScalar norm = sqrt(1.0/pow((PetscReal)sites,(PetscReal)filling)); //normalisation
	
	PetscScalar momentum_phase = 1.0; // defautls to one if not using translation invariance
	
	if ( spin_momentum != -1 )  //if using translation invariance
	{
	    for ( int i = 1 ; i < parameters->AB_sites ; i++ )
	    {
		momentum_phase += phase_array[1][(i*(M+spin_momentum))%parameters->AB_sites];
	    }		    
	}
	norm *= sqrt(momentum_phase);
	
	//now start the real optimisation attempt. We have a tree of ceil(log_2(filling)) levels
	//level zero will have same number of collections as filling.
	//levle two will have this many divided by two and rounded up and so on.
	
	int **BinomialCoefs;
	ierr = PetscMalloc((sites+1)*sizeof(int*), &BinomialCoefs);CHKERRQ(ierr);
	for ( int n = 0 ; n <= sites ; n++ )
	  {
	    ierr = PetscMalloc((sites+1)*sizeof(int*), &BinomialCoefs[n]);CHKERRQ(ierr);
	    for (int c = 0 ; c <= sites ; c++ )
	      {
	        BinomialCoefs[n][c] = (int)NchooseC(n,c);
	      }
	  }
	
	int NumberLevels = (int)ceil(log2(filling)) + 1; 
	//declare and allocate arrays for number of levels.
	int *NumberCollectionsPerLevel; ierr = PetscMalloc(NumberLevels*sizeof(int),&NumberCollectionsPerLevel);CHKERRQ(ierr);
	int **NumberConfigurationsPerCollection; ierr = PetscMalloc(NumberLevels*sizeof(int*),&NumberConfigurationsPerCollection);CHKERRQ(ierr);
	int **NumberParticlesPerCollection; ierr = PetscMalloc(NumberLevels*sizeof(int*),&NumberParticlesPerCollection);CHKERRQ(ierr);
	int ****CombinedConfigurationConstituentIndices; ierr = PetscMalloc(NumberLevels*sizeof(int***),&CombinedConfigurationConstituentIndices);CHKERRQ(ierr);
	int ***CombinedCollectionConstituentIndices; ierr = PetscMalloc(NumberLevels*sizeof(int**),&CombinedCollectionConstituentIndices);CHKERRQ(ierr);
	uPetscInt ***Configurations; ierr = PetscMalloc(NumberLevels*sizeof(uPetscInt**),&Configurations);CHKERRQ(ierr);
	
	//now work out number of collections in each level along with number of particles and allocate further space
	NumberCollectionsPerLevel[NumberLevels-1] = 1;
	ierr = PetscMalloc(NumberCollectionsPerLevel[NumberLevels-1]*sizeof(int),&NumberConfigurationsPerCollection[NumberLevels-1]);CHKERRQ(ierr);
	ierr = PetscMalloc(NumberCollectionsPerLevel[NumberLevels-1]*sizeof(int),&NumberParticlesPerCollection[NumberLevels-1]);CHKERRQ(ierr);
	ierr = PetscMalloc(NumberCollectionsPerLevel[NumberLevels-1]*sizeof(uPetscInt*),&Configurations[NumberLevels-1]);CHKERRQ(ierr);
	ierr = PetscMalloc(NumberCollectionsPerLevel[NumberLevels-1]*sizeof(int**),&CombinedConfigurationConstituentIndices[NumberLevels-1]);CHKERRQ(ierr);
	ierr = PetscMalloc(NumberCollectionsPerLevel[NumberLevels-1]*sizeof(int*),&CombinedCollectionConstituentIndices[NumberLevels-1]);CHKERRQ(ierr);
	for ( int i = 0; i < NumberCollectionsPerLevel[NumberLevels-1]; i++ ) { ierr = PetscMalloc(2*sizeof(int),&CombinedCollectionConstituentIndices[NumberLevels-1][i]);CHKERRQ(ierr);}
	NumberParticlesPerCollection[NumberLevels-1][0] = filling;	
	for( int Level = NumberLevels-2; Level >= 0 ; Level-- )
	  {
	    int count = 0;
	    for( int i = 0; i < NumberCollectionsPerLevel[Level+1]; i++ )
	      {
		if( NumberParticlesPerCollection[Level+1][i] > 1 )
		  {
		    count+=2;
		  } 
		else
		  {
		    count+=1;
		  }
	      }
	    NumberCollectionsPerLevel[Level] = count;
	    ierr = PetscMalloc(NumberCollectionsPerLevel[Level]*sizeof(int),&NumberConfigurationsPerCollection[Level]);CHKERRQ(ierr);
	    ierr = PetscMalloc(NumberCollectionsPerLevel[Level]*sizeof(uPetscInt*),&Configurations[Level]);CHKERRQ(ierr);
	    ierr = PetscMalloc(NumberCollectionsPerLevel[Level]*sizeof(int**),&CombinedConfigurationConstituentIndices[Level]);CHKERRQ(ierr);
	    ierr = PetscMalloc(NumberCollectionsPerLevel[Level]*sizeof(int*),&CombinedCollectionConstituentIndices[Level]);CHKERRQ(ierr);
	    for(int i = 0 ; i < NumberCollectionsPerLevel[Level] ; i++ ) { ierr = PetscMalloc(2*sizeof(int),&CombinedCollectionConstituentIndices[Level][i]);CHKERRQ(ierr);}
	    ierr = PetscMalloc(NumberCollectionsPerLevel[Level]*sizeof(int),&NumberParticlesPerCollection[Level]);CHKERRQ(ierr);	    
	    
	    count = 0;
	    for( int i = 0; i < NumberCollectionsPerLevel[Level+1]; i++ )
	      {
		if( NumberParticlesPerCollection[Level+1][i] > 1 )
		  {
		    NumberParticlesPerCollection[Level][count] = (NumberParticlesPerCollection[Level+1][i]+1)/2;
		    NumberParticlesPerCollection[Level][count+1] = NumberParticlesPerCollection[Level+1][i] - NumberParticlesPerCollection[Level][count];
		    CombinedCollectionConstituentIndices[Level+1][i][0] = count;
		    CombinedCollectionConstituentIndices[Level+1][i][1] = count+1;
		    count+=2;
		  } 
		else
		  {
		    NumberParticlesPerCollection[Level][count] = 1;
		    CombinedCollectionConstituentIndices[Level+1][i][0] = count;
		    CombinedCollectionConstituentIndices[Level+1][i][1] = -1;
		    count+=1;
		  }
	      }	    	    
	  }
	  		  
	NumberConfigurationsPerCollection[NumberLevels-1][0] = orig_vec_length;
	ierr = PetscMalloc(orig_vec_length*sizeof(uPetscInt), &Configurations[NumberLevels-1][0]); CHKERRQ(ierr);
	ierr = PetscMalloc(orig_vec_length*sizeof(int*),&CombinedConfigurationConstituentIndices[NumberLevels-1][0]);CHKERRQ(ierr);
	for ( int i = 0 ; i < orig_vec_length ; i++ ) { ierr = PetscMalloc(2*sizeof(int),&CombinedConfigurationConstituentIndices[NumberLevels-1][0][i]);CHKERRQ(ierr); }	
	//set the configurations on the topmost level
	
	for ( int i = 0 ; i < orig_vec_length ; i++ ) 
	  {
	    Configurations[NumberLevels-1][0][i] = spin_basis[i];
	  }
	  	  
	//now continually decompose till others are all populated.
	for ( int Level = NumberLevels - 1 ; Level > 0 ; Level-- )
	  {
	    for ( int Collection = 0 ; Collection < NumberCollectionsPerLevel[Level] ; Collection++ )
	      {
		if (CombinedCollectionConstituentIndices[Level][Collection][1] >= 0 ) 
		  {
		    DecomposeConfiguration(Configurations[Level][Collection], NumberConfigurationsPerCollection[Level][Collection], NumberParticlesPerCollection[Level][Collection], sites, &NumberConfigurationsPerCollection[Level-1][CombinedCollectionConstituentIndices[Level][Collection][0]],  &NumberConfigurationsPerCollection[Level-1][CombinedCollectionConstituentIndices[Level][Collection][1]], CombinedConfigurationConstituentIndices[Level][Collection], Configurations[Level-1], CombinedCollectionConstituentIndices[Level][Collection][0]);   
		  }
		else
		  {
		    DecomposeConfiguration(Configurations[Level][Collection], NumberConfigurationsPerCollection[Level][Collection], NumberParticlesPerCollection[Level][Collection], sites, &NumberConfigurationsPerCollection[Level-1][CombinedCollectionConstituentIndices[Level][Collection][0]], NULL , CombinedConfigurationConstituentIndices[Level][Collection], Configurations[Level-1], CombinedCollectionConstituentIndices[Level][Collection][0]);   
		  }
	      }
	    for ( int Collection = 0 ; Collection < NumberCollectionsPerLevel[Level-1] ; Collection++ )
	      {
		ierr = PetscMalloc(NumberConfigurationsPerCollection[Level-1][Collection]*sizeof(int*),&CombinedConfigurationConstituentIndices[Level-1][Collection]);CHKERRQ(ierr);
		for (int i = 0 ; i < NumberConfigurationsPerCollection[Level-1][Collection]; i++ )
		  {
		    ierr = PetscMalloc(2*sizeof(int*),&CombinedConfigurationConstituentIndices[Level-1][Collection][i]);CHKERRQ(ierr);
		  }
	      }
	  }

	
	PetscScalar ****PartialSums; 
	ierr = PetscMalloc(NumberLevels*sizeof(PetscScalar***),&PartialSums);CHKERRQ(ierr);
	int **NumberCombinationsPerCollection;
	ierr = PetscMalloc(NumberLevels*sizeof(int*),&NumberCombinationsPerCollection);CHKERRQ(ierr);
	for ( int Level = 0 ; Level < NumberLevels ; Level++ )
	  {
	    ierr = PetscMalloc(NumberCollectionsPerLevel[Level]*sizeof(PetscScalar**),&PartialSums[Level]);CHKERRQ(ierr);
	    ierr = PetscMalloc(NumberCollectionsPerLevel[Level]*sizeof(int),&NumberCombinationsPerCollection[Level]);CHKERRQ(ierr);	  
	    for ( int Collection = 0 ; Collection < NumberCollectionsPerLevel[Level] ; Collection++ )
	      {
		NumberCombinationsPerCollection[Level][Collection] = BinomialCoefs[filling][ NumberParticlesPerCollection[Level][Collection]];
		ierr = PetscMalloc(NumberConfigurationsPerCollection[Level][Collection]*sizeof(PetscScalar*),&PartialSums[Level][Collection]);CHKERRQ(ierr);		
		for ( int Configuration = 0 ; Configuration < NumberConfigurationsPerCollection[Level][Collection] ; Configuration++ )
		  {
		    ierr = PetscMalloc(NumberCombinationsPerCollection[Level][Collection]*sizeof(PetscScalar), &PartialSums[Level][Collection][Configuration]); CHKERRQ(ierr);
		  }
	      }
	  }	
	
	//print out structure for debugging purposes.
	/*for ( int Level = 0 ; Level < NumberLevels; Level++ )
	  {
	    PetscPrintf(PETSC_COMM_WORLD,"Level: %d with %d collections.\n", Level, NumberCollectionsPerLevel[Level]);
	    for( int Collection = 0; Collection < NumberCollectionsPerLevel[Level]; Collection++ )
	      {
	        PetscPrintf(PETSC_COMM_WORLD,"Collection %d, %d particles, %d configurations.\n", Collection, NumberParticlesPerCollection[Level][Collection], NumberConfigurationsPerCollection[Level][Collection]); 
	        for ( int Configuration = 0 ; Configuration < NumberConfigurationsPerCollection[Level][Collection]; Configuration++ )
	          {
	            char buf[MAX_STRING];
	            dec2bin(Configurations[Level][Collection][Configuration], sites, buf);
	            if ( Level > 0 )
	              {
	                PetscPrintf(PETSC_COMM_WORLD,"%d: %s, left: %d, right: %d\n", Configuration, buf, CombinedConfigurationConstituentIndices[Level][Collection][Configuration][0], CombinedConfigurationConstituentIndices[Level][Collection][Configuration][1]);
		      }
		    else
		      {
		        PetscPrintf(PETSC_COMM_WORLD,"%d: %s\n", Configuration, buf);
		      }
	          }
	      }
	     PetscPrintf(PETSC_COMM_WORLD,"\n"); 
	  }*/
	  
	int *Monomial; ierr = PetscMalloc(filling*sizeof(int),&Monomial);CHKERRQ(ierr);
	int *BosonIndices; ierr = PetscMalloc(sizeof(int) * filling, &BosonIndices);CHKERRQ(ierr);
	int *FermionPositions; ierr = PetscMalloc(sizeof(int) * filling, &FermionPositions);CHKERRQ(ierr);
	  
	  
	int *****CombinationDecompositions; ierr = PetscMalloc(filling*sizeof(int****), &CombinationDecompositions); CHKERRQ(ierr);  
	//now take a look at the combinations
	for ( int c = 1 ; c <= filling; c++ )
	  {
	    //PetscPrintf(PETSC_COMM_WORLD,"Combinations for %d out of %d.\n", c, filling);
	    ierr = PetscMalloc(BinomialCoefs[filling][ c]*sizeof(int***), &CombinationDecompositions[c-1]);CHKERRQ(ierr);
	    for ( int idx = 0 ; idx < BinomialCoefs[filling][ c]; idx++ )
	      {
	        get_choosen_configuration(filling, c , idx,  BosonIndices, BinomialCoefs);
	        for ( int i = 0 ; i < c ; i++ )
	          {
	            //PetscPrintf(PETSC_COMM_WORLD,"%d",BosonIndices[i]);
	          }
	        if ( c > 1 ) 
	          {
	            ierr = PetscMalloc((c-1)*sizeof(int**), &CombinationDecompositions[c-1][idx]);CHKERRQ(ierr);
	            for ( int left = 1; left < c ; left++ )
	              {
	              	int right = c - left;
	                //PetscPrintf(PETSC_COMM_WORLD,", (%d, %d): ", left, right);
	                int leftCombination, rightCombination;
	                ierr = PetscMalloc(BinomialCoefs[c][left]*sizeof(int*), &CombinationDecompositions[c-1][idx][left-1]);CHKERRQ(ierr);
	                for ( int leftidx = 0 ; leftidx < BinomialCoefs[c][ left]; leftidx++ )
	                  {
	                    ierr = PetscMalloc(2*sizeof(int), &CombinationDecompositions[c-1][idx][left-1][leftidx]);CHKERRQ(ierr);
	                    get_combination_indices(BosonIndices, filling, c, left, right, &leftCombination, &rightCombination, leftidx, BinomialCoefs);    
	                    //PetscPrintf(PETSC_COMM_WORLD,", %d, %d + %d ", leftidx, leftCombination, rightCombination);
	                    CombinationDecompositions[c-1][idx][left-1][leftidx][0] = leftCombination;
	                    CombinationDecompositions[c-1][idx][left-1][leftidx][1] = rightCombination;
	                  }
	              }
	          }	        
	        //PetscPrintf(PETSC_COMM_WORLD,"\n");  
	      }
	  }
	
	
	
	PetscScalar *row_vals;	
	if ( use_matrix == true ) 
	  {
	    ierr = PetscMalloc(sizeof(PetscScalar)*orig_vec_length, &row_vals);CHKERRQ(ierr);			
	  }
	  
	//PetscInt fstart = get_rank_chunk_start(vec_length, parameters->rank, parameters->size);
	PetscInt fstart = starting_idx;
	//PetscInt fend = get_rank_chunk_size(vec_length, parameters->rank, parameters->size);
	PetscInt fend = end_indices[parameters->rank];
	int phase_idx;  			
	PetscReal start_time, end_time; 
	PetscReal state_normal;
	start_time = MPI_Wtime();
	int k,l,tmp;		
	if ( stop_at > 0 && stop_at < (fend - fstart) ) fend = fstart + stop_at;
	for ( PetscInt fidx = fstart ; fidx < fend ; fidx++ )    	      
	//for ( int fidx  = 0; fidx < 1 ; fidx++ )
	  {	
	    //for each bosonic configuration we start filling PartialSums in from the bottom level. 	    
	    //get_monomial_bosons(parameters, basis[fidx], sites, filling, Monomial);
	    get_monomial_bosons(parameters, basis[fidx], sites, filling, Monomial, &state_normal);	    	    
	    
	    //populate level 0 first. 
	    for ( int Collection = 0; Collection < NumberCollectionsPerLevel[0] ; Collection++ )
	      {		
		for (int Combination = 0; Combination < NumberCombinationsPerCollection[0][Collection]; Combination++ )
		  {
		    get_choosen_configuration(filling, NumberParticlesPerCollection[0][Collection], Combination,  BosonIndices, BinomialCoefs);
		    for ( int Configuration = 0; Configuration < NumberConfigurationsPerCollection[0][Collection]; Configuration++)
		      {
			get_monomial_fermions(parameters, Configurations[0][Collection][Configuration], sites, NumberParticlesPerCollection[0][Collection], FermionPositions);	
			int phase_idx2 = 0;
			for ( int i = 0 ; i < NumberParticlesPerCollection[0][Collection]; i++ )
			  {
			    phase_idx2 += FermionPositions[i] * Monomial[BosonIndices[i]];			    
			  }
			PartialSums[0][Collection][Configuration][Combination] = phase_array[1][(phase_idx2) % sites];
		      }
		  }
	      }	
	      
	    //now iterate up to the other levels.
	    for ( int Level = 1 ; Level < NumberLevels ; Level++ )
	      {
		for ( int Collection = 0; Collection < NumberCollectionsPerLevel[Level]; Collection++ )
		  {
		    for ( int Combination =0 ; Combination < NumberCombinationsPerCollection[Level][Collection]; Combination++ )
		      {
			//get_choosen_configuration(filling, NumberParticlesPerCollection[Level][Collection], Combination,  BosonIndices, BinomialCoefs);
			for( int Configuration = 0; Configuration < NumberConfigurationsPerCollection[Level][Collection]; Configuration++ )
			  {
			    if( CombinedCollectionConstituentIndices[Level][Collection][1] == -1 ) 
			      {
				PartialSums[Level][Collection][Configuration][Combination] = PartialSums[Level-1][CombinedCollectionConstituentIndices[Level][Collection][0]][Configuration][Combination];
			      }
			    else
			      {
				PartialSums[Level][Collection][Configuration][Combination]=0;
				for (int leftCombIdx = 0; leftCombIdx < BinomialCoefs[NumberParticlesPerCollection[Level][Collection]][NumberParticlesPerCollection[Level-1][CombinedCollectionConstituentIndices[Level][Collection][0]]]; leftCombIdx++ )
				  {
				    int leftCombination; int rightCombination;
				    //get_combination_indices(BosonIndices, filling, NumberParticlesPerCollection[Level][Collection], NumberParticlesPerCollection[Level-1][CombinedCollectionConstituentIndices[Level][Collection][0]], NumberParticlesPerCollection[Level-1][CombinedCollectionConstituentIndices[Level][Collection][1]], &leftCombination, &rightCombination, leftCombIdx, BinomialCoefs);
				    leftCombination = CombinationDecompositions[NumberParticlesPerCollection[Level][Collection]-1][Combination][NumberParticlesPerCollection[Level-1][CombinedCollectionConstituentIndices[Level][Collection][0]]-1][leftCombIdx][0];
				    rightCombination = CombinationDecompositions[NumberParticlesPerCollection[Level][Collection]-1][Combination][NumberParticlesPerCollection[Level-1][CombinedCollectionConstituentIndices[Level][Collection][0]]-1][leftCombIdx][1];
				    if ( (leftCombination < NumberCombinationsPerCollection[Level-1][CombinedCollectionConstituentIndices[Level][Collection][0]] )
				       && (rightCombination < NumberCombinationsPerCollection[Level-1][CombinedCollectionConstituentIndices[Level][Collection][1]] ) )
				       { 				          
				          PartialSums[Level][Collection][Configuration][Combination] += PartialSums[Level-1][CombinedCollectionConstituentIndices[Level][Collection][0]][CombinedConfigurationConstituentIndices[Level][Collection][Configuration][0]][leftCombination]
				          *PartialSums[Level-1][CombinedCollectionConstituentIndices[Level][Collection][1]][CombinedConfigurationConstituentIndices[Level][Collection][Configuration][1]][rightCombination];
				       }
				       
				  }				
			      }
			  }
		      }
		  }
		
	      }
	    
	    /*if ( fidx == 0 ) 
	      {
	        PetscPrintf(PETSC_COMM_WORLD,"Bosonic config:");
	        for ( int i = 0 ; i < filling ; i++ )
	          {
	            PetscPrintf(PETSC_COMM_WORLD, "%d", Monomial[i]);
	          }
	        PetscPrintf(PETSC_COMM_WORLD,"\n");  
	        //now iterate up to the other levels.
  	        for ( int Level = 0 ; Level < NumberLevels ; Level++ )
	          {
		    for ( int Collection = 0; Collection < NumberCollectionsPerLevel[Level]; Collection++ )
		      {
		        for( int Configuration = 0; Configuration < NumberConfigurationsPerCollection[Level][Collection]; Configuration++ )
		          {
			    for ( int Combination =0 ; Combination < NumberCombinationsPerCollection[Level][Collection]; Combination++ )
			      {
			    	PetscPrintf(PETSC_COMM_WORLD,"[%d,%d,%d,%d]: %.14g + i %.14g (%.14g + i %.14g)\n", Level, Collection, Configuration, Combination, PetscRealPart(PartialSums[Level][Collection][Configuration][Combination]), PetscImaginaryPart(PartialSums[Level][Collection][Configuration][Combination]),PetscRealPart(PartialSums[Level][Collection][Configuration][Combination]*norm/state_normal), PetscImaginaryPart(PartialSums[Level][Collection][Configuration][Combination]*norm/state_normal) );
			      }
		          }
		      }
		
	          }
	      }*/
	    
	    	    	    
	    //Now finall put the sums into the final array and apply normalisation.
	    if ( use_matrix == true ) 
	      {
		for ( PetscInt oidx = 0 ; oidx < orig_vec_length ; oidx++ )
		  {
		    row_vals[oidx] = PartialSums[NumberLevels-1][0][oidx][0] * norm / state_normal;
		    //ierr = MatSetValue(*A, fidx, oidx, row_vals[oidx], INSERT_VALUES);CHKERRQ(ierr);
		  }		  
		ierr = MatSetValuesRowLocal(*FTMat, fidx, row_vals); CHKERRQ(ierr);  
	      }
	    else
	      {
		fourier_array[fidx - fstart] = 0;
		for ( PetscInt oidx = 0 ; oidx < orig_vec_length ; oidx++ )
		  {
		    fourier_array[fidx - fstart] += original_array[oidx] * PartialSums[NumberLevels-1][0][oidx][0] * norm / state_normal;		    
		  }		  
	      }
	  }	
	
	if ( use_matrix == true ) 
	  {
	    MatAssemblyBegin(*FTMat,MAT_FINAL_ASSEMBLY);
	    MatAssemblyEnd(*FTMat,MAT_FINAL_ASSEMBLY);
	    ierr = PetscFree(row_vals);CHKERRQ(ierr);
	  }
	else
	 {
	    VecRestoreArray(*A, &original_array);
	    VecRestoreArray(*B, &fourier_array);
	 }
		
	end_time = MPI_Wtime();
	
	if ( use_matrix ) 
	  {
	    PetscPrintf(PETSC_COMM_SELF,"Rank: %d, Matrix construction took: %lf seconds.\n", parameters->rank, end_time - start_time);
	  }
	else
	  {
	    PetscPrintf(PETSC_COMM_SELF,"Rank: %d, Fourier transform took: %lf seconds.\n", parameters->rank, end_time - start_time);
	  }
		
	//need to now deallocate all the space.
	for( int i = 1  ; i <= sites; i++ ) 
	{
	    ierr = PetscFree(phase_array[i]);CHKERRQ(ierr);	
	}
	ierr = PetscFree(phase_array);CHKERRQ(ierr);	
	

	//now take a look at the combinations
	for ( int c = 1 ; c <= filling; c++ )
	  {    
	    for ( int idx = 0 ; idx < BinomialCoefs[filling][ c]; idx++ )
	      {
	        if ( c > 1 ) 
	          {
	            for ( int left = 1; left < c ; left++ )
	              {
	                for ( int leftidx = 0 ; leftidx < BinomialCoefs[c][ left]; leftidx++ )
	                  {
	                    ierr = PetscFree(CombinationDecompositions[c-1][idx][left-1][leftidx]);CHKERRQ(ierr);	                    
	                  }
	                ierr = PetscFree(CombinationDecompositions[c-1][idx][left-1]);CHKERRQ(ierr);  
	              }
	            ierr = PetscFree(CombinationDecompositions[c-1][idx]);CHKERRQ(ierr);  
	          }	        
	      }
	    ierr = PetscMalloc(BinomialCoefs[filling][ c]*sizeof(int***), &CombinationDecompositions[c-1]);CHKERRQ(ierr);  
	  }
	ierr = PetscFree(CombinationDecompositions);CHKERRQ(ierr);
	
	
	for ( int n = 0; n <= sites ; n++ )
	  {
	    ierr = PetscFree(BinomialCoefs[n]);CHKERRQ(ierr);
	  }
	ierr = PetscFree(BinomialCoefs);CHKERRQ(ierr);  
	
	
	for(int Level = 0 ; Level < NumberLevels ; Level++) 
	  {
	    for ( int Collection = 0 ; Collection < NumberCollectionsPerLevel[Level] ; Collection++ )
	      {
	        for ( int Configuration = 0 ; Configuration < NumberConfigurationsPerCollection[Level][Collection] ; Configuration++ )
	          {	            
	            ierr = PetscFree(CombinedConfigurationConstituentIndices[Level][Collection][Configuration]);CHKERRQ(ierr);
	            ierr = PetscFree(CombinedConfigurationConstituentIndices[Level][Collection][Configuration]);CHKERRQ(ierr);
	            ierr = PetscFree(PartialSums[Level][Collection][Configuration]);CHKERRQ(ierr);
	          }
	        ierr = PetscFree(Configurations[Level][Collection]);CHKERRQ(ierr);    
	        ierr = PetscFree(CombinedConfigurationConstituentIndices[Level][Collection]);CHKERRQ(ierr);
	        ierr = PetscFree(CombinedCollectionConstituentIndices[Level][Collection]);CHKERRQ(ierr);
	        ierr = PetscFree(PartialSums[Level][Collection]);CHKERRQ(ierr);
	      }
	    ierr = PetscFree(NumberConfigurationsPerCollection[Level]);CHKERRQ(ierr);
	    ierr = PetscFree(NumberParticlesPerCollection[Level]);CHKERRQ(ierr);  
	    ierr = PetscFree(NumberCombinationsPerCollection[Level]);CHKERRQ(ierr);  
  	    ierr = PetscFree(Configurations[Level]);CHKERRQ(ierr);    
  	    ierr = PetscFree(CombinedConfigurationConstituentIndices[Level]);CHKERRQ(ierr);
  	    ierr = PetscFree(CombinedCollectionConstituentIndices[Level]);CHKERRQ(ierr);
  	    ierr = PetscFree(PartialSums[Level]);CHKERRQ(ierr);
	  }
	
	ierr = PetscFree(NumberCollectionsPerLevel);CHKERRQ(ierr);
	ierr = PetscFree(NumberConfigurationsPerCollection);CHKERRQ(ierr);
	ierr = PetscFree(NumberParticlesPerCollection);CHKERRQ(ierr);
	ierr = PetscFree(CombinedConfigurationConstituentIndices);CHKERRQ(ierr);
	ierr = PetscFree(CombinedCollectionConstituentIndices);CHKERRQ(ierr);
	ierr = PetscFree(PartialSums);CHKERRQ(ierr);
	ierr = PetscFree(Monomial);CHKERRQ(ierr);
	ierr = PetscFree(FermionPositions);CHKERRQ(ierr);
	ierr = PetscFree(BosonIndices);CHKERRQ(ierr);
	
	return 0;
}


int main( int argc, char **argv ) {
    PetscErrorCode ierr; 
    string adj_file, vec_file, a_sites, out_prefix, basis_file, name, matrix_file;
    struct parameters params; 
    char buf[MAX_STRING];
    PetscTruth flg;
    stringstream ss;   
    bool verbose, basis_given = false, make_real = false, truncate = true, use_complex, calc_variance_ordered = false;
    PetscInt stop_at, spin_momentum ;    
    
    
    
    /*------------------------------------------------------------------------------------------------
    Initialise some values.
    ------------------------------------------------------------------------------------------------*/
    char help[] = "Program to generate calculate entanglement entropy and spectrum in parallel for cut in fourier space.\n";
    params.A_filling = -1;	params.B_filling = -1; params.filling = -1; //by default any fillings can be used.
    params.nn_exclusion = PETSC_FALSE;
    
    SlepcInitialize(&argc,&argv,(char*)0,help);	
    
    //get rank of this process and number of processes 
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&(params.rank));CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&(params.size));CHKERRQ(ierr);
    
    /*------------------------------------------------------------------------------------------------
    Read in input details
    ------------------------------------------------------------------------------------------------*/	
    ierr = PetscOptionsGetInt(NULL,"-num_sites",&params.AB_sites,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	PetscPrintf(PETSC_COMM_WORLD,"Number of sites not specified.\n");
	params.AB_sites = -1;
	return 1;
      } 
    else 
      {
	PetscPrintf(PETSC_COMM_WORLD,"Number of sites: %d\n",params.AB_sites);
      }
    
    ierr = PetscOptionsGetInt(NULL,"-filling",&params.filling,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	PetscPrintf(PETSC_COMM_WORLD,"No filling specified.\n");
	params.filling = -1;
	return 1;
      } 
    else 
      {
	PetscPrintf(PETSC_COMM_WORLD,"Filling specified: %d\n",params.filling);
      }
      
    ierr = PetscOptionsGetInt(NULL,"-momentum",&params.total_momentum,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	params.total_momentum = -1;
	PetscPrintf(PETSC_COMM_WORLD,"No momentum specified. Will scan all possible");		
      } 
    else 
      {
	PetscPrintf(PETSC_COMM_WORLD,"Momentum specified: %d\n",params.total_momentum);
      }
    
    ierr = PetscOptionsGetInt(NULL,"-A_momentum",&params.A_total_momentum,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	params.A_total_momentum = -1;
	PetscPrintf(PETSC_COMM_WORLD,"No A momentum specified. Will not use specific value.\n");
	
      } 
    else 
      {
	PetscPrintf(PETSC_COMM_WORLD,"A momentum specified: %d\n",params.A_total_momentum);
      }	
    
    ierr = PetscOptionsGetInt(NULL,"-A_filling",&params.A_filling,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	PetscPrintf(PETSC_COMM_WORLD,"No A filling specified.\n");
	params.A_filling = -1;
	params.B_filling = -1;
	//return 1;
      } 
    else 
      {
	PetscPrintf(PETSC_COMM_WORLD,"A filling specified: %d\n",params.A_filling);
	if ( params.filling != -1 ) 
	  {
	    params.B_filling = params.filling - params.A_filling ;
	  }
      }
    
    ierr = PetscOptionsGetString(NULL,"-vec",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	PetscPrintf(PETSC_COMM_WORLD,"Must provide vector on the command line with switch -vec.\n");
	return 1;
      }
    PetscPrintf(PETSC_COMM_WORLD,"Vec file: %s\n",buf);
    vec_file.assign(buf);
    
    ierr = PetscOptionsGetString(NULL,"-name",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	PetscPrintf(PETSC_COMM_WORLD,"Using name \"default\".\n");	
	name.assign("default");
      }
    else 
      {
	PetscPrintf(PETSC_COMM_WORLD,"Using name: %s\n",buf);
	name.assign(buf);
      }
    
    ierr = PetscOptionsGetString(NULL,"-ascii",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	params.ascii = PETSC_FALSE;
      } 
    else 
      {
	if ( strcmp(buf,"true") == 0 ) 
	  {
	    params.ascii = PETSC_TRUE;
	  } 
	else 
	  {	
	    params.ascii = PETSC_FALSE;
	  }
      }
        
    bool normalise = true;
    ierr = PetscOptionsGetString(NULL,"-normalise",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE) 
      {
	if ( strcmp(buf,"false") == 0 ) 
	  {
	    normalise = false;
	  } 
	else if ( strcmp(buf,"true") == 0 ) 
	  {	
	    normalise = true;
	  }
      }
    
    verbose = false;
    ierr = PetscOptionsGetString(NULL,"-verbose",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE) 
      {
	if ( strcmp(buf,"false") == 0 ) 
	  {
	    verbose = false;
	  } 
	else if ( strcmp(buf,"true") == 0 ) 
	  {	
	    verbose = true;
	  }
      }
      
    //this option is for the case where the output vector is fully imaginary to make it real by multiplying by -i. this is false by default.
    make_real = false;
    ierr = PetscOptionsGetString(NULL,"-make_real",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE) 
      {
	if ( strcmp(buf,"false") == 0 ) 
	  {
	    make_real = false;
	  } 
	else if ( strcmp(buf,"true") == 0 ) 
	  {	
	    make_real = true;
	  }
      }
      
    // this truncate option tell whether or not should write output in truncated basis where highest momentum state is removed. By default this is the case.
    truncate = true;
    ierr = PetscOptionsGetString(NULL,"-truncate",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE) 
      {
	if ( strcmp(buf,"false") == 0 ) 
	  {
	    truncate = false;
	  } 
	else if ( strcmp(buf,"true") == 0 ) 
	  {	
	    truncate = true;
	  }
      }
    
    ierr = PetscOptionsGetInt(NULL,"-stopat",&stop_at,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE) 
      {	
	PetscPrintf(PETSC_COMM_WORLD,"Will stop after %ld elements.\n", stop_at);		
      } 
    else 
      {
	stop_at = -1;
      }
      
    
    ierr = PetscOptionsGetString(NULL,"-basis",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE) 
      {
	PetscPrintf(PETSC_COMM_WORLD,"Basis file specified: %s.\n", buf);
	basis_file.assign(buf);
	basis_given = true;
      }
      
    ierr = PetscOptionsGetInt(NULL,"-spin_momentum",&spin_momentum,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE ) 
      {	
	PetscPrintf(PETSC_COMM_WORLD,"Using spin momentum %ld.\n", spin_momentum);		
      } 
    else 
      {
	spin_momentum = -1;
      }
      
    params.spin_one = false;
    ierr = PetscOptionsGetString(NULL,"-spin_one",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE) 
      {
	if ( strcmp(buf,"false") == 0 ) 
	  {
	    params.spin_one = false;
	  } 
	else if ( strcmp(buf,"true") == 0 ) 
	  {	
	    params.spin_one = true;
	    if ( basis_given == false ) 
	      {
		PetscPrintf(PETSC_COMM_WORLD,"For spin one systems the basis must be provided.\n");
		return -1;
	      }
	  }
      }      
      
    use_complex = false;
    ierr = PetscOptionsGetString(NULL,"-complex",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE) 
      {
	if ( strcmp(buf,"false") == 0 ) 
	  {
	    use_complex = false;
	  } 
	else if ( strcmp(buf,"true") == 0 ) 
	  {	
	    use_complex = true;	    
	    PetscPrintf(PETSC_COMM_WORLD,"Using complex part of wavefunction.\n");		
	  }
      }    
    
    ierr = PetscOptionsGetInt(NULL,"-step",&params.step,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE) 
      {	
	PetscPrintf(PETSC_COMM_WORLD,"Translation step is: %d.\n", params.step);		
      } 
    else 
      {
	params.step = 1;
      }
      
    bool use_matrix = false;    
    matrix_file.assign("");
    ierr = PetscOptionsGetString(NULL,"-matrix",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE) 
      {
	if ( strcmp(buf,"false") == 0 ) 
	  {
	    use_matrix = false;
	  } 
	else if ( strcmp(buf,"true") == 0 ) 
	  {	
	    use_matrix = true;	    
	    PetscPrintf(PETSC_COMM_WORLD,"Calculating matrix.\n");		
	  }
	else  //otherwise if a file is specified then read this as the matrix
	  {
	    matrix_file.assign(buf);
	    use_matrix = true;
	    PetscPrintf(PETSC_COMM_WORLD,"Using matrix from file: %s\n",buf);
	  }
      }    
          
    ierr = PetscOptionsGetString(NULL,"-variance",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE) 
      {
	if ( strcmp(buf,"false") == 0 ) 
	  {
	    calc_variance_ordered = false;
	  } 
	else if ( strcmp(buf,"true") == 0 ) 
	  {	
	    calc_variance_ordered = true;	    
	    PetscPrintf(PETSC_COMM_WORLD,"Calculating variance ordered matrix too.\n");		
	  }
      }    
      
    bool fqft = false;  
    ierr = PetscOptionsGetString(NULL,"-fqft",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE ) 
      {
	if ( strcmp(buf,"false") == 0 ) 
	  {
	    fqft = false;
	  } 
	else if ( strcmp(buf,"true") == 0 ) 
	  {	
	    fqft = true;
	    PetscPrintf(PETSC_COMM_WORLD,"Using fast quantum fourier transform method.\n");		
	  }
      }    
          
    /*------------------------------------------------------------------------------------------------
    Finished reading input so now go to work.
    ------------------------------------------------------------------------------------------------*/	
    PetscPrintf(PETSC_COMM_WORLD,"Number of sites: %d\n",params.AB_sites);
    
    params.number_particles = params.AB_sites;	
    
    //load in the vector as local vector on each processor.
    PetscViewer viewer;
    PetscViewerCreate(PETSC_COMM_SELF,&viewer);
    if ( params.ascii ) 
      {
	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, vec_file.c_str() ,&viewer);CHKERRQ(ierr); //ascii not working well for large sizes at the moment.
      } 
    else 
      {
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, vec_file.c_str() ,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      }
    ierr = VecLoad(params.psi, viewer);CHKERRQ(ierr);
    PetscViewerDestroy(&viewer);
    
    //do check to make sure the length of the vector matches the number of AB configs.
    PetscInt vec_size;
    VecGetSize(params.psi,&vec_size);
    
    PetscPrintf(PETSC_COMM_WORLD,"Vector size: %d\n",vec_size);
    uPetscInt full_index;
    int *monomial;
    ierr = PetscMalloc(params.filling*sizeof(int),&monomial);CHKERRQ(ierr);
    
    uPetscInt *spin_basis;
    PetscInt spin_basis_dim;
    if (basis_given)
    {
      spin_basis = read_spin_basis(basis_file, &spin_basis_dim);
      if ( spin_basis == NULL || vec_size != spin_basis_dim )
      {
	 PetscPrintf(PETSC_COMM_WORLD, "Basis dimension does not match input vector or problem reading basis file.\n");
	 return 1;
      }
    } else 
    {
      spin_basis_dim = vec_size;
      ierr = PetscMalloc(spin_basis_dim * sizeof(uPetscInt), &spin_basis);
      for ( PetscInt idx = 0 ; idx < spin_basis_dim ; idx++ )
      {
	  spin_basis[idx] = get_full_basis_index(&params, params.AB_sites, params.filling, idx);
      }
    }        
    
    if ( verbose ) 
      {
	PetscScalar *vals;
	VecGetArray(params.psi,&vals);
	for ( PetscInt i = 0 ; i < vec_size ; i++ ) 
	  {
	    //full_index = get_full_basis_index(&params, params.AB_sites, params.filling, i);
	    full_index = spin_basis[i];
	    if ( params.spin_one ) {
		get_monomial_spin_one(&params, full_index, params.AB_sites, params.filling, monomial);
		dec2bin(full_index, 2*params.AB_sites, buf );
	    } else {
		get_monomial_fermions(&params, full_index, params.AB_sites, params.filling, monomial);
		dec2bin(full_index, params.AB_sites, buf );
	    }	    
	    PetscPrintf(PETSC_COMM_WORLD,"%d: %u: %s, ",i,full_index, buf);
	    for ( int j = 0 ; j < params.filling ; j++ )
	      {
		PetscPrintf(PETSC_COMM_WORLD, "%d :", monomial[j]);
	      }
	    PetscPrintf(PETSC_COMM_WORLD,"%lg + i %lg", PetscRealPart(vals[i]), PetscImaginaryPart(vals[i])); 
	    PetscPrintf(PETSC_COMM_WORLD,"\n");
	  }
	VecRestoreArray(params.psi, &vals);
      } 
    
      if ( use_matrix == false )
	{
	  PetscScalar *vals;
	  PetscInt nnz = 0;
	  VecGetArray(params.psi,&vals);
	  for ( PetscInt i = 0 ; i < vec_size ; i++ ) 
	    {
	      if ( PetscAbsScalar(vals[i]) > INPUT_VEC_TOLERANCE ) 
		{
		  nnz++;
		}	     
	    }
	  VecRestoreArray(params.psi, &vals);
	  if ( nnz < vec_size )
	    {
	      Vec tmp_vec;
	      VecCreateSeq(PETSC_COMM_SELF, nnz, &tmp_vec);	
	      VecGetArray(params.psi,&vals);
	      nnz = 0; 
	      for ( PetscInt i = 0 ; i < vec_size ; i++ ) 
		{
		  if ( PetscAbsScalar(vals[i]) > INPUT_VEC_TOLERANCE ) 
		    {
		      VecSetValue(tmp_vec, nnz, vals[i], INSERT_VALUES);
		      spin_basis[nnz] = spin_basis[i];
		      nnz++;
		    }	     
		}
	      VecRestoreArray(params.psi, &vals);	
	      VecAssemblyBegin(tmp_vec);
	      VecAssemblyEnd(tmp_vec);
	      PetscPrintf(PETSC_COMM_WORLD,"Not using %d zero coefficients.\n", vec_size - nnz);
	      spin_basis_dim = nnz;
	      vec_size = nnz;
	      VecDestroy(&params.psi);
	      VecCreateSeq(PETSC_COMM_SELF, nnz, &params.psi);
	      VecCopy(tmp_vec, params.psi);
	      VecDestroy(&tmp_vec);
	  }
      }
    
    //now the boson basis

    /*for ( PetscInt i = 0 ; i < NchooseC(params.AB_sites+params.filling-1,params.filling) ; i++ ) {
      full_index = get_full_basis_index(&params, params.AB_sites + params.filling, params.filling, i);
      get_monomial_bosons(&params, full_index, params.AB_sites, params.filling, monomial);
      dec2bin(full_index, params.AB_sites + params.filling, buf );
      PetscPrintf(PETSC_COMM_WORLD,"%d: %s, ",i,buf);
      for ( int j = 0 ; j < params.filling ; j++ ){
	  PetscPrintf(PETSC_COMM_WORLD, "%d ", monomial[j]);
      }
      PetscPrintf(PETSC_COMM_WORLD,"\n");
    }*/
    
    PetscReal FTMagnitude = 0.0;
    int momentum_start = params.filling;
    int momentum_end = params.AB_sites * params.filling;
    if ( params.total_momentum != -1 )
      {
	momentum_start = params.total_momentum ;
	momentum_end = params.total_momentum;
      }
      
    uPetscInt mask = (1ul << params.AB_sites) - 1ul; 
    PetscScalar *original_array;
    ierr = VecGetArray(params.psi, &original_array);CHKERRQ(ierr);
    if ( spin_momentum != -1 )
      {
	if ( params.spin_one == false )
	  {
	    for ( PetscInt oidx = 0 ; oidx < vec_size ; oidx++ ) 
	      {
		original_array[oidx] *= configuration_norm(spin_basis[oidx], params.AB_sites, mask, params.step);
	      }
	  }
	else
	  {
	    for ( PetscInt oidx = 0 ; oidx < vec_size ; oidx++ ) 
	      {
		original_array[oidx] *= configuration_norm_spin_one(spin_basis[oidx], params.AB_sites, mask);
	      }
	  }	  
      }
    ierr = VecRestoreArray(params.psi, &original_array);CHKERRQ(ierr);
      
    //for ( int M = params.filling ; M <= params.AB_sites * params.filling ; M++ )
    for ( int M = momentum_start ; M <= momentum_end ; M++ )
      {	    
	ss.str(""); ss << vec_file ;
	if ( params.filling != -1 ) ss << "_F" << params.filling ;
	if ( params.A_filling != -1 ) ss << "_AF" << params.A_filling ;
	//ss << "_As" << params.A_sites ;
	ss << "_M" << M ;
	out_prefix = ss.str();
	
	//now get the momentum basis. Will store full momentum basis on each node. 
	uPetscInt *basis;
	PetscInt dim = boson_basis_dimension(params.AB_sites, params.filling, M);
	ierr = PetscMalloc(dim*sizeof(uPetscInt), &basis);
	
	PetscInt* end_indices;
	ierr = PetscMalloc(params.size * sizeof(PetscInt), &end_indices);CHKERRQ(ierr);
		
	boson_basis_generate(params.AB_sites, params.filling, M, 0, basis, 0ul);
	PetscPrintf(PETSC_COMM_WORLD,"Total M: %d, Dimension: %d\n", M, dim);

	

	if ( ((use_matrix == false) || (strcmp(matrix_file.c_str(),"") == 0)) && (fqft == false) )	  
	  {
	    PetscInt* permutations;
	    ierr = PetscMalloc(dim*sizeof(PetscInt), &permutations);CHKERRQ(ierr);	
	    uPetscInt total_permutations = 0;
	    for ( int j = 0 ; j < dim ; j++ ) 
	      {
		permutations[j] = get_number_permutations(&params, basis[j], params.AB_sites, params.filling);
		total_permutations += permutations[j];
	      }	    	    
	    
	    PetscInt proc_idx = 0;
	    PetscInt running_permutation_count = 0;
	    PetscInt running_idx_count = 0;
	    while ( running_idx_count < dim )
	      {	
		running_permutation_count += permutations[running_idx_count];
		if ( running_permutation_count >= get_rank_chunk_start(total_permutations, proc_idx, params.size) + get_rank_chunk_size(total_permutations, proc_idx, params.size) )
		  {
		    end_indices[proc_idx++] = running_idx_count + 1;
		    PetscPrintf(PETSC_COMM_WORLD,"Ending idx for %d: %d\n",proc_idx-1, end_indices[proc_idx-1]);
		  }
		running_idx_count++;
	      }
	    
	    if ( proc_idx < params.size ) 
	      {
		 while ( proc_idx < params.size ) 
		  {
		    end_indices[proc_idx++] = dim;
		  }
	      }
	    PetscFree(permutations);
	  }
	else
	  {
	    PetscInt proc_idx = 0;
	    while ( proc_idx < (params.size-1) ) 
	      {
		end_indices[proc_idx] = get_rank_chunk_start(dim, proc_idx+1, params.size);
		proc_idx++;
	      }
	    end_indices[params.size-1] = dim;
	  }
		
	Vec ft;
	Mat A;
    	
	if ( params.spin_one == true) 
	  {
	    fourier_transform_vector_spin_one_opt(&params, &params.psi, &ft, M, spin_basis, spin_basis_dim, spin_momentum, basis, dim, normalise, stop_at, end_indices);
	  }
	else 
	  {	    
	    /*if ( verbose ) 
	     {
		fourier_transform_vector_debug(&params, &params.psi, &ft, M, spin_basis, spin_basis_dim, spin_momentum, basis, dim, normalise, stop_at);
	     }
	    else 
	     {*/
	     if ( use_matrix ) 
	       {     
		  if ( strcmp(matrix_file.c_str(), "" ) != 0 ) //then there is a matrix specified
		    {		      
		      PetscViewer viewer;
		      PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
		      //PetscViewerASCIIOpen(PETSC_COMM_WORLD,ss.str().c_str(),&viewer);
		      //PetscViewerMatlabOpen(PETSC_COMM_WORLD, ss.str().c_str(),FILE_MODE_WRITE,&viewer);
		      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,matrix_file.c_str(),FILE_MODE_READ,&viewer);CHKERRQ(ierr);
		      ierr = MatLoad(A, viewer);CHKERRQ(ierr);
		      PetscViewerDestroy(&viewer);
		    }
		  else 
		    {
		      if ( fqft == false ) 
			{
			  fourier_transform_matrix_opt(&params, &A, M, spin_basis, spin_basis_dim, spin_momentum, basis, dim, normalise, stop_at, end_indices);
			}
		      else
			{
			  fourier_transform_fast_method(&params, &A, &params.psi, &ft, M, spin_basis, spin_basis_dim, spin_momentum, basis, dim, normalise, stop_at, end_indices, use_matrix);			  
			}
		    }		    
		  PetscInt starting_idx = 0;
		  if ( params.rank != 0 )
		    {
		      starting_idx = end_indices[params.rank-1];
		    }		   		  
		  ierr = VecCreateMPI(MPI_COMM_WORLD, end_indices[params.rank] - starting_idx ,dim,&ft);CHKERRQ(ierr);			  
		  		  
		  //create global version of input vector to match matrix distribution
		  Vec C;  		  		  		  
		  ierr = VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE ,spin_basis_dim,&C);CHKERRQ(ierr);	
		  PetscScalar *Cvals, *vec_vals;
		  VecGetArray(params.psi, &vec_vals);
		  VecGetArray(C, &Cvals);
		  for ( int i = 0 ; i < get_rank_chunk_size(spin_basis_dim, params.rank, params.size); i++ )
		    {
		      Cvals[i] = vec_vals[get_rank_chunk_start(spin_basis_dim, params.rank, params.size) + i];
		    }
		  VecRestoreArray(params.psi, &vec_vals);
		  VecRestoreArray(C, &Cvals);
		  ierr = MatMult(A, C, ft);CHKERRQ(ierr);
		  VecDestroy(&C);
	       }
	     else
	       {
		 if ( fqft == true )
		   {
		     fourier_transform_fast_method(&params, &A, &params.psi, &ft, M, spin_basis, spin_basis_dim, spin_momentum, basis, dim, normalise, stop_at, end_indices, use_matrix);			  
		   }
		 else
		   { 
		     fourier_transform_vector_opt(&params, &params.psi, &ft, M, spin_basis, spin_basis_dim, spin_momentum, basis, dim, normalise, stop_at, end_indices);
		   }
	       }	     	     
	  }
		
	if ( use_matrix  ) 
	  {
	  	ss.str(""); ss << out_prefix << "_FT.mat" ;
		PetscViewer viewer;
		PetscViewerCreate(PETSC_COMM_WORLD,&viewer);		
		//PetscViewerASCIIOpen(PETSC_COMM_WORLD,ss.str().c_str(),&viewer);
		//PetscViewerMatlabOpen(PETSC_COMM_WORLD, ss.str().c_str(),FILE_MODE_WRITE,&viewer);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD,ss.str().c_str(),FILE_MODE_WRITE,&viewer);
		PetscViewerSetFormat(viewer, PETSC_VIEWER_NATIVE);
		MatView(A,viewer);
		PetscViewerDestroy(&viewer);	    
		MatDestroy(&A);
	  }
	
	PetscScalar prod;
	VecDot(params.psi,params.psi,&prod);
	PetscPrintf(PETSC_COMM_WORLD,"Original magnitude: %lg\n",PetscRealPart(prod));
	VecDot(ft,ft,&prod);
	PetscPrintf(PETSC_COMM_WORLD,"FT magnitude: %lg\n",PetscRealPart(prod));
	FTMagnitude += PetscRealPart(prod);
	
	if ( PetscAbsScalar(prod) > 1e-12 ) 
	  {		    	  	
	    if ( make_real ) 	     
	      {
		ierr = VecScale(ft,-PETSC_i);CHKERRQ(ierr);
	      }
	    
	    ss.str(""); ss << out_prefix << "_FT.vec" ;
	    PetscViewer viewer;
	    PetscViewerCreate(PETSC_COMM_WORLD,&viewer);	    
	    PetscViewerASCIIOpen(PETSC_COMM_WORLD,ss.str().c_str(),&viewer);
	    //PetscViewerBinaryOpen(PETSC_COMM_WORLD,ss.str().c_str(),FILE_MODE_WRITE,&viewer);
	    VecView(ft,viewer);
	    PetscViewerDestroy(&viewer);	    
	    
	    PetscInt TotalLz = 0;
	    if ( params.spin_one == PETSC_TRUE || truncate == false ) 
	      {	    
		PetscPrintf(PETSC_COMM_WORLD,"Writing text output files for spin one or non truncated case.\n");
				
		TotalLz = 2*(M - params.filling) - (params.AB_sites -1)*params.filling; 
		
// 		if ( params.size == 1 ) 
// 		  {
		ss.str(""); ss << "bosons_" << name << "_n_" << params.filling << "_2s_" << (params.AB_sites-1) << "_lz_" << TotalLz << ".0_M_" << M << ".vec.txt" ;
// 		  }
// 		else
// 		  { 
// 		    ss.str(""); ss << "bosons_" << name << "_n_" << params.filling << "_2s_" << (params.AB_sites-1) << "_lz_" << TotalLz << ".0_M_" << M << ".vec.txt_rank_" << params.rank  ;
// 		  }
		  		  
		for ( int rank = params.size-1 ; rank >= 0 ; rank-- ) 
		  {
		    MPI_Barrier(PETSC_COMM_WORLD);
		    if ( params.rank == rank ) 
		      {
			ofstream myfile;
			if ( rank == (params.size - 1) ) 
			  {
			    myfile.open(ss.str().c_str());
			  }
			else
			  {
			    myfile.open(ss.str().c_str(), ios::app | ios::out );
			  }
			myfile.precision(14);
			PetscScalar *ft_vals;
				
			PetscInt fstart = 0;
			if ( params.rank != 0 ) fstart = end_indices[params.rank-1];
				
			VecGetArray(ft,&ft_vals);
			for ( PetscInt i = end_indices[params.rank] - 1 ; i >= fstart ; i-- ) 
			  {	      
			    //full_index = get_full_basis_index(&params, params.AB_sites + params.filling, params.filling, i);		    		    
			    if ( use_complex == false )   
			      {
				if ( PetscAbsScalar(PetscImaginaryPart(ft_vals[i-fstart])) > 1e-11 )
				  {
				    PetscPrintf(PETSC_COMM_SELF, "Zeroing non zero complex part of vector %.14g + i%.14g\n",PetscRealPart(ft_vals[i-fstart]), PetscImaginaryPart(ft_vals[i-fstart]));
				  }
				myfile << scientific << PetscRealPart(ft_vals[i-fstart]) << endl;		    
			      }		    
			    else
			      {
				myfile << scientific << "(" << PetscRealPart(ft_vals[i-fstart]) << ","<< PetscImaginaryPart(ft_vals[i-fstart]) << ")" << endl;	
			      }
			    if ( verbose && params.spin_one != PETSC_TRUE ) 
			      {			
				PetscReal state_normal;
				get_monomial_bosons(&params, basis[i], params.AB_sites, params.filling, monomial, &state_normal);
				PetscPrintf(PETSC_COMM_SELF,"[");
				for ( int j = params.filling-1 ; j >= 0  ; j-- )
				  {
				    if (j != (params.filling-1) && (monomial[j]-1) > 0 ) PetscPrintf(PETSC_COMM_SELF,",");
				    if ( (monomial[j]-1) > 0 ) PetscPrintf(PETSC_COMM_SELF,"%d",monomial[j]-1);		  
				  }
				  PetscPrintf(PETSC_COMM_SELF,"] : %.14g + i%.14g\n",PetscRealPart(ft_vals[i-fstart]), PetscImaginaryPart(ft_vals[i-fstart]));
			      }
			  }
			VecRestoreArray(ft,&ft_vals);
			myfile.close();
		      }
		  }
	      }
	    else 
	      {
		PetscPrintf(PETSC_COMM_WORLD,"Writing text output files for spin 1/2 case.\n");
		TotalLz = ( M - params.filling*params.filling ) * 2; 
		
// 		if ( params.size == 1 ) 
// 		  {
		    ss.str(""); ss << "bosons_" << name << "_n_" << params.filling << "_2s_" << (params.AB_sites - 2) << "_lz_" << TotalLz << ".0_M_" << M << ".vec.txt" ;
// 		  }
// 		else
// 		  { 
// 		    ss.str(""); ss << "bosons_" << name << "_n_" << params.filling << "_2s_" << (params.AB_sites - 2) << "_lz_" << TotalLz << ".0_M_" << M << ".vec.txt_rank_" << params.rank  ;
// 		  }
		ofstream myfile;
		
		for ( int rank = params.size-1 ; rank >= 0 ; rank-- ) 
		  {
		    MPI_Barrier(PETSC_COMM_WORLD);
		    if ( params.rank == rank ) 
		      {
			ofstream myfile;
			if ( rank == (params.size - 1) ) 
			  {
			    myfile.open(ss.str().c_str());
			  }
			else
			  {
			    myfile.open(ss.str().c_str(), ios::app | ios::out );
			  }
		
			myfile.precision(14);
			PetscScalar *ft_vals;
				
			PetscInt fstart = 0;
			if ( params.rank != 0 ) fstart = end_indices[params.rank-1];
				
			VecGetArray(ft,&ft_vals);
			for ( PetscInt i = end_indices[params.rank] - 1 ; i >= fstart ; i-- ) 
			  {	      
			    //full_index = get_full_basis_index(&params, params.AB_sites + params.filling, params.filling, i);
			    PetscReal state_normal;
			    get_monomial_bosons(&params, basis[i], params.AB_sites, params.filling, monomial, &state_normal);
			    if ( monomial[params.filling-1] != params.AB_sites ) 
			      {			
				if ( use_complex == false )   
				  {
				    if ( PetscAbsScalar(PetscImaginaryPart(ft_vals[i-fstart])) > 1e-11 )
				      {
					PetscPrintf(PETSC_COMM_SELF, "Zeroing non zero complex part of vector %.14g + i%.14g\n",PetscRealPart(ft_vals[i-fstart]), PetscImaginaryPart(ft_vals[i-fstart]));
				      }
				    myfile << scientific << PetscRealPart(ft_vals[i-fstart]) << endl;		    
				  }		    
				else
				  {
				    myfile << scientific << "(" << PetscRealPart(ft_vals[i-fstart]) << ","<< PetscImaginaryPart(ft_vals[i-fstart]) << ")" << endl;	
				  }			
				//myfile << scientific << "(" << PetscRealPart(ft_vals[i-fstart]) << "," << PetscImaginaryPart(ft_vals[i-fstart]) << ")" << endl ;		    
			      }
			    else if ( PetscAbsScalar(ft_vals[i-fstart]) > 1e-11 ) 
			      {
				PetscPrintf(PETSC_COMM_SELF, "Zeroing non zero element with value %.14g + i%.14g\n",PetscRealPart(ft_vals[i-fstart]), PetscImaginaryPart(ft_vals[i-fstart]));
			      }
			    /*for ( int j = params.filling-1 ; j >= 0  ; j-- )
			      {
				if (j != (params.filling-1) &&  (monomial[j]-1) > 0 ) myfile << "," ;
				if ( (monomial[j]-1) > 0 ) myfile << monomial[j]-1;		  
			      }
			      myfile << "]" << endl;*/
			  }
			VecRestoreArray(ft,&ft_vals);
			myfile.close();
		      }
		  }
	      }
		
	    if ( calc_variance_ordered == true )
	      {		
		if ( params.size != 1 )
		  {
		    PetscPrintf(PETSC_COMM_WORLD,"Not supported in parallel yet.\n");		    
		  }
		else 
		  {
		    PetscInt local_dim = end_indices[params.rank]; 
		    if ( params.rank > 0 ) local_dim -= end_indices[params.rank-1];
		    //my_vec_element* basis_elements = new my_vec_element[local_dim];
		    my_vec_element* basis_elements;
		    ierr = PetscMalloc(local_dim * sizeof(my_vec_element), &basis_elements);CHKERRQ(ierr);
		    PetscScalar *vals;
		    VecGetArray(ft, &vals);			
		    uPetscInt variance;
		    int prev, count;
		    for ( int idx = 0; idx < local_dim ; idx++ ) 
		      {			    
			basis_elements[idx].basis_element = basis[idx];
			basis_elements[idx].val = vals[idx];
			get_monomial_bosons(&params, basis[idx], params.AB_sites, params.filling, monomial);
			variance = 0;
			prev = -1;
			count = 0;
			for ( int i = 0 ; i < params.filling ; i++ ) 
			  {							    
			    if ( monomial[i] != prev && prev != -1 )
			      {
				variance += (uPetscInt)count * (uPetscInt)(((prev)%params.AB_sites) - params.AB_sites/2)*(((prev)%params.AB_sites) - params.AB_sites/2);
				count = 0;
			      }
			    prev = monomial[i];		
			    count++;
			  }
			variance += (uPetscInt)count * (uPetscInt)(((prev)%params.AB_sites) - params.AB_sites/2)*(((prev)%params.AB_sites) - params.AB_sites/2);
			basis_elements[idx].sorting_token = variance; 
		      }		    
		    VecRestoreArray(ft, &vals);		
		    
		    sort(basis_elements, basis_elements+dim); 
		    //parallel_qsort(&(PETSC_COMM_WORLD), basis_elements, local_dim); 
		    
		    ss.str(""); ss << "bosons_" << name << "_n_" << params.filling << "_2s_" << (params.AB_sites-1) << "_lz_" << TotalLz << ".0_M_" << M << "_variance_sorted.vec"  ;			  
		    
		    for ( int rank = 0 ; rank < params.size ; rank++ ) 
		      {
			MPI_Barrier(PETSC_COMM_WORLD);
			if ( rank == params.rank ) 
			  {						
			    ofstream myfile;
			    if ( rank == 0 ) 
			      {
				myfile.open(ss.str().c_str(), ios::out);
			      }
			    else
			      {
				myfile.open(ss.str().c_str(), ios::out | ios::app );
			      }
			    myfile.precision(14);
			    PetscScalar *ft_vals;
			    int occupation[params.AB_sites];
			    for ( int idx = 0; idx < local_dim ; idx++ ) 
			      {		
				if ( PetscAbsScalar(basis_elements[idx].val) > 1e-10) 
				  {
				  for ( int i = 0 ; i < params.AB_sites ; i++ ) occupation[i] = 0;
				  get_monomial_bosons(&params, basis_elements[idx].basis_element, params.AB_sites, params.filling, monomial);						
				  int prev = -1;
				  int count = 0;
				  for ( int i = 0 ; i < params.filling ; i++ ) 
				    {							    
				      if ( monomial[i] != prev && prev != -1 )
					{
					  occupation[prev-1] = count;				
					  count = 0;
					}
				      prev = monomial[i];	
				      count++;
				    }
				  occupation[prev-1] = count;      
				  for ( int i = 0 ; i < params.AB_sites ; i++ ) 
				    {
				      //PetscPrintf(PETSC_COMM_SELF,"%d",monomial[i]);
				      myfile << occupation[i];
				    }
				  //PetscPrintf(PETSC_COMM_SELF,"\t");
				  //PetscPrintf(PETSC_COMM_SELF,"(%.14g,%.14g)",PetscRealPart(basis_elements[idx].val),PetscImaginaryPart(basis_elements[idx].val));
				  //PetscPrintf(PETSC_COMM_SELF,"\t%d",(int)basis_elements[idx].sorting_token);
				  myfile << "\t\t" << basis_elements[idx].val << "\t\t" << basis_elements[idx].sorting_token << endl;			 
				}
			      }
			    myfile.close();	
			  }
		      }
		    PetscFree(basis_elements);
		  }
	      }
	  }				
		
	PetscFree(basis);	
	PetscFree(end_indices);
	VecDestroy(&ft);	
      }	//end for over momentum	   
    PetscPrintf(PETSC_COMM_WORLD, "Total FT magnitude: %lf\n", FTMagnitude);
    
    PetscFree(spin_basis);
    
    ierr = SlepcFinalize();CHKERRQ(ierr);
    return 0;
}
