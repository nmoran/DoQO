#include "slepceps.h"
#include "slepcsvd.h"
#include "tinyxml.h"

#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h> 

#if defined(PETSC_USE_64BIT_INDICES)
typedef 	unsigned long long uPetscInt;
#define 	PETSC_MPI_INT MPI_LONG_LONG_INT
#define 	PETSC_UNSIGNED_MPI_INT MPI_UNSIGNED_LONG_LONG
#else
typedef 	unsigned int uPetscInt;
#define 	PETSC_MPI_INT MPI_INT
#define 	PETSC_UNSIGNED_MPI_INT MPI_UNSIGNED
#endif


#define MAX_STRING 200
#define TOLERANCE  1e-100
#define NSV_LIMIT	   200


using namespace std;

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

struct adjacency_information {
	int 		number_masks;
	uPetscInt 	*masks;
};

//info necessary for generating all configs in a list
struct conf_list_info {
	struct 		adjacency_information adj_info;
	PetscInt 	num_sites, filling, max_filling;
};

struct parameters {
	PetscMPIInt rank, size;
	uPetscInt *AB_confs,*A_confs,*B_confs;
	PetscInt sites, filling, A_filling, B_filling, l_A_confs, g_A_confs, l_B_confs, g_B_confs, l_AB_confs, g_AB_confs;
	PetscInt AB_sites, A_sites, B_sites; 
	uPetscInt A_mask, B_mask, AB_mask; 
	struct adjacency_information adj_info;
	Vec  psi; 
	Mat  W;
	PetscTruth nn_exclusion,ascii;
	class my_vec_entry *my_vec_entries;
};

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
	PetscErrorCode ierr;
	
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

int find_valid_configs_with_this_subset_list(struct parameters *params, struct conf_list_info *c_info, int set_size, int set_filling, uPetscInt set , list<uPetscInt> *ll) {
	//first check the stop condition if set size is already the size of the full graph or there is not enough space left to reach required filling.
	if ( set_size == c_info->num_sites ) {
		if ( c_info->filling == -1 || set_filling == c_info->filling ) { //if there is a filling limit and we are at it just add the configuration.
			ll->push_back(set);
			return 0;
		} else {
			return 0;
		}
	} else if ( c_info->filling != -1 ) {
		if ((c_info->filling - set_filling) > (c_info->num_sites - set_size) )  { //if filling limit and will not be able to make it just return.
			return 0 ;
		} else if ( set_filling == c_info->filling ) {
			ll->push_back(set);
			return 0;
		}
	} else if ( c_info->filling == -1 ) { //we check that everything is alright with respect to max filling.
		if ( set_filling == c_info->max_filling ) {
			ll->push_back(set);
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
			set2 = set;
			basis_file->write((char*)&set2,sizeof(uPetscInt));	
			//cout << "Found " << set << endl;
			(*local_reps)++;
			return 0;
		} else {
			return 0;
		}
	} else if ( c_info->filling != -1 ) {
		if ((c_info->filling - set_filling) > (c_info->num_sites - set_size)) {
			return 0 ;
		} else if ( set_filling == c_info->filling ) {
			set2 = set;
			basis_file->write((char*)&set2,sizeof(uPetscInt));	
			//cout << "Found " << set << endl;
			(*local_reps)++;
			return 0;
		}
	} else if ( c_info->filling == -1 ) {
		if ( set_filling == c_info->max_filling ) {
			//ll->push_back(set);
			set2 = set;
			basis_file->write((char*)&set2,sizeof(uPetscInt));	
			//cout << "Found " << set << endl;
			(*local_reps)++;
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
}

int find_valid_configs_list(struct parameters *params, struct conf_list_info *c_info, list<uPetscInt> *ll){
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
	
	find_valid_configs_with_this_subset_list(params,c_info, set_size,(int)(*current_fillings)[0],(uPetscInt)(*current_configs)[0],ll);	
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
			PetscInt start_of_upper = llength - (total_above - l); //the local index of the start of the upper index.
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
				uPetscInt *recv_ptr = &(array[end_of_lower]);
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
	c_info.filling = filling;		  //filling to use or -1 if filling not defined on the subsection.
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

int main( int argc, char **argv ) {
	PetscErrorCode ierr; 
	string adj_file, vec_list_file, a_sites_file;
	list<string> vec_list, a_sites_list;
	PetscInt A_filling_from = -1, A_filling_to = -1; 
	struct parameters params; 
	char buf[MAX_STRING];
	PetscTruth flg;
	
	/*------------------------------------------------------------------------------------------------
	Initialise some values.
	------------------------------------------------------------------------------------------------*/
	char help[] = "Program to generate calculate entanglement entropy and spectrum in parallel.\n";
	params.A_filling = -1;	params.B_filling = -1; params.filling = -1; //by default any fillings can be used.
	params.nn_exclusion = PETSC_FALSE;
	
	SlepcInitialize(&argc,&argv,(char*)0,help);	
	
	//get rank of this process and number of processes 
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&(params.rank));CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&(params.size));CHKERRQ(ierr);
	
	/*------------------------------------------------------------------------------------------------
	Read in input details
	------------------------------------------------------------------------------------------------*/
	ierr = PetscOptionsGetString(NULL,"-adj_file",buf,MAX_STRING,&flg);CHKERRQ(ierr);
	if ( flg == PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD,"No adjacency file provided so assuming no NN exclusion condition.\n");
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"Adj file: %s\n",buf);
		params.nn_exclusion = PETSC_TRUE;
		adj_file.assign(buf);
	}
	
	ierr = PetscOptionsGetInt(NULL,"-num_sites",&params.AB_sites,&flg);CHKERRQ(ierr);
	if ( flg == PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD,"Number of sites not specified.\n");
		params.AB_sites = -1;
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"Number of sites: %d\n",params.AB_sites);
	}
	
	ierr = PetscOptionsGetInt(NULL,"-filling",&params.filling,&flg);CHKERRQ(ierr);
	if ( flg == PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD,"No filling specified.\n");
		params.filling = -1;
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"Filling specified: %d\n",params.filling);
	}
	
	ierr = PetscOptionsGetInt(NULL,"-A_filling_from",&A_filling_from,&flg);CHKERRQ(ierr);
	if ( flg == PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD,"No A filling from specified.\n");
//		params.A_filling = -1;
//		params.B_filling = -1;
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"A filling from specified: %d\n",params.A_filling);
/*		if ( params.filling != -1 ) {
			params.B_filling = params.filling - params.A_filling ;
		}*/
	}
	
	ierr = PetscOptionsGetInt(NULL,"-A_filling_to",&A_filling_to,&flg);CHKERRQ(ierr);
	if ( flg == PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD,"No A filling to specified.\n");
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"A filling to specified: %d\n",params.A_filling);
	}
	
	
	ierr = PetscOptionsGetString(NULL,"-vec_list",buf,MAX_STRING,&flg);CHKERRQ(ierr);
	if ( flg == PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD,"Must provide vector on the command line with switch -vec list.\n");
		return 1;
	}
	PetscPrintf(PETSC_COMM_WORLD,"Vec list file: %s\n",buf);
	vec_list_file.assign(buf);
	
	ierr = PetscOptionsGetString(NULL,"-ascii",buf,MAX_STRING,&flg);CHKERRQ(ierr);
	if ( flg == PETSC_FALSE) {
		params.ascii = PETSC_FALSE;
	} else {
		if ( strcmp(buf,"true") == 0 ) {
			params.ascii = PETSC_FALSE;
		} else {
			params.ascii = PETSC_FALSE;
		}
	}
	
	ierr = PetscOptionsGetString(NULL,"-A_sites",buf,MAX_STRING,&flg);CHKERRQ(ierr);
	if ( flg == PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD,"Must provide sites of A on command line with switch -A_sites.\n");
		return 1;
	}
	PetscPrintf(PETSC_COMM_WORLD,"A sites file: %s\n",buf);
	a_sites_file.assign(buf);
	/*------------------------------------------------------------------------------------------------
	Finished reading input so now go to work.
	------------------------------------------------------------------------------------------------*/
	if ( params.nn_exclusion ) {
		parse_adjacency_info_file(adj_file, &(params.adj_info),&params.AB_sites); //parse the adjacency file.
	}
	PetscPrintf(PETSC_COMM_WORLD,"Number of sites: %d\n",params.AB_sites);
	
	//read in vectors and a sites files.
	{
		struct stat stFileInfo; 
		int intStat; 
	
		ifstream vecs(vec_list_file.c_str());
		string line;
		if (vecs.is_open() ) {
			while (vecs.good()) {
				getline(vecs,line);
				if ( line != "" ) {
					intStat  = stat(line.c_str(),&stFileInfo); 
					if ( intStat == 0 ) {
						vec_list.push_back(line);
					} else {
						PetscFPrintf(PETSC_COMM_WORLD,stderr,"File does not exist %d.\n",line.c_str());
					}
				}
			}
		}
		vecs.close();
	}
	
	{
		struct stat stFileInfo; 
		int intStat; 
		ifstream a_sites(a_sites_file.c_str());
		string line;
		if (a_sites.is_open() ) {
			while (a_sites.good() ) {
				getline(a_sites,line);
				if ( line != "" ) {
					intStat  = stat(line.c_str(),&stFileInfo); 
					if ( intStat == 0 ) {
						a_sites_list.push_back(line);
					} else {
						PetscFPrintf(PETSC_COMM_WORLD,stderr,"File does not exist %d.\n",line.c_str());
					}
				}
			}
		}
		a_sites.close();
	}
	
	parse_a_sites_file(*a_sites_list.begin(),&params);
	PetscInt max_filling = params.filling;
	build_config_list(&params, &params.AB_confs, PETSC_TRUE, &params.l_AB_confs, &params.g_AB_confs, params.AB_mask, params.filling,max_filling);
	PetscPrintf(PETSC_COMM_WORLD,"Num AB configs: local: %d global: %d\n",params.l_AB_confs, params.g_AB_confs);
	/*for ( int i = 0 ; i < params.l_AB_confs ; i++ ) {
	 dec2bin(params.AB_confs[i], params.AB_sites, buf );
	 PetscPrintf(PETSC_COMM_SELF,"AB: %s\n",buf);
	 }*/
	
	for ( list<string>::iterator iter = a_sites_list.begin() ; iter != a_sites_list.end() ; iter++ ) {
		parse_a_sites_file(*iter,&params);
		dec2bin(params.A_mask, params.AB_sites, buf );
		PetscPrintf(PETSC_COMM_WORLD,"A mask: %s\n",buf);
		dec2bin(params.B_mask, params.AB_sites, buf );
		PetscPrintf(PETSC_COMM_WORLD,"B mask: %s\n",buf);
		dec2bin(params.AB_mask, params.AB_sites, buf );
		PetscPrintf(PETSC_COMM_WORLD,"AB mask: %s\n",buf);
	
		
		for ( params.A_filling = A_filling_from ; params.A_filling <= A_filling_to ; params.A_filling++ ) {		
			max_filling = params.A_filling;
			if ( params.A_filling == -1 ) {
				max_filling = params.filling;
				params.B_filling = -1;
			} else {
				params.B_filling = params.filling - params.A_filling;
				//If fillings are too ridiculous then just skip them.
				if ( (params.nn_exclusion && params.A_filling > (float(params.A_sites+1.0)/2.0) ) || params.A_filling > params.A_sites ) {
					continue;
				}
				if ( (params.nn_exclusion && params.B_filling > (float(params.B_sites+1.0)/2.0) ) || params.B_filling > params.B_sites ) {
					continue;
				}
			}
			
			build_config_list(&params, &params.A_confs, PETSC_FALSE, &params.l_A_confs, &params.g_A_confs, params.A_mask, params.A_filling,max_filling);
			if ( params.B_filling == -1 ) {
				max_filling = params.filling;
			} else max_filling = params.B_filling;
			if ( params.A_confs > 0 ) {
				build_config_list(&params, &params.B_confs, PETSC_TRUE, &params.l_B_confs, &params.g_B_confs, params.B_mask, params.B_filling,max_filling);
			}
			
			PetscPrintf(PETSC_COMM_WORLD,"Num A configs: %d\n",params.g_A_confs);
			/*for ( int i = 0 ; i < params.g_A_confs ; i++ ) {
				dec2bin(params.A_confs[i], params.AB_sites, buf );
				PetscPrintf(PETSC_COMM_WORLD,"A: %s\n",buf);
			}*/
			PetscPrintf(PETSC_COMM_WORLD,"Num B configs: local: %d global: %d\n",params.l_B_confs, params.g_B_confs);
			/*for ( int i = 0 ; i < params.l_B_confs ; i++ ) {
				dec2bin(params.B_confs[i], params.AB_sites, buf );
				PetscPrintf(PETSC_COMM_SELF,"B: %s\n",buf);
			}*/
			
		
			if ( params.g_A_confs > 0 && params.g_B_confs > 0 ) {
				//PetscPrintf(PETSC_COMM_SELF,"Rank: %d, m: %d, n: %d, M: %d, N: %d\n",params.rank,params.l_B_confs, params.g_A_confs, params.g_B_confs ,params.g_A_confs);
				//create matrix and zero entries
				ierr = MatCreate(PETSC_COMM_WORLD,&params.W);CHKERRQ(ierr);
				ierr = MatSetSizes(params.W,params.l_B_confs, PETSC_DECIDE, params.g_B_confs ,params.g_A_confs);CHKERRQ(ierr);
				ierr = MatSetType(params.W, MATMPIDENSE);CHKERRQ(ierr);
				PetscPrintf(PETSC_COMM_WORLD,"Matrix dimensions %d x %d.\n",(int)params.g_B_confs, (int)params.g_A_confs); 
				
				//find the limits of the B configs on each process.
				uPetscInt *B_conf_limits;
				ierr = PetscMalloc((params.size)*sizeof(uPetscInt),&B_conf_limits);CHKERRQ(ierr);
				if ( params.l_B_confs > 0 ) {
					MPI_Allgather(&params.B_confs[params.l_B_confs-1],1,PETSC_UNSIGNED_MPI_INT,B_conf_limits,1,PETSC_UNSIGNED_MPI_INT,PETSC_COMM_WORLD);
				} else {
					MPI_Allgather(&params.B_mask,1,PETSC_UNSIGNED_MPI_INT,B_conf_limits,1,PETSC_UNSIGNED_MPI_INT,PETSC_COMM_WORLD);
				}
				
				for ( list<string>::iterator vec_iter = vec_list.begin() ; vec_iter != vec_list.end() ; vec_iter++) {	
					//load in the vector
					PetscFPrintf(PETSC_COMM_WORLD,stderr,"Vector file <%s>.\n",vec_iter->c_str());
					PetscViewer viewer;
					PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
					if ( params.ascii ) {
						ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, vec_iter->c_str() ,&viewer);CHKERRQ(ierr); //ascii not working well for large sizes at the moment.
					} else {
						ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, vec_iter->c_str() ,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
					}
					ierr = VecLoad(viewer,VECMPI,&params.psi);CHKERRQ(ierr);
					PetscViewerDestroy(viewer);
					
					//do check to make sure the length of the vector matches the number of AB configs.
					PetscInt vec_size;
					VecGetSize(params.psi,&vec_size);
					if ( vec_size != params.g_AB_confs) {
						PetscPrintf(PETSC_COMM_WORLD,"Mismatch in vec length and number of configs. Vec: %d and confs: %d.\n",(int)vec_size, (int)params.g_AB_confs);
						return 1;
					}
					
					MatZeroEntries(params.W);					

					//now create array of data structures which contain both vec value and conf value.
					PetscScalar *vec_elements;
					VecGetArray(params.psi,&vec_elements);
					PetscInt non_zero_count = 0;
					for ( PetscInt idx = 0 ; idx < params.l_AB_confs ; idx++ ) {
						if ( PetscAbsScalar(vec_elements[idx]) > TOLERANCE ) {
							non_zero_count++;
						}
					}
					ierr = PetscMalloc(non_zero_count*sizeof(class my_vec_entry),&params.my_vec_entries);CHKERRQ(ierr);
					PetscInt tmp_idx = 0;
					for ( PetscInt idx = 0 ; idx < params.l_AB_confs ; idx++ ) {
						if ( PetscAbsScalar(vec_elements[idx]) > TOLERANCE ) {	
							params.my_vec_entries[tmp_idx].conf = params.AB_confs[idx];
							params.my_vec_entries[tmp_idx].val = vec_elements[idx];
							params.my_vec_entries[tmp_idx++].mask = params.B_mask;
						}
					}
					VecRestoreArray(params.psi,&vec_elements);
					sort(params.my_vec_entries, params.my_vec_entries + non_zero_count); //sort according to the b mask.

					
					//now we start the peer to peer communication. 
					
					uPetscInt lower_limit, upper_limit;
					PetscInt chunk_start, chunk_size, row, col, recv_chunk_size; 
					int dest_proc, source_proc;
					dest_proc = params.rank;
					if ( dest_proc > 0 ) {
						lower_limit = B_conf_limits[dest_proc-1];
					} else lower_limit = 0;
					upper_limit = B_conf_limits[dest_proc];
					chunk_start = find_start_of_range(params.my_vec_entries,non_zero_count,lower_limit,upper_limit,params.B_mask,&chunk_size); 
					for ( PetscInt idx = chunk_start ; idx < chunk_start + chunk_size ; idx++ ) {
						row = find_in_sorted_array(params.B_confs,params.l_B_confs, params.my_vec_entries[idx].conf & params.B_mask);
						col = find_in_sorted_array(params.A_confs,params.l_A_confs, params.my_vec_entries[idx].conf & params.A_mask);
						if ( row >= 0 && col >= 0 ) MatSetValue(params.W,row+get_rank_chunk_start(params.g_B_confs,params.rank,params.size),col,params.my_vec_entries[idx].val,ADD_VALUES);
					}
					
					MPI_Request send_req;
					MPI_Status send_stat;
					class my_vec_entry *recv_buf;
					int proc_offset = 1; 
					while ( proc_offset < params.size ) { //loop over an iterator so that each process communicates with each other process. 
						//find out the off set in the local my_vec_entries array of the proc with id (rank + offset)
						dest_proc = (params.rank + proc_offset) % params.size;
						source_proc = (params.rank + params.size - proc_offset) % params.size;
						if ( dest_proc > 0 ) {
							lower_limit = B_conf_limits[dest_proc-1];
						} else lower_limit = 0;
						upper_limit = B_conf_limits[dest_proc];
						chunk_start = find_start_of_range(params.my_vec_entries,non_zero_count,lower_limit,upper_limit,params.B_mask,&chunk_size); 
						MPI_Isend(&chunk_size,1,PETSC_MPI_INT,dest_proc,1,PETSC_COMM_WORLD,&send_req);//send using non blocking communication.
						MPI_Recv(&recv_chunk_size,1,PETSC_MPI_INT,source_proc,1,PETSC_COMM_WORLD,&send_stat); //receive the size of the chunk to receive
						MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
						if ( chunk_size > 0 ) {
							MPI_Isend((char*)&params.my_vec_entries[chunk_start],chunk_size*sizeof(class my_vec_entry),MPI_BYTE,dest_proc,1,PETSC_COMM_WORLD,&send_req);//send using non blocking communication.
						}
						if ( recv_chunk_size > 0 ) {
							ierr = PetscMalloc(recv_chunk_size*sizeof(class my_vec_entry),&recv_buf);CHKERRQ(ierr);			
							MPI_Recv((char*)recv_buf,recv_chunk_size*sizeof(class my_vec_entry),MPI_BYTE,source_proc,1,PETSC_COMM_WORLD,&send_stat); //receive the size of the chunk to receive			
							for ( PetscInt idx = 0 ; idx < recv_chunk_size ; idx++ ) {
								row = find_in_sorted_array(params.B_confs,params.l_B_confs, recv_buf[idx].conf & params.B_mask);
								col = find_in_sorted_array(params.A_confs,params.l_A_confs, recv_buf[idx].conf & params.A_mask);
								if ( row >= 0 && col >= 0 ) MatSetValue(params.W,row+get_rank_chunk_start(params.g_B_confs,params.rank,params.size),col,recv_buf[idx].val,ADD_VALUES);
							}
							PetscFree(recv_buf);
						}
						if ( chunk_size > 0 ) {
							MPI_Wait(&send_req,&send_stat); //wait for non blocking send to complete.
						}
						proc_offset++;
					}
					
					MatAssemblyBegin(params.W,MAT_FINAL_ASSEMBLY);
					MatAssemblyEnd(params.W,MAT_FINAL_ASSEMBLY);
					
					
					PetscPrintf(PETSC_COMM_WORLD,"Matrix constructed.\n"); 
					
					stringstream ss;
					
					/*ss.str("");
					ss << vec_file << ".W2.m" ;
					PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
					PetscViewerASCIIOpen(PETSC_COMM_WORLD,ss.str().c_str(),&viewer);
					PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
					MatView(params.W,viewer);
					PetscViewerDestroy(viewer);*/
					
					
					SVD svd;
					ierr = SVDCreate(PETSC_COMM_WORLD,&svd);CHKERRQ(ierr);
					ierr = SVDSetOperator(svd,params.W);CHKERRQ(ierr);
					
					PetscInt min_dimension = min(params.g_A_confs,params.g_B_confs);
					if ( min_dimension <  NSV_LIMIT ) {
						ierr = SVDSetDimensions(svd,min_dimension,min_dimension,min_dimension);CHKERRQ(ierr);
					} else {
						ierr = SVDSetDimensions(svd,NSV_LIMIT,((NSV_LIMIT*3/2)>min_dimension)?min_dimension:(NSV_LIMIT*3/2),NSV_LIMIT);CHKERRQ(ierr);
					}
					
					/*
					 Use thick-restart Lanczos as default solver
					 */
					ierr = SVDSetType(svd,SVDTRLANCZOS);CHKERRQ(ierr);
					
					/*
					 Set solver parameters at runtime
					 */
					ierr = SVDSetFromOptions(svd);CHKERRQ(ierr);
					PetscPrintf(PETSC_COMM_WORLD,"Starting solver.\n");
					ierr = SVDSolve(svd);CHKERRQ(ierr);
					PetscInt its, nconv;
					ierr = SVDGetIterationNumber(svd, &its);CHKERRQ(ierr);
					ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of iterations of the method: %d\n",its);CHKERRQ(ierr);

					ierr = SVDGetConverged(svd,&nconv);CHKERRQ(ierr);
					ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of converged approximate singular triplets: %d\n\n",nconv);CHKERRQ(ierr);
					
					PetscReal tol;
					PetscInt maxits;
					ierr = SVDGetTolerances(svd,&tol,&maxits);CHKERRQ(ierr);
					ierr = PetscPrintf(PETSC_COMM_WORLD,"Tolerance: %g and maximum number of iterations %d\n\n",tol,maxits);CHKERRQ(ierr);
					
					
					FILE *fp,*fp_ee;
					if ( params.rank == 0 ) {
						stringstream ss;
						ss << *vec_iter ;
						if ( params.filling != -1 ) ss << "_f" << params.filling ;
						if ( params.A_filling != -1 ) ss << "_Af" << params.A_filling ;
						ss << "_As" << params.A_sites << ".es" ;
						string output_filename = ss.str();
						fp = fopen(output_filename.c_str(),"w");
						ss.str("");ss << *vec_iter ;
						if ( params.filling != -1 ) ss << "_f" << params.filling ;
						if ( params.A_filling != -1 ) ss << "_Af" << params.A_filling ;
						ss << "_As" << params.A_sites << ".ee" ;
						output_filename = ss.str();
						fp_ee = fopen(output_filename.c_str(),"w");
					}
					
					PetscReal sigma, error, ent_ent = 0.0,xi;
					if (nconv>0) {
						
						Vec    u,v; 
						ierr = MatGetVecs(params.W,&v,&u);CHKERRQ(ierr);
						/*
						 Display singular values and relative errors
						 */
						ierr = PetscPrintf(PETSC_COMM_WORLD,
										   "Entanglement spectrum and errors\n"
										   "  ---------------------------------------\n" );CHKERRQ(ierr);
						for( PetscInt i=0; i<nconv; i++ ) {
							/* 
							 Get converged singular triplets: i-th singular value is stored in sigma
							 */
							ierr = SVDGetSingularTriplet(svd,i,&sigma,u,v);CHKERRQ(ierr);
							
							/*
							 Compute the error associated to each singular triplet 
							 */
							ierr = SVDComputeRelativeError(svd,i,&error);CHKERRQ(ierr);
							
							xi = -PetscLogScalar(sigma)*2.0;
							ent_ent += -sigma*sigma * PetscLogScalar(sigma*sigma);
							
							//ierr = PetscPrintf(PETSC_COMM_WORLD,"       % 6f      ",sigma); CHKERRQ(ierr);
							//ierr = PetscPrintf(PETSC_COMM_WORLD,"       % 6f      ",-PetscLogScalar(sigma)*2.0); CHKERRQ(ierr);
							//ierr = PetscPrintf(PETSC_COMM_WORLD," % 12g\n",error);CHKERRQ(ierr);
							ierr = PetscPrintf(PETSC_COMM_WORLD,"%12lg\t%12lg\t%12lg\n",sigma,-PetscLogScalar(sigma)*2.0,error); CHKERRQ(ierr);
							if ( params.rank == 0 ) {
								ierr = PetscFPrintf(PETSC_COMM_WORLD,fp,"%12lg\t%12lg\t%12lg\n",sigma,-PetscLogScalar(sigma)*2.0,error); CHKERRQ(ierr);
							}
						}
						ierr = PetscPrintf(PETSC_COMM_WORLD,"\n" );CHKERRQ(ierr);
						ierr = PetscPrintf(PETSC_COMM_WORLD,"Entanglement entropy: %12lg\n",ent_ent);CHKERRQ(ierr);
						ierr = PetscFPrintf(PETSC_COMM_WORLD,fp_ee,"Entanglement entropy: %12lg\n",ent_ent);CHKERRQ(ierr);
						ierr = VecDestroy(u);CHKERRQ(ierr);
						ierr = VecDestroy(v);CHKERRQ(ierr);
					}
					
					if ( params.rank == 0 ) {
						fclose(fp);
						fclose(fp_ee);
					}
					
					ierr = SVDDestroy(svd);CHKERRQ(ierr);
				} //end loop over vectors
				ierr = MatDestroy(params.W);CHKERRQ(ierr);
				PetscFree(B_conf_limits);	
				PetscFree(params.my_vec_entries);
			} else { //end check dimensions
				PetscPrintf(PETSC_COMM_WORLD,"One of the dimensions is 0.\n");
			}
			if ( params.l_A_confs >  0 ) PetscFree(params.A_confs);
			if ( params.l_B_confs >  0 ) PetscFree(params.B_confs);
		} //end loop over A fillings
	} //end loop over a sites files.
				
				
	if ( params.l_AB_confs > 0 ) PetscFree(params.AB_confs);
	
	ierr = SlepcFinalize();CHKERRQ(ierr);
	return 0;
}