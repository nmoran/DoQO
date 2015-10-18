#include "slepceps.h"
#include "tinyxml.h"

#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <algorithm>


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


typedef PetscBool PetscTruth;

using namespace std;

class my_uPetscInt {
	public: 
		uPetscInt num;
		uPetscInt mask;
	
	 bool operator==(const my_uPetscInt &other) const {
    	return (num & mask) == (other.num & other.mask);  
	 }
	 
	 bool operator<(const my_uPetscInt &other) const {
    	return (num & mask) < (other.num & other.mask);  
	 }
	 
	 bool operator>(const my_uPetscInt &other) const {
    	return (num & mask) > (other.num & other.mask);  
	 }
};



struct adjacency_information {
	int 		number_masks;
	uPetscInt 	*masks;
};

struct parameters {
	PetscMPIInt rank, size;
	//vector<uPetscInt> ab_configs_vec;
	list<uPetscInt> ab_configs_vec;
	uPetscInt *ab_configs_array;
	PetscInt sites, filling;
	struct adjacency_information adj_info;
	Vec  psi; 
	Mat  rho_a;
	uPetscInt A_mask, B_mask; 
	PetscInt ab_configs, a_configs, b_configs, ab_sites, a_sites, b_sites;
	PetscInt a_filling; 
	uPetscInt *A_configs,*B_configs; 
};


int parse_adjacency_info_file(string filename, struct adjacency_information *adj_info,PetscInt *sites){
	//string filename(adjacency_file);
	PetscErrorCode ierr;
	
	*sites = 0 ;
	
	TiXmlDocument doc(filename);
	
	bool loadOkay = doc.LoadFile();
	
	if ( !loadOkay )
	{
		PetscPrintf(PETSC_COMM_WORLD,"Could not load adjacency file '%s'. Error='%s'. Exiting.\n", filename.c_str(), doc.ErrorDesc() );
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
			PetscPrintf(PETSC_COMM_WORLD,"Invalid number of edges %d.\n",num_sites);
			exit(1);
		}
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
		for ( int i = 0 ; i < params->ab_sites ; i++ ) {
			params->B_mask += (uPetscInt)pow(2.0,double(i));
		}
		params->B_mask ^= params->A_mask; 
	}
	
	return 0;
}


//this function checks the adjacency information for the given index and returns false if it contains adjacent particles.
PetscTruth check_adjacency(struct parameters *params, uPetscInt index) {
	for ( int i = 0 ; i < params->adj_info.number_masks; i++ ) {
		if ( (index & params->adj_info.masks[i] ) == params->adj_info.masks[i] ) {
			return PETSC_FALSE;
		}
	}
	return PETSC_TRUE;
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

int find_valid_configs_with_this_subset(struct parameters *params, int set_size, int set_filling, uPetscInt set ) {
	//first check the stop condition if set size is already the size of the full graph or there is not enough space left to reach required filling.
	if ( set_size == params->sites ) {
		if (set_filling == params->filling ) {
			params->ab_configs_vec.push_back(set);
			return 0;
		} else {
			return 0;
		}
	} else if ( (params->filling - set_filling) > (params->sites - set_size) ) {
		return 0 ;
	} else if ( set_filling == params->filling ) {
		params->ab_configs_vec.push_back(set);
		return 0;
	}
	
	//have two branches one with 0 at site set_size and one with 1 at set_size
	//now add a 1 to position set_size, check if valid and increment filling.
	uPetscInt new_set = set + (((uPetscInt)1) << set_size) ; 
	if ( check_adjacency(params,new_set) ) {
		find_valid_configs_with_this_subset(params, set_size+1 , set_filling+1, new_set);
	}
	
	
	//do not need to check the addition of a 0 as will already be full.
	find_valid_configs_with_this_subset(params, set_size+1 , set_filling, set);
}

int find_valid_configs(struct parameters *params){
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
	
	//while there are less valid configs that processes 
	while ( count < params->size && set_size < params->sites ) {
		next_configs = new vector<uPetscInt>(); 
		next_fillings = new vector<int>();
		int i = 0 ;
		int new_count = 0; 
		while ( i < count ) {		
			if ( (params->filling - (*current_fillings)[i]) <= (params->sites - set_size) ) {
				next_configs->push_back((*current_configs)[i]);
				next_fillings->push_back((*current_fillings)[i]);
				new_count ++;
			}
				
			if ( (*current_fillings)[i] <= params->filling ) {
				uPetscInt tmp;
				tmp = (*current_configs)[i] + (((uPetscInt)1) << set_size) ;
				if ( check_adjacency(params,tmp)) {
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
	
	if ( count >= params->size && set_size < params->sites ) {
		cout << "Count " << count << " rank " << params->rank << " size " << get_rank_chunk_size(count,params->rank,params->size) << " and start " <<  get_rank_chunk_start(count,params->rank,params->size)  << endl;
		for ( PetscInt idx = 0; idx <  get_rank_chunk_size(count, params->rank,params->size) ; idx++ ) {
			cout << "Rank " << params->rank << " getting started with " << (*current_configs)[get_rank_chunk_start(count, params->rank,params->size)+idx] << endl;
			find_valid_configs_with_this_subset(params,set_size,(*current_fillings)[get_rank_chunk_start(count, params->rank,params->size)+idx],  (*current_configs)[get_rank_chunk_start(count, params->rank,params->size)+idx]);	
		}		
	} 
	
	count = 0;
	for ( list<uPetscInt>::iterator iter = params->ab_configs_vec.begin();  iter != params->ab_configs_vec.end() ; iter++ ) {
		//dec2bin(*iter, params.sites, buf );
		//cout << buf << endl ; 
		count++;
	}
	cout << "Rank : " << params->rank << " number : " << count << endl;
	
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



int main( int argc, char **argv ) {
	PetscErrorCode ierr; 
	string adj_file, vec_file, a_sites;
	struct parameters params; 
	char buf[MAX_STRING];
	PetscTruth flg;
	PetscTruth calc_weight_matrix = PETSC_FALSE;
	
	char help[] = "Program to generate valid configs.\n";
	SlepcInitialize(&argc,&argv,(char*)0,help);
	//get rank of this process and number of processes 
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&(params.rank));CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&(params.size));CHKERRQ(ierr);
	
	/*------------------------------------------------------------------------------------------------
	Read in input details
	------------------------------------------------------------------------------------------------*/
	ierr = PetscOptionsGetString(NULL,"-adj_file",buf,MAX_STRING,&flg);CHKERRQ(ierr);
	if ( flg == PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD,"Must provide an adjacency file on the command line with switch -adj_file.\n");
		return 1;
	}
	PetscPrintf(PETSC_COMM_WORLD,"Adj file: %s\n",buf);
	adj_file.assign(buf);
	
	ierr = PetscOptionsGetInt(NULL,"-filling",&params.filling,&flg);CHKERRQ(ierr);
	if ( flg == PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD,"Must provide filling on the command line with switch -filling.\n");
		return 1;
	}
	PetscPrintf(PETSC_COMM_WORLD,"Filling: %d\n",params.filling);
	
	ierr = PetscOptionsGetString(NULL,"-vec",buf,MAX_STRING,&flg);CHKERRQ(ierr);
	if ( flg == PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD,"Must provide vector on the command line with switch -vec.\n");
		return 1;
	}
	PetscPrintf(PETSC_COMM_WORLD,"Vec file: %s\n",buf);
	vec_file.assign(buf);
	
	
	ierr = PetscOptionsGetString(NULL,"-A_sites",buf,MAX_STRING,&flg);CHKERRQ(ierr);
	if ( flg == PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD,"Must provide sites of A on command line with switch -A_sites.\n");
		return 1;
	}
	PetscPrintf(PETSC_COMM_WORLD,"A sites file: %s\n",buf);
	a_sites.assign(buf);
	
	ierr = PetscOptionsGetString(NULL,"-w_matrix",buf,MAX_STRING,&flg);CHKERRQ(ierr);
	if ( flg ) {
		if (strcmp(buf,"true") == 0 ) {
			calc_weight_matrix = PETSC_TRUE;
		}else {
			calc_weight_matrix = PETSC_FALSE;
		}
	}
	/*------------------------------------------------------------------------------------------------
	Finished reading input so now go to work.
	------------------------------------------------------------------------------------------------*/
	parse_adjacency_info_file(adj_file, &(params.adj_info),&params.ab_sites); //parse the adjacency file.
	PetscPrintf(PETSC_COMM_WORLD,"Number of sites: %d\n",params.ab_sites);
	params.sites = params.ab_sites;
	parse_a_sites_file(a_sites,&params);
	dec2bin(params.A_mask, params.sites, buf );
	PetscPrintf(PETSC_COMM_WORLD,"A mask: %s\n",buf);
	dec2bin(params.B_mask, params.sites, buf );
	PetscPrintf(PETSC_COMM_WORLD,"B mask: %s\n",buf);
	
	find_valid_configs(&params); //find valid configurations and add to a list.
	
	//count valid configs
	params.ab_configs = 0;
	for ( list<uPetscInt>::iterator iter = params.ab_configs_vec.begin();  iter != params.ab_configs_vec.end() ; iter++ ) {
		//dec2bin(*iter, params.sites, buf );
		//cout << buf << endl ; 
		params.ab_configs++;
	}
	cout << "Number of configs over AB : " << params.ab_configs << endl;
	
	//transfer them to an array
	ierr = PetscMalloc(params.ab_configs*sizeof(uPetscInt),&params.ab_configs_array);CHKERRQ(ierr);
	int idx = 0 ; 
	for ( list<uPetscInt>::iterator iter = params.ab_configs_vec.begin();  iter != params.ab_configs_vec.end() ; iter++ ) {
		params.ab_configs_array[idx++] = *iter;
	}
	params.ab_configs_vec.clear();
	
	//sort the array
	sort(params.ab_configs_array,params.ab_configs_array + params.ab_configs);
	/*for ( idx = 0 ; idx < params.ab_configs ; idx++ ) {
		dec2bin(params.ab_configs_array[idx], params.sites, buf );
		cout << buf << endl ; 
	}*/
	
	uPetscInt *tmp_A_configs, *tmp_B_configs;
	ierr = PetscMalloc(params.ab_configs*sizeof(uPetscInt),&tmp_A_configs);CHKERRQ(ierr);
	ierr = PetscMalloc(params.ab_configs*sizeof(uPetscInt),&tmp_B_configs);CHKERRQ(ierr);
	for ( idx = 0 ; idx < params.ab_configs ; idx++ ) {
		tmp_A_configs[idx] = params.ab_configs_array[idx] & params.A_mask ;
		tmp_B_configs[idx] = params.ab_configs_array[idx] & params.B_mask ;
	}
	sort(tmp_A_configs, tmp_A_configs + params.ab_configs); //sort array so that identical configs are beside each other.
	sort(tmp_B_configs, tmp_B_configs + params.ab_configs); //sort array so that identical configs are beside each other.
	//now count the number of unique configs.
	params.a_configs = 1 ; 
	uPetscInt current = tmp_A_configs[0]; 
	for ( idx = 1 ; idx < params.ab_configs ; idx++ ) {
		if ( tmp_A_configs[idx] != current ){
			params.a_configs++;
			current = tmp_A_configs[idx];
		}
	}
	
	params.b_configs = 1 ; 
	current = tmp_B_configs[0]; 
	for ( idx = 1 ; idx < params.ab_configs ; idx++ ) {
		if ( tmp_B_configs[idx] != current ){
			params.b_configs++;
			current = tmp_B_configs[idx];
		}
	}
	PetscPrintf(PETSC_COMM_WORLD,"Number of configs on A: %d\n",params.a_configs);
	PetscPrintf(PETSC_COMM_WORLD,"Number of configs on B: %d\n",params.b_configs);

	ierr = PetscMalloc(params.a_configs*sizeof(uPetscInt),&params.A_configs);CHKERRQ(ierr);
	params.A_configs[0] = tmp_A_configs[0];
	current = tmp_A_configs[0];
	PetscInt current_idx = 1;
	for ( idx = 1 ; idx < params.ab_configs ; idx++ ) {
		if ( tmp_A_configs[idx] != current ){
			params.A_configs[current_idx++] = tmp_A_configs[idx];
			current = tmp_A_configs[idx];
		}
	}
	PetscFree(tmp_A_configs);
	
	ierr = PetscMalloc(params.b_configs*sizeof(uPetscInt),&params.B_configs);CHKERRQ(ierr);
	params.B_configs[0] = tmp_B_configs[0];
	current = tmp_B_configs[0];
	current_idx = 1;
	for ( idx = 1 ; idx < params.ab_configs ; idx++ ) {
		if ( tmp_B_configs[idx] != current ){
			params.B_configs[current_idx++] = tmp_B_configs[idx];
			current = tmp_B_configs[idx];
		}
	}
	PetscFree(tmp_B_configs);
	
	/*cout << "A Configs " << endl << "--------------" << endl ; 
	for ( idx = 0 ; idx < params.a_configs ; idx++ ) {
		dec2bin(params.A_configs[idx], params.sites, buf );
		cout << buf << endl ; 
	}
	cout << "B Configs " << endl << "--------------" << endl ; 
	for ( idx = 0 ; idx < params.b_configs ; idx++ ) {
		dec2bin(params.B_configs[idx], params.sites, buf );
		cout << buf << endl ; 
	}*/
	
	//create dense matrix
	MatCreateSeqDense(PETSC_COMM_WORLD,params.a_configs,params.a_configs,PETSC_NULL,&params.rho_a);
	MatZeroEntries(params.rho_a);
	
	//read in vector that we are calculating entanglement entropy of
	PetscViewer viewer;
	PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
//    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, vec_file.c_str() ,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, vec_file.c_str() ,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecLoad(params.psi, viewer);CHKERRQ(ierr);
    PetscViewerDestroy(&viewer);
	
	
	PetscScalar *vec_elements;
	VecGetArray(params.psi,&vec_elements);
	
	/* commenting out tests
	my_uPetscInt a,b ;
	
	a.num = 5;
	a.mask = 1;
	b.num = 3;
	b.mask = 1;
	
	if ( a == b ) PetscPrintf(PETSC_COMM_WORLD,"A == B is true.\n");
	else if ( a > b ) PetscPrintf(PETSC_COMM_WORLD,"A > B is true.\n");
	else if ( a < b ) PetscPrintf(PETSC_COMM_WORLD,"A < B is true.\n");

	a.num = 4;
	if ( a == b ) PetscPrintf(PETSC_COMM_WORLD,"A == B is true.\n");
	else if ( a > b ) PetscPrintf(PETSC_COMM_WORLD,"A > B is true.\n");
	else if ( a < b ) PetscPrintf(PETSC_COMM_WORLD,"A < B is true.\n");

	a.num = 5;
	a.mask = 6;
	b.mask = 6;

	if ( a == b ) PetscPrintf(PETSC_COMM_WORLD,"A == B is true.\n");
	else if ( a > b ) PetscPrintf(PETSC_COMM_WORLD,"A > B is true.\n");
	else if ( a < b ) PetscPrintf(PETSC_COMM_WORLD,"A < B is true.\n");
	*/
	
	my_uPetscInt *my_array;
	ierr = PetscMalloc(params.ab_configs*sizeof(my_uPetscInt),&my_array);
	for ( PetscInt i = 0 ; i < params.ab_configs; i++ ) {
		my_array[i].num = params.ab_configs_array[i];
		my_array[i].mask = params.B_mask;
	}
	sort(my_array, my_array + params.ab_configs); //sort the array according to the configs on sites of A 

	/*for ( PetscInt i = 0 ; i < params.ab_configs; i++ ) {
		dec2bin(my_array[i].num, params.sites, buf );
		cout << buf << endl ; 
	}*/
	
	//loop over the A configs
	PetscInt ab_idx = 0 ;
	while ( ab_idx < params.ab_configs ) {
		//count the number of consecutive patterns that are the same. 
		int count = 0;
		uPetscInt b_pattern = my_array[ab_idx].num & my_array[ab_idx].mask;
		while ( ( (my_array[ab_idx + count].num & params.B_mask) == b_pattern) && ab_idx + count < params.ab_configs ) {
			count++; 
		}
		
		PetscInt i_ab_idx,k_ab_idx,i_a_idx,k_a_idx;
		
		//go through all the ways of pairing up these
		for ( PetscInt i_off = 0 ; i_off < count ; i_off++ ) {
			i_ab_idx = find_in_sorted_array(params.ab_configs_array, params.ab_configs,my_array[ab_idx+i_off].num);
			i_a_idx = find_in_sorted_array(params.A_configs, params.a_configs, my_array[ab_idx+i_off].num & params.A_mask);   
			for ( PetscInt k_off = 0 ; k_off < count ; k_off++ ) {
				k_ab_idx = find_in_sorted_array(params.ab_configs_array, params.ab_configs,my_array[ab_idx+k_off].num);
				k_a_idx = find_in_sorted_array(params.A_configs, params.a_configs, my_array[ab_idx+k_off].num & params.A_mask);   
				MatSetValue(params.rho_a,i_a_idx,k_a_idx,vec_elements[k_ab_idx]*PetscConj(vec_elements[i_ab_idx]),ADD_VALUES);
			}
		}		
		ab_idx += count;
		//PetscPrintf(PETSC_COMM_WORLD,"AB idx: %d\n",ab_idx);
	}
	
	/*PetscInt i_idx, k_idx;
	for ( PetscInt i = 0 ; i < params.ab_configs; i++ ) {
		i_idx = find_in_sorted_array(params.A_configs, params.a_configs, params.ab_configs_array[i] & params.A_mask); 
		for ( PetscInt k = 0 ; k < params.ab_configs; k++ ) {
			if ( (params.ab_configs_array[k] & params.B_mask) == (params.ab_configs_array[i]&params.B_mask) ) {
				k_idx = find_in_sorted_array(params.A_configs, params.a_configs, params.ab_configs_array[k] &params.A_mask); 
				MatSetValue(params.rho_a,i_idx,k_idx,vec_elements[k]*PetscConj(vec_elements[i]),ADD_VALUES);
			}
		}
	}*/
	
	VecRestoreArray(params.psi,&vec_elements);
	MatAssemblyBegin(params.rho_a,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(params.rho_a,MAT_FINAL_ASSEMBLY);
	
	//write the matrix to disk
	stringstream ss;
	ss << vec_file << ".rho_a.m" ;
	PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,ss.str().c_str(),&viewer);
	PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	MatView(params.rho_a,viewer);
	PetscViewerDestroy(&viewer);
	MatDestroy(&params.rho_a);	
	

	
	if ( calc_weight_matrix ) {
		//calcualte the weight matrix
		//loop over all the elements in the full basis. For each find the index of the a part and the b part. 
		MatCreateSeqDense(PETSC_COMM_WORLD,params.a_configs,params.b_configs,PETSC_NULL,&params.rho_a);
		MatZeroEntries(params.rho_a);
			
		VecGetArray(params.psi,&vec_elements);	
		PetscInt row, col ;
		for ( PetscInt i = 0 ; i < params.ab_configs; i++ ) {
			row = find_in_sorted_array(params.A_configs, params.a_configs, params.ab_configs_array[i] & params.A_mask);
			col = find_in_sorted_array(params.B_configs, params.b_configs, params.ab_configs_array[i] & params.B_mask);
			//dec2bin(params.ab_configs_array[i],params.sites,buf);
			//cout << buf << " with val " << vec_elements[i] << " and coords :" << row << ", " << col << endl ;
			MatSetValue(params.rho_a,row,col,vec_elements[i],ADD_VALUES);
		}
		VecRestoreArray(params.psi,&vec_elements);
		MatAssemblyBegin(params.rho_a,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(params.rho_a,MAT_FINAL_ASSEMBLY);
		
		//write weights matrix to disk
		ss.str("");
		ss << vec_file << ".W.m" ;
		PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
		PetscViewerASCIIOpen(PETSC_COMM_WORLD,ss.str().c_str(),&viewer);
		PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
		MatView(params.rho_a,viewer);
		PetscViewerDestroy(&viewer);
		
		MatDestroy(&params.rho_a);
	}
	
	
	PetscFree(params.B_configs);
	PetscFree(params.A_configs);
	PetscFree(params.ab_configs_array);
	PetscFree(my_array);

	ierr = SlepcFinalize();CHKERRQ(ierr);
	return 0;
}
