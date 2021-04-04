/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 * 
 * File containing functions called from the main file solver.c  
 * 
 * Version: $Id$ 
 * 
*/

#include "utils.h" 
#include "expectation_value_funcs.h"


/*! \brief This function prints the DoQO notice including version information upon DoQO starting.
	Printing is performed from the rank zero process and only if the verbosity is above zero.
*/
void 	write_doqo_heading(struct parameter_struct *parameters) {
	if ( parameters->rank == 0 ) {
		 PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------\n");
		 PetscPrintf(PETSC_COMM_WORLD,"| ______           _______  _______   |\n");
		 PetscPrintf(PETSC_COMM_WORLD,"| |   _  \\  .-----.|   _   ||   _   | |\n");
		 PetscPrintf(PETSC_COMM_WORLD,"| |.  |   \\ |  _  ||.  |   ||.  |   | |\n");
		 PetscPrintf(PETSC_COMM_WORLD,"| |.  |    \\|_____||.  |   ||.  |   | |\n");
		 PetscPrintf(PETSC_COMM_WORLD,"| |:  1    /       |:  1   ||:  1   | |\n");
		 PetscPrintf(PETSC_COMM_WORLD,"| |::.. . /        |::..   ||::.. . | |\n");
		 PetscPrintf(PETSC_COMM_WORLD,"| `------'         `----|:.|`-------' |\n");
		 PetscPrintf(PETSC_COMM_WORLD,"|                       `--'          |\n");
		 PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------\n");
		 PetscPrintf(PETSC_COMM_WORLD,"Diagonalisation of Quantum Operators\n",DOQO_VERSION);
		 PetscPrintf(PETSC_COMM_WORLD,"Version: %s\n",DOQO_VERSION);
 		 PetscPrintf(PETSC_COMM_WORLD,"Copyright 2010\n",DOQO_VERSION);
 		 PetscPrintf(PETSC_COMM_WORLD,"Free for non commercial use.\n",DOQO_VERSION);
 		 PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------\n");
	}								
}


#undef __FUNCT__
#define __FUNCT__ "Initialise"
/*! \brief Function that initialises SLEPc which in turn initialises PETSc and MPI. Default values for parameters also set.*/
int initialise(int argc, char **argv , struct parameter_struct *parameters ){
	char help[] = "Program to generate Hamiltonian for selected system and compute eigenvalues.\n";
	PetscErrorCode ierr;

	//Initialise Slepc
	SlepcInitialize(&argc,&argv,(char*)0,help);
	//get rank of this process and number of processes 
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&(parameters->rank));CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&(parameters->size));CHKERRQ(ierr);


	//set default values these will be overwritten with values from the input files if they exist
	strcpy(parameters->method,"krylovschur");
	parameters->model_type = SPIN_HALF;
	parameters->eigPairs = 2 ;
	parameters->prepare = PETSC_TRUE;
	parameters->distribute = PETSC_TRUE;
	parameters->solver = PETSC_TRUE;
	parameters->use_sectors = 0; //this mean use all the sectors
	parameters->benchmark = PETSC_FALSE;
	parameters->save_states = PETSC_FALSE;
	parameters->save_states_ascii = PETSC_FALSE;
	parameters->save_states_matlab = PETSC_FALSE;
	parameters->save_states_real = PETSC_FALSE;
	parameters->multiply_test = PETSC_FALSE;
	parameters->max_its = 500;
	parameters->verbosity = 1;
	parameters->save_matrix = PETSC_FALSE;
	parameters->save_matrix_ascii = PETSC_FALSE;
	parameters->save_matrix_matlab = PETSC_FALSE;
	parameters->use_hermitian = PETSC_FALSE;
	parameters->use_disk = PETSC_FALSE;
	parameters->use_bst = PETSC_FALSE;
	parameters->other_wf_overlap = PETSC_FALSE;
	parameters->staggered_subspace = PETSC_FALSE; // this is a flag which when true restricts to space of states where each staggered site is blocked by an occupied neighbour.
	parameters->tasks.number_parameters = 0 ;
	parameters->tasks.start_task = 0;
	parameters->save_basis = PETSC_FALSE;

	//set options to defaults before reading command line options
	parameters->tol = 1.0e-13;
	parameters->phase_tol = 1.0e-10;
	parameters->deg_tol = 1.0e-10;

	parameters->use_parity_sectors = PETSC_FALSE;
	parameters->parity_info.num_relevant_sectors=0;
	parameters->use_conservation_sectors = PETSC_FALSE;
	parameters->conservation_filling = -1;
	parameters->use_momentum_sectors = PETSC_FALSE;
	parameters->momentum_info.num_relevant_sectors=0;
	parameters->use_spectral_flow = PETSC_FALSE;
	parameters->spectral_info.number_relevant_points[0] = 1;
	parameters->spectral_info.number_relevant_points[1] = 1;
	parameters->use_rotation_invariance = PETSC_FALSE;
	parameters->rotation_number_relevant_sectors = 1; 
	
	parameters->nn_exclusion = PETSC_FALSE;
	parameters->nn_recursive_algorithm = PETSC_FALSE;
	parameters->calculate_correlations = PETSC_FALSE;
	parameters->prod_wf_overlap = PETSC_FALSE;
	parameters->calculate_expectation_values = PETSC_FALSE;
	
	//set the default function pointers so that the full general basis is used. 
	parameters->get_full_basis_index = &get_full_basis_index_passthru;
	parameters->get_reduced_basis_index = &get_reduced_basis_index_passthru;
	parameters->simple_basis_check = &dummy_basis_check;

	parameters->representative_check = &dummy_representative_check; 
	parameters->group_representative = &dummy_group_representative;
	parameters->group_normal = &dummy_group_normal; 
	parameters->translate_using_masks = &translate_using_masks_spins;
	parameters->rotate_using_masks = &apply_rotation;
	parameters->getrow = &getrow_spin_half;

	return 0;
} 

#undef __FUNCT__
#define __FUNCT__ "Parse_input"
/*! \brief This function parses the master input file which set describes sets parameters and other input files. */
int parse_main_input_file(struct parameter_struct *parameters) { 
	PetscErrorCode ierr;
	PetscTruth flg;

	//Moving to new input mechanism with just one xml input file specifying all parameters
	ierr = PetscOptionsGetString(NULL, NULL, "-input",parameters->input_file,MAX_STRING,&flg);CHKERRQ(ierr);
	if ( flg == PETSC_FALSE) {
	  PetscPrintf(PETSC_COMM_WORLD, "No input file specified exiting.\n");
	  return 1;
	}

	string input_filename(parameters->input_file);
	string buf =  input_filename.substr(0, input_filename.find_last_of(".") );
	strcpy(parameters->output_prefix,buf.c_str());
	
	
	//now we attempt to read the input xml file and set the relevant parameters.
	//TiXmlDocument doc(parameters->input_file);
	//TiXmlDocument doc("sample_input.xml");
	string filename(parameters->input_file);
	TiXmlDocument doc(filename);
	
	bool loadOkay = doc.LoadFile();
	
	if ( !loadOkay )
	{
		PetscPrintf(PETSC_COMM_WORLD,"Could not load the input file '%s'. Error='%s'. Exiting.\n", parameters->input_file, doc.ErrorDesc() );
		return 1;
	}
	
	TiXmlHandle docHandle( &doc );
	
	TiXmlElement* element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "VERBOSITY").Element();
	if ( element)
	{	
		parameters->verbosity = atoi(element->GetText());
		stringstream ss;
		ss << "Verbosity level: " << parameters->verbosity << endl;
		if ( parameters->verbosity >= 1) { PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());}
	}
	
	if ( parameters->verbosity >= 2) { 
		PetscPrintf(PETSC_COMM_WORLD,"Input file: '%s'\n", parameters->input_file);
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "MODEL_TYPE").Element();
	if ( element )
	{	
		if ( strcmp(element->GetText(),"SPIN_HALF") == 0 ) {
			parameters->model_type = SPIN_HALF;
			parameters->getrow = &getrow_spin_half;
			parameters->getrowlocal = &getrowlocal_spin_half;
			parameters->translate_using_masks = &translate_using_masks_spins;
		} else if ( strcmp(element->GetText(),"FERMIONIC") == 0 ) {
			parameters->model_type = FERMIONIC;
			parameters->getrow = &getrow_fermionic;
			parameters->getrowlocal = &getrowlocal_fermionic;
			parameters->translate_using_masks = &translate_using_masks_fermions;
		} else if ( strcmp(element->GetText(),"SPIN_ONE") == 0 ) {
			parameters->model_type = SPIN_ONE;
			parameters->getrow = &getrow_spin_one;
			parameters->getrowlocal = &getrowlocal_spin_one;
			parameters->translate_using_masks = &translate_using_masks_spins;
		}		
	}
	if ( parameters->verbosity >= 2) {
		if ( parameters->model_type == SPIN_HALF ) {
			PetscPrintf(PETSC_COMM_WORLD, "Model type: Spin half\n");
		} else if ( parameters->model_type == FERMIONIC ) {
			PetscPrintf(PETSC_COMM_WORLD, "Model type: Fermionic\n");
		} else if ( parameters->model_type == SPIN_ONE ) {
			PetscPrintf(PETSC_COMM_WORLD, "Model type: Spin one\n");
		}
	}
	
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "MODEL_FILE").Element();
	if ( element)
	{	
		strcpy(parameters->hamiltonian.filename,element->GetText());
		if ( parameters->verbosity >= 2 ) PetscPrintf(PETSC_COMM_WORLD,"Model file: '%s'\n", parameters->hamiltonian.filename);
	} else  { //if there was no hamiltonian file specified
		PetscPrintf(PETSC_COMM_WORLD, "Error no hamiltonian file specified.\n");
		return 1;
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "TASK_LIST").Element();
	if ( element)
	{	
		strcpy(parameters->tasks.filename,element->GetText());
		if ( parameters->verbosity >= 2 ) PetscPrintf(PETSC_COMM_WORLD,"Tasks file: '%s'\n", parameters->tasks.filename);
		if ( element->Attribute("start_task") != NULL ) {
			parameters->tasks.start_task = atoi(element->Attribute("start_task"));
		}
	} else  { //if there was no task file specified
		PetscPrintf(PETSC_COMM_WORLD, "Error no task file specified.\n");
		return 1;
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "OUTPUT_PREFIX").Element();
	if ( element)
	{	
		strcpy(parameters->output_prefix,element->GetText());
		if ( parameters->verbosity >= 2 ) PetscPrintf(PETSC_COMM_WORLD,"Output prefix: '%s'\n", parameters->output_prefix);
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "EIGENVALUES").Element();
	if ( element)
	{	
		parameters->eigPairs = atoi(element->GetText());
		if ( parameters->verbosity >= 2 ) PetscPrintf(PETSC_COMM_WORLD,"Number of eigenpairs: '%d'\n", (int)parameters->eigPairs);
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "PREPARE_STAGE").Element();
	if ( element)
	{	
		if ( strcmp(element->GetText(),"true") == 0 ) parameters->prepare = PETSC_TRUE;
		else parameters->prepare = PETSC_FALSE;
	}

	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "NN_EXCLUSION").Element();
	if ( element)
	{	
		if ( strcmp(element->GetText(),"true") == 0 ) {
			parameters->nn_exclusion = PETSC_TRUE;
			char buf2[MAX_STRING];
			if ( element->Attribute("adjacency_file") == NULL ) {
				PetscPrintf(PETSC_COMM_WORLD,"Adjacency file not specified.\n");
				return 1;
			}
			strcpy(buf2,element->Attribute("adjacency_file"));
			if ( parse_adjacency_info_file(parameters,buf2) != 0 ) {
			PetscFPrintf(PETSC_COMM_WORLD,stderr,"Error parsing adjacency file %s.\n", buf2);
				return 1 ;
			}
			if ( element->Attribute("recursive") != NULL ) {
				strcpy(buf2,element->Attribute("recursive"));
				if (strcmp(buf2,"true") == 0 ) {
					parameters->nn_recursive_algorithm = PETSC_TRUE;
				}
			}
			parameters->simple_basis_check = &check_adjacency;
			if ( parameters->verbosity >= 2 ) PetscPrintf(PETSC_COMM_WORLD,"Nearest neighbour exclusion: enabled\n");
		}  else parameters->nn_exclusion = PETSC_FALSE;
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "DISTRIBUTE_STAGE").Element();
	if ( element)
	{	
		if ( strcmp(element->GetText(),"true") == 0 ) 
			parameters->distribute = PETSC_TRUE;
 	 	else parameters->distribute = PETSC_FALSE;
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SOLVE_STAGE").Element();
	if ( element)
	{	
		if ( strcmp(element->GetText(),"true") == 0 ) {
			parameters->solver = PETSC_TRUE;
		} else {parameters->solver = PETSC_FALSE;}
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "MULTIPLY_TEST").Element();
	if ( element)
	{	
		if ( strcmp(element->GetText(),"true") == 0 ) {
			parameters->multiply_test = PETSC_TRUE;
			if ( parameters->verbosity >= 2 ) PetscPrintf(PETSC_COMM_WORLD,"Multiplication test: enabled\n");
		} else parameters->multiply_test = PETSC_FALSE;	
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SAVE_STATES").Element();
	if ( element)
	{	
		if ( strcmp(element->GetText(),"true") == 0 ) {	
			parameters->save_states = PETSC_TRUE;
			if ( element->Attribute("format") != NULL ) {
				if (strcmp("ascii",element->Attribute("format")) == 0 ) {
					parameters->save_states_ascii = PETSC_TRUE;
				} else if (strcmp("matlab",element->Attribute("format")) == 0 ) {
					parameters->save_states_ascii = PETSC_TRUE;
					parameters->save_states_matlab = PETSC_TRUE;
				} 
			}
			if ( element->Attribute("real") != NULL ) {
				if ( strcmp("true",element->Attribute("real")) == 0 ) {
					parameters->save_states_real = PETSC_TRUE;
				} else parameters->save_states_real = PETSC_FALSE;
			}
			if ( parameters->verbosity >= 2) { 
				if (parameters->save_states_ascii ) {
					PetscPrintf(PETSC_COMM_WORLD,"Saving states: enabled format ascii\n");
				} else PetscPrintf(PETSC_COMM_WORLD,"Saving states: enabled format binary\n");
				if (parameters->save_states_real ) {
					PetscPrintf(PETSC_COMM_WORLD,"Saving states as reals if possible.\n");
				}
			}
			
		} else parameters->save_states = PETSC_FALSE;	

	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "USE_HERMITIAN").Element();
	if ( element)
	{	
		if ( strcmp(element->GetText(),"true") == 0 ) parameters->use_hermitian = PETSC_TRUE;
		else  parameters->use_hermitian = PETSC_FALSE;	
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "BENCHMARK").Element();
	if ( element)
	{	
		if ( strcmp(element->GetText(),"true") == 0 ) {
			parameters->benchmark = PETSC_TRUE;
			if ( parameters->verbosity >= 2) { PetscPrintf(PETSC_COMM_WORLD,"Benchmark enabled\n");}
		} else parameters->benchmark = PETSC_FALSE;
		
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SOLVER_TYPE").Element();
	if ( element)
	{	
		strcpy(parameters->method,element->GetText());
		if ( parameters->verbosity >= 2) { PetscPrintf(PETSC_COMM_WORLD,"Solver type: %s\n",parameters->method);}
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "MAX_ITERATIONS").Element();
	if ( element)
	{	
		parameters->max_its = (PetscInt)atoi(element->GetText());
		stringstream ss;
		ss << "Maximum number of iterations: " << parameters->max_its << endl;
		if ( parameters->verbosity >= 2) PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
	}

	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "DEGENERACY_TOLERANCE").Element();
	if ( element)
	{	
		parameters->deg_tol = atof(element->GetText());
		if ( parameters->verbosity >= 2) {PetscPrintf(PETSC_COMM_WORLD,"Degeneracy tolerance %le\n",parameters->deg_tol);}
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SAVE_MATRIX").Element();
	if ( element)
	{	
		if ( strcmp(element->GetText(),"true") == 0 ) {
			parameters->save_matrix = PETSC_TRUE;
			if ( parameters->verbosity >= 2) PetscPrintf(PETSC_COMM_WORLD,"Saving matrces: enabled\n");
			if ( element->Attribute("format") != NULL ) {
				if (strcmp("ascii",element->Attribute("format")) == 0 ) {
					parameters->save_matrix_ascii = PETSC_TRUE;
				} else if (strcmp("matlab",element->Attribute("format")) == 0 ) {
					parameters->save_matrix_ascii = PETSC_TRUE;
					parameters->save_matrix_matlab = PETSC_TRUE;
				} 
			}
		}	else parameters->save_matrix = PETSC_FALSE;
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "USE_DISK").Element();
	if ( element)
	{	
		if ( strcmp(element->GetText(),"true") == 0 ) {
			parameters->use_disk = PETSC_TRUE;
			if ( parameters->verbosity >= 2) PetscPrintf(PETSC_COMM_WORLD,"Using disk: enabled\n");
		} else parameters->use_disk = PETSC_FALSE;
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "USE_BST").Element();
	if ( element)
	{	
		if ( strcmp(element->GetText(),"true") == 0 ) {
			#ifdef _BOOST_BST_
			parameters->use_bst = PETSC_TRUE;
			if ( parameters->verbosity >= 2) PetscPrintf(PETSC_COMM_WORLD,"Using binary sorting trees: enabled\n");
			#else
			if ( parameters->verbosity > 0) PetscPrintf(PETSC_COMM_WORLD,"Binary source tree option specified but DoQO is not compiled with support for BST. See readme.\n");
			#endif
		} else parameters->use_bst = PETSC_FALSE;
	}
	
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SOLVER_TOLERANCE").Element();
	if ( element)
	{	
		parameters->tol = atof(element->GetText());
		if ( parameters->verbosity >= 2) {PetscPrintf(PETSC_COMM_WORLD,"Solver tolerance: %le\n",parameters->tol);}
	}
	
	
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "PHASE_TOLERANCE").Element();
	if ( element)
	{	
		parameters->phase_tol = atof(element->GetText());
		if ( parameters->verbosity >= 2) {PetscPrintf(PETSC_COMM_WORLD,"Phase tolerance: %le\n",parameters->phase_tol);}
	}


	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SYMMETRIES").FirstChild("MOMENTUM").Element();
	if ( element )
	{	
		parameters->use_momentum_sectors = PETSC_TRUE;
		parameters->representative_check = &momentum_representative_check; 
    	parameters->group_representative = &momentum_group_representative;
    	parameters->group_normal = &momentum_group_normal; 
		strcpy(parameters->momentum_info_file,element->Attribute("file"));
		if ( parameters->verbosity >= 2) {PetscPrintf(PETSC_COMM_WORLD,"Translation invariance: enabled, lattice details read from '%s'\n",parameters->momentum_info_file);}
		parameters->momentum_info.num_relevant_sectors = 0 ; //if zero all sectors are used
		element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SYMMETRIES").FirstChild("MOMENTUM").FirstChild("RELEVANT_SECTORS").Element();
		if ( element ) {
			char buf2[MAX_STRING];
			if ( element->Attribute("number") == NULL ) {
				PetscPrintf(PETSC_COMM_WORLD,"Relevant sectors must have attribute \"number\".\n");
				return 1;
			}
			strcpy(buf2,element->Attribute("number"));
			parameters->momentum_info.num_relevant_sectors = atoi(buf2);
			ierr = PetscMalloc(parameters->momentum_info.num_relevant_sectors*sizeof(int),&parameters->momentum_info.relevant_sectors);CHKERRQ(ierr);
			for ( int i = 0 ; i < parameters->momentum_info.num_relevant_sectors ; i++ ) {
				element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SYMMETRIES").FirstChild("MOMENTUM").FirstChild("RELEVANT_SECTORS").Child("SECTOR",i).Element();
				if ( element ) {
					parameters->momentum_info.relevant_sectors[i] = atoi(element->GetText());
				} else {
					PetscPrintf(PETSC_COMM_WORLD,"Number of relevant momentum sectors incorrect.\n");
					return 1;
				}
			}
			if ( parameters->verbosity >= 3) {
				PetscPrintf(PETSC_COMM_WORLD,"\tNumber of sectors to use: %d (",parameters->momentum_info.num_relevant_sectors);
				for ( int i = 0 ; i < parameters->momentum_info.num_relevant_sectors ; i++ ) {
					PetscPrintf(PETSC_COMM_WORLD,"%d ", parameters->momentum_info.relevant_sectors[i]);
				}
				PetscPrintf(PETSC_COMM_WORLD,")\n");
			}
		}
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SYMMETRIES").FirstChild("PARITY").Element();
	if ( element )
	{	
		if ( parameters->nn_recursive_algorithm ==  PETSC_FALSE || parameters->nn_exclusion == PETSC_FALSE || parameters->use_disk == PETSC_FALSE)  {
			parameters->get_full_basis_index = &get_full_parity_basis_index;
			parameters->get_reduced_basis_index = &get_reduced_parity_basis_index;
		}
		parameters->use_parity_sectors = PETSC_TRUE;
		strcpy(parameters->parity_masks,element->Attribute("file"));
		if ( parameters->verbosity >= 2) { PetscPrintf(PETSC_COMM_WORLD,"Conservation of parity: enabled, sites affected read from '%s'\n",parameters->parity_masks);}
		parameters->parity_info.num_relevant_sectors = 0; //if zero all sectors are used
		element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SYMMETRIES").FirstChild("PARITY").FirstChild("RELEVANT_SECTORS").Element();
		if ( element ) {
			char buf2[MAX_STRING];
			if ( element->Attribute("number") == NULL ) {
				PetscPrintf(PETSC_COMM_WORLD,"Relevant sectors must have attribute \"number\".\n");
				return 1;
			}
			strcpy(buf2,element->Attribute("number"));
			parameters->parity_info.num_relevant_sectors = atoi(buf2);
			ierr = PetscMalloc(parameters->parity_info.num_relevant_sectors*sizeof(int),&parameters->parity_info.relevant_sectors);CHKERRQ(ierr);
			for ( int i = 0 ; i < parameters->parity_info.num_relevant_sectors ; i++ ) {
				element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SYMMETRIES").FirstChild("PARITY").FirstChild("RELEVANT_SECTORS").Child("SECTOR",i).Element();
				if ( element ) {
					parameters->parity_info.relevant_sectors[i] = atoi(element->GetText());
				} else {
					PetscPrintf(PETSC_COMM_WORLD,"Number of relevant parity sectors incorrect.\n");
					return 1;
				}
			}
			if ( parameters->verbosity >= 3) {
				PetscPrintf(PETSC_COMM_WORLD,"\tNumber of sectors to use: %d (",parameters->parity_info.num_relevant_sectors);
				for ( int i = 0 ; i < parameters->parity_info.num_relevant_sectors ; i++ ) {
					PetscPrintf(PETSC_COMM_WORLD,"%d ", parameters->parity_info.relevant_sectors[i]);
				}
				PetscPrintf(PETSC_COMM_WORLD,")\n");
			}
		}
	}
	
	//read in details of conservation symmetries if they are specified in the input file.
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SYMMETRIES").FirstChild("FILLING").Element();
	if ( element )
	{	
		if ( parameters->use_parity_sectors == PETSC_TRUE ) {
			PetscPrintf(PETSC_COMM_WORLD,"Cannot use both conservation of filling and parity at the same time.\n");
			return 1;
		}
		
		if ( parameters->nn_recursive_algorithm ==  PETSC_FALSE || parameters->nn_exclusion == PETSC_FALSE || parameters->use_disk == PETSC_FALSE)  {
			parameters->get_full_basis_index = &get_full_conservation_basis_index;
			parameters->get_reduced_basis_index = &get_reduced_conservation_basis_index;
		}
		
		parameters->use_conservation_sectors = PETSC_TRUE;
		strcpy(parameters->conservation_masks,element->Attribute("file"));
		if ( parameters->verbosity >= 2) { PetscPrintf(PETSC_COMM_WORLD,"Conservation of filling: enabled, sites affected read from '%s'\n",parameters->conservation_masks);}
		parameters->parity_info.num_relevant_sectors = 0; //if its zero all sectors are used
		element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SYMMETRIES").FirstChild("FILLING").FirstChild("RELEVANT_SECTORS").Element();
		if ( element ) {
			parameters->parity_info.num_relevant_sectors = atoi(element->Attribute("number"));
			ierr = PetscMalloc(parameters->parity_info.num_relevant_sectors*sizeof(int),&parameters->parity_info.relevant_sectors);CHKERRQ(ierr);
			for ( int i = 0 ; i < parameters->parity_info.num_relevant_sectors ; i++ ) {
				element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SYMMETRIES").FirstChild("FILLING").FirstChild("RELEVANT_SECTORS").Child("SECTOR",i).Element();
				if ( element ) {
					parameters->parity_info.relevant_sectors[i] = atoi(element->GetText());
				} else {
					PetscPrintf(PETSC_COMM_WORLD,"Number of relevant conservation sectors incorrect.\n");
					return 1;
				}
			}
			if ( parameters->verbosity >= 3) {
				stringstream ss ;
				ss << "\tNumber of sectors to use: " << parameters->parity_info.num_relevant_sectors << " (" ; 
				PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
				for ( int i = 0 ; i < parameters->parity_info.num_relevant_sectors ; i++ ) {
					PetscPrintf(PETSC_COMM_WORLD,"%d ", (int)parameters->parity_info.relevant_sectors[i]);
				}
				PetscPrintf(PETSC_COMM_WORLD,")\n");
			}
		}
	}
	
	//need to setup operators before setting up correlators so we know the number of particles.
	if ( setup_operators(parameters) > 0 ) { //setup the operators so they are ready to be distributed
		PetscFinalize();
		return 1 ;
    }
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "CALCULATE_CORRELATIONS").Element();
	if ( element ) {	
		parameters->number_ops = atoi(element->Attribute("number"));
		if ( parameters->number_ops > 0 ) {
			parameters->calculate_correlations = PETSC_TRUE;
			ierr = PetscMalloc(parameters->number_ops*sizeof(struct matrix_meta_data), &parameters->ops_data);CHKERRQ(ierr);
			for ( int i = 0 ; i < parameters->number_ops ; i++ ) {
				parameters->ops_data[i].number_terms = 1;
				parameters->ops_data[i].number_parameters = 0 ;
				allocate_matrix_data(&parameters->ops_data[i]);
				//strcpy(parameters->ops_data[i].parameter_names[0],"corr"); 
				element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "CALCULATE_CORRELATIONS").Child("CORRELATOR",i).Element();
				if ( element ) { 
					int points = atoi(element->Attribute("number")); //the number of sites for the correlator.
					TiXmlElement *element2;
					stringstream ss;
					for ( int j = 0 ; j < points ; j++ ) {
						element2 = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "CALCULATE_CORRELATIONS").Child("CORRELATOR",i).Child("SITE",j).Element();
						if ( element2 ) {
							int site = atoi(element2->GetText());
							if ( j > 0 ) ss << "," ;
							if ( parameters->model_type == SPIN_HALF ) {
								ss << site << " Z" ;
							} else if ( parameters->model_type == FERMIONIC ) {
								ss << site << " C," << site << " A" ;
							}
						} else {
							PetscPrintf(PETSC_COMM_WORLD,"Number of sites specified for correlator does not match number of points specified.\n");
							return 1;
						}
					} //end loop over points
					ss << "*" ;
					strcpy(parameters->ops_data[i].terms[0], ss.str().c_str());
					prepare_term_data(parameters,&parameters->ops_data[i],0);
					parameters->ops_data[i].term_multipliers[0] = 1.0;
					//parameters->ops_data[i].imag_term_multipliers[0] = 0.0;
					//deallocate raw terms to save memory
					ierr = PetscFree(parameters->ops_data[i].terms[0]);CHKERRQ(ierr);
					ierr = PetscFree(parameters->ops_data[i].terms);CHKERRQ(ierr);
				} else {
					PetscPrintf(PETSC_COMM_WORLD,"Number of correlators specified does not match number given.\n");
					return 1;
				}
			} //end loop over ops 
		} //end if number_ops > 0 
	} //end if element


	//Now read in data that is used to calculation spectral flow information.
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("SPECTRAL_FLOW").Element();
	if ( element ) {	
		parameters->use_spectral_flow = PETSC_TRUE;
		parameters->translate_using_masks = &translate_using_masks_fermions_spectral;
		parameters->group_normal = &momentum_group_normal_spectral;
		
#ifndef PETSC_USE_COMPLEX
		PetscPrintf(PETSC_COMM_WORLD,"Spectral flow feature will not work without complex support.\n");
		return 1;
#endif
		
		element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("SPECTRAL_FLOW").FirstChild("BORDER_PARAM_0").Element();
		if ( element ) {
			strcpy(parameters->spectral_info.alpha_param[0],element->GetText());
		}
		element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("SPECTRAL_FLOW").FirstChild("BORDER_CONJUGATE_PARAM_0").Element();
		if ( element ) {
			strcpy(parameters->spectral_info.alpha_conj_param[0],element->GetText());
		}
		
		element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("SPECTRAL_FLOW").FirstChild("BORDER_PARAM_1").Element();
		if ( element ) {
			strcpy(parameters->spectral_info.alpha_param[1],element->GetText());
		}
		element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("SPECTRAL_FLOW").FirstChild("BORDER_CONJUGATE_PARAM_1").Element();
		if ( element ) {
			strcpy(parameters->spectral_info.alpha_conj_param[1],element->GetText());
		}
		
		if ( parameters->verbosity >= 2) { 
			PetscPrintf(PETSC_COMM_WORLD,"Edge 0 parameter:'%s'\n",parameters->spectral_info.alpha_param[0]);
			PetscPrintf(PETSC_COMM_WORLD,"Conjugate edge 0 parameter:'%s'\n",parameters->spectral_info.alpha_conj_param[0]);
			PetscPrintf(PETSC_COMM_WORLD,"Edge 1 parameter:'%s'\n",parameters->spectral_info.alpha_param[1]);
			PetscPrintf(PETSC_COMM_WORLD,"Conjugate edge 1 parameter:'%s'\n",parameters->spectral_info.alpha_conj_param[1]);
		}
		
		element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("SPECTRAL_FLOW").FirstChild("POINTS_0").Element();
		if ( element ) {
			parameters->spectral_info.number_points[0] = atoi(element->GetText());
			parameters->spectral_info.number_relevant_points[0] = parameters->spectral_info.number_points[0];
			if ( parameters->verbosity >= 2) PetscPrintf(PETSC_COMM_WORLD,"Number of points for spectral flow 0: %d\n",(int)parameters->spectral_info.number_points[0]);
		} else {
			PetscPrintf(PETSC_COMM_WORLD,"Number of points to use for spectral flow calculation not specified in 0 direction.\n");
			return 1;
		}
		
		element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("SPECTRAL_FLOW").FirstChild("POINTS_1").Element();
		if ( element ) {
			parameters->spectral_info.number_points[1] = atoi(element->GetText());
			parameters->spectral_info.number_relevant_points[1] = parameters->spectral_info.number_points[1];
			if ( parameters->verbosity >= 2) PetscPrintf(PETSC_COMM_WORLD,"Number of points for spectral flow 1: %d\n",(int)parameters->spectral_info.number_points[1]);
		} else {
			PetscPrintf(PETSC_COMM_WORLD,"Number of points to use for spectral flow calculation not specified in 1 direction.\n");
			return 1;
		}
		
		//parameters->spectral_info.number_relevant_points = 0 ; //if zero all sectors are used
		element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("SPECTRAL_FLOW").FirstChild("RELEVANT_SECTORS_0").Element();
		if ( element ) {
			char buf2[MAX_STRING];
			if ( element->Attribute("number") == NULL ) {
				PetscPrintf(PETSC_COMM_WORLD,"Relevant sectors must have attribute \"number\".\n");
				return 1;
			}
			strcpy(buf2,element->Attribute("number"));
			parameters->spectral_info.number_relevant_points[0] = atoi(buf2);
			ierr = PetscMalloc(parameters->spectral_info.number_relevant_points[0]*sizeof(int),&parameters->spectral_info.relevant_points[0]);CHKERRQ(ierr);
			for ( int i = 0 ; i < parameters->spectral_info.number_relevant_points[0] ; i++ ) {
				element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SPECTRAL_FLOW").FirstChild("RELEVANT_SECTORS_0").Child("SECTOR",i).Element();
				if ( element ) {
					parameters->spectral_info.relevant_points[0][i] = atoi(element->GetText());
				} else {
					PetscPrintf(PETSC_COMM_WORLD,"Number of relevant points in direction 0 for spectral flow is incorrect.\n");
					return 1;
				}
			}
			if ( parameters->verbosity >= 3) {
				PetscPrintf(PETSC_COMM_WORLD,"\tNumber of points to use for spectral flow calculations in 0 direction: %d (",(int)parameters->spectral_info.number_relevant_points[0]);
				for ( int i = 0 ; i < parameters->spectral_info.number_relevant_points[0] ; i++ ) {
					PetscPrintf(PETSC_COMM_WORLD,"%d ", (int)parameters->spectral_info.relevant_points[0][i]);
				}
				PetscPrintf(PETSC_COMM_WORLD,")\n");
			}
		}
		//parameters->spectral_info.number_relevant_points = 0 ; //if zero all sectors are used
		element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("SPECTRAL_FLOW").FirstChild("RELEVANT_SECTORS_1").Element();
		if ( element ) {
			char buf2[MAX_STRING];
			if ( element->Attribute("number") == NULL ) {
				PetscPrintf(PETSC_COMM_WORLD,"Relevant sectors must have attribute \"number\".\n");
				return 1;
			}
			strcpy(buf2,element->Attribute("number"));
			parameters->spectral_info.number_relevant_points[1] = atoi(buf2);
			ierr = PetscMalloc(parameters->spectral_info.number_relevant_points[1]*sizeof(int),&parameters->spectral_info.relevant_points[1]);CHKERRQ(ierr);
			for ( int i = 0 ; i < parameters->spectral_info.number_relevant_points[1] ; i++ ) {
				element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SPECTRAL_FLOW").FirstChild("RELEVANT_SECTORS_1").Child("SECTOR",i).Element();
				if ( element ) {
					parameters->spectral_info.relevant_points[1][i] = atoi(element->GetText());
				} else {
					PetscPrintf(PETSC_COMM_WORLD,"Number of relevant points in direction 1 for spectral flow is incorrect.\n");
					return 1;
				}
			}
			if ( parameters->verbosity >= 3) {
				PetscPrintf(PETSC_COMM_WORLD,"\tNumber of points to use for spectral flow calculations in 1 direction: %d (",(int)parameters->spectral_info.number_relevant_points[1]);
				for ( int i = 0 ; i < parameters->spectral_info.number_relevant_points[1] ; i++ ) {
					PetscPrintf(PETSC_COMM_WORLD,"%d ", (int)parameters->spectral_info.relevant_points[1][i]);
				}
				PetscPrintf(PETSC_COMM_WORLD,")\n");
			}
		}
		
	} //end if <SPECTRAL_FLOW> section list exists.
	
	//section to read in product wavefunction over lap details
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "PROD_WF_OVERLAP").Element();
	if ( element)
	{	
		parameters->prod_wf_overlap = PETSC_TRUE;
		char buf2[MAX_STRING];
		if ( element->Attribute("file") == NULL ) {
			PetscPrintf(PETSC_COMM_WORLD,"Product wavefunction overlap file not specified.\n");
			return 1;
		}
		strcpy(buf2,element->Attribute("file"));
		if ( parse_prod_wf_overlap_file(parameters,buf2) != 0 ) {
		PetscFPrintf(PETSC_COMM_WORLD,stderr,"Error parsing product wavefunction overlap file %s.\n", buf2);
			return 1 ;
		}
		if ( parameters->verbosity >= 2 ) PetscPrintf(PETSC_COMM_WORLD,"Product wavefunction overlap: enabled\n");
	}

	//section to read in product wavefunction over lap details
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "OTHER_WF_OVERLAP").Element();
	if ( element)
	{	
		parameters->other_wf_overlap = PETSC_TRUE;
		//char buf2[MAX_STRING];
		if ( element->Attribute("file") == NULL ) {
			PetscPrintf(PETSC_COMM_WORLD,"Other wavefunction overlap file not specified.\n");
			return 1;
		}
		strcpy(parameters->other_wf_file,element->Attribute("file"));

		if ( parameters->verbosity >= 2 ) PetscPrintf(PETSC_COMM_WORLD,"Other wavefunction overlap enabled with wf from file %s.\n",parameters->other_wf_file);
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "STAGGERED_SUBSPACE").Element();
	if ( element)
	{	
		if ( strcmp(element->GetText(),"true") == 0 ) { 
			parameters->staggered_subspace = PETSC_TRUE;
			parameters->simple_basis_check = &susy_staggered_subspace_check;
		} 
	}
	
	int return_val;
	if ( (return_val = parse_rotation_info(parameters,&docHandle)) > 0 ) {
		return return_val;
	}
	
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild( "SAVE_BASIS").Element();
	if ( element)
	{	
		if ( strcmp(element->GetText(),"true") == 0 ) {	
		  parameters->save_basis = PETSC_TRUE;
		}
	}
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("EXPECTATION_VALUES").Element();
	if ( element)
	{	
		if ( element->Attribute("file") == NULL ) {
			PetscPrintf(PETSC_COMM_WORLD,"Must specify file with information about vectors and sectors.\n");
			return 1;
		}
		if ( element->Attribute("overlap_vec") != NULL ) 
		{
			parameters->exp_vals_info.overlap_vec.assign(element->Attribute("overlap_vec"));
		} 
		else
		{
			parameters->exp_vals_info.overlap_vec.assign("");
		}
		
		strcpy(parameters->expectation_values_file,element->Attribute("file"));
		parameters->calculate_expectation_values = PETSC_TRUE;
	}		

	return 0;	
}

#undef __FUNCT__
#define __FUNCT__ "read_input_files"
/*! \brief Function to read additional input files. */
int read_input_files(struct parameter_struct *parameters) {
	if ( parse_task_file(&parameters->tasks, parameters, &parameters->hamiltonian) > 0 ) { // parse the task file 
		PetscPrintf(PETSC_COMM_WORLD,"Error parsing tasks file. Make sure this file exists and has the correct format.\n");
		return 1 ; 
	}
	
	if ( parameters->use_momentum_sectors == PETSC_TRUE ) {
		if ( (parse_momentum_file(parameters) ) != 0 ) {
			PetscPrintf(PETSC_COMM_WORLD,"Error parsing momentum information file.\n");
			return 1;
		}
		//print_momentum_information(parameters);
	}
	
	struct matrix_meta_data *md;
	md = &parameters->hamiltonian;
	if ( parameters->use_parity_sectors == PETSC_TRUE ) {
		if ( read_parity_masks(parameters->parity_masks, &(parameters->parity_info.masks), &(parameters->parity_info.number_masks)) != 0 ) { 
			PetscPrintf(PETSC_COMM_WORLD,"Error parsing parity information file.\n");
			return 1;
		}
		parameters->parity_info.number_parity_sectors = (int)pow((double)2,(double)parameters->parity_info.number_masks);
		//print_parity_information(parameters);
	}
	
	if ( parameters->use_conservation_sectors == PETSC_TRUE ) {
		if ( read_parity_masks(parameters->conservation_masks, &(parameters->parity_info.masks), &(parameters->parity_info.number_masks)) != 0 ) { 
			PetscPrintf(PETSC_COMM_WORLD,"Error parsing conservation information file.\n");
			return 1;
		}
		//print_conservation_information(parameters);
	}
	
	if ( parameters->use_parity_sectors == PETSC_TRUE ) {
		search_for_parity_sectors(&parameters->hamiltonian, parameters);
	} else if ( parameters->use_conservation_sectors == PETSC_TRUE ) {
		search_for_conservation_sectors(&parameters->hamiltonian, parameters);
	} 
	
	if ( parameters->use_rotation_invariance == PETSC_TRUE ) {
		read_rotation_information(parameters);
	}
	
	read_expectation_values_file(parameters);
	
	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "get_rank_chunk_size"
/*! \brief This function calculates the size of the portion of a distributed data structure that the process with rank number 'rank' out of 'size' processes owns.*/ 
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

#undef __FUNCT__
#define __FUNCT__ "get_max_chunk_size"
/*! \brief This function calculates the maximum size of a distributed data structure that any process out of 'size' processes owns.*/ 
PetscInt get_max_chunk_size(PetscInt n, PetscInt size) {
	return (PetscInt)ceil((double)n / (double)size);
}


#undef __FUNCT__
#define __FUNCT__ "get_rank_chunk_start"
/*! \brief This function calculates the starting index of the portion of a distributed data structure that the process with rank number 'rank' out of 'size' processes owns.*/ 
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

#undef __FUNCT__
#define __FUNCT__ "get_range_with_item"
/*! \brief This function calculates the rank of the process out of 'size' processes that the item with index 'item' has in a distributed data structure of size 'n'.*/ 
PetscInt get_range_with_item(PetscInt n, PetscInt item, PetscInt size){
	//we will use a divide an conquer approach here
	int left = 0 ;
	int right = size ;
	
	while ( right - left > 1 ) {
		if ( item >= get_rank_chunk_start(n,left + ((right-left)/2),size) ) {
			left = left + ((right-left)/2);
		} else {
			right = left + ((right-left)/2) ;
		}
	}
	return left;
}


#undef __FUNCT__
#define __FUNCT__ "dec2bin"
/*! \brief This function converts the number 'a' to the character array of length 'len' containing the binary representation of 'a'.*/ 
char *dec2bin(uPetscInt a, int len, char *num) {
	//char *num ;
	int i;
	uPetscInt tmp;

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

#undef __FUNCT__
#define __FUNCT__ "bin2dec"
/*! \brief This function converts the binary representation stored in 'bnum' to its decimal representation.*/ 
uPetscInt bin2dec(char *bnum) {
	int len, i;
	uPetscInt num;
	
	num = 0; 
	len = strlen(bnum);
	
	for ( i = 0; i < len ; i++ ){
		if ( bnum[len-1-i] == '1' ){
			num += (uPetscInt)pow((double)2,(double)i);
		}
	}
	return num;
}

#undef __FUNCT__
#define __FUNCT__ "output_memory_usage"
/*! \brief This function outputs the memory usage tracked by PETSc.*/
void output_memory_usage(){
  PetscLogDouble mem, tot;
 
  PetscMallocGetMaximumUsage(&mem);

  MPI_Allreduce(&mem,&tot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  PetscFPrintf(PETSC_COMM_WORLD,stderr,"Maximum memory used by process: %f (%fMB)\n",mem,mem/(PetscLogDouble)(1024*1024));
  PetscFPrintf(PETSC_COMM_WORLD,stderr,"Total maximum memory used by all process: %f (%fMB)\n",tot,tot/(PetscLogDouble)(1024*1024));

  PetscMallocGetCurrentUsage(&mem);

  MPI_Allreduce(&mem,&tot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  PetscFPrintf(PETSC_COMM_WORLD,stderr,"Current memory used by process: %f (%fMB)\n",mem,mem/(PetscLogDouble)(1024*1024));
  PetscFPrintf(PETSC_COMM_WORLD,stderr,"Total current memory used by all process: %f (%fMB)\n",tot,tot/(PetscLogDouble)(1024*1024));

}

#undef __FUNCT__
#define __FUNCT__ "find_in_sorted_array"
/*! \brief This function finds the index of the item 'item' in the sorted array 'array' which is of length 'len'. Returns -1 if not found.*/
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


#undef __FUNCT__
#define __FUNCT__ "insert_into_sorted_list"
/*! \brief This function inserts a new item into a sorted list of PetscInt. It returns 1 if the item is successfully inserted and 0 if the item already exists in the list. */
int insert_into_sorted_list(list<uPetscInt> *sorted_list, uPetscInt item) {
	//PetscPrintf(PETSC_COMM_WORLD,"empty %d, size %d, true %d\n", sorted_list->empty(), sorted_list->size(),true);
	if ( sorted_list->empty() ) { 
		//PetscPrintf(PETSC_COMM_WORLD,"List address %d, item %d\n",sorted_list,item);
		sorted_list->push_back(item);
		//PetscPrintf(PETSC_COMM_WORLD,"List size %d, empty %d, first %d\n", sorted_list->size(), sorted_list->empty(), sorted_list->front());
		return 1;
	}
	list<uPetscInt>::iterator iter = sorted_list->begin();
	while ( iter != sorted_list->end() && item > *iter ) {
		iter++;
	}
	if (iter != sorted_list->end() && item == *iter) {
		return 0;
	} else {
		//PetscPrintf(PETSC_COMM_WORLD,"Inserting %d\n",item);
		sorted_list->insert(iter,item);
		return 1;
	}
}

#undef __FUNCT__
#define __FUNCT__ "trim"
/*! \brief This function trims whitespace from the start and end of a string.*/
int trim(const string str, string *trimmed) {
	size_t start = str.find_first_not_of(" \t\n\r");
	if(start == string::npos) return 1;
	trimmed->assign(str,start,str.find_last_not_of(" \t\n\r") - start + 1);
	return 0; 
}

#undef __FUNCT__
#define __FUNCT__ "parse_adjacency_info_file"
/*! \brief This function parses the xml file containing the adjacency information for the lattice.*/
int parse_adjacency_info_file(struct parameter_struct *parameters, char *adjacency_file){
	string filename(adjacency_file);
	PetscErrorCode ierr;
	
	TiXmlDocument doc(filename);
	
	bool loadOkay = doc.LoadFile();
	
	if ( !loadOkay )
	{
		PetscPrintf(PETSC_COMM_WORLD,"Could not load adjacency file '%s'. Error='%s'. Exiting.\n", adjacency_file, doc.ErrorDesc() );
		return 1;
	}
	
	TiXmlHandle docHandle( &doc );
	
	TiXmlElement* element = docHandle.FirstChild("EDGES").Element();
	TiXmlElement* element2;
	if ( element)
	{	
		parameters->adjacency_info.number_masks = atoi(element->Attribute("number"));
		if ( parameters->adjacency_info.number_masks <= 0 ){
			PetscPrintf(PETSC_COMM_WORLD,"Invalid number of edges %d.\n", parameters->adjacency_info.number_masks);
			return 1;
		}
		ierr = PetscMalloc( parameters->adjacency_info.number_masks * sizeof(uPetscInt), &parameters->adjacency_info.masks);CHKERRQ(ierr);
		for ( int i = 0 ; i < parameters->adjacency_info.number_masks ; i++ ) {
			parameters->adjacency_info.masks[i] = 0;
			element2 = docHandle.FirstChild("EDGES").Child("EDGE",i).Element();
			if ( element2 ) {
				parameters->adjacency_info.masks[i] += ((uPetscInt)1) << (atoi(element2->Attribute("from"))-1);
				parameters->adjacency_info.masks[i] += ((uPetscInt)1) << (atoi(element2->Attribute("to"))-1);
			} else {
				PetscPrintf(PETSC_COMM_WORLD,"Edge %d does not exist.\n",i);
				return 1;
			}
		}
	}
	
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "deallocate_adjacency_info"
/*! \brief This function deallocates space allocated while parsing the adjacency file.*/
int deallocate_adjacency_info(struct parameter_struct *parameters) {
	if ( parameters->nn_exclusion == PETSC_TRUE) {
		PetscErrorCode ierr = PetscFree(parameters->adjacency_info.masks);CHKERRQ(ierr);
	}
	return 0 ;
}

#undef __FUNCT__
#define __FUNCT__ "check_adjacency"
/*! \brief This function checks the adjacency information for the given index and returns false if it contains adjacent particles.*/
PetscTruth check_adjacency(struct parameter_struct *parameters, uPetscInt index) {
	for ( int i = 0 ; i < parameters->adjacency_info.number_masks; i++ ) {
		if ( (index & parameters->adjacency_info.masks[i] ) == parameters->adjacency_info.masks[i] ) {
			return PETSC_FALSE;
		}
	}
	return PETSC_TRUE;
}

#undef __FUNCT__
#define __FUNCT__ "susy_staggered_subspace_check"
/*! \brief This function checks the adjacency information for the given index and returns false if it contains adjacent particles.*/
PetscTruth susy_staggered_subspace_check(struct parameter_struct *parameters, uPetscInt index) {
	if ( check_adjacency(parameters,index) == PETSC_FALSE ) {
		return PETSC_FALSE;
	}

	/*int count = 0; 
	for ( int i = 0 ; i < parameters->number_particles / 3 ; i++ ) {
		if ( (index & (((uPetscInt)1) << 3*i) ) == 0 && (index & (((uPetscInt)1) << 3*i + 2) ) == 0 ) {
			count ++; 
			if ( count > 1 ) return PETSC_FALSE; 
		}
	}*/
	
	int count = 0; 
	for ( int i = 0 ; i < parameters->number_particles / 3 ; i++ ) {
		if ( (index & (((uPetscInt)1) << (3*i+1)) ) != 0 ) {
			count ++; 
			if ( count > 2 ) return PETSC_FALSE; 
		}
	}
	
	return PETSC_TRUE;
}


#undef __FUNCT__
#define __FUNCT__ "dummy_basis_check"
/*! \brief This function is a dummy function which will always return true. When the standard basis is being used the relvant function reference points to this function.*/
PetscTruth dummy_basis_check(struct parameter_struct *parameters, uPetscInt index){
	return PETSC_TRUE;
}	

#undef __FUNCT__
#define __FUNCT__ "save_matrix"
/*! \brief This function writes the matrices in binary format to disk.*/
int save_matrix(struct parameter_struct *parameters,Mat *A, const char *filename){
	PetscViewer viewer;
	PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
	if ( parameters->save_matrix_ascii ) {
		PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
		if ( parameters->save_matrix_matlab ) {
			PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
		}
	} else {
		PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);
	}
	MatView(*A,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "tokenize_string"
/*! \brief This function tokenizes a string according to the given delimination characters.*/
int tokenize_string(string line, string delim, vector<string> *parts) {
	string::size_type comma_pos;

	while ( (comma_pos = line.find_first_of(delim)) != line.npos) { 
		if ( comma_pos > 0 ) {
			parts->push_back(line.substr(0,comma_pos));	
		}
		line = line.substr(comma_pos + 1);
	}
	if ( line.length() > 0 ) {
		parts->push_back(line);
	}
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "str2scalar"
/*! \brief This function converts a string to a scalar. */
PetscScalar str2scalar(string str) {
	//first thing to do is break it into each of the parts by dividing it at the plus or minus signes.
	//Must make sure that there is not an "e" before the plus or minus or that it is not at the start of 
	//a number. 

	size_t pos, start=0;
	PetscScalar value = 0.0;
	PetscScalar minus = 1.0;
	string buf_str;
	trim(str,&buf_str);
	str.assign(buf_str);
	
	while ( (pos = str.find_first_of("+-",start) ) != str.npos ) {
		if ( pos == 0 || ( str[pos-1] == 'e' ) ) { //if the sign is at the start or comes after an e then we set the start position past this position and go on.
			start = pos+1;
		} else {  //if neither of the above conditions were satisfied then the sign is separting terms.
			int imag = 0;
			char buf[MAX_STRING];
			char *ptr = buf;
			for ( size_t i = 0 ; i < pos ; i++ ) {
				if ( str[i] == 'i' ) {
					imag = 1; 
				} else {
					*ptr = str[i];
					ptr++;
				}
			}
			*ptr = '\0';

		    #if defined(PETSC_USE_COMPLEX)
			if (imag) {
				value += (minus * PETSC_i * (PetscScalar)atof(buf));
			} 
            #endif

            if (!imag) {
				value += (minus * (PetscScalar)atof(buf));
			}
			
			if ( str[pos] == '-' ) {
				minus = -1.0;
			} else {
				minus = 1.0;
			}
			str = str.substr(pos+1);
			trim(str,&buf_str);
			str.assign(buf_str);
		}
	}
	
	if ( str.length() > 0 ) {
		int imag = 0;
		char buf[MAX_STRING];
		char *ptr = buf;
		for ( size_t i = 0 ; i < str.length() ; i++ ) {
			if ( str[i] == 'i' ) {
				imag = 1; 
			} else {
				*ptr = str[i];
				ptr++;
			}
		}
		*ptr = '\0';
		#if defined(PETSC_USE_COMPLEX)
		if (imag) {
			value += (minus * PETSC_i * (PetscScalar)atof(buf));
		} 
		#endif
		
		if (!imag) {
			value += (minus * (PetscScalar)atof(buf));
		}	
	}

	return value;
}




#undef __FUNCT__
#define __FUNCT__ "find_valid_configs"
/*! \brief This function recursively builds up a list of valid configurations that satisfy the nearest neightbour exclusion in parallel.*/
int find_valid_configs(struct parameter_struct *parameters, ofstream *basis_file, PetscInt *local_reps){
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
	while ( count < parameters->size && set_size < parameters->number_particles ) {
		next_configs = new vector<uPetscInt>(); 
		next_fillings = new vector<int>();
		int i = 0 ;
		int new_count = 0; 
		while ( i < count ) {		
				
			if ( (parameters->use_conservation_sectors == PETSC_FALSE) || ((parameters->real_conservation_sector - (*current_fillings)[i]) <= (parameters->number_particles - set_size))) {
				uPetscInt tmp;
				tmp = (*current_configs)[i] + (((uPetscInt)1) << set_size) ;
				if ( (*parameters->simple_basis_check)(parameters,tmp)) { //the simple basis check should check the adjacency NN_EXLUSION IS SPECIFIED.
					next_configs->push_back(tmp);
					next_fillings->push_back((*current_fillings)[i]+1);
					new_count++;
				}	
			}
			
			if ( (parameters->use_conservation_sectors == PETSC_FALSE) || ((*current_fillings)[i] <= parameters->real_conservation_sector ) ) {
				next_configs->push_back((*current_configs)[i]);
				next_fillings->push_back((*current_fillings)[i]);
				new_count ++;
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
	
	if ( count >= parameters->size && set_size < parameters->number_particles ) {
		//cout << "Count " << count << " rank " << parameters->rank << " size " << get_rank_chunk_size((PetscInt)count,parameters->rank,parameters->size) << " and start " <<  get_rank_chunk_start((PetscInt)count,parameters->rank,parameters->size)  << endl;
		for ( PetscInt idx = 0; idx <  get_rank_chunk_size((PetscInt)count, parameters->rank,parameters->size) ; idx++ ) {
			//cout << "Rank " << parameters->rank << " getting started with " << (*current_configs)[get_rank_chunk_start((PetscInt)count, parameters->rank,parameters->size)+idx] << endl;
			find_valid_configs_with_this_subset(parameters,set_size,(*current_fillings)[get_rank_chunk_start((PetscInt)count, parameters->rank,parameters->size)+idx],  (*current_configs)[get_rank_chunk_start((PetscInt)count, parameters->rank,parameters->size)+idx],basis_file,local_reps);	
		}		
	} 
	
	count = 0;

	//cout << "Rank : " << parameters->rank << " number : " << *local_reps << endl;
	
	return 0;
}



#undef __FUNCT__
#define __FUNCT__ "find_valid_configs_with_this_subset"
/*! \brief This function adds a site and calls the same function for each possible value of the additional site if it is valid. */
int find_valid_configs_with_this_subset(struct parameter_struct *parameters, int set_size, int set_filling, uPetscInt set, ofstream *basis_file, PetscInt *local_reps ) {
	//first check the stop condition if set size is already the size of the full graph or there is not enough space left to reach required filling.
	uPetscInt set2;
	if ( set_size == parameters->number_particles ) { //if we have reach the full size of the lattice
		if ( (parameters->use_conservation_sectors == PETSC_FALSE || set_filling == parameters->real_conservation_sector) && ((*parameters->representative_check)(parameters,set) == PETSC_TRUE) ) {
			
			set2 = set;
			basis_file->write((char*)&set2,sizeof(uPetscInt));	
			//cout << "Found " << set << endl;
			(*local_reps)++;
			return 0;
		} else {
			return 0;
		}
	} else if ( parameters->use_conservation_sectors && (parameters->real_conservation_sector - set_filling) > (parameters->number_particles - set_size) ) {
		return 0 ;
	} else if (  parameters->use_conservation_sectors && set_filling == parameters->real_conservation_sector && ((*parameters->representative_check)(parameters,set) == PETSC_TRUE) ) {
		set2 = set;
		basis_file->write((char*)&set2,sizeof(uPetscInt));	
		(*local_reps)++;
		//cout << "Found " << set << endl;
		return 0;
	}
	
	//have two branches one with 0 at site set_size and one with 1 at set_size
	//now add a 1 to position set_size, check if valid and increment filling.
	uPetscInt new_set = set + (((uPetscInt)1) << set_size) ; 
	if ( (*parameters->simple_basis_check)(parameters,new_set) ) {
		find_valid_configs_with_this_subset(parameters, set_size+1 , set_filling+1, new_set,basis_file,local_reps);
	}
	
	
	//do not need to check the addition of a 0 as will already be full.
	find_valid_configs_with_this_subset(parameters, set_size+1 , set_filling, set,basis_file,local_reps);
	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "str_to_mask"
/*! \brief This function parses the xml file containing the adjacency information for the lattice.*/
uPetscInt str_to_mask(const char *buf){
	uPetscInt mask = 0;
	int len = strlen(buf);
	for ( size_t i = 0 ; i < strlen(buf) ; i++ ) {
		if ( buf[len-1-i] == '1' ) {
			mask += two_pow(i);
		}
	}
	
	return mask;
}

#undef __FUNCT__
#define __FUNCT__ "parse_prod_wf_overlap_file"
/*! \brief This function parses the xml file containing the adjacency information for the lattice.*/
int parse_prod_wf_overlap_file(struct parameter_struct *parameters, char *prod_wf_file){
	string filename(prod_wf_file);
	PetscErrorCode ierr;
	
	TiXmlDocument doc(filename);
	
	bool loadOkay = doc.LoadFile();
	
	if ( !loadOkay )
	{
		PetscPrintf(PETSC_COMM_WORLD,"Could not load product wavefunction file '%s'. Error='%s'. Exiting.\n", prod_wf_file, doc.ErrorDesc() );
		return 1;
	}
	
	TiXmlHandle docHandle( &doc );
	
	TiXmlElement* element = docHandle.FirstChild("PRODUCT_TERMS").Element();
	TiXmlElement *element2, *element3, *element4, *element5;
	if ( element )
	{	
		parameters->prod_wf_info.number_terms = atoi(element->Attribute("number_terms"));
		if ( parameters->prod_wf_info.number_terms <= 0 ){
			PetscPrintf(PETSC_COMM_WORLD,"Invalid number of terms %d.\n", parameters->prod_wf_info.number_terms);
			return 1;
		}
		ierr = PetscMalloc( parameters->prod_wf_info.number_terms * sizeof(struct wf_term), &parameters->prod_wf_info.terms);CHKERRQ(ierr);
		for ( int i = 0 ; i < parameters->prod_wf_info.number_terms ; i++ ) { 
			element2 = docHandle.FirstChild("PRODUCT_TERMS").Child("PRODUCT_TERM",i).Element();	
			if ( element2 ) {
				struct wf_term *tmp_term = &parameters->prod_wf_info.terms[i];
				tmp_term->number_subterms = atoi(element2->Attribute("number_subterms"));
				tmp_term->coefficient = atof(element2->Attribute("coefficient"));
				ierr = PetscMalloc( tmp_term->number_subterms * sizeof(struct wf_subterm), &tmp_term->subterms);CHKERRQ(ierr);
				for ( int j = 0 ; j <  tmp_term->number_subterms ; j++ ) {
					element3 = docHandle.FirstChild("PRODUCT_TERMS").Child("PRODUCT_TERM",i).Child("PRODUCT_SUBTERM",j).Element();
					if ( element3 ) {
						struct wf_subterm *tmp_subterm = &tmp_term->subterms[j];
						tmp_subterm->coefficient = atof(element3->Attribute("coefficient"));
						tmp_subterm->number_groups = atoi(element3->Attribute("number_groups"));
						ierr = PetscMalloc( tmp_subterm->number_groups * sizeof(struct wf_group), &tmp_subterm->groups);CHKERRQ(ierr);	
						for ( int k = 0 ; k < tmp_subterm->number_groups ; k++ ) {
							struct wf_group *tmp_group = &tmp_subterm->groups[k];
							element4 = docHandle.FirstChild("PRODUCT_TERMS").Child("PRODUCT_TERM",i).Child("PRODUCT_SUBTERM",j).Child("GROUP",k).Element();
							if ( element4 ){ 
								tmp_group->number_confs = atoi(element4->Attribute("number_confs"));	
								tmp_group->mask = str_to_mask(element4->Attribute("mask"));
								ierr = PetscMalloc( tmp_group->number_confs * sizeof(struct wf_conf), &tmp_group->confs);CHKERRQ(ierr);	
								for ( int l = 0 ; l < tmp_group->number_confs ; l++ ) {
									//struct wf_conf *tmp_conf = &tmp_group->confs[l];
									element5 = docHandle.FirstChild("PRODUCT_TERMS").Child("PRODUCT_TERM",i).Child("PRODUCT_SUBTERM",j).Child("GROUP",k).Child("CONF",l).Element();
									if ( element5 ) {
										tmp_group->confs[l].coefficient = atof(element5->Attribute("coefficient"));
										//tmp_group->confs[l].mask = str_to_mask(element5->Attribute("mask"));
										tmp_group->confs[l].mask = str_to_mask(element5->GetText());
									} else { PetscPrintf(PETSC_COMM_WORLD,"WF Conf %d not found.\n",l); return 1; }
								}
							} else { PetscPrintf(PETSC_COMM_WORLD,"WF Group %d not found.\n",k); return 1; }
						}
					} else { PetscPrintf(PETSC_COMM_WORLD,"WF Prod subterm %d not found.\n",j); return 1; }
				} //end loop over subterms
			} else { PetscPrintf(PETSC_COMM_WORLD,"WF term %d not found.\n",i); return 1; } //end if terms
		} //end loop over terms
	} //end if product terms
				
	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "deallocate_prod_wf_info"
/*! \brief This function deallocates space allocated while parsing the prod wf file.*/
int deallocate_prod_wf_info(struct parameter_struct *parameters) {
	if ( parameters->prod_wf_overlap == PETSC_TRUE) {
		PetscErrorCode ierr;
		for ( int i = 0 ; i < parameters->prod_wf_info.number_terms ; i++ ) {
			for ( int j = 0 ; j < parameters->prod_wf_info.terms[i].number_subterms ; j++ ) {
				for ( int k = 0 ; k < parameters->prod_wf_info.terms[i].subterms[j].number_groups ; k++ ) {
					ierr = PetscFree(parameters->prod_wf_info.terms[i].subterms[j].groups[k].confs);CHKERRQ(ierr);
				} //end loop over groups
				ierr = PetscFree(parameters->prod_wf_info.terms[i].subterms[j].groups);CHKERRQ(ierr);
			} //end loop over subterms
			ierr = PetscFree(parameters->prod_wf_info.terms[i].subterms);CHKERRQ(ierr);
		} //end loop over terms
		ierr = PetscFree(parameters->prod_wf_info.terms);CHKERRQ(ierr);
	}
	return 0 ;
}


#undef __FUNCT__
#define __FUNCT__ "create_prod_wf_state"
/*! \brief This function creates a product wavefunction according to details specified in the input file.*/
int create_prod_wf_state(struct parameter_struct *parameters, Vec *x) {
	//create parallel vector for the wave function.
	Vec term_vec;
	PetscErrorCode ierr = VecCreateMPI(parameters->mat_solver_comm,PETSC_DECIDE,parameters->current_basis_size,x);CHKERRQ(ierr);
	ierr = VecDuplicate(*x,&term_vec);
	VecZeroEntries(*x); 
	uPetscInt conf;
	PetscReal value;
	struct prod_wf_information *PPS = &parameters->prod_wf_info;
	PetscTruth valid;
	uPetscInt masked_conf;
	PetscInt number_matching_confs;
	PetscScalar acc_coefficient;
	
	for ( PetscInt term = 0 ; term < PPS->number_terms; term++) {
		VecZeroEntries(term_vec); // zero all the entries of the vector for this term.
		acc_coefficient = PPS->terms[term].coefficient; // set the coefficient.
		for ( PetscInt subterm = 0 ; subterm < PPS->terms[term].number_subterms; subterm++) {
			acc_coefficient *= PPS->terms[term].subterms[subterm].coefficient;
			//loop over local configurations and check which ones are product states.
			if ( parameters->nn_exclusion ) { //if using nearest neighbour exclusion then will need to use stored array. Otherwise can work out the configurations on the fly.
				for ( PetscInt i = 0  ; i < get_rank_chunk_size(parameters->current_basis_size, parameters->rank, parameters->size); i++ ) {
					conf = parameters->get_full_basis_index(parameters,parameters->basis_info.local_reps_array[i]);
					value = 1.0;
					//now we loop over the product groups 
					valid = PETSC_TRUE;
					number_matching_confs = 0; 
					for ( PetscInt group = 0 ; group < PPS->terms[term].subterms[subterm].number_groups; group++ ) { 
						masked_conf = conf & PPS->terms[term].subterms[subterm].groups[group].mask ;
						for ( PetscInt conf_idx = 0 ; conf_idx < PPS->terms[term].subterms[subterm].groups[group].number_confs; conf_idx++ ) {
							struct wf_conf *tmp_conf = &PPS->terms[term].subterms[subterm].groups[group].confs[conf_idx];
							if ( masked_conf == tmp_conf->mask ) {
								value = value * PetscRealPart(tmp_conf->coefficient);
								number_matching_confs++;
								break; //break out of the loop over confs.
							}
						} //end the loop over confs in the group.
					} //end loop over groups.
					
					if ( number_matching_confs == PPS->terms[term].subterms[subterm].number_groups ) {
						//VecSetValue(*x,get_rank_chunk_start(parameters->current_basis_size,parameters->rank,parameters->size)+i,acc_coefficient*value,ADD_VALUES);
						VecSetValue(term_vec,get_rank_chunk_start(parameters->current_basis_size,parameters->rank,parameters->size)+i,value,ADD_VALUES);
					}
				}
			} else {
				PetscPrintf(PETSC_COMM_WORLD,"Feature to use PPS wavefunctions for models without NN exclusion not implmented.\n");
				return 0;
			}
		} //end the loop over subterms
		VecAssemblyBegin(term_vec);
		VecAssemblyEnd(term_vec);
		VecNormalize(term_vec,&value); // normalise.
		VecScale(term_vec,acc_coefficient); // after normalising scale using coefficent.
		
		PetscScalar *from_vec, *to_vec;
		ierr = VecGetArray(term_vec,&from_vec);CHKERRQ(ierr);
		ierr = VecGetArray(*x,&to_vec);CHKERRQ(ierr);
		PetscInt local_size;
		ierr = VecGetLocalSize(*x,&local_size);
		for ( PetscInt i = 0 ; i < local_size; i++ ) {
			to_vec[i] += from_vec[i];
		}
		ierr = VecRestoreArray(term_vec,&from_vec);CHKERRQ(ierr);
		ierr = VecRestoreArray(*x,&to_vec);CHKERRQ(ierr);
		VecNormalize(*x,&value); //normalise after each term.
	} //end the loop over terms. 
	
	/*int *local_counts;
	int max_group = 0;
	for ( int i = 0 ; i < parameters->prod_wf_info.number_groups ; i++ ){
		if (parameters->prod_wf_info.group_sizes[i] > max_group) {
			max_group = parameters->prod_wf_info.group_sizes[i];
		}
	}
	max_group++; //need an extra space for fully populated.
	ierr = PetscMalloc(max_group * sizeof(int),&local_counts);CHKERRQ(ierr);	
	
	
	for ( int term = 0 ; term < number_terms ; term++ ) { 
		PetscScalar coefficient = 1.0;
		if ( parameters->prod_wf_info.extra_terms  ) {
			coefficient = parameters->prod_wf_info.coefficients[term];
		}
		
		//loop over local configurations and check which ones are product states.
		if ( parameters->nn_exclusion ) { //if using nearest neighbour exclusion then will need to use stored array. Otherwise can work out the configurations on the fly.
			for ( PetscInt i = 0  ; i < get_rank_chunk_size(parameters->current_basis_size, parameters->rank, parameters->size); i++ ) {
				conf = parameters->get_full_basis_index(parameters,parameters->basis_info.local_reps_array[i]);
				value = 1.0;
				int group_count = 0;
				//initialise counters to 0.
				for (int k = 0 ; k < max_group ; k++ ) {
					local_counts[k] = 0;
				}
				int previous_count = -1;
				PetscTruth still_good = PETSC_TRUE;
				for ( int j = 0 ; j < parameters->prod_wf_info.number_groups ; j++ ) { //loop over the groups of sites.
					int count = 0;
					for ( int k = 0 ; k < parameters->prod_wf_info.group_sizes[j] ; k++ ) { //the loop over the sites in each group.
						if ( (conf & parameters->prod_wf_info.masks[j][k]) != 0 ) {  //if the configuration has a fermion at that site then increment the counter and take the sign into account.
							count++;
							value *= (PetscReal)parameters->prod_wf_info.signs[j][k];
						}
					} //end loop over sites in group
					local_counts[count]++; //increment the relevant counter. 
					if ( count == 1 && previous_count == -1) {
						still_good = PETSC_FALSE;	
					}
					if ( count == 0 && ( previous_count != -1 && previous_count != 1 ) ) {
						still_good = PETSC_FALSE;	
					}
					if ( count == 2 && ( previous_count != -1 && previous_count != 1 ) ) {
						still_good = PETSC_FALSE;	
					}
					previous_count = count;
					if (count == 1 ) { 
						group_count++; //if the wavefunction contained only one fermion from the group the increment the group counter. 
					} //else break;
				} //end loop over groups
				if ( still_good ) value *= -1.0;
				if ( parameters->prod_wf_info.extra_terms == PETSC_FALSE ) {
					if ( group_count == parameters->prod_wf_info.number_groups ) { //if the group count is right then add this configuration. 
						VecSetValue(*x,get_rank_chunk_start(parameters->current_basis_size,parameters->rank,parameters->size)+i,coefficient*value,INSERT_VALUES);
					}
				} else {
					PetscTruth valid = PETSC_TRUE;
					for ( int n = 0 ; n < parameters->prod_wf_info.number_term_counts; n++ ) {
						if ( local_counts[n] != parameters->prod_wf_info.term_counts[term][n] ) {
							valid = PETSC_FALSE;
						}	
					}
					if ( valid && parameters->prod_wf_info.exclusions ) {
						for ( int h = 0 ; h < parameters->prod_wf_info.number_exclusions ; h++ ) {
							int new_count = 0;
							for ( int g = 0 ; g < parameters->prod_wf_info.number_exclusion_sites[h] ; g++ ) {
								if ( (conf & parameters->prod_wf_info.exclusion_masks[h][g]) != 0 ) {
									new_count++;
								}								
							}
							if ( new_count == parameters->prod_wf_info.number_exclusion_sites[h]){
								value *= (PetscReal)parameters->prod_wf_info.exclusion_coeffs[h];
							}
						}
					}
					if (valid){
						VecSetValue(*x,get_rank_chunk_start(parameters->current_basis_size,parameters->rank,parameters->size)+i,coefficient*value,INSERT_VALUES);
					}
				}
			} // end loop over configs
			VecAssemblyBegin(*x);
			VecAssemblyEnd(*x);
		} else {
			for ( PetscInt i = 0  ; i < get_rank_chunk_size(parameters->current_basis_size, parameters->rank, parameters->size); i++ ) {
				conf = parameters->get_full_basis_index(parameters,i);
				value = 1.0;
				int group_count = 0;
				for ( int j = 0 ; j < parameters->prod_wf_info.number_groups ; j++ ) { //loop over the groups of sites.
					int count = 0;
					for ( int k = 0 ; k < parameters->prod_wf_info.group_sizes[j] ; k++ ) {
						if ( conf & parameters->prod_wf_info.masks[j][k] != 0 ) { 
							count++;
							value *= (PetscReal)parameters->prod_wf_info.signs[j][k];
						}
					} //end loop over sites in group
					if (count == 1 ) {
						group_count++;
					} else break;
				} //end loop over groups
				if ( group_count == parameters->prod_wf_info.number_groups ) {
					VecSetValue(*x,get_rank_chunk_start(parameters->current_basis_size,parameters->rank,parameters->size)+i,value,INSERT_VALUES);
				}
			} // end loop over configs
			VecAssemblyBegin(*x);
			VecAssemblyEnd(*x);
		}	
		VecNormalize(*x,&value);
	} //end the loop over terms. 
	PetscFree(local_counts);*/
	
	//VecView(*x,PETSC_VIEWER_STDOUT_WORLD);
	//VecView(*x,PETSC_VIEWER_STDOUT_WORLD);
	
	return 0 ;
}
