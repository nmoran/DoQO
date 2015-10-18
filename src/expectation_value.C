/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala.
 * 
 * Main file of C version of Eigen Solver
 * 
 * Version: $Id$ 
 * 
*/

#include "solver.h" 
#include "hamiltonian.h"
#include "utils.h"
#include "solver_funcs.h" 
#include "tasks.h"
#include "symmetry.h"
//#include <iostream>

//using namespace std;



#undef __FUNCT__
#define __FUNCT__ "main"
int main( int argc, char **argv ) {
    Mat            	A, WP;           /* operator matrix */
    EPS            	eps;         /* eigenproblem solver context */
    PetscReal		*eigenvalues, *vortices_num; 
    PetscReal 		*error;
    PetscErrorCode 	ierr;
    int				number_eigenvalues, i, maxits,task;
    PetscReal		starttime, endtime ,t1, t2, t3 ,t4;
    struct parameter_struct parameters;
    PetscInt 		num_sectors, sector; 
    Vec 		starting_vec,y,y2;
    char                vec_file[MAX_STRING];
    PetscViewer 	viewer;
    PetscScalar 	result, *tmp;
    PetscTruth	        flg;
    //Vec order, lookup ; //two vectors used for sz symmetry sectors

    //---------------------------------------------------------------------------------------------------------------------
    // Initialisation section
    //---------------------------------------------------------------------------------------------------------------------

    
    //Initialise SLEPC and PETSC and set rank and size.
    initialise(argc, argv, &parameters); 


    parameters.use_shell = PETSC_FALSE; 
	sprintf(parameters.op_meta_filename,"expectation_output.meta");

#ifdef DEBUG
    PetscFPrintf(PETSC_COMM_WORLD,stderr,"Finised initialisation.\n");
#endif

	//default maximum number iterations
    maxits = 5000; 

    //get starting time
    starttime = MPI_Wtime();
   	
    strcpy(parameters.hamiltonian.filename,"");
    strcpy(vec_file,"");
    strcpy(parameters.tasks.filename,"");
    PetscOptionsGetString(PETSC_NULL,"-op",parameters.hamiltonian.filename,MAX_STRING,&flg);
    PetscOptionsGetString(PETSC_NULL,"-vec",vec_file,MAX_STRING,&flg);
    PetscOptionsGetString(PETSC_NULL,"-taskfile",parameters.tasks.filename, MAX_STRING,&flg);
 
    if ( strcmp(parameters.hamiltonian.filename,"") == 0 || strcmp(vec_file,"") == 0  || strcmp(parameters.tasks.filename,"") == 0 ) { 
      PetscPrintf(PETSC_COMM_WORLD,"You must specify an operator, vector and taskfile.\n");
      ierr = SlepcFinalize();CHKERRQ(ierr);
      return -1;
    }

    //PetscPrintf(PETSC_COMM_SELF,"Vec file <%s>\n", vec_file);


#ifdef DEBUG
     PetscFPrintf(PETSC_COMM_WORLD,stderr,"Finised parsing arguments.\n");
#endif
	
    if ( setup_operators(&parameters) > 0 ) { //setup the operators so they are ready to be distributed
      PetscFinalize();
      return 1 ;
    }

#ifdef DEBUG
	PetscFPrintf(PETSC_COMM_WORLD,stderr,"Finised setting up operators.\n");
#endif


    if ( parse_task_file(&parameters.tasks, &parameters, &parameters.hamiltonian) > 0 ) { // parse the task file 
      PetscFinalize();
      return 1 ; 
    }
    
    //    write_meta_data_start(&parameters, argc, argv);
    //write_data_start(&parameters);
    
    //print_task_info(&parameters.tasks, &parameters.hamiltonian); // print out information about the different tasks
    //print_matrix_info(&parameters.hamiltonian); // print information about the matrix in question
    
    write_task_info(&parameters);
    write_matrix_info(&parameters,&parameters.hamiltonian);
    

    //ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,parameters.hamiltonian.numberBasisStates,&y);CHKERRQ(ierr);
    //VecDuplicate(y,&y2);
	
    PetscPrintf(PETSC_COMM_WORLD,"Attempting to load vector from %s with size %d.\n", vec_file , parameters.numberBasisStates);

    //read in vector
    PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, vec_file ,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecLoad(viewer,VECMPI,&y);CHKERRQ(ierr);
    PetscViewerDestroy(viewer);
    
    VecDuplicate(y,&y2);
  
    /*VecGetArray(y,&tmp);
	
    PetscPrintf(PETSC_COMM_WORLD,"Array elements %lf, %lf ... \n", tmp[0],tmp[1]);
	
    VecRestoreArray(y,&tmp);*/

    
    PetscFPrintf(PETSC_COMM_WORLD,stderr,"Starting main loop.\n");
#ifdef DEBUG
    PetscFPrintf(PETSC_COMM_WORLD,stderr,"Starting main loop.\n");
#endif
    
    //---------------------------------------------------------------------------------------------------------------------
    // Main loop 
    //---------------------------------------------------------------------------------------------------------------------
    //loop through each task in the task file 
    for ( task = 0 ; task < parameters.tasks.number_tasks ; task++ ) {
      PetscFPrintf(PETSC_COMM_WORLD,stderr,"Task %d.\n",task);
      parameters.current_task = task;
      apply_task_parameters(&parameters.hamiltonian, &parameters.tasks, task ) ; //apply the parameters relevent to each task
	
	  num_sectors = 1;
	
	
		
		//print_matrix_info(&parameters.hamiltonian);
		
		if ( parameters.distribute == PETSC_TRUE ) {
		
			for ( sector = 0; sector < num_sectors ; sector++ ) {
			
				t3 = MPI_Wtime();
				
				if ( parameters.use_sz_symmetries == PETSC_TRUE) { 
					MatrixDistributeSector(&A,&parameters,&parameters.hamiltonian,sector); // distribute the matrix amongst all the processors
				} else {
					MatrixDistribute(&A, &parameters, &parameters.hamiltonian); // distribute the matrix amongst all the processors
				}
				t4 = MPI_Wtime();
				PetscFPrintf(PETSC_COMM_WORLD,stderr, "Distribute matrix time: %lf\n", t4 - t3 ) ;
				//write_meta_time(&parameters, t4 - t3 , "Matrix Distribution time");
				
				if ( parameters.solver == PETSC_TRUE ) { // if not doing multiply test
					MatMult(A,y,y2);
					VecDot(y,y2,&result);					
				} //end if solve flag is true 
			} // loop over sectors 
		} // end if distribute flag is true
		
		
		PetscPrintf(PETSC_COMM_WORLD,"Expectation value: %lf\n", result);

	} //end for loop over tasks
	
	//---------------------------------------------------------------------------------------------------------------------
	// Cleanup stage
	//---------------------------------------------------------------------------------------------------------------------
	deallocate_task_space(&parameters.tasks,&parameters.hamiltonian);

	//deallocate_matrix_data(&parameters.hamiltonian); 
	
    	endtime = MPI_Wtime();
	
	
    
    //PetscFPrintf(PETSC_COMM_WORLD,stderr, "Time taken: %lf\n" , endtime - starttime); 
    //write_meta_time(&parameters, endtime - starttime , "Overall Time Taken");
    //Finalise SLEPC 
    ierr = SlepcFinalize();CHKERRQ(ierr);
    
    return 0;
}
