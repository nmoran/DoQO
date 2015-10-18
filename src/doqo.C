/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 * 
 * Main file of C version of Eigen Solver
 * 
 * Version: $Id$ 
 * 
*/
/*! \file */ 

/*! \mainpage DoQO Code Documentation
 *
 * \section intro_sec Introduction
 *
 * DoQO is a code that performs diagonalisation of quantum operators. See the README file 
 * for details of the installation procedure. See the files tab for a list of the source files 
 * with links to documentation about each.
 */

#include "doqo.h" 
#include "hamiltonian.h"
#include "utils.h"
#include "solver_funcs.h" 
#include "tasks.h"
#include "symmetry.h"
#include "output.h"
#include "momentum.h"
#include "expectation_value_funcs.h"
#include <iostream>


using namespace std;


#undef __FUNCT__
#define __FUNCT__ "main"

/*! The main function where DoQO starts execution.
*/
int main( int argc, char **argv ) {
	//Mat            	A;           /* operator matrix */
	//EPS            	eps;         /* eigenproblem solver context */
	PetscErrorCode 	ierr;
	PetscReal		starttime, endtime ,t1, t2, t3 ,t4;
	struct 			parameter_struct parameters;
	PetscInt 		num_parity_sectors, num_momentum_sectors, num_conservation_sectors, number_rotation_sectors; 
	Vec 			starting_vec,y;
	struct 			solver_results       results;

	//---------------------------------------------------------------------------------------------------------------------
	// Initialisation section
	//---------------------------------------------------------------------------------------------------------------------
	
	//Initialise SLEPC and PETSC and set rank and size.
	initialise(argc, argv, &parameters); 

	//Get starting time.
	starttime = MPI_Wtime();

	//print DoQO banner and version details.
	write_doqo_heading(&parameters);

	//Parse the input file and assign parameters. 
	if ( parse_main_input_file(&parameters) != 0 )
	  {
	    //Finalise SLEPC 
	    ierr = SlepcFinalize();CHKERRQ(ierr);
	    return -1;
	  }
	  		  
	//This file will parse the other input files referenced in the main input file.
	if ( read_input_files(&parameters) != 0 )
	  {
		ierr = SlepcFinalize();CHKERRQ(ierr);
		return -1;
	  }
	  
	if ( parameters.calculate_expectation_values == PETSC_FALSE ) 
	  {
	    //this writes an xml output file with the parameters of the run and lists the output files that will be written.
	    write_main_xml_output(&parameters);
			    
	    //work out how many sectors there are taking into consideration parity, conservation and momentum
	    if ( parameters.use_parity_sectors == PETSC_TRUE ) {
		    num_parity_sectors = parameters.parity_info.number_parity_sectors;
		    if ( parameters.parity_info.num_relevant_sectors > 0 ) num_parity_sectors = parameters.parity_info.num_relevant_sectors;
	    } else { 
		    num_parity_sectors = 1;
	    }
	    
	    if ( parameters.use_conservation_sectors == PETSC_TRUE ) {
		    num_conservation_sectors = parameters.parity_info.number_parity_sectors;
		    if ( parameters.parity_info.num_relevant_sectors > 0 ) num_conservation_sectors = parameters.parity_info.num_relevant_sectors;
	    } else { 
		    num_conservation_sectors = 1;
	    }
	    
	    if ( parameters.use_momentum_sectors == PETSC_TRUE ) { 
		    num_momentum_sectors = parameters.momentum_info.num_sectors[0] * parameters.momentum_info.num_sectors[1];
		    if ( parameters.momentum_info.num_relevant_sectors > 0 ) num_momentum_sectors = parameters.momentum_info.num_relevant_sectors ; 
	    } else {
		    num_momentum_sectors = 1;
	    }
	    
	    if ( parameters.use_rotation_invariance == PETSC_TRUE ) { 
		    if ( parameters.rotation_number_relevant_sectors > 0 ) 
		    {
			    number_rotation_sectors = parameters.rotation_number_relevant_sectors;
		    } else number_rotation_sectors = parameters.number_rotation_sectors;
	    } else {
		    number_rotation_sectors = 1;
	    }
			    
	    if ( parameters.verbosity >= 1 && parameters.rank == 0) {
		    PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------\n"); 
		    stringstream ss;
		    ss.str("");ss << "Number of tasks           :\t\t" << parameters.tasks.number_tasks << endl; 
		    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
		    if (parameters.use_parity_sectors ) {
			    ss.str("");ss << "Number of parity sectors  :\t\t" << num_parity_sectors << endl; 
			    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
		    }
		    if (parameters.use_conservation_sectors ) {
			    ss.str("");ss << "Number of filling sectors :\t\t" << num_conservation_sectors << endl; 
			    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
		    }
		    if (parameters.use_momentum_sectors ) {
			    ss.str("");ss << "Number of momentum sectors:\t\t" << num_momentum_sectors << endl; 
			    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
		    }
		    if (parameters.use_spectral_flow ) {
			    ss.str("");ss << "Number of spectral flow points:\t\t" << parameters.spectral_info.number_relevant_points[0]*parameters.spectral_info.number_relevant_points[1]<< endl; 
			    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
		    }
		    if (parameters.use_rotation_invariance ) 
		    {		
			    ss.str("");ss << "Number of rotation sectors:\t\t" << number_rotation_sectors << endl; 
			    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
		    }
	    ss.str("");ss << "Number of diagonalisations:\t\t" << parameters.tasks.number_tasks * num_parity_sectors * num_conservation_sectors * num_momentum_sectors * parameters.spectral_info.number_relevant_points[0]*parameters.spectral_info.number_relevant_points[1]*number_rotation_sectors<< endl; 
	    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
	    }
	    
	    PetscInt total_diags = parameters.tasks.number_tasks * num_parity_sectors * num_conservation_sectors * num_momentum_sectors * parameters.spectral_info.number_relevant_points[0]*parameters.spectral_info.number_relevant_points[1] * number_rotation_sectors ;
	    
	    //---------------------------------------------------------------------------------------------------------------------
	    // Main loop 
	    //---------------------------------------------------------------------------------------------------------------------
	    //loop through each task in the task file 
	    for ( parameters.current_task = parameters.tasks.start_task ; parameters.current_task < parameters.tasks.number_tasks ; parameters.current_task++ ) {
	    
		    apply_task_parameters(&parameters.hamiltonian, &parameters.tasks, (int)parameters.current_task ) ; //apply the parameters relevent to each task
		    
		    //loop through sectors and perform diagonalisation
		    for ( parameters.parity_sector = 0 ; parameters.parity_sector < num_parity_sectors ; parameters.parity_sector++ ) { 
			    for ( parameters.conservation_sector = 0 ; parameters.conservation_sector < num_conservation_sectors ; parameters.conservation_sector++ ) { 
				    for ( parameters.spectral_info.current_point[1] = 0 ; parameters.spectral_info.current_point[1] < parameters.spectral_info.number_relevant_points[1]; parameters.spectral_info.current_point[1]++) { 
					    for ( parameters.spectral_info.current_point[0] = 0 ; parameters.spectral_info.current_point[0] < parameters.spectral_info.number_relevant_points[0]; parameters.spectral_info.current_point[0]++) { 
						    for ( parameters.momentum_sector = 0 ; parameters.momentum_sector < num_momentum_sectors ; parameters.momentum_sector++ ) { 
							    for ( parameters.rotation_current_sector = 0 ; parameters.rotation_current_sector < number_rotation_sectors;  parameters.rotation_current_sector++ ) {
								    PetscInt diag_num = (parameters.current_task) * (num_parity_sectors*num_conservation_sectors*num_momentum_sectors*parameters.spectral_info.number_relevant_points[0]*parameters.spectral_info.number_relevant_points[1]*number_rotation_sectors) 
													    + (parameters.parity_sector) * ( num_conservation_sectors*num_momentum_sectors*parameters.spectral_info.number_relevant_points[0]*parameters.spectral_info.number_relevant_points[1]*number_rotation_sectors) 
													    + (parameters.conservation_sector) * (num_momentum_sectors*parameters.spectral_info.number_relevant_points[0]*parameters.spectral_info.number_relevant_points[1]*number_rotation_sectors)
													    + (parameters.spectral_info.current_point[1]) * (num_momentum_sectors*parameters.spectral_info.number_relevant_points[0]*number_rotation_sectors)
													    + (parameters.spectral_info.current_point[0]) * (num_momentum_sectors*number_rotation_sectors)
													    + parameters.momentum_sector*(number_rotation_sectors)
													    + parameters.rotation_current_sector+1;
								    if ( parameters.verbosity >= 1 && parameters.rank == 0) {
									    PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------\n"); 
									    stringstream ss;
									    ss.str("");ss << "Diagonalisation: " << diag_num << " of " <<  total_diags << endl;  
									    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
									    PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------\n"); 
								    }
															    
								    Mat *A;
								    PetscMalloc(sizeof(Mat),&A);
								    if ( parameters.distribute == PETSC_TRUE ) {						
									    t3 = MPI_Wtime();
									    parameters.current_basis_size = MatrixDistribute(A,&parameters,&parameters.hamiltonian);
									    t4 = MPI_Wtime();
									    if ( parameters.verbosity >= 3 ) PetscPrintf(PETSC_COMM_WORLD, "Time taken to build matrix: %lf\n", t4 - t3 ) ;
									    
				    
									    if ( parameters.save_matrix == PETSC_TRUE ) { 
										    char mat_filename[MAX_STRING];
										    sprintf(mat_filename, "%s_%s.mat" , parameters.output_prefix, parameters.current_output_additional_prefix);
										    save_matrix(&parameters,A,mat_filename); 
									    }
									    
									    if ( parameters.multiply_test == PETSC_TRUE && parameters.in_communicator == PETSC_TRUE ) { //the multiply test performs a simple multiply operation with a vector of ones
										    ierr = VecCreateMPI(parameters.mat_solver_comm, PETSC_DECIDE, parameters.current_basis_size,&starting_vec);CHKERRQ(ierr);
										    VecDuplicate(starting_vec,&y);
										    VecSet(starting_vec,1.0);
										    t1 = MPI_Wtime();
										    MatMult(*A,starting_vec,y);
										    t2 = MPI_Wtime();
										    if ( parameters.verbosity >= 3 ) PetscPrintf(parameters.mat_solver_comm,"Time taken for matrix multiplication test %f.\n", t2 - t1 ); 
										    VecDestroy(&starting_vec);
										    VecView(y,PETSC_VIEWER_STDOUT_WORLD );
										    VecDestroy(&y);
									    } else if ( parameters.solver == PETSC_TRUE ) { //&& parameters.in_communicator == PETSC_TRUE) {
										    results.time_taken = 0;
										    t1 = MPI_Wtime();
										    EPS  *eps;
										    PetscMalloc(sizeof(EPS),&eps);
										    if ( get_rank_chunk_size(parameters.current_basis_size,parameters.rank,parameters.size)  > 0 ) {
											    if  ( parameters.current_basis_size > 1 )  {
												    PetscInt neigs = parameters.eigPairs;
												    if (neigs > (parameters.current_basis_size )) {
													    neigs = parameters.current_basis_size; 
												    }
												    
												    createEigenSolver(&parameters,A, eps, neigs, parameters.tol, parameters.max_its); // create an eigensolver
							    
												    if ( parameters.benchmark == PETSC_TRUE ) {
													    ierr = VecCreateMPI(parameters.mat_solver_comm, PETSC_DECIDE, parameters.current_basis_size,&starting_vec);CHKERRQ(ierr);
													    VecSet(starting_vec,1.0);
													    #if SLEPC_VERSION_MAJOR == 3 &&  SLEPC_VERSION_MINOR > 0 
													    EPSSetInitialSpace(*eps,1,&starting_vec);
													    #else
													    EPSSetInitialVector(*eps,starting_vec);
													    #endif

												    }
										    
												    solveEigenProblem(eps,&parameters,&results); //solve the system
												    if ( parameters.verbosity >= 1 ) PetscPrintf(PETSC_COMM_WORLD,"Diagonalsiation finished, %d values converged.\n", (int)results.eigenvalues_converged);
											    } else if ( parameters.current_basis_size == 1 ) {
												    results.eigenvalues_converged = 1 ;
												    results.iterations_taken = 0 ; 
											    } else if ( parameters.current_basis_size == 0 ) {
												    results.eigenvalues_converged = 0 ; 
												    results.degenerate_bands = 0 ;
											    }
											    if ( results.eigenvalues_converged == 0 ) {
												    results.degenerate_bands = 0 ;
											    }
					    
										    
											    //ierr = MatDestroy(A);CHKERRQ(ierr);
											    if ( results.eigenvalues_converged > 0 ) {
												    //assign space to store eigen values and vortex numbers
												    ierr = PetscMalloc(results.eigenvalues_converged * sizeof(PetscReal),&results.eigenvalues);CHKERRQ(ierr);
												    ierr = PetscMalloc(results.eigenvalues_converged * sizeof(PetscReal),&results.errors);CHKERRQ(ierr);
												    if ( parameters.prod_wf_overlap ) {
													    ierr = PetscMalloc(results.eigenvalues_converged * sizeof(PetscReal),&results.overlaps);CHKERRQ(ierr);
												    }
												    if ( parameters.other_wf_overlap ) {
													    ierr = PetscMalloc(results.eigenvalues_converged * sizeof(PetscReal),&results.other_overlaps);CHKERRQ(ierr);
												    }
												    //retrieve the results and calculate the vortex numbers 
												    if ( parameters.current_basis_size > 1 && results.eigenvalues_converged > 0 ){
													    getSolverResults(eps,&parameters,&results);
												    } else if ( parameters.current_basis_size == 1  ) {
													    if ( parameters.rank == 0 ){
														    PetscInt row = 0;
														    PetscInt ncols;
														    const PetscInt *cols;
														    const PetscScalar *values;
														    MatGetRow(*A,row,&ncols,&cols,&values);
														    if ( ncols > 0 ) {
															    results.eigenvalues[0] = PetscRealPart(values[0]);
															    results.errors[0] = 0.0 ;
															    PetscPrintf(PETSC_COMM_WORLD, "%.16lf\t%.16lf\n" ,results.eigenvalues[0], results.errors[0]);
														    }
														    MatRestoreRow(*A,row,&ncols,&cols,&values);
													    }
													    ierr = PetscMalloc(2*sizeof(PetscInt),&results.degeneracies);CHKERRQ(ierr);
													    results.degeneracies[0] = 0;
													    results.degeneracies[1] = 1; 
													    results.degenerate_bands = 1 ;
													    if (  parameters.calculate_correlations == PETSC_TRUE ) {
														    ierr = PetscMalloc(sizeof(PetscScalar**),&results.correlations);CHKERRQ(ierr);
														    ierr = PetscMalloc(parameters.number_ops*sizeof(PetscScalar*),&results.correlations[0]);CHKERRQ(ierr);
														    for ( int i = 0 ; i < parameters.number_ops ; i++ ) {
															    ierr = PetscMalloc(sizeof(PetscScalar),&results.correlations[0][i]);CHKERRQ(ierr);
															    Mat Op;	
															    ierr = MatCreate(parameters.mat_solver_comm,&Op);CHKERRQ(ierr);
															    ierr = MatSetSizes(Op, get_rank_chunk_size(parameters.current_basis_size,parameters.rank,parameters.size), PETSC_DECIDE, parameters.current_basis_size ,parameters.current_basis_size);CHKERRQ(ierr);
															    if ( parameters.use_hermitian == PETSC_TRUE ) {
																    ierr = MatSetType(Op,MATMPISBAIJ);CHKERRQ(ierr);
																    ierr = MatSetOption(Op,MAT_GETROW_UPPERTRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
															    } else { 
																    ierr = MatSetType(Op,MATMPIAIJ);CHKERRQ(ierr);
															    }
															    MatSetFromOptions(Op);CHKERRQ(ierr);
															    MatMPIAIJSetPreallocation(Op,1,PETSC_NULL,0,PETSC_NULL);CHKERRQ(ierr);
															    set_matrix_elements(&Op,&parameters,&parameters.ops_data[i],parameters.basis_info.local_reps_array,parameters.basis_info.ending_idx_for_proc,NULL,NULL,0);
															    PetscInt row = 0;
															    PetscInt ncols;
															    const PetscInt *cols;
															    const PetscScalar *values;
															    MatGetRow(Op,row,&ncols,&cols,&values);
															    if ( ncols > 0 ) {
																    results.correlations[0][i][0] = values[0];
															    }
															    MatRestoreRow(Op,row,&ncols,&cols,&values);
															    MatDestroy(&Op);
														    }
													    }
												    } // end basis size == 1 
											    } // end if eigenvalues_converged > 0 
											    
											    write_xml_data(&parameters,&results);
											    
											    if ( results.eigenvalues_converged > 0 ) {
												    ierr = PetscFree(results.eigenvalues);CHKERRQ(ierr);
												    ierr = PetscFree(results.errors);CHKERRQ(ierr);
												    ierr = PetscFree(results.degeneracies);CHKERRQ(ierr);
												    if ( parameters.prod_wf_overlap ) {
													    ierr = PetscFree(results.overlaps);CHKERRQ(ierr);
												    }
												    if ( parameters.other_wf_overlap ) {
													    ierr = PetscFree(results.other_overlaps);CHKERRQ(ierr);
												    }
												    if ( parameters.calculate_correlations == PETSC_TRUE ) {
													    for ( int i = 0 ; i < results.degenerate_bands ; i++ ) {
														    for ( int j = 0 ; j < parameters.number_ops ; j++ ) {
															    ierr = PetscFree(results.correlations[i][j]);CHKERRQ(ierr);
														    }
														    ierr = PetscFree(results.correlations[i]);CHKERRQ(ierr);
													    }
													    ierr = PetscFree(results.correlations);CHKERRQ(ierr);
												    }
											    }
											    
											    t2 = MPI_Wtime();
											    if ( parameters.verbosity >= 3 ) PetscPrintf(parameters.mat_solver_comm,"Time taken for diagonalisation: %f\n", t2 - t1 ); 
											    
											    //Destroy the instance of the eigen solver and matrix.
											    if ( parameters.in_communicator == PETSC_TRUE && parameters.current_basis_size > 1 ) {
												    ierr = EPSDestroy(eps);CHKERRQ(ierr);
												    //ierr = EPSDestroy(*eps);CHKERRQ(ierr);
											    }
											    PetscFree(eps);
											    
											    //looks like the EPSDestroy takes care of destroying the initial vector too. 
											    /*if ( parameters.benchmark == PETSC_TRUE ) {
												    VecDestroy(starting_vec);
											    }*/ 
						    
											    if ( parameters.in_communicator == PETSC_TRUE && parameters.current_basis_size > 0 ) { //if a matrix was created 
												    ierr = MatDestroy(A);CHKERRQ(ierr);
												    //ierr = MatDestroy(*A);CHKERRQ(ierr);
												    if ( parameters.calculate_correlations == PETSC_TRUE || parameters.prod_wf_overlap ) {
													    PetscFree(parameters.basis_info.local_reps_array);
													    PetscFree(parameters.basis_info.ending_idx_for_proc);
												    }
												    
												    if ( MPI_Comm_free(&(parameters.mat_solver_comm)) != MPI_SUCCESS ) {
													    PetscPrintf(PETSC_COMM_SELF,"Rank %d: Freeing communicator not successful.\n", parameters.rank);
												    }									
											    }
										    } else { //else for if current_basis_size == 0
											    results.eigenvalues_converged = 0 ; 
											    results.degenerate_bands = 0 ;
											    write_xml_data(&parameters,&results);
										    } //end if this process is part of the solving group.
									    } //end if solve flag is true 
								    } // end if distribute flag is true
								    PetscFree(A);
								    MPI_Barrier(PETSC_COMM_WORLD);
							    } //end loop over rotation sectors
						    } //end loop over momentum sectors 
					    } //end loop over spectral flow points in direction 0
				    } //end loop over spectral flow points in direction 1 
			    } //end loop over conservation sectors 
		    } // end loop over parity sectors 
	    } //end for loop over tasks
	  }
	else 
	  {
	    PetscPrintf(PETSC_COMM_WORLD, "Calculating expection values for %d vectors in filling sector %d and rotation sector %d.\n", (int)parameters.exp_vals_info.number_vectors, (int)parameters.exp_vals_info.filling_sector, (int)parameters.exp_vals_info.rotation_sector);
	    for ( int i = 0 ; i < (int)parameters.exp_vals_info.number_spaces ; i++ )
	      {
		for ( int j = 0 ; j < (int)parameters.exp_vals_info.number_vectors_per_space[i] ; j++ )
		  {
		    PetscPrintf(PETSC_COMM_WORLD, "%s\n", parameters.exp_vals_info.vectors[i][j].c_str());
		  }
	      }
	    parameters.parity_sector = 0;
	    parameters.spectral_info.current_point[0] = 0;
	    parameters.spectral_info.current_point[1] = 0;
	    parameters.momentum_sector = 0;
	    
	    if ( parameters.exp_vals_info.filling_sector != -1 ) 
	      {
		parameters.conservation_sector = parameters.exp_vals_info.filling_sector ;
		parameters.real_conservation_sector = parameters.exp_vals_info.filling_sector;
		parameters.parity_info.num_relevant_sectors = 0 ;
	      }
	    else
	      {
		parameters.conservation_sector = 0;
	      }
	    
	    if ( parameters.exp_vals_info.rotation_sector != -1 ) 
		parameters.rotation_current_sector = parameters.exp_vals_info.rotation_sector ;
	    else
		parameters.rotation_current_sector = 0;
	    	    	    
	     //loop through each task in the task file 
	    for ( parameters.current_task = parameters.tasks.start_task ; parameters.current_task < parameters.tasks.number_tasks ; parameters.current_task++ ) 
	      {	    
		apply_task_parameters(&parameters.hamiltonian, &parameters.tasks, (int)parameters.current_task ) ; //apply the parameters relevent to each task 
		Mat *A;
		PetscMalloc(sizeof(Mat),&A);
		if ( parameters.distribute == PETSC_TRUE )
		  {						
		    t3 = MPI_Wtime();
		    parameters.current_basis_size = MatrixDistribute(A,&parameters,&parameters.hamiltonian);
		    t4 = MPI_Wtime();
		    if ( parameters.verbosity >= 3 ) PetscPrintf(PETSC_COMM_WORLD, "Time taken to build matrix: %lf\n", t4 - t3 ) ;		    
		    calculate_expectation_values(&parameters, A);		    		    
		  }	
		PetscFree(A);
		MPI_Barrier(PETSC_COMM_WORLD);
	      }
	  }
      	
	//---------------------------------------------------------------------------------------------------------------------
	// Cleanup stage
	//---------------------------------------------------------------------------------------------------------------------
	deallocate_task_space(&parameters.tasks,&parameters.hamiltonian);

	deallocate_matrix_data(&parameters.hamiltonian); 
	
	if ( parameters.calculate_correlations == PETSC_TRUE ) {
		for ( int i = 0 ; i < parameters.number_ops ; i++ ) {
			deallocate_matrix_data(&parameters.ops_data[i]);
		}
		PetscFree(parameters.ops_data);
	}
	
	if ( parameters.prod_wf_overlap ) {
		deallocate_prod_wf_info(&parameters);
	}
	
	
	if ( parameters.use_momentum_sectors == PETSC_TRUE && parameters.momentum_info.num_relevant_sectors > 0 ) {
		ierr = PetscFree(parameters.momentum_info.relevant_sectors);CHKERRQ(ierr);
	}
	
	if ( (parameters.use_parity_sectors == PETSC_TRUE || parameters.use_conservation_sectors == PETSC_TRUE ) && parameters.parity_info.num_relevant_sectors > 0 ) {
		ierr = PetscFree(parameters.parity_info.relevant_sectors);CHKERRQ(ierr);
	}
	if ( parameters.use_parity_sectors == PETSC_TRUE  || parameters.use_conservation_sectors == PETSC_TRUE ) { //free up space used by parity structure 
			free_parity_info(&parameters);
	}
	
	if ( parameters.use_spectral_flow ) {
		if ( parameters.spectral_info.number_relevant_points[0] <  parameters.spectral_info.number_points[0] ) {
			PetscFree(parameters.spectral_info.relevant_points[0]);
		}
		if ( parameters.spectral_info.number_relevant_points[1] <  parameters.spectral_info.number_points[1] ) {
			PetscFree(parameters.spectral_info.relevant_points[1]);
		}
	}
	
	if ( parameters.use_rotation_invariance ) {
		deallocate_rotation_memory(&parameters);
	}
	
	deallocate_adjacency_info(&parameters);
	deallocate_expectaction_value_space(&parameters);   
	    
	endtime = MPI_Wtime();
	
	if ( parameters.verbosity >= 3 ) {
	    PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------\n"); 
	    PetscPrintf(PETSC_COMM_WORLD, "Overall time taken: %lf\n" , endtime - starttime); 
	}
	
	// output_memory_usage();
	
	//Finalise SLEPC
	ierr = SlepcFinalize();CHKERRQ(ierr);
	
	return 0;
}


