/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 * 
 * Solver functions.
 * 
 * Version: $Id$ 
 * 
*/

#include "solver_funcs.h"


#undef __FUNCT__
#define __FUNCT__ "createEigenSolver"
/*! \brief Function to create an eigen solver for the operator A and setup parameters.*/
int createEigenSolver(struct parameter_struct *parameters,Mat *A, EPS *eps, PetscInt eigenPairs, PetscReal tol, int maxits){
    PetscErrorCode ierr;
    //Create eigensolver context
    ierr = EPSCreate(parameters->mat_solver_comm,eps);CHKERRQ(ierr);

    //ierr = EPSSetType(*eps,EPSKRYLOVSCHUR);CHKERRQ(ierr);
    
    //Set operators. In this case, it is a standard eigenvalue problem
    ierr = EPSSetOperators(*eps,*A,PETSC_NULL);CHKERRQ(ierr);
    
    //set the problem type to a hermitian problem
    ierr = EPSSetProblemType(*eps, EPS_HEP);CHKERRQ(ierr);

    //Set the number of eigenPairs to try and get. Let Petsc Decide on the max subspace size. 
#if SLEPC_VERSION_MAJOR == 3 
    ierr = EPSSetDimensions(*eps,eigenPairs,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
#else
    ierr = EPSSetDimensions(*eps,eigenPairs,PETSC_DECIDE);CHKERRQ(ierr);
#endif

    //Set solver to retrieve smalled eigen values 
    ierr = EPSSetWhichEigenpairs(*eps, EPS_SMALLEST_REAL);CHKERRQ(ierr);

    //Set the tollerance to use and the maximum number of iterations
    ierr = EPSSetTolerances(*eps,tol, maxits); CHKERRQ(ierr); 
    
    //Set solver parameters at runtime
    ierr = EPSSetFromOptions(*eps);CHKERRQ(ierr);
    
    return 0 ;
}


#undef __FUNCT__
#define __FUNCT__ "solverEigenProblem"
/*! \brief Function that calls the solver to solve the eigen problem. */ 
PetscInt solveEigenProblem(EPS *eps, struct parameter_struct *parameters, struct solver_results *results){
//#if SLEPC_VERSION_MAJOR == 3 
//    const EPSType type;
//#else
    EPSType type;
//#endif

    PetscErrorCode ierr;

	PetscReal start_time = MPI_Wtime();
    //Solve the eigensystem
    ierr = EPSSolve(*eps);CHKERRQ(ierr);
	PetscReal end_time = MPI_Wtime();
	results->time_taken = end_time - start_time;

    stringstream ss ;
	if ( parameters->verbosity >= 1 ) {
		// Get number of converged approximate eigenpairs
		ierr = EPSGetConverged(*eps,&results->eigenvalues_converged);CHKERRQ(ierr);
		ss.str("");ss << "Number of converged eigenpairs: " << results->eigenvalues_converged << endl  ; 
		ierr = PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());CHKERRQ(ierr);
	}

	if ( parameters->verbosity >= 2 ) {
    //Retrieve the number of iterations it took and print it 
		ierr = EPSGetIterationNumber(*eps, &results->iterations_taken);CHKERRQ(ierr);
		ss.str(""); ss << "Number of iterations: " << results->iterations_taken << endl ; 
		ierr = PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());CHKERRQ(ierr);
	}

	if ( parameters->verbosity >= 2 ) {
		//Optional: Get some information from the solver and display it
		ierr = EPSGetType(*eps,&type);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Method used: %s\n",type);CHKERRQ(ierr);
	}
    
    if ( parameters->verbosity >= 3 ) {
		EPSConvergedReason reason;
		ierr = EPSGetConvergedReason(*eps,&reason);CHKERRQ(ierr);
		ss.str("");ss << "Convergence reason: " << reason << endl  ; 
		ierr = PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());CHKERRQ(ierr);
   	}

    //ierr = EPSGetDimensions(*eps,&nev,PETSC_NULL);CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %d\n",nev);CHKERRQ(ierr);
    //ierr = EPSGetTolerances(*eps,&tol,&maxit);CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",tol,maxit);CHKERRQ(ierr);

    return results->eigenvalues_converged;
}

#undef __FUNCT__
#define __FUNCT__ "getEigenValues"
/*! \brief Function retrieves the eigenvalues that converged. */
int getEigenValues(EPS *eps, PetscReal *eigenvalues, PetscInt nconv){
    PetscScalar kr, ki;
    int i; 
    PetscErrorCode ierr;

    //if any values converged 
    if (nconv>0 ) {

		for ( i = 0 ; i < nconv ; i++ ) {
			ierr = EPSGetEigenpair(*eps,i,&kr,&ki,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
			eigenvalues[i] = PetscRealPart(kr); 
#else
			eigenvalues[i] = kr; 
#endif
			//PetscPrintf(PETSC_COMM_WORLD,"%lf\n",eigenvalues[i]);
		}
    }
    return 0 ; 
}

#undef __FUNCT__
#define __FUNCT__ "getEigenVector"
/*! \brief Function retrieves the nth eigenvector. */
int getEigenVector(EPS *eps, Vec eigenvector_real, Vec eigenvector_imag, PetscInt n){
    PetscScalar kr, ki;
    PetscErrorCode ierr;

    ierr = EPSGetEigenpair(*eps,n,&kr,&ki,eigenvector_real,eigenvector_imag);CHKERRQ(ierr);
	
    return 0 ; 
}


#undef __FUNCT__
#define __FUNCT__ "cleanupEigenSolver"
/*! \brief Function cleans up after the eigensolver by calling the proper routines to deallocate space. */
int cleanupEigenSolver(struct parameter_struct *parameters, EPS *eps, Mat *A){
    PetscErrorCode ierr;
    //struct shell_matrix_context *context;
    //void *ctx; 
    
    ierr = EPSDestroy(eps);CHKERRQ(ierr);
    
    /*if ( parameters->use_shell == PETSC_TRUE ) {
    	ierr = MatShellGetContext( *A, &ctx ); CHKERRQ(ierr);
    	context = (struct shell_matrix_context*)ctx;
    	cleanup_matrix_free(context);
    }
    
    ierr = MatDestroy(*A);CHKERRQ(ierr);*/
    
    return 0;
}

#undef __FUNCT__
#define __FUNCT__ "degeneracy_analysis"
/*! \brief Function that works out how many degenerate bands there are and how many eigenvalues are in each. */
int degeneracy_analysis(struct solver_results *results, PetscReal tol){
    int i, band ;
    PetscErrorCode ierr;

    results->degenerate_bands = 1;
    
    //first count the degenerate bands so we can allocate how much space to give
    for ( i = 0 ; i < results->eigenvalues_converged-1  ; i++ ) {
		if ( fabs(results->eigenvalues[i] - results->eigenvalues[i+1]) > tol ) {
			 results->degenerate_bands++; 
		}
    }
    
    //allocate the space
    ierr = PetscMalloc((results->degenerate_bands + 1 ) * sizeof(PetscInt), &results->degeneracies);CHKERRQ(ierr);

    //now loop through eigen values again setting each element in the degeneracies array to the first index of each band.
    results->degeneracies[0] = 0; 
    band = 1 ;
    for ( i = 0 ; i < results->eigenvalues_converged-1 ; i++ ) { 
      if ( fabs(results->eigenvalues[i] - results->eigenvalues[i+1]) > tol ) {
			results->degeneracies[band++] = i+1;
      } 
    }
    results->degeneracies[band] = results->eigenvalues_converged;

    return 0;
}

#undef __FUNCT__
#define __FUNCT__ "getSolverResults"
/*! \brief Function that gets results back from completed solver and sets the relevant fields in the results structure. */
int getSolverResults(EPS *eps, struct parameter_struct *parameters, struct solver_results *results) {
	PetscScalar     kr, ki, val, vali;
	PetscErrorCode ierr;
	PetscViewer  	viewer;
	//Mat 	B; 
	Vec x,xi ;//, y , yi,vr,vi ; 
	//EPS eps2;
	char state_filename[MAX_STRING];
	
	
	//retrieve an array of these eigen values
	getEigenValues(eps,results->eigenvalues,results->eigenvalues_converged);
	
	for ( PetscInt i = 0 ; i < results->eigenvalues_converged ; i++ ) {
		EPSGetErrorEstimate(*eps,i,&results->errors[i]);
	}
	
	degeneracy_analysis(results, parameters->deg_tol);
	
	if ( parameters->verbosity >= 1 ) { 
		//print out eigenvalues
		for ( PetscInt i = 0 ; i < results->eigenvalues_converged ; i++ ) {
			PetscPrintf(PETSC_COMM_WORLD, "%.16lf\t%.16lf\n" ,results->eigenvalues[i], results->errors[i]);
		}
	}
	
	// write states to disk
	if ( parameters->save_states == PETSC_TRUE ) {
		ierr = VecCreateMPI(parameters->mat_solver_comm,PETSC_DECIDE,parameters->current_basis_size,&x);CHKERRQ(ierr);
		VecDuplicate(x,&xi);
		for ( PetscInt i = 0; i < results->eigenvalues_converged ; i++ ) {
			//sprintf(state_filename, "%s_TASK_%d_STATE_%d.vec" , parameters->output_prefix, parameters->current_task ,(int) i ) ;
			sprintf(state_filename, "%s.%s_state_%d.vec" , parameters->output_prefix, parameters->current_output_additional_prefix, (int)i ) ;
			PetscViewerCreate(parameters->mat_solver_comm,&viewer);
			if ( parameters->save_states_ascii ) {
				PetscViewerASCIIOpen(parameters->mat_solver_comm, state_filename,&viewer);
				if ( parameters->save_states_matlab ) {
					PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
				}
			} else {
				PetscViewerBinaryOpen(parameters->mat_solver_comm, state_filename,FILE_MODE_WRITE,&viewer);
			}
			getEigenVector(eps, x, xi, i);
			if ( parameters->save_states_real ) 
			  {
			    PetscScalar *vals;
			    PetscReal l_max = 0.0; 
			    PetscInt max_idx = 0; 
			    VecGetArray(x,&vals);
			    for ( int idx = 0 ; idx < get_rank_chunk_size(parameters->current_basis_size, parameters->rank, parameters->size); idx++ )
			      {
				if ( PetscAbsScalar(vals[idx]) > l_max ) 
				  {
				      l_max = PetscAbsScalar(vals[idx]);
				      max_idx = idx;
				  }
			      }
			    PetscReal *send_buf; 
			    ierr = PetscMalloc(sizeof(PetscReal)*3, &send_buf); CHKERRQ(ierr);
			    send_buf[0] = l_max;
			    send_buf[1] = PetscRealPart(vals[max_idx]);
			    send_buf[2] = PetscImaginaryPart(vals[max_idx]);
			    VecRestoreArray(x, &vals);
			    PetscReal *all_vals;			    
			    ierr = PetscMalloc(parameters->restricted_size*sizeof(PetscReal)*3, &all_vals); CHKERRQ(ierr);			    
			    
			    MPI_Allgather(send_buf, 3, MPI_DOUBLE, all_vals, 3, MPI_DOUBLE, parameters->mat_solver_comm);			    
			    l_max = 0.0;
			    for ( int idx = 0 ; idx < parameters->restricted_size ; idx++ ) 
			      {
				if ( all_vals[idx*3] > l_max )
				  {
				    l_max = all_vals[idx*3];
				    max_idx = idx*3;
				  }
			      }	 
			    PetscScalar max_val;
                #if defined(PETSC_USE_COMPLEX)
			    max_val = all_vals[max_idx +1] + PETSC_i * all_vals[max_idx +2];
                #else
			    max_val = all_vals[max_idx +1];
                #endif
			    VecScale(x, 1.0/max_val);
			    ierr = PetscFree(all_vals); CHKERRQ(ierr);
			    ierr = PetscFree(send_buf); CHKERRQ(ierr);			    
			    VecNormalize(x, &l_max);			    
			  }
			VecView(x,viewer);
			PetscViewerDestroy(&viewer);
		}
		ierr = VecDestroy(&x);CHKERRQ(ierr);
		ierr = VecDestroy(&xi);CHKERRQ(ierr);
	}

	Vec trial;	

	if ( parameters->other_wf_overlap ) {
		PetscViewer viewer;
		PetscViewerCreate(parameters->mat_solver_comm,&viewer);
		ierr = PetscViewerBinaryOpen(parameters->mat_solver_comm,parameters->other_wf_file,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
		//ierr = VecLoad(viewer,VECMPI,&trial);CHKERRQ(ierr);
		ierr = VecLoad(trial, viewer);CHKERRQ(ierr);
		PetscViewerDestroy(&viewer);
		PetscInt vec_size;
		VecGetSize(trial,&vec_size);
		
		Vec eigenstate,tmp; 
		ierr = VecCreateMPI(parameters->mat_solver_comm,PETSC_DECIDE,parameters->current_basis_size,&eigenstate);CHKERRQ(ierr);
		VecDuplicate(eigenstate,&tmp);
		
		if ( vec_size == parameters->current_basis_size ) {
			for ( PetscInt i = 0 ; i < results->eigenvalues_converged ; i++ ) {
				getEigenVector(eps, eigenstate, tmp, i);
				PetscScalar tmp_scalar;
				VecDot(trial,eigenstate,&tmp_scalar);
				results->other_overlaps[i] = PetscRealPart(tmp_scalar);
			}
		} else {
			if ( parameters->verbosity > 0 ) PetscPrintf(parameters->mat_solver_comm,"Size of vectors do not match up.\n"); 
		}
		
		VecDestroy(&tmp);
		VecDestroy(&eigenstate);
		VecDestroy(&trial);
		//will deallocate vector space later after correlations are calculated. 
	}


	//now calculate overlaps
	if ( parameters->prod_wf_overlap ) {

		create_prod_wf_state(parameters,&trial); //create the trial wavefunction
		Vec eigenstate,tmp; 
		ierr = VecCreateMPI(parameters->mat_solver_comm,PETSC_DECIDE,parameters->current_basis_size,&eigenstate);CHKERRQ(ierr);
		VecDuplicate(eigenstate,&tmp);
		for ( PetscInt i = 0 ; i < results->eigenvalues_converged ; i++ ) {
			getEigenVector(eps, eigenstate, tmp, i);
			PetscScalar tmp_scalar;
			VecDot(trial,eigenstate,&tmp_scalar);
			//results->overlaps[i] = PetscRealPart(tmp_scalar);
			results->overlaps[i] = PetscRealPart(PetscPowScalar(tmp_scalar,2.0));
		}
		VecDestroy(&tmp);
		VecDestroy(&eigenstate);
		
		if ( parameters->save_states ) {
			PetscViewer viewer;
			sprintf(state_filename, "%s.prod_wf.vec" , parameters->output_prefix) ;
			PetscViewerCreate(parameters->mat_solver_comm,&viewer);
			if ( parameters->save_states_ascii ) {
				PetscViewerASCIIOpen(parameters->mat_solver_comm, state_filename,&viewer);
				if ( parameters->save_states_matlab ) {
					PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
				}
			} else {
				PetscViewerBinaryOpen(parameters->mat_solver_comm, state_filename,FILE_MODE_WRITE,&viewer);
			}
			VecView(trial,viewer);
			PetscViewerDestroy(&viewer);
		}
		//will deallocate vector space later after correlations are calculated. 
	}
	
	
	//calculate the correlations which were specified in the input file.
	if ( parameters->calculate_correlations == PETSC_TRUE ) {
		//first need to allocate space for the correlations.
		ierr = PetscMalloc(results->degenerate_bands * sizeof(PetscScalar**),&results->correlations);CHKERRQ(ierr); //allocate the array of pointers to matrices.
		for ( PetscInt i = 0 ; i < results->degenerate_bands ; i++ ) {
			int dim = results->degeneracies[i+1] - results->degeneracies[i];
			ierr = PetscMalloc(parameters->number_ops * sizeof(PetscScalar*),&results->correlations[i]);CHKERRQ(ierr);
			for ( PetscInt j = 0 ; j < parameters->number_ops ; j++ ) {
				ierr = PetscMalloc(dim*dim*sizeof(PetscScalar),&results->correlations[i][j]);CHKERRQ(ierr);
			}
		}
		if ( parameters->prod_wf_overlap ) {
			ierr = PetscMalloc(parameters->number_ops * sizeof(PetscScalar),&results->trial_correlations);CHKERRQ(ierr);
		}
		
		
		for ( PetscInt j = 0 ; j < parameters->number_ops ; j++ ) {		
			Mat Op;	
			ierr = MatCreate(parameters->mat_solver_comm,&Op);CHKERRQ(ierr);
			ierr = MatSetSizes(Op, get_rank_chunk_size(parameters->current_basis_size,parameters->rank,parameters->size), PETSC_DECIDE, parameters->current_basis_size ,parameters->current_basis_size);CHKERRQ(ierr);
			if ( parameters->use_hermitian == PETSC_TRUE ) {
				ierr = MatSetType(Op,MATMPISBAIJ);CHKERRQ(ierr);
				ierr = MatSetOption(Op,MAT_GETROW_UPPERTRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
			} else { 
				ierr = MatSetType(Op,MATMPIAIJ);CHKERRQ(ierr);
			}
			MatSetFromOptions(Op);CHKERRQ(ierr);
			MatMPIAIJSetPreallocation(Op,1,PETSC_NULL,0,PETSC_NULL);CHKERRQ(ierr);
			set_matrix_elements(&Op,parameters,&parameters->ops_data[j],parameters->basis_info.local_reps_array,parameters->basis_info.ending_idx_for_proc,NULL,NULL,0);

			Vec tmp; //this is just a dummy vector to pass to the geteigenvector function. 
			ierr = VecCreateMPI(parameters->mat_solver_comm,PETSC_DECIDE,parameters->current_basis_size,&tmp);CHKERRQ(ierr);
			
			for ( PetscInt i = 0 ; i < results->degenerate_bands ; i++ ) {
				Vec *vecs ;

				//When complex is enabled all the imaginary parts will be in the one vector.
				int dim = results->degeneracies[i+1] - results->degeneracies[i];
				ierr = PetscMalloc(dim*sizeof(Vec),&vecs);CHKERRQ(ierr);
				
				for ( int k = 0 ; k < dim ; k++ ) {
					VecDuplicate(tmp,&vecs[k]);
					getEigenVector(eps, vecs[k], tmp, results->degeneracies[i]+k);
				}
				
				
				for ( int bra = 0 ; bra < dim ; bra++ ) {
					for ( int ket = 0 ; ket < dim ; ket++ ) {
						ierr = MatMult(Op,vecs[ket],tmp) ; CHKERRQ(ierr);
						ierr = VecDot(vecs[bra],tmp,&results->correlations[i][j][bra*dim+ket]); CHKERRQ(ierr);					
					}
				}		
				
			
				for ( int k = 0 ; k < dim ; k++ ) {
					VecDestroy(&vecs[k]);
				}
				ierr = PetscFree(vecs);CHKERRQ(ierr);

			}//end for over degenerate bands.	
			
			//if we have a trial wavefunction calculate the correlations for this too
			if ( parameters->prod_wf_overlap ) {
				ierr = MatMult(Op,trial,tmp) ; CHKERRQ(ierr);
				ierr = VecDot(trial,tmp,&results->trial_correlations[j]); CHKERRQ(ierr);					
			}
	
			VecDestroy(&tmp);					
			MatDestroy(&Op);
		} //end for over correlation operators.
	}// end if using correlations
	
	if ( parameters->prod_wf_overlap ) {
		VecDestroy(&trial);
	}
	
	return 0;
}

