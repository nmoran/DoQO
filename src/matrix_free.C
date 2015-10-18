/* 
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala.
 */
 
#include "matrix_free.h"


#undef __FUNCT__
#define __FUNCT__ "MatSpinSys_Mult"
PetscErrorCode MatSpinSys_Mult( Mat A, Vec x, Vec y ) {
	struct matrix_meta_data *matrix_data;
//	struct matrix_row_values row_vals;
	struct shell_matrix_context *context;
	struct parameter_struct *parameters; 
	void *ctx; 
	PetscInt Istart, Iend,i,j,Istart_other, Iend_other,dest,source; 
	PetscScalar *px,*py;
	PetscErrorCode ierr;;
	PetscScalar *buf1, *buf2;
	MPI_Request send_req;
	MPI_Status recv_stat,send_stat;
	int num_received;
	
	
	
	PetscFunctionBegin;
	ierr = MatShellGetContext( A, &ctx ); CHKERRQ(ierr);
	context = (struct shell_matrix_context*)ctx;
	matrix_data = context->matrix_data;
	parameters = context->parameters; 
	buf1 = context->buf1;
	buf2 = context->buf2;
	ierr = VecGetOwnershipRange(x, &Istart, &Iend );
	VecGetArray(x, &px );
	VecGetArray(y, &py );
	
	for ( i = Istart ; i < Iend ; i++ ) {
		py[i - Istart] = 0.0;
	}
	
//	allocate_row_space(&row_vals, matrix_data);
	
	/*for ( i = Istart ; i < Iend ; i++ ) {
			row_vals.row = i; //set which row to retrieve
			getrowlocal(&row_vals, matrix_data, parameters,Istart,Iend); //get the values from that row	
			
			for ( j = 0 ; j < row_vals.val_count ; j++ ) {
#if defined(PETSC_USE_COMPLEX)			
				//row_vals.values[j] *= px[i-Istart];
				if ( (row_vals.cols[j]-Istart) < 0 || row_vals.cols[j] >= Iend ) PetscFPrintf(PETSC_COMM_WORLD,stderr,"Out of range %d, range %d-%d, col %d.\n" ,row_vals.cols[j]-Istart, Istart, Iend, row_vals.cols[j]);
				py[row_vals.cols[j]-Istart] += (px[i-Istart] * row_vals.values[j]);
#else
				//row_vals.real_vals[j] *= px[i-Istart];
				//row_vals.imag_vals[j] = row_vals.imag_vals[j] * px[i-Istart];
				py[row_vals.cols[j]-Istart] += (px[i-Istart] * row_vals.real_vals[j]);
#endif
			}
	}*/
/*#if defined(PETSC_USE_COMPLEX)
			VecSetValues(y, row_vals.val_count , row_vals.cols, row_vals.values , ADD_VALUES );
#else
			VecSetValues(y, row_vals.val_count , row_vals.cols, row_vals.real_vals , ADD_VALUES );
#endif*/

	
	getrow_matrix_free(context,Istart,Iend,Istart,Iend,px,py);
	
	//VecAssemblyBegin(y);
	//VecAssemblyEnd(y);

	//ierr=PetscMalloc(sizeof(PetscScalar) * get_max_chunk_size(parameters->hamiltonian.numberBasisStates,parameters->size), &buf1);
	//ierr=PetscMalloc(sizeof(PetscScalar) * get_max_chunk_size(parameters->hamiltonian.numberBasisStates,parameters->size), &buf2);
	//	buf1 = (PetscScalar *) malloc( sizeof(PetscScalar) * get_max_chunk_size(parameters->hamiltonian.numberBasisStates,parameters->size));
	//buf2 = (PetscScalar *) malloc( sizeof(PetscScalar) * get_max_chunk_size(parameters->hamiltonian.numberBasisStates,parameters->size));
	
	//loop over each block in the column starting from block after current block
	for ( i = 1 ; i < parameters->size  ; i++ ) { 
		dest = (parameters->rank + i ) % parameters->size;//destination rank 
		source = (parameters->rank - i ) % parameters->size;//source rank  
		//if ( source == parameters->rank ) source--; 
		if ( source < 0 ) source = parameters->size + source;      
		
		Istart_other = get_rank_chunk_start(parameters->numberBasisStates, dest, parameters->size); //start of destination chunk
		Iend_other = Istart_other + get_rank_chunk_size(parameters->numberBasisStates, dest, parameters->size); //end of destination chunk
		
		
		
		if ( i != 1 ){
		  MPI_Wait(&send_req, &send_stat); //if not the first chunk
		}
		
		//zero the buffer
		for ( j = 0 ; j < Iend_other - Istart_other; j++ ) {
			buf1[j] = 0.0 ;
		}
		
		getrow_matrix_free(context,Istart,Iend,Istart_other,Iend_other,px,buf1);
		
		//loop over each column and get local relevant block
		/*for ( j = 0 ; j < Iend - Istart; j++ ) {
			row_vals.row = j + Istart ; 
			getrowlocal(&row_vals, matrix_data, parameters, Istart_other,Iend_other); //get the values from that row	
			for ( k = 0 ; k < row_vals.val_count; k++ ) {
#if defined(PETSC_USE_COMPLEX)			
				//row_vals.values[j] *= px[i-Istart];
				if ( (row_vals.cols[k]-Istart_other) < 0 || row_vals.cols[k] >= Iend_other ) PetscFPrintf(PETSC_COMM_WORLD,stderr,"Out of range %d, range %d-%d, col %d.\n" ,row_vals.cols[k]-Istart_other, Istart_other, Iend_other, row_vals.cols[k]);
				buf1[row_vals.cols[k]-Istart_other] += (px[j] * row_vals.values[k]);
#else
				//row_vals.real_vals[j] *= px[i-Istart];
				//row_vals.imag_vals[j] = row_vals.imag_vals[j] * px[i-Istart];
				buf1[row_vals.cols[k]-Istart_other] += (px[j] * row_vals.real_vals[k]);
#endif
			}
		}*/
		
		
		//use a non blocking send to send chunk of final matrix to destination process.
		MPI_Isend(buf1,Iend_other-Istart_other,MPIU_SCALAR,dest,0,PETSC_COMM_WORLD,&send_req);
		MPI_Recv(buf2,Iend - Istart,MPIU_SCALAR,source,0,PETSC_COMM_WORLD,&recv_stat);

		MPI_Get_count(&recv_stat,MPIU_SCALAR,&num_received);	
		for ( j = 0 ; j < Iend - Istart; j++ ) {
			py[j] += (buf2[j]); 
		}
	
	} //end loop over chunks

	//MPI_Win_wait(mpi_window);
	
	if (parameters->size > 1 ) {    
		MPI_Wait(&send_req, &send_stat); //if we are not on a singe processor	
	}
	//PetscFree(buf1);
	//PetscFree(buf2);
	//free arrays
	//free(buf1);
	//free(buf2);
	
	//restore vectors
	VecRestoreArray(x,&px);
	VecRestoreArray(y,&py);

	//deallocate_row_space(&row_vals);

	//MPI_Barrier(MPI_COMM_WORLD);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "prepare_matrix_free"
int prepare_matrix_free(struct shell_matrix_context *shell_context){
  struct matrix_meta_data *matrix_data;
//  struct matrix_row_values row_vals;
//  struct shell_matrix_context *context;
  struct parameter_struct *parameters;
  int i, term_idx, j;
  PetscErrorCode ierr;
  
  parameters = shell_context->parameters;
  matrix_data = shell_context->matrix_data;
  
  //allocate space for chunk tranfers
  ierr=PetscMalloc(sizeof(PetscScalar) * get_max_chunk_size(parameters->numberBasisStates,parameters->size), &(shell_context->buf1));
  ierr=PetscMalloc(sizeof(PetscScalar) * get_max_chunk_size(parameters->numberBasisStates,parameters->size), &(shell_context->buf2));

  //initialise the number of relevant terms
  shell_context->number_relevant_terms = 0 ;

  //count terms with non zero coefficient
  for ( i = 0 ; i < matrix_data->number_terms ; i++ ) {
    if ( PetscAbsScalar(matrix_data->term_multipliers[i]) > 0.0 ) shell_context->number_relevant_terms++;
  }

  //allocate space for relevant terms
  //shell_context->flips = (uPetscInt*)malloc(shell_context->number_relevant_terms * sizeof(uPetscInt));
  ierr = PetscMalloc(shell_context->number_relevant_terms * sizeof(uPetscInt),&shell_context->flips);CHKERRQ(ierr);
  //shell_context->sub_term_numbers = (int *) malloc(sizeof(int) * shell_context->number_relevant_terms);
  ierr = PetscMalloc(sizeof(int) * shell_context->number_relevant_terms,&shell_context->sub_term_numbers);CHKERRQ(ierr);
  for ( i = 0; i < shell_context->number_relevant_terms; i++ ) shell_context->sub_term_numbers[i] = 0;
  //shell_context->subterms = (int**)malloc(sizeof(int*) * shell_context->number_relevant_terms);
  ierr = PetscMalloc(sizeof(int*) * shell_context->number_relevant_terms,&shell_context->subterms);CHKERRQ(ierr);

  //initialise the final number of terms that will be used.
  shell_context->final_number_terms = 0;
  
  //loop through terms with non zero coefficient
  for ( i = 0 ; i < matrix_data->number_terms ; i++ ) {
    if ( PetscAbsScalar(matrix_data->term_multipliers[i]) > 0.0 ) {
      term_idx = -1; //set it to negative value to flag it has not been found.
      for ( j = 0 ; j < shell_context->final_number_terms ; j++ ) { //loop over already established terms
		if ( shell_context->flips[j] == matrix_data->flips[i] ) {
	  		term_idx = j; //term with the same flip mask already exists so just append the details of this one.
	  		break;
		} // end check to see if it already exists
      }
      if (term_idx == -1 ) { // if a term with this flip mask does not exist then create a new one.
		term_idx = shell_context->final_number_terms++;
		shell_context->flips[term_idx] = matrix_data->flips[i];
      }
      shell_context->sub_term_numbers[term_idx]++; 
    }
  }

  for ( i = 0 ; i < shell_context->final_number_terms; i++ ) {
    //shell_context->subterms[i] = (int *) malloc(sizeof(int) * shell_context->sub_term_numbers[i]);
    ierr = PetscMalloc(sizeof(int) * shell_context->sub_term_numbers[i],&shell_context->subterms[i]);CHKERRQ(ierr);
  }

  for ( i = 0; i < shell_context->number_relevant_terms; i++ ) shell_context->sub_term_numbers[i] = 0;//reset to zero so we can use them to set values
 
  //loop through terms with non zero coefficient
  for ( i = 0 ; i < matrix_data->number_terms ; i++ ) {
    if ( PetscAbsScalar(matrix_data->term_multipliers[i]) > 0.0 ) {
      term_idx = -1; //set it to negative value to flag it has not been found.
      for ( j = 0 ; j < shell_context->final_number_terms ; j++ ) { //loop over already established terms
		if ( shell_context->flips[j] == matrix_data->flips[i] ) {
	  		term_idx = j; //term with the same flip mask already exists so just append the details of this one.
	 	 	break;
		} // end check to see if it already exists
      }
      if (term_idx == -1 ) { // if a term with this flip mask does not exist then create a new one.
	//should not get here 
		PetscFPrintf(PETSC_COMM_SELF,stderr,"Term not accounted for.\n");
      } else {
		shell_context->subterms[term_idx][shell_context->sub_term_numbers[term_idx]] = i ; 
		shell_context->sub_term_numbers[term_idx]++;
      }
    }
  }
	return 0;
      
}

#undef __FUNCT__
#define __FUNCT__ "cleanup_matrix_free"
int cleanup_matrix_free(struct shell_matrix_context *shell_context){
  int i ;
  PetscErrorCode ierr;

  //free up space now
  PetscFree(shell_context->buf1);
  PetscFree(shell_context->buf2);

  //free up space used by term details
  //free(shell_context->flips);
  //free(shell_context->sub_term_numbers);
  ierr = PetscFree(shell_context->flips);CHKERRQ(ierr);
  ierr = PetscFree(shell_context->sub_term_numbers);CHKERRQ(ierr);
  //for ( i = 0 ; i < shell_context->final_number_terms ; i++ ) free(shell_context->subterms[i]);
  for ( i = 0 ; i < shell_context->final_number_terms ; i++ ) {
  	ierr = PetscFree(shell_context->subterms[i]);CHKERRQ(ierr);
  }
  //free(shell_context->subterms);
  ierr = PetscFree(shell_context->subterms);

	return 0 ;
}


#undef __FUNCT__
#define __FUNCT__ "getrow_matrix_free"
int getrow_matrix_free(struct shell_matrix_context *shell_context, PetscInt Istart, PetscInt Iend , PetscInt Istart_other, PetscInt Iend_other, PetscScalar *px, PetscScalar *buf) {
	PetscInt i;//iterate over columns
	int j,k, zi, yi,term; //iterate ove terms and subterms
	uPetscInt col ; //store current column index
	uPetscInt row , zero_const;; //store current row index
	struct matrix_meta_data *matrix_data;
	PetscScalar tmp;
#if defined(PETSC_USE_COMPLEX)
	PetscScalar val_tmp ; 
#else
	PetscReal real_tmp, imag_tmp;
#endif
	
	zero_const = (uPetscInt)0;
	
	matrix_data = shell_context->matrix_data; 
	
	//loop over each column 
	for ( i = 0; i < Iend - Istart ; i++ ) {
		if ( px[i] != 0.0 ) { //if the element of the x vector is not zero
			col = i + Istart; //get the global column index
			//loop over terms 
			for ( j = 0 ; j <  shell_context->final_number_terms ; j++ ) {
				row = col ^ shell_context->flips[j]; // appply flip mask to get the row 
				if ( row >= (uPetscInt)Istart_other && row < (uPetscInt)Iend_other  ) { //if it is within the correct range
					for ( k = 0 ; k < shell_context->sub_term_numbers[j] ; k++ ) {
						term = shell_context->subterms[j][k];
#if defined(PETSC_USE_COMPLEX)
						val_tmp =  shell_context->matrix_data->term_multipliers[term];
#else
						real_tmp =  shell_context->matrix_data->term_multipliers[term];
						imag_tmp = 0.0;
#endif

						for ( zi = 0 ; zi < matrix_data->z_count[term] ; zi++ ) {
							if ( (col  & matrix_data->zs[term][zi]) != zero_const ) { //if there is a '1' at that position
#if defined(PETSC_USE_COMPLEX)				
								val_tmp = -1.0 * val_tmp;
#else				
								real_tmp = -real_tmp;
								//imag_tmp = -imag_tmp // no point in this as it has to be 0 at this point
#endif				
							}
						}	 //	end for over zs	
						
						//look at the Ys on the after flip which is the column number
						for ( yi = 0 ; yi < matrix_data->y_count[term] ; yi++ ) {
							if ( (row & matrix_data->ys[term][yi]) > zero_const ) { //if there is a '1' at that position
							//multiply by -i
#if defined(PETSC_USE_COMPLEX)				
								val_tmp =  -1.0 * val_tmp * PETSC_i;
#else				
								tmp = imag_tmp; 
								imag_tmp = -real_tmp;
								real_tmp = tmp;
#endif				
					
							} else {
						//multiply by i 
#if defined(PETSC_USE_COMPLEX)				
								val_tmp =  val_tmp * PETSC_i;
#else
								tmp = imag_tmp; 
								imag_tmp = real_tmp;
								real_tmp = -tmp;	
#endif				
							}
						} // end for over ys 	
						
						//apply complex congugate operation
#if defined(PETSC_USE_COMPLEX)	
						val_tmp = PetscConj(val_tmp);
#else
						imag_tmp = - imag_tmp ; 
#endif

#if defined(PETSC_USE_COMPLEX)	
						buf[row - Istart_other] += val_tmp * px[i]; 	
#else
						buf[row - Istart_other] += real_tmp * px[i]; 
#endif
						
					}// end loop over the subterms
				}
			}
		}
	}
	
	
	
	
	
	
/*	int i, j, not_found; 
	uPetscInt col_tmp, zero_const; 
	PetscReal tmp;
	//char *num;
#if defined(PETSC_USE_COMPLEX)
	PetscScalar val_tmp ; 
#else
	PetscReal real_tmp, imag_tmp;
#endif
	//initialise the value counter
	row->val_count = 0 ;
	zero_const = (uPetscInt)0;

	//PetscPrintf(PETSC_COMM_WORLD, "row\n");
	//For each row of the structure matrix
	for ( i = 0; i < matrix_data->number_terms ; i++ ) {
		if ( fabs(matrix_data->term_multipliers[i]) > 0.0 ) {
	
			col_tmp = row->row ^ matrix_data->flips[i]; //this expresion will flip all the spins due to X and Y terms.
			//PetscPrintf(PETSC_COMM_WORLD, "Row: %d, Part: %d, Col: %d\n" , row->row , i, col_tmp ) ;


			//if ( col_tmp >= (uPetscInt)parameters->Istart && col_tmp <= (uPetscInt)parameters->Iend ) { //if it is local the use it but otherwise done
			if ( col_tmp >= (uPetscInt)Istart && col_tmp < (uPetscInt)Iend ) { //if it is local the use it but otherwise done

#if defined(PETSC_USE_COMPLEX)
			val_tmp = 1.0;  
#else
			real_tmp = 1.0;
			imag_tmp = 0.0;
#endif
	    
			//PetscPrintf(PETSC_COMM_WORLD, "%lf\n", tmp_val.real()); 
	    
			//look at the Zs before the flip so we look at 'row' rather than 'col_tmp'
			for ( j = 0 ; j < matrix_data->z_count[i] ; j++ ) {
				if ( (row->row  & matrix_data->zs[i][j]) != zero_const ) { //if there is a '1' at that position
#ifdef SUPERDEBUG
		num = (char *) malloc((matrix_data->number_spins + 1) * sizeof(char));
		PetscPrintf(PETSC_COMM_WORLD, "%s ", dec2bin(row->row, matrix_data->number_spins,num) ); 
		PetscPrintf(PETSC_COMM_WORLD, "& %s ", dec2bin(matrix_data->zs[i][j] , matrix_data->number_spins,num) );
		PetscPrintf(PETSC_COMM_WORLD, "= %s ",dec2bin(row->row  & matrix_data->zs[i][j], sizeof(uPetscInt) * 8 , num )); 
		PetscPrintf(PETSC_COMM_WORLD, "!= %s\n",dec2bin( zero_const, sizeof(uPetscInt) * 8, num )) ;
		free(num);
		ierr = PetscFree(num);CHKERRQ(ierr);
#endif
				
#if defined(PETSC_USE_COMPLEX)				
					val_tmp = -val_tmp;
#else				
					real_tmp = -real_tmp;
				//imag_tmp = -imag_tmp // no point in this as it has to be 0 at this point
#endif				
					}
			}	 //	end for over zs	
	    
			//look at the Ys on the after flip which is the column number
			for ( j = 0 ; j < matrix_data->y_count[i] ; j++ ) {
				if ( (col_tmp & matrix_data->ys[i][j]) > zero_const ) { //if there is a '1' at that position
					//multiply by -i
#if defined(PETSC_USE_COMPLEX)				
					val_tmp =  val_tmp *   -PETSC_i;
#else				
					tmp = imag_tmp; 
					imag_tmp = -real_tmp;
					real_tmp = tmp;
#endif				
					
				} else {
				//multiply by i 
#if defined(PETSC_USE_COMPLEX)				
					val_tmp =  val_tmp * PETSC_i;
#else
					tmp = imag_tmp; 
					imag_tmp = real_tmp;
					real_tmp = -tmp;	
#endif				
				}
			} // end for over ys 		
	    
			//apply multiplier for each term
#if defined(PETSC_USE_COMPLEX)	
			val_tmp *= matrix_data->term_multipliers[i];
#else
			real_tmp *= matrix_data->term_multipliers[i];
			imag_tmp *= matrix_data->term_multipliers[i];
#endif
	 
			not_found = 1; //indicates not found
			//see if that column already exists. This can be sped up by keeping array sorted by column number  
			for ( j = 0 ; j < row->val_count ; j++ ) {
			//if this column is already there then add the values 
				if (row->cols[j] == col_tmp ) {
					//PetscPrintf(PETSC_COMM_WORLD, "Cols the same on part %d: col %d == col %d\n" , i, row->cols[j] , col_tmp );
#if defined(PETSC_USE_COMPLEX)
					row->values[j] += val_tmp;
#else
					row->real_vals[j] += real_tmp; 
					row->imag_vals[j] += imag_tmp; 
#endif
					not_found = 0; //indicate that same column was found
					break; 
				}
			}
	    
			//PetscPrintf(PETSC_COMM_WORLD, "Found %d\n", not_found ); 
			//if column does not already exist then add new value 
			if ( not_found ) {
				//PetscPrintf(PETSC_COMM_WORLD, "Adding new nz value %lf at col %d, this val num %d, multiplier %lf\n", real_tmp, col_tmp, row->val_count, matrix_data->term_multipliers[i] ); 
				row->cols[row->val_count] = col_tmp; 
#if defined(PETSC_USE_COMPLEX)
				row->values[row->val_count] = val_tmp;
#else
				row->real_vals[row->val_count] = real_tmp;
				row->imag_vals[row->val_count] = imag_tmp;
#endif
				row->val_count++; 	
				
				} //end if in local range
			} // end if fabs(term_multipliers) > 0
		} //end loop over terms
		
	} //end for loop over parts*/
	    	
	return 0;
}
