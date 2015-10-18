/* 
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala.
 *
 * Functions for creating and manipulating hamiltonians 
 * 
 * Version: $Id$
 * 
 * 
*/


#include "hamiltonian.h"

#define STR_LEN 50


#undef __FUNCT__
#define __FUNCT__ "setup_operators"
/*! \brief This function calls the setup_operator function on each operator. */
int setup_operators(struct parameter_struct *parameters) {
	if ( setup_operator( parameters, &parameters->hamiltonian) > 0 ) {
		return 1; 
	}

#ifdef DEBUG
	PetscFPrintf(PETSC_COMM_WORLD, stderr, "Setup first operator.\n");
#endif
	
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "setup_operator"
/*! \brief Sets up the operator. This includes reading parameters and terms form the input files and calculating the flips. */
int setup_operator(struct parameter_struct *parameters, struct matrix_meta_data *matrix_data) {
	FILE *fp;
	char  *ptr ;
	int c;
	char buffer[MAX_STRING];
	int count;
//	Vec tmp;
//	PetscErrorCode ierr;
	

	//read data from the input file
	fp = fopen(matrix_data->filename, "r");
	if ( fp != NULL ) {
		//read the number of spins in the model
		fscanf(fp,"SITES %d\n",&(parameters->number_particles));
		fscanf(fp,"PARAMETERS\n");
		
		
		parameters->numberBasisStates = ((PetscInt)1) << parameters->number_particles ; //work out the number of basis states from the number of spins
		if ( parameters->model_type == SPIN_ONE ) {
			parameters->numberBasisStates = (PetscInt)pow(3.0,(double)parameters->number_particles); //work out the number of basis states for spin one systems.
		}
		
		//matrix_data->mloc = get_rank_chunk_size(matrix_data->numberBasisStates, parameters->rank, parameters->size);
		//destroy again.
		
		//now we count the amount of paramters for this operator
		matrix_data->number_parameters = 0;
		strcpy(buffer, "");
		ptr = buffer ;
		while ( (c = fgetc(fp)) != EOF) {  // While there are still characters left in the file and we have not come across the string "TERMS"
			if ( c == '\n' ) { //when we get to a line break check the parameter to see it is not terms
				*ptr = '\0';
				if ( strcmp(buffer, "TERMS") != 0 ) {
					matrix_data->number_parameters++; //increment the parameter count 
					ptr = buffer; //change the ptr so it points to the start of the buffer again.
				} else {
					break; //we have come to the list of terms so break
				}
			} else if ( c != ' ') { // otherwise if the character is not a space we add it to the buffer
				*ptr++ = (char)c;  
			}
		}
		
#ifdef DEBUG
		PetscFPrintf(PETSC_COMM_WORLD, stderr, "Counted parameters. %d.\n", matrix_data->number_parameters);
#endif

		//now we count the amount of terms for this operator
		matrix_data->number_terms = 0;
		
		int char_count = 0; 
		while ( (c = fgetc(fp)) != EOF) {
			if ( c == '\n' ) { 
				if (char_count > 0 ) matrix_data->number_terms++; //increment the terms count 
				char_count = 0;
			} else if ( c != ' ') {
				char_count++; 
			}
		}
		if (char_count > 0 ) matrix_data->number_terms++; //increment the terms count 
		
#ifdef DEBUG
		PetscFPrintf(PETSC_COMM_WORLD, stderr, "Counted terms. %d.\n", matrix_data->number_terms);
#endif	
	
		fclose(fp);
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"Error opening file %s.\n", matrix_data->filename);
		return 1;
	}	
	
	//allocate space for parameters and term
	allocate_matrix_data(matrix_data);

#ifdef DEBUG
	PetscFPrintf(PETSC_COMM_WORLD, stderr, "Allocated space.\n");
#endif

	//read in the parameter names and the terms
	fp = fopen(matrix_data->filename, "r") ;
	
	if ( fp != NULL ) {
		fscanf(fp,"SITES %d\n", &(parameters->number_particles));
		fscanf(fp,"PARAMETERS\n");
		
		//read parameter names 
		ptr = buffer; 
		count = 0;
		while ( (c = fgetc(fp)) != EOF ) {
			if ( c != '\n' ) {
				if ( c != ' ') *ptr++ = c ;
			} else {
				*ptr = '\0';
				if ( strcmp(buffer, "TERMS") == 0 ){
					if ( count != matrix_data->number_parameters ) {
						PetscPrintf(PETSC_COMM_WORLD,"Parameter count missmatch.\n");
					}
					break;
				}
				int found = 0 ;
				for ( int i = 0 ; i < parameters->tasks.number_parameters; i++ ) {
					if ( strcmp(parameters->tasks.parameter_names[i].c_str(), buffer) == 0 ) { // if the parameter already exists.
						found = 1;
						matrix_data->parameter_indices[count] = i ;
					}
				}
				if ( !found ) {
					string tmp(buffer);
					parameters->tasks.parameter_names.push_back(tmp);
					matrix_data->parameter_indices[count] = parameters->tasks.number_parameters ;
					parameters->tasks.number_parameters++; 
				}
				count++; 
				//strcpy(matrix_data->parameter_names[count++], buffer );
				ptr = buffer;
			}
		}

#ifdef DEBUG
		PetscFPrintf(PETSC_COMM_WORLD, stderr, "Read parameters.\n");
#endif
		
		//read in the terms
		ptr = buffer; 
		count = 0;
		int ptr_len = 0;
		while ( (c = fgetc(fp)) != EOF ) {
			if ( c != '\n' ) {
				if ( c != ' ' ){ 
					*ptr++ = c ;
					ptr_len++;
				}
			} else {
				*ptr = '\0';
				if ( ptr_len > 0 ) {
					strcpy(matrix_data->terms[count++], buffer );
				}
				ptr_len = 0 ;
				ptr = buffer;
			}
		}	

#ifdef DEBUG
		PetscFPrintf(PETSC_COMM_WORLD, stderr, "Read terms.\n");
#endif

		fclose(fp);
	} else {
		return 1 ;
	}
	
	prepare_matrix_data(parameters,matrix_data);

	return 0; 
}




#undef __FUNCT__
#define __FUNCT__ "prepare_matrix_data"
/*! \brief Prepare the data necessary for using the getrow function to retrieve a row for this matrix.*/
int prepare_matrix_data(struct parameter_struct *parameters,struct matrix_meta_data *matrix_data) {
	int i;
	PetscErrorCode ierr;
		
	PetscTruth all_z; //bool used to indicate if all the operators are z 

	matrix_data->nz = 0 ;
	all_z = PETSC_FALSE ; //setting to non zero value to indicate no line of all Z found yet

	//first prepare each of the flip masks
	//loop through each part getting flip masks from Xs and Ys and countings Zs and Ys 
	for ( i = 0 ; i < matrix_data->number_terms ; i++ ) {
		prepare_term_data(parameters,matrix_data, i ) ;
	}
	
	//deallocate raw terms to save memory
	for ( i = 0 ; i < matrix_data->number_terms ; i++ ){
		ierr = PetscFree(matrix_data->terms[i]);CHKERRQ(ierr);
	}
	ierr = PetscFree(matrix_data->terms);CHKERRQ(ierr);
	
	return 0 ;
}  // end prepare matrix data function 	


#undef __FUNCT__
#define __FUNCT__ "perpare_term_data"
/*! \brief Parse a single term and extract information needed for using it to construct the matrix. */
int prepare_term_data(struct parameter_struct *parameters, struct matrix_meta_data *matrix_data, int term ) {
		int spin, pos, ws, idx;
		char type,*ptr,*buf_ptr;
		char buf[MAX_PARM];
		int  yidx, zidx;
		PetscErrorCode ierr;

		//initialise the flip masks and counts to 0s 
		matrix_data->flips[term] = 0; 
		matrix_data->y_count[term] = 0; 
		matrix_data->z_count[term] = 0; 
		matrix_data->x_count[term] = 0; 
		matrix_data->term_multipliers[term] = 0.0 ;
		
		//go through each item in the term
		pos = 0;
		while ( matrix_data->terms[term][pos] != '*' && matrix_data->terms[term][pos] != '\0' ) { //loop through string with pos as index while not at the end of the string or at a '*' character
			get_term(matrix_data->terms[term], &spin, &type, &pos); // read a term and return the spin number, type and the position in the string one is left at
			//PetscFPrintf(PETSC_COMM_WORLD, stderr,"Spin %d, type %c " , spin, type );
			//PetscFPrintf(PETSC_COMM_WORLD, stderr,"\n" );
			//start from the last element in the row 
			switch(type){
			case('X'): 
				//matrix_data->flips[term] ^= (uPetscInt) pow((double)2,(double)(matrix_data->number_spins - spin)); //difference between ordering to match python code
				//matrix_data->flips[term] ^= (uPetscInt) pow((double)2,(double)(spin-1));
				matrix_data->flips[term] ^= ((uPetscInt) 1)<< (spin-1);
				matrix_data->x_count[term]++; //increment the number of xs 
				break;
			case('Y'):
				//matrix_data->flips[term] ^= (uPetscInt) pow((double)2,(double)(matrix_data->number_spins - spin));
				//matrix_data->flips[term] ^= (uPetscInt) pow((double)2,(double)(spin-1));
				matrix_data->flips[term] ^= ((uPetscInt) 1)<< (spin-1);
				matrix_data->y_count[term]++; //increment the number of ys
				break;
			case('Z'):
				matrix_data->z_count[term]++; //increment the number of zs 
				break;
			//going to store the masks for both annihilation and creation operators in the ys and will save an indicated in the zs to say whether it 
			//is an annihilation or creation operator. Bit messy but should work fine.
			case('C'): //in the case where we are dealing with fermionic systems will use the sy data structures to store positions of creation operators.
			case('A'):  //in the case where we are dealing with fermionic systems will use the sz data structures to store positions of annihilation operators.
			case('S'): //will use this for the sz operator for spin one systems.
			case('R'):
			case('L'):  
				//matrix_data->flips[term] ^= (uPetscInt) pow((double)2,(double)(spin-1)); //this might need to be changed.
				matrix_data->y_count[term]++; //increment the number of ys
				matrix_data->z_count[term]++; //increment the number of zs
				break;

			} // end switch
		} //end while: we must have come to the end or encountered a '*' character
		
		//allocate space for the y and z masks
		ierr = PetscMalloc(matrix_data->y_count[term] * sizeof(uPetscInt),&matrix_data->ys[term]);CHKERRQ(ierr);
		ierr = PetscMalloc(matrix_data->z_count[term] * sizeof(uPetscInt),&matrix_data->zs[term]);CHKERRQ(ierr);
		
		//loop through again and populate the y and z masks
		yidx = 0;
		zidx = 0;
		pos = 0;
		while ( matrix_data->terms[term][pos] != '*' && matrix_data->terms[term][pos] != '\0' ) {
			get_term(matrix_data->terms[term], &spin, &type, &pos);
			
			switch(type){
				case('Y'):
					//matrix_data->ys[term][yidx++] = (uPetscInt)pow((double)2,(double)(matrix_data->number_spins - spin)); //place a 1 in the position of the Y 
					matrix_data->ys[term][yidx++] = ((uPetscInt) 1)<< (spin-1); //place a 1 in the position of the Y 
					break;
				case('Z'):
					//matrix_data->zs[term][zidx++] = (uPetscInt)pow((double)2,(double)(matrix_data->number_spins - spin)); //place a 1 in the position if the Z 
					matrix_data->zs[term][zidx++] = ((uPetscInt) 1)<< (spin-1); //place a 1 in the position if the Z 					
					break;
				case('C'):
					matrix_data->ys[term][yidx] = ((uPetscInt) 1)<< (spin-1); //place a 1 in the position of the Y 
					matrix_data->zs[term][yidx++] = CREATION_OPERATOR;
					break;
				case('A'):
					matrix_data->ys[term][yidx] = ((uPetscInt) 1)<< (spin-1); //place a 1 in the position of the Y 
					matrix_data->zs[term][yidx++] = ANNIHILATION_OPERATOR;
					break;
				case('S'):
					matrix_data->ys[term][yidx] = ((uPetscInt) 1)<< (spin-1); //place a 1 in the position of the Y 
					matrix_data->zs[term][yidx++] = SPIN_ONE_Z_OPERATOR;
					break;
				case('R'):
					matrix_data->ys[term][yidx] = ((uPetscInt) 1)<< (spin-1); //place a 1 in the position of the Y 
					matrix_data->zs[term][yidx++] = SPIN_ONE_RAISING;
					break;
				case('L'):
					matrix_data->ys[term][yidx] = ((uPetscInt) 1)<< (spin-1); //place a 1 in the position of the Y 
					matrix_data->zs[term][yidx++] = SPIN_ONE_LOWERING;
					break;
			default:
				break;
			} // end switch
		} //end while: we must have encountered the end of string or '*' character again
		
		matrix_data->number_term_params[term] = 0 ;
		if ( matrix_data->terms[term][pos] == '*' ) { // if there are parameters			
			//now we count the parameters for the term
			ptr = &matrix_data->terms[term][pos+1] ;
			ws = 1;  //set to 1 to indicate that we are in white space
			while ( *ptr != '\0' ) {  //loop till end of line
				if ( (*ptr == ' ' || *ptr == ',') && ws != 1 ) {
					ws = 1;
				} else if ( *ptr != ' ' && ws == 1 ) {
					matrix_data->number_term_params[term]++; 
					ws = 0; 
				}
				ptr++;
			}
			//allocate space to store the array of term parameters
			//matrix_data->term_parameters[term] = (int *)malloc( matrix_data->number_term_params[term] * sizeof(int));
			ierr = PetscMalloc(matrix_data->number_term_params[term] * sizeof(int),&matrix_data->term_parameters[term] );CHKERRQ(ierr);
			
			//populate list of term parameters
			ptr = &matrix_data->terms[term][pos+1] ;
			ws = 1; 
			idx = 0; 
			buf_ptr = buf;
			while ( *ptr != '\0' ) { 
				if ( (*ptr == ' ' || *ptr == ',') && ws != 1 ) {
					ws = 1;
					*buf_ptr = '\0';
					//strcpy(matrix_data->term_parameters[idx++] , buf ) ;
					//matrix_data->term_parameters[term][idx++] = get_parameter_name_idx( matrix_data , buf ) ;
					
					for ( int i = 0 ; i < parameters->tasks.number_parameters ; i++ ) {
						if ( strcmp( parameters->tasks.parameter_names[i].c_str(), buf) == 0 ) {
							matrix_data->term_parameters[term][idx++] = i;
						}
					}
					
					buf_ptr = buf;
				} else if ( *ptr != ' ' && *ptr != ',' ) {
					if ( ws == 1 ) {
						ws = 0;
						buf_ptr = buf; 
					}
					*buf_ptr++ = *ptr;
	
				}
				ptr++;
			}
			
			if ( ws == 0 ) {
				*buf_ptr = '\0';
				//strcpy(matrix_data->term_parameters[idx++] , buf ) ;
				//matrix_data->term_parameters[term][idx++] = get_parameter_name_idx( matrix_data , buf ) ;
				for ( int i = 0 ; i < parameters->tasks.number_parameters ; i++ ) {
					if ( strcmp( parameters->tasks.parameter_names[i].c_str(), buf) == 0 ) {
						matrix_data->term_parameters[term][idx++] = i;
					}
				}
			}
		}
	return 0 ;
}


#undef __FUNCT__
#define __FUNCT__ "get_term"
/*! \brief When given a pointer to a position in a term it will return the type X, Y, Z, A or C and the spin number 1 - <number spins>' */ 
int get_term(char *str, int *spin, char *type, int *pos){
	char buf[MAX_PARM], *buf_ptr, *ptr;
	
	strcpy(buf,"");
	ptr = &str[*pos];
	buf_ptr = buf;
	*spin = 0; 
	while ( *ptr != '\0' && *ptr != '*' && *ptr != ',' ) { //while not at the end of the line or at a '*'
		if ( *ptr <= '9' && *ptr >= '0' ) { //if the character is numeric
			*buf_ptr++ = *ptr;
		} else { //must be finished the character.
			if ( *spin == 0 ) {
				*buf_ptr = '\0'; // add termination character
				*spin = atoi(buf); // convert string to int
				buf_ptr = buf; //reset pointer to the beginning
			}
		}
		
		if ( *ptr == 'X' || *ptr == 'Y' || *ptr == 'Z' || *ptr == 'C' || *ptr == 'A' || *ptr == 'R' || *ptr == 'L' || *ptr == 'S' ) { // the A and the C stand for annihilation operator and creation operator on that fermion site if dealing with a fermionic system.
			*type = *ptr; 
		}
		ptr++;
		(*pos)++;
	}
	if (*ptr == ',' ) (*pos)++ ; //if we are at the comma move it on one
	return 0 ;
} 


#undef __FUNCT__
#define __FUNCT__ "allocate_row_space"
/*! \brief Allocate space to store the column indices and values of the non zero elements on a single row of the sparse matrix. */
int allocate_row_space(struct matrix_row_values *row, struct matrix_meta_data *matrix_data){
	PetscErrorCode ierr;
	//row->cols = (PetscInt*)malloc (matrix_data->number_terms * sizeof(PetscInt));
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(PetscInt),&row->cols);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
	//row->values = (PetscScalar *)malloc (matrix_data->number_terms * sizeof(PetscScalar));
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(PetscScalar),&row->values);CHKERRQ(ierr);
#else
	//row->real_vals = (PetscReal *)malloc (matrix_data->number_terms * sizeof(PetscReal));
	//row->imag_vals = (PetscReal *)malloc (matrix_data->number_terms * sizeof(PetscReal));
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(PetscReal),&row->real_vals);CHKERRQ(ierr);
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(PetscReal),&row->imag_vals);CHKERRQ(ierr);
	
#endif

	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "deallocate_row_space"
/*! \brief Deallocate space used to store the column indices and values of the non zero elements on a single row of the sparse matrix. */
int deallocate_row_space(struct matrix_row_values *row){
	PetscErrorCode ierr;
	//free(row->cols);
	ierr = PetscFree(row->cols);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
	//free(row->values);
	ierr = PetscFree(row->values);CHKERRQ(ierr);
#else
	//free(row->real_vals);
	//free(row->imag_vals);
	ierr = PetscFree(row->real_vals);CHKERRQ(ierr);
	ierr = PetscFree(row->imag_vals);CHKERRQ(ierr);
#endif
	return 0 ;
}

#undef __FUNCT__
#define __FUNCT__ "getrow_spin_half"
/*! \brief Given a row index calculate the columns and values of the non zero elements in that row for spin half systems. */
int getrow_spin_half(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters, int nz) {
	int i, j, not_found, idx; 
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
//	PetscErrorCode ierr;

	//PetscPrintf(PETSC_COMM_WORLD, "row\n");
	//For each row of the structure matrix
	for ( i = 0; i < matrix_data->number_terms ; i++ ) {
		if ( PetscAbsScalar(matrix_data->term_multipliers[i]) > 0.0 ) {
	
			col_tmp = row->row ^ matrix_data->flips[i]; //this expresion will flip all the spins due to X and Y terms.
			//PetscPrintf(PETSC_COMM_WORLD, "Row: %d, Part: %d, Col: %d\n" , row->row , i, col_tmp ) ;


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
		//num = (char *) malloc((matrix_data->number_spins + 1) * sizeof(char));
		ierr = PetscMalloc((matrix_data->number_spins + 1) * sizeof(char),&num);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD, "%s ", dec2bin(row->row, matrix_data->number_spins,num) ); 
		PetscPrintf(PETSC_COMM_WORLD, "& %s ", dec2bin(matrix_data->zs[i][j] , matrix_data->number_spins,num) );
		PetscPrintf(PETSC_COMM_WORLD, "= %s ",dec2bin(row->row  & matrix_data->zs[i][j], sizeof(uPetscInt) * 8 , num )); 
		PetscPrintf(PETSC_COMM_WORLD, "!= %s\n",dec2bin( zero_const, sizeof(uPetscInt) * 8, num )) ;
		//free(num);
		ierr = PetscFree(num);CHKERRQ(ierr);
#endif
				
#if defined(PETSC_USE_COMPLEX)				
					val_tmp = -1.0 * val_tmp;
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
					val_tmp =  -1.0 * val_tmp *  PETSC_i;
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
			val_tmp *= matrix_data->term_multipliers[i] ;
#else
			real_tmp *= matrix_data->term_multipliers[i];
			imag_tmp *= matrix_data->term_multipliers[i];
#endif
	 
			not_found = 1; //indicates not found
			//see if that column already exists. This can be sped up by keeping array sorted by column number  
			for ( j = 0 ; j < row->val_count ; j++ ) {
			//if this column is already there then add the values 
					if (row->cols[j] == (PetscInt)col_tmp ) {
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
			} // end if fabs(term_multipliers) > 0
		} //end loop over terms

		if ( nz == 1 ) { 
			idx = 0;
			for( j = 0 ; j < row->val_count ; j++ ) {

#if defined(PETSC_USE_COMPLEX)
				if ( PetscAbsScalar(row->values[j]) > parameters->tol ) {
#else
				if ( fabs(row->real_vals[j]) > parameters->tol  || row->imag_vals[j] > parameters->tol) {
#endif
					if ( idx != j ) {
						row->cols[idx] = row->cols[j];
#if defined(PETSC_USE_COMPLEX)
						row->values[idx] = row->values[j];
#else
						row->real_vals[idx] = row->real_vals[j];					
						row->imag_vals[idx] = row->imag_vals[j];
#endif
					} //end if idx != j
					idx++;
#if defined(PETSC_USE_COMPLEX)					
				} //end if abs(value) > 0
#else
				} //end if abs(value) > 0
#endif
			}//end for j
			row->val_count = idx;
		} // end if nz == 1 
		
	} //end for loop over parts	    	
		
	return 0;
}



#undef __FUNCT__
#define __FUNCT__ "getrow_fermionic"
/*! \brief Given a row index calculate the columns and values of the non zero elements in that row for fermionic systems. */
int getrow_fermionic(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters, int nz) {
	int i, j, not_found, idx; 
	uPetscInt col_tmp, zero_const; 
	
	//char *num;
#if defined(PETSC_USE_COMPLEX)
	PetscScalar val_tmp ; 
#else
	PetscReal real_tmp, imag_tmp;
#endif
	//initialise the value counter
	row->val_count = 0 ;
	zero_const = (uPetscInt)0;
//	PetscErrorCode ierr;

	//PetscPrintf(PETSC_COMM_WORLD, "row\n");
	//For each row of the structure matrix
	for ( i = 0; i < matrix_data->number_terms ; i++ ) {
		if ( PetscAbsScalar(matrix_data->term_multipliers[i]) > 0.0 ) {
	
			//col_tmp = row->row ^ matrix_data->flips[i]; //this expresion will flip all the spins due to X and Y terms.
			//PetscPrintf(PETSC_COMM_WORLD, "Row: %d, Part: %d, Col: %d\n" , row->row , i, col_tmp ) ;
			col_tmp = row->row;


#if defined(PETSC_USE_COMPLEX)
			val_tmp = 1.0;  
#else
			real_tmp = 1.0;
			imag_tmp = 0.0;
#endif
	    
			//look at the Zs before the flip so we look at 'row' rather than 'col_tmp'
			//The Zs here will be the annihilation operators so the behavior we want is that for all sites occupied up to this site the parity will be 
			//changed and if the site upon which this operator is being applied is not occupied the value should be set to 0.
			for ( j = matrix_data->z_count[i]-1 ; j >=0  ; j--) {	
				if ( matrix_data->zs[i][j] == ANNIHILATION_OPERATOR) { 
					if ( (col_tmp  & matrix_data->ys[i][j]) == zero_const ) { //if there is no fermion at this site to be annilated then we set the value to zero.
					
	#if defined(PETSC_USE_COMPLEX)				
						val_tmp = 0.0;
	#else				
						real_tmp = 0.0;
					//imag_tmp = -imag_tmp // no point in this as it has to be 0 at this point
	#endif				
					} else { //otherwise if there was a fermion to be annihilated we annihilate it here and also we account for the parity here.
						//first we have to find the position of the annihilation operator or we can just keep going till we come to it.
						col_tmp -= matrix_data->ys[i][j] ;
						int pos = 0 ;
						while ( ( ( ((uPetscInt)1) << pos) & matrix_data->ys[i][j]) == zero_const ) {
							if ( ( ( ((uPetscInt)1) << pos) & col_tmp) > zero_const ) {
	#if defined(PETSC_USE_COMPLEX)				
						val_tmp *= -1.0;
	#else				
						real_tmp *= -1.0;
					//imag_tmp = -imag_tmp // no point in this as it has to be 0 at this point
	#endif
							}
							pos++;
						}
					}
				} else if ( matrix_data->zs[i][j] == CREATION_OPERATOR) {
					if ( (col_tmp  & matrix_data->ys[i][j]) > zero_const ) { //if there a fermion at this site then we set the value to zero.
					
	#if defined(PETSC_USE_COMPLEX)				
						val_tmp = 0.0;
	#else				
						real_tmp = 0.0;
					//imag_tmp = -imag_tmp // no point in this as it has to be 0 at this point
	#endif				
					} else { //otherwise if there was a fermion to be annihilated we annihilate it here and also we account for the parity here.
						//first we have to find the position of the annihilation operator or we can just keep going till we come to it.
						col_tmp += matrix_data->ys[i][j] ;
						int pos = 0 ;
						while ( (( ((uPetscInt)1) << pos) & matrix_data->ys[i][j]) == zero_const ) {
							if ( (( ((uPetscInt)1) << pos) & col_tmp) > zero_const ) {
	#if defined(PETSC_USE_COMPLEX)				
						val_tmp *= -1.0;
	#else				
						real_tmp *= -1.0;
					//imag_tmp = -imag_tmp // no point in this as it has to be 0 at this point
	#endif
							}
							pos++;
						}
					}
				}
			}	 
			
			
			//apply multiplier for each term
#if defined(PETSC_USE_COMPLEX)	
			val_tmp *=  matrix_data->term_multipliers[i] ;
#else
			real_tmp *= matrix_data->term_multipliers[i] ;
			imag_tmp *= matrix_data->term_multipliers[i];
#endif
	 
			not_found = 1; //indicates not found
			//see if that column already exists. This can be sped up by keeping array sorted by column number  
			for ( j = 0 ; j < row->val_count ; j++ ) {
			//if this column is already there then add the values 
					if (row->cols[j] == (PetscInt)col_tmp ) { 
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
				//row->values[row->val_count] += (PetscScalar)((double)i * PETSC_i);
#else
				row->real_vals[row->val_count] = real_tmp;
				row->imag_vals[row->val_count] = imag_tmp;
#endif
				row->val_count++; 	
			} // end if fabs(term_multipliers) > 0
		} //end loop over terms

		if ( nz == 1 ) { 
			idx = 0;
			for( j = 0 ; j < row->val_count ; j++ ) {

#if defined(PETSC_USE_COMPLEX)
				if ( PetscAbsScalar(row->values[j]) > parameters->tol ) {
#else
				if ( fabs(row->real_vals[j]) > parameters->tol  || row->imag_vals[j] > parameters->tol) {
#endif
					if ( idx != j ) {
						row->cols[idx] = row->cols[j];
#if defined(PETSC_USE_COMPLEX)
						row->values[idx] = row->values[j];
#else
						row->real_vals[idx] = row->real_vals[j];					
						row->imag_vals[idx] = row->imag_vals[j];
#endif
					} //end if idx != j
					idx++;
#if defined(PETSC_USE_COMPLEX)	
				} //end if abs(value) > 0
#else
				} //end if abs(value) > 0
#endif
			}//end for j
			row->val_count = idx;
		} // end if nz == 1 
		
	} //end for loop over parts

	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "getrow_spin_one"
/*! \brief Given a row index calculate the columns and values of the non zero elements in that row for spin one systems. */
int getrow_spin_one(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters, int nz) {
	int i, j, not_found, idx; 
	uPetscInt col_tmp, zero_const; 
	
	//char *num;
#if defined(PETSC_USE_COMPLEX)
	PetscScalar val_tmp ; 
#else
	PetscReal real_tmp, imag_tmp;
#endif
	//initialise the value counter
	row->val_count = 0 ;
	zero_const = (uPetscInt)0;
//	PetscErrorCode ierr;

	//PetscPrintf(PETSC_COMM_WORLD, "row\n");
	//For each row of the structure matrix
	for ( i = 0; i < matrix_data->number_terms ; i++ ) {
		if ( PetscAbsScalar(matrix_data->term_multipliers[i]) > 0.0 ) {	
			//col_tmp = row->row ^ matrix_data->flips[i]; //this expresion will flip all the spins due to X and Y terms.
			//PetscPrintf(PETSC_COMM_WORLD, "Row: %d, Part: %d, Col: %d\n" , row->row , i, col_tmp ) ;
			col_tmp = row->row;

#if defined(PETSC_USE_COMPLEX)
			val_tmp = 1.0;  
#else
			real_tmp = 1.0;
			imag_tmp = 0.0;
#endif
	    
			//look at the Zs before the flip so we look at 'row' rather than 'col_tmp'
			// If the bit position corresponding to the spin position is occupied then it is in the sz=-1 state. If the bit position offset from this by the number of spins is occupied then it is in the 
			// sz=1 state. If both are unoccupied then it is in the sz=0 state. 			
			for ( j = matrix_data->z_count[i]-1 ; j >=0  ; j--) {	
				if ( matrix_data->zs[i][j] == SPIN_ONE_LOWERING) { 
					//if this is a lowering operator if in sz=-1 state then will give zero. check this now.
					if ( (col_tmp  & matrix_data->ys[i][j]) != zero_const ) { //if this spin is in sz=-1 state.					
	#if defined(PETSC_USE_COMPLEX)				
						val_tmp = 0.0;
	#else				
						real_tmp = 0.0;
					//imag_tmp = -imag_tmp // no point in this as it has to be 0 at this point
	#endif				
					} else if ( (col_tmp  & (matrix_data->ys[i][j] << parameters->number_particles)) == zero_const ) { //if it is in the sz = 0 state then want it to go to -1 state.
						col_tmp += matrix_data->ys[i][j];
					} else { //otherwise it must be in the sz = 1 state. and need to send it to 0 state.  
						col_tmp -= (matrix_data->ys[i][j] << parameters->number_particles);						
					}
				} else if ( matrix_data->zs[i][j] == SPIN_ONE_RAISING) {
					//if this is a lowering operator if in sz=-1 state then will give zero. check this now.
					if ( (col_tmp  & (matrix_data->ys[i][j] << parameters->number_particles))  != zero_const ) { //if this spin is in sz=1 state.					
	#if defined(PETSC_USE_COMPLEX)				
						val_tmp = 0.0;
	#else				
						real_tmp = 0.0;
					//imag_tmp = -imag_tmp // no point in this as it has to be 0 at this point
	#endif				
					} else if ( (col_tmp  & matrix_data->ys[i][j] ) == zero_const ) { //if it is in the sz = 0 state then want it to go to 1 state.
						col_tmp += (matrix_data->ys[i][j] << parameters->number_particles);
					} else { //otherwise it must be in the sz = -1 state. and need to send it to 0 state.  
						col_tmp -= matrix_data->ys[i][j] ;						
					}
				} else if ( matrix_data->zs[i][j] == SPIN_ONE_Z_OPERATOR ) {
					//does not change the column changes value depending on state. 
					if ( (col_tmp  & matrix_data->ys[i][j])  != zero_const ) { //if this spin is in sz=-1 state.					
						#if defined(PETSC_USE_COMPLEX)				
						val_tmp *= -1.0;
						#else				
						real_tmp *= -1.0;					
						#endif				
					} else if ( (col_tmp  & (matrix_data->ys[i][j] << parameters->number_particles))  == zero_const ) { //if this spin is in sz=0 state. No need to do anything about sz=1 state as this multiplies by 1.
						#if defined(PETSC_USE_COMPLEX)				
						val_tmp = 0.0;
						#else				
						real_tmp = 0.0;					
						#endif				
					} 				  
				}
			}	 
			
			
			//apply multiplier for each term
#if defined(PETSC_USE_COMPLEX)	
			val_tmp *=  matrix_data->term_multipliers[i] ;
#else
			real_tmp *= matrix_data->term_multipliers[i] ;
			imag_tmp *= matrix_data->term_multipliers[i];
#endif
	 
			not_found = 1; //indicates not found
			//see if that column already exists. This can be sped up by keeping array sorted by column number  
			for ( j = 0 ; j < row->val_count ; j++ ) {
			//if this column is already there then add the values 
					if (row->cols[j] == (PetscInt)col_tmp ) { 
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
				//row->values[row->val_count] += (PetscScalar)((double)i * PETSC_i);
#else
				row->real_vals[row->val_count] = real_tmp;
				row->imag_vals[row->val_count] = imag_tmp;
#endif
				row->val_count++; 	
			} // end if fabs(term_multipliers) > 0
		} //end loop over terms

		if ( nz == 1 ) { 
			idx = 0;
			for( j = 0 ; j < row->val_count ; j++ ) {

#if defined(PETSC_USE_COMPLEX)
				if ( PetscAbsScalar(row->values[j]) > parameters->tol ) {
#else
				if ( fabs(row->real_vals[j]) > parameters->tol  || row->imag_vals[j] > parameters->tol) {
#endif
					if ( idx != j ) {
						row->cols[idx] = row->cols[j];
#if defined(PETSC_USE_COMPLEX)
						row->values[idx] = row->values[j];
#else
						row->real_vals[idx] = row->real_vals[j];					
						row->imag_vals[idx] = row->imag_vals[j];
#endif
					} //end if idx != j
					idx++;
#if defined(PETSC_USE_COMPLEX)	
				} //end if abs(value) > 0
#else
				} //end if abs(value) > 0
#endif
			}//end for j
			row->val_count = idx;
		} // end if nz == 1 
		
	} //end for loop over parts

	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "getrowlocal_spin_half"
/*! \brief Given a row index calculate the columns and values of the non zero elements in that row for spin half systems for a given range of columns. */
int getrowlocal_spin_half(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters,int nz,PetscInt Istart, PetscInt Iend ) {
	int i, j, not_found; 
	uPetscInt col_tmp, zero_const; 
	PetscReal tmp;
	//char *num;
#if defined(PETSC_USE_COMPLEX)
	PetscScalar val_tmp ; 
#else
	PetscReal real_tmp, imag_tmp;
#endif
//	PetscErrorCode ierr;
	
	//initialise the value counter
	row->val_count = 0 ;
	zero_const = (uPetscInt)0;

	//PetscPrintf(PETSC_COMM_WORLD, "row\n");
	//For each row of the structure matrix
	for ( i = 0; i < matrix_data->number_terms ; i++ ) {
		if ( PetscAbsScalar(matrix_data->term_multipliers[i]) > 0.0 ) {
	
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
		//num = (char *) malloc((matrix_data->number_spins + 1) * sizeof(char));
		ierr = PetscMalloc((matrix_data->number_spins + 1) * sizeof(char),&num);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD, "%s ", dec2bin(row->row, matrix_data->number_spins,num) ); 
		PetscPrintf(PETSC_COMM_WORLD, "& %s ", dec2bin(matrix_data->zs[i][j] , matrix_data->number_spins,num) );
		PetscPrintf(PETSC_COMM_WORLD, "= %s ",dec2bin(row->row  & matrix_data->zs[i][j], sizeof(uPetscInt) * 8 , num )); 
		PetscPrintf(PETSC_COMM_WORLD, "!= %s\n",dec2bin( zero_const, sizeof(uPetscInt) * 8, num )) ;
		//free(num);
		ierr = PetscFree(num);CHKERRQ(ierr);
#endif
				
#if defined(PETSC_USE_COMPLEX)				
					val_tmp = -1.0 * val_tmp;
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
					val_tmp =  -1.0 * val_tmp *  PETSC_i;
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
				if (row->cols[j] == (PetscInt)col_tmp ) {
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
		
		if ( nz == 1 ) { 
			int idx = 0;
			for( j = 0 ; j < row->val_count ; j++ ) {

#if defined(PETSC_USE_COMPLEX)
				if ( PetscAbsScalar(row->values[j]) > 0.0 ) {
#else
				if ( fabs(row->real_vals[j]) > 0.0  || row->imag_vals[j] > 0.0 ) {
#endif
					if ( idx != j ) {
						row->cols[idx] = row->cols[j];
#if defined(PETSC_USE_COMPLEX)
						row->values[idx] = row->values[j];
#else
						row->real_vals[idx] = row->real_vals[j];					
						row->imag_vals[idx] = row->imag_vals[j];
#endif
					} //end if idx != j
					idx++;
#if defined(PETSC_USE_COMPLEX)	
				} //end if abs(value) > 0
#else
				} //end if abs(value) > 0
#endif
			}//end for j
			row->val_count = idx;
		}
		
		
	} //end for loop over parts
	    	
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "getrowlocal_fermionic"
/*! \brief Given a row index calculate the columns and values of the non zero elements in that row for fermionic systems for a given range of columns. */
int getrowlocal_fermionic(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters,int nz,PetscInt Istart, PetscInt Iend ) {
	int i, j, not_found; 
	uPetscInt col_tmp, zero_const; 
	PetscReal tmp;
	//char *num;
#if defined(PETSC_USE_COMPLEX)
	PetscScalar val_tmp ; 
#else
	PetscReal real_tmp, imag_tmp;
#endif
//	PetscErrorCode ierr;
	
	//initialise the value counter
	row->val_count = 0 ;
	zero_const = (uPetscInt)0;

	//PetscPrintf(PETSC_COMM_WORLD, "row\n");
	//For each row of the structure matrix
	for ( i = 0; i < matrix_data->number_terms ; i++ ) {
		if ( PetscAbsScalar(matrix_data->term_multipliers[i]) > 0.0 ) {
	
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
		//num = (char *) malloc((matrix_data->number_spins + 1) * sizeof(char));
		ierr = PetscMalloc((matrix_data->number_spins + 1) * sizeof(char),&num);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD, "%s ", dec2bin(row->row, matrix_data->number_spins,num) ); 
		PetscPrintf(PETSC_COMM_WORLD, "& %s ", dec2bin(matrix_data->zs[i][j] , matrix_data->number_spins,num) );
		PetscPrintf(PETSC_COMM_WORLD, "= %s ",dec2bin(row->row  & matrix_data->zs[i][j], sizeof(uPetscInt) * 8 , num )); 
		PetscPrintf(PETSC_COMM_WORLD, "!= %s\n",dec2bin( zero_const, sizeof(uPetscInt) * 8, num )) ;
		//free(num);
		ierr = PetscFree(num);CHKERRQ(ierr);
#endif
				
#if defined(PETSC_USE_COMPLEX)				
					val_tmp = -1.0 * val_tmp;
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
					val_tmp =  -1.0 * val_tmp *  PETSC_i;
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
				if (row->cols[j] == (PetscInt)col_tmp ) {
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
		
		if ( nz == 1 ) { 
			int idx = 0;
			for( j = 0 ; j < row->val_count ; j++ ) {

#if defined(PETSC_USE_COMPLEX)
				if ( PetscAbsScalar(row->values[j]) > 0.0 ) {
#else
				if ( fabs(row->real_vals[j]) > 0.0  || row->imag_vals[j] > 0.0 ) {
#endif
					if ( idx != j ) {
						row->cols[idx] = row->cols[j];
#if defined(PETSC_USE_COMPLEX)
						row->values[idx] = row->values[j];
#else
						row->real_vals[idx] = row->real_vals[j];					
						row->imag_vals[idx] = row->imag_vals[j];
#endif
					} //end if idx != j
					idx++;
#if defined(PETSC_USE_COMPLEX)	
				} //end if abs(value) > 0
#else
				} //end if abs(value) > 0
#endif
			}//end for j
			row->val_count = idx;
		}
		
		
	} //end for loop over parts
	    	
	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "getrowlocal_spin_one"
/*! \brief Given a row index calculate the columns and values of the non zero elements in that row for spin one systems. */
int getrowlocal_spin_one(struct matrix_row_values *row, struct matrix_meta_data *matrix_data, struct parameter_struct *parameters, int nz, PetscInt Istart, PetscInt Iend ) {
	int i, j, not_found, idx; 
	uPetscInt col_tmp, zero_const; 
	
	//char *num;
#if defined(PETSC_USE_COMPLEX)
	PetscScalar val_tmp ; 
#else
	PetscReal real_tmp, imag_tmp;
#endif
	//initialise the value counter
	row->val_count = 0 ;
	zero_const = (uPetscInt)0;
//	PetscErrorCode ierr;

	//PetscPrintf(PETSC_COMM_WORLD, "row\n");
	//For each row of the structure matrix
	for ( i = 0; i < matrix_data->number_terms ; i++ ) {
		if ( PetscAbsScalar(matrix_data->term_multipliers[i]) > 0.0 ) {	
			//col_tmp = row->row ^ matrix_data->flips[i]; //this expresion will flip all the spins due to X and Y terms.
			//PetscPrintf(PETSC_COMM_WORLD, "Row: %d, Part: %d, Col: %d\n" , row->row , i, col_tmp ) ;
			col_tmp = row->row;


#if defined(PETSC_USE_COMPLEX)
			val_tmp = 1.0;  
#else
			real_tmp = 1.0;
			imag_tmp = 0.0;
#endif
	    
			//look at the Zs before the flip so we look at 'row' rather than 'col_tmp'
			// If the bit position corresponding to the spin position is occupied then it is in the sz=-1 state. If the bit position offset from this by the number of spins is occupied then it is in the 
			// sz=1 state. If both are unoccupied then it is in the sz=0 state. 			
			for ( j = matrix_data->z_count[i]-1 ; j >=0  ; j--) {	
				if ( matrix_data->zs[i][j] == SPIN_ONE_LOWERING) { 
					//if this is a lowering operator if in sz=-1 state then will give zero. check this now.
					if ( (col_tmp  & matrix_data->ys[i][j]) != zero_const ) { //if this spin is in sz=-1 state.					
	#if defined(PETSC_USE_COMPLEX)				
						val_tmp = 0.0;
	#else				
						real_tmp = 0.0;
					//imag_tmp = -imag_tmp // no point in this as it has to be 0 at this point
	#endif				
					} else if ( (col_tmp  & (matrix_data->ys[i][j] << parameters->number_particles)) == zero_const ) { //if it is in the sz = 0 state then want it to go to -1 state.
						col_tmp += matrix_data->ys[i][j];
					} else { //otherwise it must be in the sz = 1 state. and need to send it to 0 state.  
						col_tmp -= (matrix_data->ys[i][j] << parameters->number_particles);						
					}
				} else if ( matrix_data->zs[i][j] == SPIN_ONE_RAISING) {
					//if this is a lowering operator if in sz=-1 state then will give zero. check this now.
					if ( (col_tmp  & (matrix_data->ys[i][j] << parameters->number_particles))  != zero_const ) { //if this spin is in sz=1 state.					
	#if defined(PETSC_USE_COMPLEX)				
						val_tmp = 0.0;
	#else				
						real_tmp = 0.0;
					//imag_tmp = -imag_tmp // no point in this as it has to be 0 at this point
	#endif				
					} else if ( (col_tmp  & (matrix_data->ys[i][j] << parameters->number_particles)) == zero_const ) { //if it is in the sz = 0 state then want it to go to 1 state.
						col_tmp += (matrix_data->ys[i][j] << parameters->number_particles);
					} else { //otherwise it must be in the sz = -1 state. and need to send it to 0 state.  
						col_tmp -= matrix_data->ys[i][j] ;						
					}
				} else if ( matrix_data->zs[i][j] == SPIN_ONE_Z_OPERATOR ) {
					//does not change the column changes value depending on state. 
					if ( (col_tmp  & matrix_data->ys[i][j])  != zero_const ) { //if this spin is in sz=-1 state.					
						#if defined(PETSC_USE_COMPLEX)				
						val_tmp *= -1.0;
						#else				
						real_tmp *= -1.0;					
						#endif				
					} else if ( (col_tmp  & (matrix_data->ys[i][j] << parameters->number_particles))  == zero_const ) { //if this spin is in sz=0 state. No need to do anything about sz=1 state as this multiplies by 1.
						#if defined(PETSC_USE_COMPLEX)				
						val_tmp = 0.0;
						#else				
						real_tmp = 0.0;					
						#endif				
					} 				  
				}
			}	 
			
			
			//apply multiplier for each term
#if defined(PETSC_USE_COMPLEX)	
			val_tmp *=  matrix_data->term_multipliers[i] ;
#else
			real_tmp *= matrix_data->term_multipliers[i] ;
			imag_tmp *= matrix_data->term_multipliers[i];
#endif
	 
			not_found = 1; //indicates not found
			//see if that column already exists. This can be sped up by keeping array sorted by column number  
			for ( j = 0 ; j < row->val_count ; j++ ) {
			//if this column is already there then add the values 
					if (row->cols[j] == (PetscInt)col_tmp ) { 
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
			if ( not_found && (PetscInt)col_tmp >= Istart && (PetscInt)col_tmp < Iend ) {
				//PetscPrintf(PETSC_COMM_WORLD, "Adding new nz value %lf at col %d, this val num %d, multiplier %lf\n", real_tmp, col_tmp, row->val_count, matrix_data->term_multipliers[i] ); 
				row->cols[row->val_count] = col_tmp; 
#if defined(PETSC_USE_COMPLEX)
				row->values[row->val_count] = val_tmp;
				//row->values[row->val_count] += (PetscScalar)((double)i * PETSC_i);
#else
				row->real_vals[row->val_count] = real_tmp;
				row->imag_vals[row->val_count] = imag_tmp;
#endif
				row->val_count++; 	
			} // end if fabs(term_multipliers) > 0
		} //end loop over terms

		if ( nz == 1 ) { 
			idx = 0;
			for( j = 0 ; j < row->val_count ; j++ ) {

#if defined(PETSC_USE_COMPLEX)
				if ( PetscAbsScalar(row->values[j]) > parameters->tol ) {
#else
				if ( fabs(row->real_vals[j]) > parameters->tol  || row->imag_vals[j] > parameters->tol) {
#endif
					if ( idx != j ) {
						row->cols[idx] = row->cols[j];
#if defined(PETSC_USE_COMPLEX)
						row->values[idx] = row->values[j];
#else
						row->real_vals[idx] = row->real_vals[j];					
						row->imag_vals[idx] = row->imag_vals[j];
#endif
					} //end if idx != j
					idx++;
#if defined(PETSC_USE_COMPLEX)	
				} //end if abs(value) > 0
#else
				} //end if abs(value) > 0
#endif
			}//end for j
			row->val_count = idx;
		} // end if nz == 1 
		
	} //end for loop over parts

	return 0;
}



#undef __FUNCT__
#define __FUNCT__ "allocate_matrix_data"
/*! \brief Allocate space to store data on how to construct the matrix.*/
int allocate_matrix_data(struct matrix_meta_data *matrix_data) {
	int i; 
	PetscErrorCode ierr;
	

	//allocate storage for masks
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(uPetscInt),&matrix_data->flips);CHKERRQ(ierr);
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(int),&matrix_data->y_count);CHKERRQ(ierr);
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(int),&matrix_data->z_count);CHKERRQ(ierr);
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(int),&matrix_data->x_count);CHKERRQ(ierr);
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(uPetscInt*),&matrix_data->ys);CHKERRQ(ierr);
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(uPetscInt*),&matrix_data->zs);CHKERRQ(ierr);
//	ierr = PetscMalloc(matrix_data->number_terms * sizeof(double),&matrix_data->term_multipliers);CHKERRQ(ierr);
//	ierr = PetscMalloc(matrix_data->number_terms * sizeof(double),&matrix_data->imag_term_multipliers);CHKERRQ(ierr);
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(PetscScalar),&matrix_data->term_multipliers);CHKERRQ(ierr);
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(int),&matrix_data->number_term_params);CHKERRQ(ierr);
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(int*),&matrix_data->term_parameters);CHKERRQ(ierr);
	ierr = PetscMalloc(matrix_data->number_parameters * sizeof(int),&matrix_data->parameter_indices);CHKERRQ(ierr);

	//allocate storage for parameters and terms
	//ierr = PetscMalloc(matrix_data->number_parameters * sizeof(char *),&matrix_data->parameter_names);CHKERRQ(ierr);
	//for (i = 0 ; i < matrix_data->number_parameters ; i++ ) {
	//	ierr = PetscMalloc(MAX_PARM * sizeof(char),&matrix_data->parameter_names[i]);CHKERRQ(ierr);
	//}
	
	//matrix_data->terms = (char **)malloc(matrix_data->number_terms * sizeof(char *) ) ;
	ierr = PetscMalloc(matrix_data->number_terms * sizeof(char*),&matrix_data->terms);CHKERRQ(ierr);
	for (i = 0 ; i < matrix_data->number_terms ; i++) {
		ierr = PetscMalloc(MAX_TERM * sizeof(char),&matrix_data->terms[i]);CHKERRQ(ierr);
	}
	
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "deallocate_matrix_data"
/*! \brief Deallocate space used to store matrix details.*/
int deallocate_matrix_data(struct matrix_meta_data *matrix_data) {
	int i; 
	PetscErrorCode ierr;

	if (matrix_data->number_terms >= 1 ) {
		for ( i = 0 ; i < matrix_data->number_terms ; i++ ) {
			ierr = PetscFree(matrix_data->ys[i]);CHKERRQ(ierr);
			ierr = PetscFree(matrix_data->zs[i]);CHKERRQ(ierr);
			if ( matrix_data->number_parameters > 0 ) {
			  ierr = PetscFree(matrix_data->term_parameters[i]);CHKERRQ(ierr);
			}
		}
	}
	
	//if ( matrix_data->number_parameters >= 1 ) {
	//	for ( i = 0 ; i < matrix_data->number_parameters ; i++ ){
	//		ierr = PetscFree(matrix_data->parameter_names[i]);CHKERRQ(ierr);
	//	}
	//	ierr = PetscFree(matrix_data->parameter_names);CHKERRQ(ierr);
	//}
	

	//ierr = PetscFree(matrix_data->terms);CHKERRQ(ierr); //already deallocated after used. 
	ierr = PetscFree(matrix_data->flips);CHKERRQ(ierr);
	ierr = PetscFree(matrix_data->term_multipliers);CHKERRQ(ierr);
	//ierr = PetscFree(matrix_data->imag_term_multipliers);CHKERRQ(ierr);
	ierr = PetscFree(matrix_data->ys);CHKERRQ(ierr);
	ierr = PetscFree(matrix_data->zs);CHKERRQ(ierr);
	ierr = PetscFree(matrix_data->x_count);CHKERRQ(ierr);
	ierr = PetscFree(matrix_data->y_count);CHKERRQ(ierr);
	ierr = PetscFree(matrix_data->z_count);CHKERRQ(ierr);
	ierr = PetscFree(matrix_data->term_parameters);CHKERRQ(ierr);
	ierr = PetscFree(matrix_data->number_term_params);CHKERRQ(ierr);
	ierr = PetscFree(matrix_data->parameter_indices);CHKERRQ(ierr);

	
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "MatrixDistribute"
/*! \brief This is to be a general purpose matrix distribution function that is capable of using all the symmetries and features.
	All details relating to the operator terms will be gotten from the matrix_data structure and all other parameters from the parameters structure. */
PetscInt MatrixDistribute(Mat *A,struct parameter_struct *parameters, struct matrix_meta_data *matrix_data) {
	
	//PetscInt perfect_mapping_sector ;
	PetscErrorCode ierr;

	//set up which sector we are in.
	parameters->real_parity_sector  = parameters->parity_sector;
	parameters->real_momentum_sector  = parameters->momentum_sector;
	parameters->real_conservation_sector = parameters->conservation_sector;
	if ( parameters->use_rotation_invariance)
	{
		if ( parameters->rotation_number_relevant_sectors > 0 ) 
		{		  		 
		      set_rotation_sector(parameters, parameters->rotation_relevant_sectors[parameters->rotation_current_sector]);
		} 
		else 
		{
		      set_rotation_sector(parameters, parameters->rotation_current_sector);
		}
	}

	sprintf(parameters->current_output_additional_prefix,"task_%d",(int)parameters->current_task);
	if (parameters->use_parity_sectors == PETSC_TRUE ) {
		if ( parameters->parity_info.num_relevant_sectors > 0 ) { 
			parameters->real_parity_sector = parameters->parity_info.relevant_sectors[parameters->parity_sector];
		}
		char tmp_buf[MAX_STRING];
		strcpy(tmp_buf,parameters->current_output_additional_prefix);
		sprintf(parameters->current_output_additional_prefix,"%s_parity_%d",tmp_buf,(int)parameters->real_parity_sector);
	}
	if (parameters->use_conservation_sectors == PETSC_TRUE ) {
		if ( parameters->parity_info.num_relevant_sectors > 0 ) { 
			parameters->real_conservation_sector = parameters->parity_info.relevant_sectors[parameters->conservation_sector];
		}
		char tmp_buf[MAX_STRING];
		strcpy(tmp_buf,parameters->current_output_additional_prefix);
		sprintf(parameters->current_output_additional_prefix,"%s_filling_%d",tmp_buf,(int)parameters->real_conservation_sector);
	}

	if ( parameters->use_momentum_sectors == PETSC_TRUE )  {
		if ( parameters->momentum_info.num_relevant_sectors > 0 ) { 
			parameters->real_momentum_sector = parameters->momentum_info.relevant_sectors[parameters->momentum_sector];
		}
		parameters->momentum_info.current_sector[0] = parameters->real_momentum_sector / parameters->momentum_info.num_sectors[1];
		parameters->momentum_info.current_sector[1] = parameters->real_momentum_sector % parameters->momentum_info.num_sectors[1];
	}
	
	//this needs to be after the momentum sectors have been decided but before the phases have been calculated.
	if ( parameters->use_spectral_flow ) {
		setup_spectral_flow_info(parameters,matrix_data);
		char tmp_buf[MAX_STRING];
		strcpy(tmp_buf,parameters->current_output_additional_prefix);
		sprintf(parameters->current_output_additional_prefix,"%s_spectral_%d",tmp_buf,(int)(parameters->spectral_info.real_current_point[1]*parameters->spectral_info.number_points[0]+parameters->spectral_info.real_current_point[0]));
	}
	
	if ( parameters->use_momentum_sectors == PETSC_TRUE )  {
		char tmp_buf[MAX_STRING];
		strcpy(tmp_buf,parameters->current_output_additional_prefix);
		sprintf(parameters->current_output_additional_prefix,"%s_momentum_%d",tmp_buf,(int)parameters->real_momentum_sector);
		parameters->momentum_info.phases = calculate_translation_phases(parameters);
		calculate_translation_masks(parameters);
	}
	
	if ( parameters->use_rotation_invariance == PETSC_TRUE ) {
		char tmp_buf[MAX_STRING];
		/*if ( parameters->rotation_info.num_relevant_sectors < parameters->rotation_info.num_sectors ) {
			parameters->rotation_info.real_sector = parameters->rotation_info.relevant_sectors[parameters->rotation_info.current_sector] ; 
		} else {
			parameters->rotation_info.real_sector = parameters->rotation_info.current_sector; 
		}*/
		strcpy(tmp_buf,parameters->current_output_additional_prefix);
		if ( parameters->rotation_number_relevant_sectors > parameters->rotation_current_sector ) 
		{
			sprintf(parameters->current_output_additional_prefix,"%s_rotation_%d",tmp_buf,(int)parameters->rotation_relevant_sectors[parameters->rotation_current_sector]);
		}
		else
		{
		  	sprintf(parameters->current_output_additional_prefix,"%s_rotation_%d",tmp_buf,(int)parameters->rotation_current_sector);
		}
		for ( int rot_op_idx = 0 ; rot_op_idx < parameters->number_rotation_ops ; rot_op_idx++ )
		{
			calculate_rotation_phases(parameters, rot_op_idx);
		}
	}
		
	uPetscInt *local_reps_array, *ending_idx_for_proc;	
	
	if ( parameters->use_disk == PETSC_TRUE ) {
		parameters->current_basis_size = build_local_basis_lists_using_disk(parameters, &local_reps_array, &ending_idx_for_proc); 
	} else {
		//Now if we have any thing that restrict the basis we need to build local lists
		parameters->current_basis_size = build_local_basis_lists(parameters, &local_reps_array, &ending_idx_for_proc); 
	}	
	
	if ( parameters->verbosity >= 15 ) {
		for ( int i = 0 ; i < get_rank_chunk_size(parameters->current_basis_size,parameters->rank,parameters->size); i++ ) {
			stringstream ss;
			char num[MAX_STRING];
			dec2bin((*parameters->get_full_basis_index)(parameters,local_reps_array[i]),parameters->number_particles,num);
			ss << i+get_rank_chunk_start(parameters->current_basis_size,parameters->rank,parameters->size) << ": " << num << endl ;
			PetscPrintf(PETSC_COMM_SELF,"%s",ss.str().c_str());
		}
	}

	if ( parameters->save_basis ) {
		stringstream basis_ss;
		basis_ss << parameters->output_prefix << "." << parameters->current_output_additional_prefix << ".basis" ;
		for ( int working_rank = 0 ; working_rank < parameters->size ; working_rank ++ ) 
		{	 
			stringstream ss;
			if ( parameters->rank == working_rank)
			{
			    FILE *fp;
			    if ( parameters->rank == 0 )		      
			    {
				fp = fopen(basis_ss.str().c_str(),"w");
				ss.str("");ss << "BASIS_SIZE = " << parameters->current_basis_size << endl;
				PetscFPrintf(PETSC_COMM_SELF, fp, "%s", ss.str().c_str());
			    }
			    else
			    {
			      fp = fopen(basis_ss.str().c_str(),"a");		      		      		      
			    }
			    
			    if ( parameters->use_momentum_sectors  || parameters->nn_exclusion || parameters->use_rotation_invariance )
			    {
				for ( int i = 0 ; i < get_rank_chunk_size(parameters->current_basis_size,parameters->rank,parameters->size); i++ )
				{
				  ss.str("");ss << (*parameters->get_full_basis_index)(parameters,local_reps_array[i]) << endl;
				  PetscFPrintf(PETSC_COMM_SELF, fp, "%s", ss.str().c_str());		      		      
				}
			    } 
			    else 
			    {
				for ( int i = 0 ; i < get_rank_chunk_size(parameters->current_basis_size,parameters->rank,parameters->size); i++ )
				{
				  ss.str("");ss << (*parameters->get_full_basis_index)(parameters,i) << endl;
				  PetscFPrintf(PETSC_COMM_SELF, fp, "%s", ss.str().c_str());		      		      
				}
			    }
			    fclose(fp);		    		  
			}				
			MPI_Barrier(PETSC_COMM_WORLD);
		}
	}
	
	
	if ( parameters->verbosity >= 1 ) { 
			stringstream ss;
			ss << "Number of basis elements: " << parameters->current_basis_size << endl ;
			PetscPrintf(PETSC_COMM_WORLD,"%s", ss.str().c_str());
	}
	
	PetscInt *d_nnz, *o_nnz;
	uPetscInt *reps,*indices;
	if ( parameters->current_basis_size > 0 && get_rank_chunk_size(parameters->current_basis_size,parameters->rank,parameters->size) > 0 ) { 
		uPetscInt num_unique_cols;
		//this function builds the lists of required basis elements for each process and then retrieves the indices from the relevant process.
		exchange_required_basis_indices(parameters,matrix_data,local_reps_array,ending_idx_for_proc,&d_nnz,&o_nnz,&reps,&indices,&num_unique_cols);
		
		ierr = MatCreate(parameters->mat_solver_comm,A);CHKERRQ(ierr);
		ierr = MatSetSizes(*A, get_rank_chunk_size(parameters->current_basis_size,parameters->rank,parameters->size), PETSC_DECIDE, parameters->current_basis_size ,parameters->current_basis_size);CHKERRQ(ierr);
		if ( parameters->use_hermitian == PETSC_TRUE ) {
			ierr = MatSetType(*A,MATMPISBAIJ);CHKERRQ(ierr);
			ierr = MatSetOption(*A,MAT_GETROW_UPPERTRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
		} else { 
			ierr = MatSetType(*A,MATMPIAIJ);CHKERRQ(ierr);
		}
		MatSetFromOptions(*A);CHKERRQ(ierr);
		MatMPIAIJSetPreallocation(*A,1,d_nnz ,0, o_nnz );CHKERRQ(ierr);
		ierr = PetscFree(d_nnz);CHKERRQ(ierr);
		ierr = PetscFree(o_nnz);CHKERRQ(ierr);
		//done with the preallocation now need to set the actual value
		
		set_matrix_elements(A,parameters,matrix_data,local_reps_array,ending_idx_for_proc,reps,indices,num_unique_cols);
		
		//Create and set the matrix for the correlators
		if ( parameters->calculate_correlations == PETSC_TRUE || parameters->prod_wf_overlap == PETSC_TRUE) {
			parameters->basis_info.local_reps_array = local_reps_array;  //save pointers to the basis information so that it can be used later.
			parameters->basis_info.ending_idx_for_proc = ending_idx_for_proc; 		 
		}
		
		if ( parameters->use_momentum_sectors == PETSC_TRUE || parameters->nn_exclusion == PETSC_TRUE || parameters->use_rotation_invariance == PETSC_TRUE ) {
			if ( parameters->calculate_correlations == PETSC_FALSE && parameters->prod_wf_overlap == PETSC_FALSE) { //if not calculating correlations then basis information is not required later so can be deallocated.
				ierr = PetscFree(local_reps_array);CHKERRQ(ierr);
				ierr = PetscFree(ending_idx_for_proc);CHKERRQ(ierr);
			}
			ierr = PetscFree(reps);CHKERRQ(ierr);
			ierr = PetscFree(indices);CHKERRQ(ierr);
		}
		
		if ( parameters->use_momentum_sectors == PETSC_TRUE )  {
			ierr = PetscFree(parameters->momentum_info.phases);CHKERRQ(ierr);
			for ( int i = 0 ; i < parameters->momentum_info.num_sectors[0]; i++ ) {
				for ( int j = 0 ; j < parameters->momentum_info.num_sectors[1]; j++ ) {
					ierr = PetscFree(parameters->momentum_info.translation_masks[i][j]);CHKERRQ(ierr);	
				}
				ierr = PetscFree(parameters->momentum_info.translation_masks[i]);CHKERRQ(ierr);
			}
			ierr = PetscFree(parameters->momentum_info.translation_masks);CHKERRQ(ierr);
		}
	} //end if basis > 0
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	return parameters->current_basis_size;
}

#undef __FUNCT__
#define __FUNCT__ "build_local_basis_lists_using_disk"
/*! \brief This function will do the same thing that build local basis lists does except it will use the local disk to keep the list temporarily. */
PetscInt build_local_basis_lists_using_disk(struct parameter_struct *parameters, uPetscInt **local_reps_array, uPetscInt **ending_idx_for_proc){
	//Determing the size of the possible basis set.
	PetscErrorCode ierr;
	PetscInt i;
	uPetscInt possible_basis_size = (uPetscInt)parameters->numberBasisStates; //if no restriction is present where basis size can be calculated.
	if ( parameters->use_conservation_sectors == PETSC_TRUE ) {
		possible_basis_size = get_size_of_conservation_sector(parameters,(uPetscInt)parameters->real_conservation_sector);
	} else if ( parameters->use_parity_sectors == PETSC_TRUE ) {
		possible_basis_size = get_size_of_parity_sector(parameters,(uPetscInt)parameters->real_parity_sector);
	}
	parameters->current_basis_size = possible_basis_size;
	
	PetscInt n;
	PetscInt *all_local_reps;
	PetscReal starttime1;
	uPetscInt mapped_idx;
	PetscInt local_reps ;
	PetscReal additionaltime;
	string basis_file_prefix("basis_list_portion_");
	if ( parameters->use_momentum_sectors  || parameters->nn_exclusion || parameters->use_rotation_invariance  ) {
		local_reps = 0 ;
		starttime1  = MPI_Wtime();
		
		stringstream ss;
		ss << basis_file_prefix << parameters->rank << ".bas" ;
		string basis_filename = ss.str();
		ofstream basis_file ( basis_filename.c_str(),ios::out | ios::binary);
		
		if ( parameters->nn_exclusion  && parameters->use_disk  && parameters->nn_recursive_algorithm ) {
			find_valid_configs(parameters, &basis_file, &local_reps);
		} else {
			//loop over local chunk of full  basis set and count local reps
			for ( uPetscInt idx = (uPetscInt)get_rank_chunk_start(possible_basis_size, parameters->rank, parameters->size) ; 
				idx < (uPetscInt)get_rank_chunk_start(possible_basis_size, parameters->rank, parameters->size)+(uPetscInt)get_rank_chunk_size(possible_basis_size, parameters->rank, parameters->size); 
				idx++ ) {
				mapped_idx = parameters->get_full_basis_index(parameters,idx);
				// Next we check to ensure that the configuration satisfies ony conditions that have been specified. For example the nearest neighbour exclusion principle.
				// Also if the momentum resolution is being used we check if the given configuration is a representative. 
				if ( (*parameters->simple_basis_check)(parameters,mapped_idx) == PETSC_TRUE && (*parameters->representative_check)(parameters,mapped_idx) == PETSC_TRUE ){
					basis_file.write((char*)&idx,sizeof(uPetscInt));
					local_reps ++;
				}
			}
		}
		
		basis_file.close();
		additionaltime = MPI_Wtime();
		
		if ( parameters->verbosity >= 10 ) {
			stringstream ss;
			ss << "Rank " << parameters->rank << ": Time taken to run through all reps in local partition (" << local_reps << ") " << additionaltime - starttime1 << endl ;
			PetscPrintf(PETSC_COMM_SELF,ss.str().c_str());
		}
		
	
		//not needed anymore since we can get this from the all gather
		//sum number of representatives on each process to get overall number of representatives
		MPI_Allreduce(&local_reps, &n, 1, PETSC_MPI_INT, MPI_SUM, PETSC_COMM_WORLD); //sum the number of reps found on each process to find new basis size.
		//now know the size of the basis set. If it is zero then we just move on and make sure we deallocate everything right. If it is one we can 
		//go ahead with the construction but not pass to the solver. 
		
		
		//adding one extra element to end.
		ierr = PetscMalloc((parameters->size+1)*sizeof(PetscInt),&all_local_reps); 
		MPI_Allgather(&local_reps,1,PETSC_MPI_INT,&all_local_reps[1],1,PETSC_MPI_INT,PETSC_COMM_WORLD); //gather the number of reps found on each process.
		
		parameters->current_basis_size = n;
	} //end if using momentum or nn exculsion.
		
	create_matrix_communicator(parameters);//creats a new communicator with only the processes that are to be used in it.
	
	if ( parameters->use_momentum_sectors  || parameters->nn_exclusion || parameters->use_rotation_invariance ) {	
		if ( get_rank_chunk_size(n,parameters->rank,parameters->size) > 0 ) { //if there are any basis elements for the current process.
			
			PetscReal additionaltime2 = MPI_Wtime();
			if ( parameters->verbosity >= 10) {
				PetscPrintf(PETSC_COMM_WORLD,"Time for all gather: %lf.\n", additionaltime2 - additionaltime);
			}
			
			//allocate space to store chunk on each process
			//PetscPrintf(PETSC_COMM_SELF,"Rank %d, local reps %d and total %d and chunk is %d.\n", parameters->rank, local_reps,n, get_rank_chunk_size(n,parameters->rank,parameters->size));	
			ierr = PetscMalloc(get_rank_chunk_size(n,parameters->rank,parameters->size)*sizeof(uPetscInt),local_reps_array);CHKERRQ(ierr);
			//ierr = PetscMalloc(get_rank_chunk_size(n,parameters->rank,parameters->size)*sizeof(PetscReal),&local_reps_normals);CHKERRQ(ierr);
			ierr = PetscMalloc(parameters->size*sizeof(uPetscInt),ending_idx_for_proc);CHKERRQ(ierr);
			
			all_local_reps[0] = 0; 
			for ( i = 1 ; i <= parameters->size ; i++ ) {
				all_local_reps[i] += all_local_reps[i-1];
			}
			
			PetscInt offset = 0;
			//need to now find which process and offset to start at.			
			PetscInt to_fill = get_rank_chunk_size(n,parameters->rank,parameters->size); //keep account of how many array places are left to fill. 
			PetscInt proc = 0; //start at the first process.
			while (to_fill > 0 ) {
					while ( (get_rank_chunk_start(n,parameters->rank,parameters->size)+offset) >= all_local_reps[proc+1] ) { //while the starting index is not within this range move on.
						proc++; 
					}
					//open the file corresponding to processes proc
					stringstream ss;
					ss << basis_file_prefix << proc << ".bas" ;
					string basis_filename = ss.str();
					ifstream basis_file ( basis_filename.c_str(),ios::in | ios::binary );
					//seek if we need to
					PetscInt available = all_local_reps[proc+1] - all_local_reps[proc];
					if ( (get_rank_chunk_start(n,parameters->rank,parameters->size)+offset) > all_local_reps[proc] ) {
						basis_file.seekg((get_rank_chunk_start(n,parameters->rank,parameters->size)+offset-all_local_reps[proc]) * sizeof(uPetscInt),ios_base::beg);
						available -= (get_rank_chunk_start(n,parameters->rank,parameters->size)+offset-all_local_reps[proc]);
					}
					//now read in as many elements as we want
					if ( to_fill < available ) {
						basis_file.read((char *) &(*local_reps_array)[offset], sizeof(uPetscInt) * to_fill );
						to_fill -= to_fill;
						offset += to_fill;
					} else {
						basis_file.read((char *) &(*local_reps_array)[offset], sizeof(uPetscInt) * available );
						offset += available;
						to_fill -= available;
					}
					basis_file.close();
			}
			
			
			//now must sort the distributed array			
			if ( parameters->nn_exclusion && parameters->use_disk && parameters->nn_recursive_algorithm ) {
				parallel_qsort_uPetscInt(parameters,&(parameters->mat_solver_comm),*local_reps_array,get_rank_chunk_size(n,parameters->rank,parameters->size));
				MPI_Barrier(parameters->mat_solver_comm);
			}
			
			if ( get_rank_chunk_size(n,parameters->rank,parameters->size) == 0 ) {
				uPetscInt tmp =0;
				MPI_Allgather(&tmp,1,PETSC_UNSIGNED_MPI_INT,*ending_idx_for_proc,1,PETSC_UNSIGNED_MPI_INT,parameters->mat_solver_comm);
			} else {
				MPI_Allgather(&(*local_reps_array)[get_rank_chunk_size(n,parameters->rank,parameters->size)-1],1,PETSC_UNSIGNED_MPI_INT,*ending_idx_for_proc,1,PETSC_UNSIGNED_MPI_INT,parameters->mat_solver_comm);
			}
			
			PetscReal endtime1 = MPI_Wtime();
			if ( parameters->verbosity >= 5 ) {
				PetscPrintf(PETSC_COMM_WORLD,"Time taken for filling local reps array: %lf\n",endtime1-starttime1);
			}
			
		}  // end if n > 0
		PetscFree(all_local_reps);
		return n;
	} else {
		return possible_basis_size;
	}
	return 0 ;	
		

}

#undef __FUNCT__
#define __FUNCT__ "build_local_basis_lists"
/*! \brief This function prepares the lists of valid basis configurations that is to be used. */
PetscInt build_local_basis_lists(struct parameter_struct *parameters, uPetscInt **local_reps_array, uPetscInt **ending_idx_for_proc){
	//Determing the size of the possible basis set.
	PetscErrorCode ierr;
	PetscInt  i, j ;
	uPetscInt possible_basis_size = (uPetscInt)parameters->numberBasisStates; //if no restriction is present where basis size can be calculated.
	if ( parameters->use_conservation_sectors == PETSC_TRUE ) {
		possible_basis_size = get_size_of_conservation_sector(parameters,(uPetscInt)parameters->real_conservation_sector);
	} else if ( parameters->use_parity_sectors == PETSC_TRUE ) {
		possible_basis_size = get_size_of_parity_sector(parameters,(uPetscInt)parameters->real_parity_sector);
	}
	parameters->current_basis_size = possible_basis_size;
	
	PetscInt n;
	PetscInt *all_local_reps;
	PetscReal starttime1;
	uPetscInt mapped_idx;
	PetscInt local_reps ;
	PetscReal additionaltime;
	if ( parameters->use_momentum_sectors == PETSC_TRUE || parameters->nn_exclusion == PETSC_TRUE || parameters->use_rotation_invariance == PETSC_TRUE ) {
		local_reps = 0 ;
		starttime1  = MPI_Wtime();
		
		//loop over local chunk of full  basis set and count local reps
		for ( uPetscInt idx = (uPetscInt)get_rank_chunk_start(possible_basis_size, parameters->rank, parameters->size) ; 
			idx < (uPetscInt)get_rank_chunk_start(possible_basis_size, parameters->rank, parameters->size)+(uPetscInt)get_rank_chunk_size(possible_basis_size, parameters->rank, parameters->size); 
			idx++ ) {
			mapped_idx = parameters->get_full_basis_index(parameters,idx);
			if ( (*parameters->simple_basis_check)(parameters,mapped_idx) == PETSC_TRUE && (*parameters->representative_check)(parameters,mapped_idx) == PETSC_TRUE ){
				local_reps ++;
			}
		}
		
		additionaltime = MPI_Wtime();
		
		if ( parameters->verbosity >= 10 ) {
			stringstream ss;
			ss << "Rank " << parameters->rank << ": Time taken to run through all reps in local partition (" << local_reps << ") " << additionaltime - starttime1 << endl ;
			PetscPrintf(PETSC_COMM_SELF,ss.str().c_str());
		}
	
		
		//not needed anymore since we can get this from the all gather
		//sum number of representatives on each process to get overall number of representatives
		MPI_Allreduce(&local_reps, &n, 1, PETSC_MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
		//now know the size of the basis set. If it is zero then we just move on and make sure we deallocate everything right. If it is one we can 
		//go ahead with the construction but not pass to the solver. 

		
		
		//adding one extra element to end.
		ierr = PetscMalloc((parameters->size+1)*sizeof(PetscInt),&all_local_reps); 
		MPI_Allgather(&local_reps,1,PETSC_MPI_INT,&all_local_reps[1],1,PETSC_MPI_INT,PETSC_COMM_WORLD);
		
		
		parameters->current_basis_size = n;
	}
		
	create_matrix_communicator(parameters);
	
	if ( parameters->use_momentum_sectors == PETSC_TRUE || parameters->nn_exclusion == PETSC_TRUE || parameters->use_rotation_invariance == PETSC_TRUE  ) {	
		if ( get_rank_chunk_size(n,parameters->rank,parameters->size) > 0 ) { //if there are any basis elements 
			//PetscInt *all_local_reps;
			//adding one extra element to end.
			//ierr = PetscMalloc((parameters->restricted_size+1)*sizeof(PetscInt),&all_local_reps); 
			//MPI_Allgather(&local_reps,1,PETSC_MPI_INT,&all_local_reps[1],1,PETSC_MPI_INT,parameters->mat_solver_comm);
			
			
			PetscReal additionaltime2 = MPI_Wtime();
			if ( parameters->verbosity >= 10) {
				PetscPrintf(PETSC_COMM_WORLD,"Time for all gather %lf.\n", additionaltime2 - additionaltime);
			}
			
			PetscInt *range_contains; //array that will keep track of the number of representatives contained in each range so we can exclude the zero ranges later.
			ierr = PetscMalloc((parameters->size+1)*sizeof(PetscInt),&range_contains);
			for ( i = 0 ; i < parameters->size ; i++ ) {
				range_contains[i] = all_local_reps[i+1];
			}
			
			all_local_reps[0] = 0; 
			for ( i = 1 ; i <= parameters->size ; i++ ) {
				all_local_reps[i] += all_local_reps[i-1];
			}
			//n = all_local_reps[parameters->size];
			
			
		
			//now we have the number of representatives in each range do a second iteration to improve load balancing. Assuming representatives are dispersed 
			//evenly within each range. 
			//find the start of the new range
			PetscInt *new_chunks;
			ierr = PetscMalloc((parameters->size+1)*sizeof(PetscInt),&new_chunks); 
			
			new_chunks[parameters->size] = possible_basis_size; // last chunk always ends on before last element.
			for ( PetscInt rank = parameters ->size - 1 ; rank >= 0 ; rank-- ) {
				PetscInt new_start = get_rank_chunk_start(n,rank,parameters->size); //get the start index of this rank in the reduced basis list.
				PetscInt new_full_start; 
				
				for ( i = 0 ; i < parameters->size ; i++ ) { 
					if ( new_start >= all_local_reps[i] && new_start < all_local_reps[i+1] ){
						new_full_start = get_rank_chunk_start(possible_basis_size,i,parameters->size) + ((double)(new_start - all_local_reps[i])*(double)(get_rank_chunk_size(possible_basis_size,i,parameters->size))/(double)(all_local_reps[i+1]-all_local_reps[i]));
						
					} else if (get_rank_chunk_size(n,rank,parameters->size) == 0){
						new_full_start = possible_basis_size;
					}
				}
				
				new_chunks[rank] = new_full_start;
				
				/*if ( rank == ( parameters->size - 1 ) ) {
					new_full_end = matrix_data->numberBasisStates;
				} else {
					for ( i = 0 ; i < parameters->size ; i++ ) { 
						if ( new_end >= all_local_reps[i] && new_end < all_local_reps[i+1] ){
							new_full_end = get_rank_chunk_start(matrix_data->numberBasisStates,i,parameters->size) + ((new_end - all_local_reps[i])*get_rank_chunk_size(matrix_data->numberBasisStates,i,parameters->size))/(all_local_reps[i+1]-all_local_reps[i]);
						}
					}
				}*/
				if ( parameters->verbosity >= 10 ) {
					stringstream ss;
					ss << "Rank " << rank << ": new range from " << new_chunks[rank] << " to " << new_chunks[rank+1] << endl;
					//PetscPrintf(PETSC_COMM_WORLD,"Rank %d: new range from %d to %d.\n", rank,new_chunks[rank],new_chunks[rank+1]);
					PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
				}
			}
			
			
			//loop over local chunk of full  basis set and count local reps
			local_reps = 0 ;  
			uPetscInt idx = (uPetscInt)new_chunks[parameters->rank] ;
			while ( idx < (uPetscInt)new_chunks[parameters->rank+1] ) {
				int range = get_range_with_item(possible_basis_size,(PetscInt)idx,parameters->size);
				//PetscPrintf(PETSC_COMM_SELF,"%d is in range %d from %d to %d.\n",j,range,get_rank_chunk_start(matrix_data->numberBasisStates,range,parameters->size),get_rank_chunk_start(matrix_data->numberBasisStates,range,parameters->size)+get_rank_chunk_size(matrix_data->numberBasisStates,range,parameters->size));
				//check if the range is fully contained in the new chunk
				if ( (uPetscInt)get_rank_chunk_start(possible_basis_size,range,parameters->size) >= idx &&  (get_rank_chunk_start(possible_basis_size,range,parameters->size) + get_rank_chunk_size(possible_basis_size,range,parameters->size)) < new_chunks[parameters->rank+1] ) {
					local_reps += range_contains[range];
					if ( range != (parameters->size - 1 ) ) {
						idx = (uPetscInt)get_rank_chunk_start(possible_basis_size,range+1,parameters->size);
					} else {
						idx = (uPetscInt)new_chunks[parameters->rank+1];
					}
				} else if ( range_contains[range] == 0 ) {
					if ( range != (parameters->size -1 ) ) {
						idx = (uPetscInt)get_rank_chunk_start(possible_basis_size,range+1,parameters->size);
					} else {
						idx = (uPetscInt)new_chunks[parameters->rank+1]; 
					}
				} else  {
					mapped_idx = parameters->get_full_basis_index(parameters,idx);
					if ( (*parameters->simple_basis_check)(parameters,mapped_idx) == PETSC_TRUE && (*parameters->representative_check)(parameters,mapped_idx) == PETSC_TRUE ){
						local_reps ++;
					}
					idx++;
				}
			}
			
			PetscReal additionaltime3 = MPI_Wtime();
			
			if ( parameters->verbosity >= 10 ) {
#if defined(PETSC_USE_64BIT_INDICES)
			PetscPrintf(PETSC_COMM_SELF,"Rank %d: Time taken to run through AGAIN all reps in local partition (%lld) %lf.\n", parameters->rank, local_reps,additionaltime3 - additionaltime2);
#else
			PetscPrintf(PETSC_COMM_SELF,"Rank %d: Time taken to run through AGAIN all reps in local partition (%d) %lf.\n", parameters->rank, local_reps,additionaltime3 - additionaltime2);
#endif
			}
		
			
			//gather all counts again
			MPI_Allgather(&local_reps,1,PETSC_MPI_INT,&all_local_reps[1],1,PETSC_MPI_INT,parameters->mat_solver_comm);
			all_local_reps[0] = 0; 
			for ( i = 1 ; i <= parameters->size ; i++ ) {
				all_local_reps[i] += all_local_reps[i-1];
				if ( parameters->verbosity >= 10 ) {
					stringstream ss;
					ss << "Local reps accummulation on count on " << i << ": " << all_local_reps[i] << " in range " << new_chunks[i-1] << " to " << new_chunks[i] << endl;
					PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
					//PetscPrintf(PETSC_COMM_WORLD,"Local reps accummulation on count on %d: %d in range %d to %d\n",i,all_local_reps[i],new_chunks[i-1],new_chunks[i]);
				}
			}
			
			//allocate space to store chunk on each process
			//PetscPrintf(PETSC_COMM_SELF,"Rank %d, local reps %d and total %d and chunk is %d.\n", parameters->rank, local_reps,n, get_rank_chunk_size(n,parameters->rank,parameters->size));	
			ierr = PetscMalloc(get_rank_chunk_size(n,parameters->rank,parameters->size)*sizeof(uPetscInt),local_reps_array);CHKERRQ(ierr);
			//ierr = PetscMalloc(get_rank_chunk_size(n,parameters->rank,parameters->size)*sizeof(PetscReal),&local_reps_normals);CHKERRQ(ierr);
			ierr = PetscMalloc(parameters->size*sizeof(uPetscInt),ending_idx_for_proc);CHKERRQ(ierr);
			
			//If our range is the very first or spread between two ranges we can start going forward first otherwise we have to wait for the first element.
			//this should get us to the correct range.
			i=0;
			while( i < parameters->size && get_rank_chunk_start(n,parameters->rank,parameters->size) >= all_local_reps[i+1]  ) {
				i++;
			}
			if ( parameters->verbosity >= 10 ) {
				stringstream ss;
				ss << "Rank " << parameters->rank << ": Start is in range " << i << "( start " << get_rank_chunk_start(n,parameters->rank,parameters->size) << ")"<< endl;
				PetscPrintf(PETSC_COMM_SELF,ss.str().c_str());
				//PetscPrintf(PETSC_COMM_SELF,"Rank %d: Start is in range %d (start %d)\n",parameters->rank,i,get_rank_chunk_start(n,parameters->rank,parameters->size));
			}
			
			PetscInt k; 
			//now we have to test if it is anchored to either side so we can start right away or if we have to wait.
			if ( get_rank_chunk_start(n,parameters->rank,parameters->size) > all_local_reps[i] &&  (get_rank_chunk_start(n,parameters->rank,parameters->size)+get_rank_chunk_size(n,parameters->rank,parameters->size)) < all_local_reps[i+1] ) {
				//if its not anchored to either side we need to wait to receive where to start from the previous process with a blocking receive.
				if ( parameters->verbosity >= 10 ) {
					stringstream ss;
					ss << "Rank " << parameters->rank << ": got to wait to receive start as " << get_rank_chunk_start(n,parameters->rank,parameters->size) << " is not anchored to " << all_local_reps[i] << " or " << all_local_reps[i+1] << endl ; 
					//PetscPrintf(PETSC_COMM_SELF,"Rank %d: got to wait to receive start as %d is not anchored to %d or %d\n",parameters->rank,get_rank_chunk_start(n,parameters->rank,parameters->size),all_local_reps[i],all_local_reps[i+1]);
					PetscPrintf(PETSC_COMM_SELF,ss.str().c_str());
				}
		//		MPI_Recv(&j,1,PETSC_MPI_INT,parameters->rank-1,0,PETSC_COMM_WORLD,&stat); //receive starting element into j
				k = 0 ;
				int found = 0;
				int tmp_rank = parameters->rank - 1 ; //go to previous rank
				while ( !found &&  tmp_rank >= 0 ) { 
					i=0;
					while( get_rank_chunk_start(n,tmp_rank,parameters->size) >= all_local_reps[i+1]  ) { //find range for this rank
						i++;
					}
					//make sure this rank is not anchored
					if ( get_rank_chunk_start(n,tmp_rank,parameters->size) == all_local_reps[i] ) { // of the start of this range coincides with the start of the reduced basis range
						found = 1;
						j = new_chunks[i];
						k = all_local_reps[i] - get_rank_chunk_start(n,parameters->rank,parameters->size);
					//else if it overlaps the border
					} else if ( get_rank_chunk_start(n,tmp_rank,parameters->size) < all_local_reps[i+1] &&  (get_rank_chunk_start(n,tmp_rank,parameters->size)+get_rank_chunk_size(n,tmp_rank,parameters->size)) > all_local_reps[i+1]) {	
						found = 1;
						j = new_chunks[i+1];
						k = all_local_reps[i+1] - get_rank_chunk_start(n,parameters->rank,parameters->size);
					} else {	//other wise try another range.
						tmp_rank--;
					}
				} //end while
		
			} else if ( get_rank_chunk_start(n,parameters->rank,parameters->size) == all_local_reps[i] ) {  //if the start coincides with the start of the range.
				//PetscPrintf(PETSC_COMM_SELF,"Rank %d matches at the beginning\n",parameters->rank);
				j = new_chunks[i];
				k = 0;
			} else if ((get_rank_chunk_start(n,parameters->rank,parameters->size)+get_rank_chunk_size(n,parameters->rank,parameters->size)) == all_local_reps[i+1] ) { //if the end corresponds with the end. 
				//PetscPrintf(PETSC_COMM_SELF,"Rank %d matches at the end\n",parameters->rank);
				j = new_chunks[i+1];
				k = get_rank_chunk_size(n,parameters->rank,parameters->size);
			} else {  //other wise it must run over the end and so can start forward at this point and then work back to get the rest.
				//PetscPrintf(PETSC_COMM_SELF,"Rank %d overlaps the start\n",parameters->rank);
				j = new_chunks[i+1];
				k = all_local_reps[i+1] - get_rank_chunk_start(n,parameters->rank,parameters->size);
			}
			
			PetscInt starter = j;
			//this is working forwards 
			while(k < get_rank_chunk_size(n,parameters->rank,parameters->size) && (uPetscInt)j < possible_basis_size ) { //condition on j not needed but to be on the safe side
				int range = get_range_with_item(possible_basis_size,j,parameters->size); //find out which range from the original count we are in and how many basis elements are in it.
				if ( range_contains[range] == 0 ) {  //if there are non we can skip this range.
					j = get_rank_chunk_start(possible_basis_size,range + 1,parameters->size);
				} 
				mapped_idx = parameters->get_full_basis_index(parameters,j); //we get the full index from the reduced index j.
				//we check that the full basis obeys the basis conditions and that it is a representative basis element (only used in momentum calculations).
				if ( (*parameters->simple_basis_check)(parameters,mapped_idx) == PETSC_TRUE && (*parameters->representative_check)(parameters,mapped_idx) == PETSC_TRUE ){
					if (k >= 0) {  //if k is positive and the checks work out we add the entry to the array.
						(*local_reps_array)[k] = (uPetscInt)j ;
					}
					k++;
				}
				j++; //check the next index.
			}
			
			
			PetscReal additionaltime5 = MPI_Wtime();
		
			
			//now work backwards and fill the rest of the list if its necessary.
			//we can work back from the last index which we know from the array 
			if ( get_rank_chunk_start(n,parameters->rank,parameters->size) < all_local_reps[i+1] && get_rank_chunk_start(n,parameters->rank,parameters->size) != all_local_reps[i] &&  (get_rank_chunk_start(n,parameters->rank,parameters->size)+get_rank_chunk_size(n,parameters->rank,parameters->size)) >= all_local_reps[i+1] ) {
				k = all_local_reps[i+1] - get_rank_chunk_start(n,parameters->rank,parameters->size) - 1;
				j = new_chunks[i+1]-1;
				while ( k >= 0 && j >= 0) { //again condition  on the j should not be necessary 
					int range = get_range_with_item(possible_basis_size,j,parameters->size);
					//PetscPrintf(PETSC_COMM_SELF,"%d is in range %d from %d to %d.\n",j,range,get_rank_chunk_start(matrix_data->numberBasisStates,range,parameters->size),get_rank_chunk_start(matrix_data->numberBasisStates,range,parameters->size)+get_rank_chunk_size(matrix_data->numberBasisStates,range,parameters->size));
					if ( range_contains[range] == 0 ) {
						j = get_rank_chunk_start(possible_basis_size,range,parameters->size) - 1;
					} 
					mapped_idx = parameters->get_full_basis_index(parameters,j);
					if ( (*parameters->simple_basis_check)(parameters,mapped_idx) == PETSC_TRUE && (*parameters->representative_check)(parameters,mapped_idx) == PETSC_TRUE ){
						(*local_reps_array)[k] = (uPetscInt)j ; 
						k--;
					}
					j--;
				}
			}
			
			PetscReal additionaltime6 = MPI_Wtime();
			
		
			
			PetscFree(new_chunks);
			PetscFree(range_contains);
			
			PetscReal additionaltime4 = MPI_Wtime();
			if ( parameters->verbosity >= 10 ) {
				stringstream ss;
				ss << "Rank " << parameters->rank << ": Time taken to fill the local arrays after initial two counts: " << additionaltime4 - additionaltime3 << endl ;
				PetscPrintf(PETSC_COMM_SELF,ss.str().c_str());
			}
			
			//now we need to find out up to which index each process is storing. Do this with an all to all gather.
			if ( get_rank_chunk_size(n,parameters->rank,parameters->size) == 0 ) {
				uPetscInt tmp =0;
				MPI_Allgather(&tmp,1,PETSC_UNSIGNED_MPI_INT,*ending_idx_for_proc,1,PETSC_UNSIGNED_MPI_INT,parameters->mat_solver_comm);
			} else {
				MPI_Allgather(&(*local_reps_array)[get_rank_chunk_size(n,parameters->rank,parameters->size)-1],1,PETSC_UNSIGNED_MPI_INT,*ending_idx_for_proc,1,PETSC_UNSIGNED_MPI_INT,parameters->mat_solver_comm);
			}
			/*for ( i = 0 ; i < parameters->size ; i ++ ){ 
				PetscPrintf(PETSC_COMM_WORLD,"Chunk %d goes up to %d(%u).\n",i,ending_idx_for_proc[i],local_reps_array[get_rank_chunk_size(n,parameters->rank,parameters->size)-1]);
			} //printing list and looks good*/
		
			PetscReal endtime1 = MPI_Wtime();
			if ( parameters->verbosity >= 5 ) {
				PetscPrintf(PETSC_COMM_WORLD,"Time taken for filling local reps array: %lf\n",endtime1-starttime1);
			}
		}  // end if n > 0
		PetscFree(all_local_reps);
		return n;
	} else {
		return possible_basis_size;
	}
	return 0 ;
}

#undef __FUNCT__
#define __FUNCT__ "exchanged_required_basis_indices"
/*! \brief This function exchanges the required basis elements between processes. */
PetscInt exchange_required_basis_indices(struct parameter_struct *parameters,struct matrix_meta_data *matrix_data, uPetscInt *local_reps_array,uPetscInt *ending_idx_for_proc,PetscInt **d_nnz, PetscInt **o_nnz,uPetscInt **reps, uPetscInt **indices,uPetscInt *num_unique_cols){
	struct matrix_row_values row_vals;
	PetscErrorCode ierr;
	PetscInt n = parameters->current_basis_size;
	
	
	PetscReal endtime1 = MPI_Wtime();
	allocate_row_space(&row_vals, matrix_data); 
	//get start and finish index of the ownership range. Note Iend is the end index + 1 
	PetscInt Istart = get_rank_chunk_start(n, parameters->rank, parameters->size);
	PetscInt Iend = get_rank_chunk_start(n, parameters->rank, parameters->size) + get_rank_chunk_size(n, parameters->rank, parameters->size);
	//allocate space for arrays containing information about non zero entries.
	ierr = PetscMalloc((Iend - Istart) * sizeof(PetscInt),d_nnz);CHKERRQ(ierr);
	ierr = PetscMalloc((Iend - Istart) * sizeof(PetscInt),o_nnz);CHKERRQ(ierr);
	
	#ifdef _BOOST_BST_
	BST **unique_cols_per_process_bst;
	#endif
	list<uPetscInt> **unique_cols_per_process;
	
	*num_unique_cols = 0;
	if ( parameters->use_momentum_sectors == PETSC_TRUE || parameters->nn_exclusion == PETSC_TRUE || parameters->use_rotation_invariance == PETSC_TRUE  ) {
		//for ( int proc = 0 ; proc < parameters->restricted_size ; proc++ ) {
			#ifdef _BOOST_BST_
			if ( parameters->use_bst == PETSC_TRUE ) {
				//BST tmp;
				//unique_cols_per_process_bst.insert(unique_cols_per_process_bst.begin(),tmp);
				ierr = PetscMalloc(sizeof(BST*)*parameters->size,&unique_cols_per_process_bst);CHKERRQ(ierr);
				//unique_cols_per_process_bst = new BST[parameters->size];
				for ( int i = 0 ; i < parameters->size ; i++ ) {
					unique_cols_per_process_bst[i] = new BST();
					//ierr = PetscMalloc(sizeof(BST),&unique_cols_per_process_bst[i]);CHKERRQ(ierr);
				}
			} else 
			#endif
			{
				//list<uPetscInt>  tmp;
				//unique_cols_per_process.insert(unique_cols_per_process.begin(),tmp);
				ierr = PetscMalloc(sizeof(list<uPetscInt>*)*parameters->size,&unique_cols_per_process);CHKERRQ(ierr);
				for ( int i = 0 ; i < parameters->restricted_size ; i++ ) {
					unique_cols_per_process[i] = new list<uPetscInt>();
					//ierr = PetscMalloc(sizeof(list<uPetscInt>),&unique_cols_per_process[i]);CHKERRQ(ierr)
				}
			}
			
		//}
		//unique_cols_per_process.reserve(parameters->size); //array of lists, one for each process to store unique column indices.	
	}
	//ierr = PetscMalloc(parameters->size * sizeof(list<uPetscInt>), &unique_cols_per_process); //allocate space for lists. TODO deallocate this array at appropriate time.
	
	for ( PetscInt j = 0 ; j < Iend - Istart ; j++ ) {
		(*d_nnz)[j] = 1; //diagonal always has to have at least one (petsc manual).
		(*o_nnz)[j] = 0; //set off diagonal to zero.
	}
	PetscTruth exists;
	PetscTruth exists2;
	PetscScalar phase;
	uPetscInt tmp_col; 
	uPetscInt group_rep;
	for ( PetscInt i = 0 ; i < Iend - Istart ; i++ ) {
		if ( parameters->use_momentum_sectors == PETSC_TRUE || parameters->nn_exclusion == PETSC_TRUE || parameters->use_rotation_invariance == PETSC_TRUE  ) { 
			row_vals.row = (*parameters->get_full_basis_index)(parameters,local_reps_array[i]); //set which row to retrieve
		} else row_vals.row = (*parameters->get_full_basis_index)(parameters,(uPetscInt)(Istart + i)); 
		(*parameters->getrow)(&row_vals, matrix_data, parameters,1); //get the values from that row	
		//list<PetscInt> col_list; //create a list that will store all non zero column indices up to this 
		//list<PetscInt>::iterator col_it; //iterator to iterate over this list
		list<uPetscInt> counted_cols; //list to store the ammount of columns that have been counted for this row. Can get duplicates when two cols have the same representative.
		for (int j = 0 ; j < row_vals.val_count ; j++ ) {
			group_rep = (*parameters->group_representative)(parameters,row_vals.cols[j],PETSC_TRUE,&phase,&exists);  //get the representative for the given basis element
			tmp_col = (*parameters->get_reduced_basis_index)(parameters,group_rep,&exists2);
			//PetscPrintf(PETSC_COMM_WORLD,"From %d we get %d (rep %d).\n",row_vals.cols[j],tmp_col,get_representative(parameters,row_vals.cols[j],PETSC_TRUE,&phase,mom_sector));
			if ( exists == PETSC_TRUE && exists2 == PETSC_TRUE && (*parameters->simple_basis_check)(parameters,group_rep) == PETSC_TRUE) {
				//find appropriate list to insert column index in
				if ( parameters->use_momentum_sectors == PETSC_TRUE || parameters->nn_exclusion == PETSC_TRUE || parameters->use_rotation_invariance == PETSC_TRUE  ) { //if basis requires lists of representatives.
					for ( int proc = 0; proc < parameters->restricted_size ; proc++ ) {
						if ( (proc == 0 && tmp_col <= (uPetscInt)ending_idx_for_proc[0] ) || ( tmp_col <= (uPetscInt)ending_idx_for_proc[proc] && tmp_col > (uPetscInt)ending_idx_for_proc[proc-1] )) { //if the representation is to be found on this processor then add it to the list
							if ( parameters->use_hermitian == PETSC_FALSE || tmp_col >= local_reps_array[i] ) { //only include elements from upper triangular if use_hermitian is true.
								if ( insert_into_sorted_list(&counted_cols,tmp_col) == 1 ) { //this is here to prevent non zero elements being counted twice for the preallocation for a particular row. 
									if ( tmp_col < (uPetscInt)local_reps_array[0] || tmp_col > (uPetscInt)local_reps_array[Iend-Istart-1] ) { //if it is in the off diagonal. On one core there is no off diagonal block 
										(*o_nnz)[i]++;
									} else { 
										if ( tmp_col != (uPetscInt)local_reps_array[i] ) { //if the element is not on the diagonal.
											(*d_nnz)[i]++ ;
										}
									}
								}	
								#ifdef _BOOST_BST_
								if ( parameters->use_bst == PETSC_TRUE ) {
									CuPetscInt *tmp;
									//ierr = PetscMalloc(sizeof(CuPetscInt),&tmp);CHKERRQ(ierr);
									tmp = new CuPetscInt();
									tmp->int_ = tmp_col;
									if ( unique_cols_per_process_bst[proc]->insert(*tmp).second == false ) {
										//ierr = PetscFree(tmp);CHKERRQ(ierr);
										delete tmp;
									} else {
										(*num_unique_cols)++;
									}
								} else 
								#endif
								{								
									if ( insert_into_sorted_list(unique_cols_per_process[proc],tmp_col) == 1 ) { 
										(*num_unique_cols)++;//this is to count the complete number so we can allocate the arrays and not have to count later.		
									}
								}  
							} // end if to check if we are in the upper half or not using the hermitian type.
						} //end if the column is within the range of the current process 
					}//end for over range of processes.
				} else { //if no list of indices are requred we just want to count numbers of non zeros for preallocation purposes. 
					if ( parameters->use_hermitian == PETSC_FALSE || tmp_col >= row_vals.row ) { //if we are using hermitian property or its in upper triangle then we count it.
						if ( tmp_col < (uPetscInt)Istart || tmp_col >= (uPetscInt)Iend ) { //if its in the off diagonal block
							(*o_nnz)[i]++;
						} else if ( tmp_col != (uPetscInt)(i + Istart) ) { //if its in the diagonal block but not on the diagonal as that is already counted.
							(*d_nnz)[i]++;
						}
					}
				}
			} //if the representative is included in this sector
		}	
		//PetscPrintf(PETSC_COMM_SELF, "Row %d, d: %d, o: %d, tot: %d\n", Istart + i, (*d_nnz)[i], (*o_nnz)[i], (*d_nnz)[i] + (*o_nnz)[i]);
	}
	
	//want to construct three sorted arrays one containing the representations, one containing the indices and one containing the corresponding normals.
	//to find the size of this array we sum the totals of the individual lists of processors.
	
	PetscReal endtime2 = MPI_Wtime();
	if ( parameters->verbosity >= 5) {
		stringstream ss;
		ss << "Number unique reps: "<< *num_unique_cols << endl ;
		PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
		ss.str(""); ss << "Time taken to fill local per process lists: "<< endtime2 - endtime1 << endl;
		PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
	}


	//if we do not need to index the basis then we are done as have counted the non zeros otherwise we exchange the needed indices.
	if ( parameters->use_momentum_sectors == PETSC_TRUE || parameters->nn_exclusion == PETSC_TRUE || parameters->use_rotation_invariance == PETSC_TRUE ) { 
		//declare and allocate space in arrays.
		
		//PetscReal *normals;
		ierr = PetscMalloc(*num_unique_cols * sizeof(uPetscInt),reps);CHKERRQ(ierr);
		ierr = PetscMalloc(*num_unique_cols * sizeof(uPetscInt),indices);CHKERRQ(ierr);
		//ierr = PetscMalloc(num_unique_cols * sizeof(PetscReal),&normals);CHKERRQ(ierr);
	
		uPetscInt *unique_cols_proc_offs; //array to store how many reps' indices are needed from each process
		ierr = PetscMalloc(parameters->size * sizeof(uPetscInt), &unique_cols_proc_offs);CHKERRQ(ierr);
		{
		int reps_idx =0;
			for ( int i = 0 ; i < parameters->restricted_size; i++ ) {
				if ( i == 0 ) unique_cols_proc_offs[i] = 0; else unique_cols_proc_offs[i] = unique_cols_proc_offs[i-1];
				#ifdef _BOOST_BST_
				if ( parameters->use_bst == PETSC_TRUE ) {
					CuPetscInt *instance;
					BST::iterator it(unique_cols_per_process_bst[i]->begin());
					while ( it != unique_cols_per_process_bst[i]->end() ) {
						instance = &(*it);
						(*reps)[reps_idx++] = instance->int_;
						it = unique_cols_per_process_bst[i]->erase(it);
						//ierr = PetscFree(instance);CHKERRQ(ierr);
						delete instance;
						unique_cols_proc_offs[i]++;
					}
				} else 
				#endif 
				{
					while( unique_cols_per_process[i]->empty() == false) {
						(*reps)[reps_idx++] = unique_cols_per_process[i]->front();
						unique_cols_per_process[i]->pop_front();
						unique_cols_proc_offs[i]++;
					}
				}
			}
		}
		
		#ifdef _BOOST_BST_
		if ( parameters->use_bst == PETSC_TRUE ) {
			//ierr = PetscFree(unique_cols_per_process_bst);CHKERRQ(ierr);
			//delete unique_cols_per_process_bst;
			for ( int i = 0 ; i < parameters->size; i++ ) {
				//PetscFree(unique_cols_per_process_bst[i]);
				delete unique_cols_per_process_bst[i];
			}
			PetscFree(unique_cols_per_process_bst);
		} else 
		#endif
		{
			//ierr = PetscFree(unique_cols_per_process);CHKERRQ(ierr);
			//delete unique_cols_per_process;
			for ( int i = 0 ; i < parameters->restricted_size; i++ ) {
				//PetscFree(unique_cols_per_process[i]);
				delete unique_cols_per_process[i];
			}
			PetscFree(unique_cols_per_process);
		}
	
		//we now have an array called reps of size num_unique_cols with all the reps and count in each stored in array unique_cols_per_process_sizes
		//now we start getting the indices and normals from the other processes starting with the current one.
		PetscReal endtime3 = MPI_Wtime();
		if ( parameters->verbosity >= 5 ) {
			PetscPrintf(PETSC_COMM_WORLD,"Time taken  to fill local unique col reps array from per process lists: %lf\n",endtime3 - endtime2);
		}
		
		//get indices and normals from current process
		for ( PetscInt i = (parameters->rank==0) ? 0 : unique_cols_proc_offs[parameters->rank-1] ; i <  (PetscInt)unique_cols_proc_offs[parameters->rank] ; i++ ) {
			(*indices)[i] = find_in_sorted_array(local_reps_array,get_rank_chunk_size(n,parameters->rank,parameters->size),(*reps)[i]);
			//normals[i] = local_reps_normals[indices[i]];
			if ( ((PetscInt)(*indices)[i]) >= 0 ) {
				(*indices)[i] += get_rank_chunk_start(n,parameters->rank,parameters->size) ;
			}
		} 
	
		PetscReal endtime4 = MPI_Wtime();
		if ( parameters->verbosity >= 5 ) {
			PetscPrintf(PETSC_COMM_WORLD,"Time taken  to fill local cols and normals from local process data: %lf\n",endtime4 - endtime3);
		}
		
		MPI_Request req;
		MPI_Status stat;
	
		//get indices and normals from the other processes.	 
		for ( int proc_off = 1 ; proc_off < parameters->restricted_size ; proc_off++ ) { //the proc_off is the offset from current rank we are getting indices and normals from.
			int source_proc = (parameters->rank + proc_off) % parameters->restricted_size;
			//int target_proc = (((parameters->rank + proc_off) % parameters->size) == 0 ) ? parameters->size-1 : ((parameters->rank + proc_off) % parameters->size);
			int target_proc = (parameters->restricted_size + parameters->rank - proc_off) % parameters->restricted_size;
			
			//PetscPrintf(PETSC_COMM_SELF,"Rank %d, Source %d and destination %d.\n", parameters->rank,source_proc,target_proc);
			
			//we hope to receive the indices and normals from the source_proc but must send the reps first
			//we in turn must receive the reps from the target_proc and set the associated indices and normals back. 
	
			//non blocking send reps to source_proc
			
			uPetscInt send_count = (source_proc==0)?unique_cols_proc_offs[0]:unique_cols_proc_offs[source_proc]-unique_cols_proc_offs[source_proc-1];
			uPetscInt recv_count;
			MPI_Isend(&send_count,1,PETSC_UNSIGNED_MPI_INT,source_proc,0,parameters->mat_solver_comm,&req); //send the size of the array we are sending
			MPI_Recv(&recv_count,1,PETSC_UNSIGNED_MPI_INT,target_proc,0,parameters->mat_solver_comm,&stat); //receive the size of the array we are receiving
			MPI_Wait(&req,&stat); //wait for first send to complete
			//send reps
			MPI_Isend(&(*reps)[(source_proc==0)?0:unique_cols_proc_offs[source_proc-1]],send_count,PETSC_UNSIGNED_MPI_INT,source_proc,0,parameters->mat_solver_comm,&req);
			uPetscInt *target_reps,*target_indices; //PetscReal *target_normals; 
			ierr = PetscMalloc(recv_count*sizeof(uPetscInt),&target_reps);CHKERRQ(ierr);
			ierr = PetscMalloc(recv_count*sizeof(uPetscInt),&target_indices);CHKERRQ(ierr);
			//ierr = PetscMalloc(recv_count*sizeof(PetscReal),&target_normals);CHKERRQ(ierr);
			MPI_Recv(target_reps,recv_count,PETSC_UNSIGNED_MPI_INT,target_proc,0,parameters->mat_solver_comm,&stat);//receive reps
	
			for (PetscInt i = 0 ; i < (PetscInt)recv_count ; i++ ) {
				target_indices[i] = find_in_sorted_array(local_reps_array,get_rank_chunk_size(n,parameters->rank,parameters->size),target_reps[i]);
				//target_normals[i] = local_reps_normals[target_indices[i]];
				if ( ((PetscInt)target_indices[i]) >= 0 ) { //ensure that it is actually found and not -1
					target_indices[i] += get_rank_chunk_start(n,parameters->rank,parameters->size); 
				}
			}
			
			//send back the indices 
			MPI_Request req2;
			MPI_Isend(target_indices, recv_count, PETSC_UNSIGNED_MPI_INT, target_proc,0, parameters->mat_solver_comm,&req2); 	
			//send back the normals
			//MPI_Request req3;
			//MPI_Isend(target_normals, recv_count, PETSC_MPI_REAL, target_proc,1, PETSC_COMM_WORLD,&req3);
			
			//receive the indices 
			MPI_Recv(&(*indices)[(source_proc==0)?0:unique_cols_proc_offs[source_proc-1]],send_count,PETSC_UNSIGNED_MPI_INT,source_proc,0,parameters->mat_solver_comm,&stat);
			//MPI_Recv(&normals[(source_proc==0)?0:unique_cols_proc_offs[source_proc-1]],send_count,PETSC_MPI_REAL,source_proc,1	,PETSC_COMM_WORLD,&stat);
			
			MPI_Wait(&req,&stat);
			MPI_Wait(&req2,&stat);
			//MPI_Wait(&req3,&stat);
			
			PetscFree(target_reps); //free the array used to store the reps received.
			PetscFree(target_indices); 
			
			//PetscFree(target_normals); 
		
		} //end loop over processors to do the staggered communication. Not performed when using only one process.
	
		PetscReal endtime5 = MPI_Wtime();
		
		if ( parameters->verbosity >= 5 ) {
			PetscPrintf(PETSC_COMM_WORLD,"Time taken  to fill local cols and normals from other process data: %lf\n",endtime5 - endtime4);
		}
		ierr = PetscFree(unique_cols_proc_offs);CHKERRQ(ierr);
	} //end if use_momentum or nn_exclusion. 
	
	deallocate_row_space(&row_vals);
	return 0; 
}

#undef __FUNCT__
#define __FUNCT__ "set_matrix_elements"
/*! \brief This function sets all the basis elements for the local portion of the matrix on each processor. */
PetscInt set_matrix_elements(Mat *A,struct parameter_struct *parameters,struct matrix_meta_data *matrix_data, uPetscInt *local_reps_array,uPetscInt *ending_idx_for_proc,uPetscInt *reps, uPetscInt *indices,uPetscInt num_unique_cols){	
	//for each procssor 
	PetscInt Istart = get_rank_chunk_start(parameters->current_basis_size, parameters->rank, parameters->size);
	PetscInt Iend = get_rank_chunk_start(parameters->current_basis_size, parameters->rank, parameters->size) + get_rank_chunk_size(parameters->current_basis_size, parameters->rank, parameters->size);
	struct matrix_row_values row_vals;
	uPetscInt tmp_col;
	
	
	PetscReal row_normal, col_normal;
	PetscScalar phase;
	PetscTruth exists, exists2;
	PetscErrorCode ierr;
	allocate_row_space(&row_vals, matrix_data); 
	matrix_data->is_diagonal = PETSC_TRUE;
	
	PetscReal endtime5 = MPI_Wtime();
	for ( PetscInt i = 0 ; i < Iend - Istart ; i++ ) {
		PetscInt row_idx = i + Istart;
		if ( parameters->use_momentum_sectors == PETSC_TRUE || parameters->nn_exclusion == PETSC_TRUE || parameters->use_rotation_invariance == PETSC_TRUE  ) { 
			row_vals.row = (*parameters->get_full_basis_index)(parameters,local_reps_array[i]); //set which row to retrieve
			(*parameters->getrow)(&row_vals, matrix_data, parameters,1); //get the values from that row
			row_normal = (*parameters->group_normal)(parameters,row_vals.row);
			for ( int j = 0 ; j < row_vals.val_count ; j++ ) { 
				
				//group_rep = (*parameters->group_representative)(parameters,row_vals.cols[j],PETSC_TRUE,&phase,&exists);
				//tmp_col = (*parameters->get_reduced_basis_index)(parameters,group_rep,&exists2); //get the representative for the given basis element
				tmp_col = (*parameters->get_reduced_basis_index)(parameters,(*parameters->group_representative)(parameters,row_vals.cols[j],PETSC_TRUE,&phase,&exists),&exists2); 
				PetscInt idx ;
				if ( (parameters->rank==0 || tmp_col > ending_idx_for_proc[parameters->rank-1] ) && tmp_col <= ending_idx_for_proc[parameters->rank] ) {
					idx = find_in_sorted_array(local_reps_array,get_rank_chunk_size(parameters->current_basis_size, parameters->rank, parameters->restricted_size),tmp_col);
				} else idx = find_in_sorted_array(reps,num_unique_cols,tmp_col);
				//PetscPrintf(PETSC_COMM_WORLD,"Row %d(%d), col %d(%d), phase %lf + i %lf.\n",row_vals.row,i,row_vals.cols[j],indices[idx],PetscRealPart(phase), PetscImaginaryPart(phase));
				if ( idx == -1 ) { //if its not found just add a 0 
					//PetscPrintf(PETSC_COMM_WORLD,"Was not found.\n");
					row_vals.cols[j] = row_idx ;
					#if defined(PETSC_USE_COMPLEX)
					row_vals.values[j] = 0;
					#else		
					row_vals.real_vals[j] = 0;
					#endif 
				} else { //else if the column representative is valid
//					if ( parameters->use_hermitian == PETSC_FALSE || indices[idx]  >= (uPetscInt)row_idx ) {
						col_normal = (*parameters->group_normal)(parameters,(uPetscInt)row_vals.cols[j]);
						//PetscPrintf(PETSC_COMM_WORLD,"Row %d(%d), col %d(%d), row_normal %lf, col_normal %lf, phase %lf + i %lf, value %lf + i%lf.\n",row_vals.row,i,row_vals.cols[j],indices[idx],row_normal,col_normal, PetscRealPart(phase), PetscImaginaryPart(phase),PetscRealPart(row_vals.values[j]), PetscImaginaryPart(row_vals.values[j]));
						//set the non zero values in  the row 
						#if defined(PETSC_USE_COMPLEX)
						row_vals.values[j] *= phase * col_normal / row_normal;
						#else		
						row_vals.real_vals[j] *= phase * col_normal / row_normal;
						#endif 
						if ( (parameters->rank==0 || tmp_col > ending_idx_for_proc[parameters->rank-1] ) && tmp_col <= ending_idx_for_proc[parameters->rank] ) {
							row_vals.cols[j] = get_rank_chunk_start(parameters->current_basis_size, parameters->rank, parameters->restricted_size) +idx;
						} else {
							row_vals.cols[j] = (PetscInt)indices[idx];
						}
/*					} else {
						//PetscPrintf(PETSC_COMM_WORLD,"Not hermitian or too small.\n");
						row_vals.cols[j] = row_idx ;
						#if defined(PETSC_USE_COMPLEX)
						row_vals.values[j] = 0;
						#else		
						row_vals.real_vals[j] = 0;
						#endif 
					}*/
				}
				//checking if there are any off diagonal elements.
				#if defined(PETSC_USE_COMPLEX)
				if ( row_vals.cols[j] != row_idx && PetscAbsScalar(row_vals.values[j]) > parameters->tol ) {
				#else
				if ( row_vals.cols[j] != row_idx && ( (PetscAbsReal(row_vals.real_vals[j]) > parameters->tol) || (PetscAbsReal(row_vals.imag_vals[j]) > parameters->tol)) ) {
				#endif
					matrix_data->is_diagonal = PETSC_FALSE;
				}
			}
		} else { // else if not using momentum sectors or nn exculsion
			row_vals.row = (*parameters->get_full_basis_index)(parameters,(uPetscInt)(Istart + i)); 
			(*parameters->getrow)(&row_vals, matrix_data, parameters,1); //get the values from that row
			for ( int j = 0 ; j < row_vals.val_count ; j++ ) { 
				tmp_col = (*parameters->get_reduced_basis_index)(parameters,row_vals.cols[j],&exists2); //get the representative for the given basis element
				if ( (exists2 == PETSC_TRUE) && ( (parameters->use_hermitian == PETSC_FALSE ) || ( tmp_col >= (uPetscInt)row_idx ))) {
					row_vals.cols[j] = (PetscInt)tmp_col;
					matrix_data->is_diagonal = PETSC_FALSE;
				} else {
					row_vals.cols[j] = row_idx ;
					#if defined(PETSC_USE_COMPLEX)
					row_vals.values[j] = 0;
					#else		
					row_vals.real_vals[j] = 0;
					#endif
				}
			}
		}
		
		if (row_vals.val_count > 0 ) {
			//set the non zero values in  the row 
#if defined(PETSC_USE_COMPLEX)
			ierr = MatSetValues(*A, 1, &row_idx, row_vals.val_count, row_vals.cols, row_vals.values, ADD_VALUES); CHKERRQ(ierr);
#else		
			ierr = MatSetValues(*A, 1, &row_idx, row_vals.val_count, row_vals.cols, row_vals.real_vals, ADD_VALUES); CHKERRQ(ierr);
#endif 
		}
	} //end loop over rows 

	ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	
	if ( parameters->use_hermitian == PETSC_TRUE ) { 
		ierr = MatSetOption(*A,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr);
		ierr = MatSetOption(*A,MAT_SYMMETRY_ETERNAL,PETSC_TRUE);CHKERRQ(ierr);
	}
	
	
	PetscInt is_diagonal = 1,check;
	if ( matrix_data->is_diagonal == PETSC_FALSE ) { //if it is not diagonal on this process set to zero so that when product is taken total will be zero.
			is_diagonal = 0; 
	}
	MPI_Allreduce(&is_diagonal,&check, 1, PETSC_MPI_INT, MPI_PROD, parameters->mat_solver_comm);
	
	if (check == 0 ) {
			matrix_data->is_diagonal = PETSC_FALSE;
	}
	
	PetscReal endtime6 = MPI_Wtime();
	if ( parameters->verbosity >= 5 ) {
		PetscFPrintf(PETSC_COMM_WORLD,stderr,"Time taken to fill local matrix elements: %lf\n",endtime6 - endtime5);
		if ( matrix_data->is_diagonal) {
			PetscFPrintf(PETSC_COMM_WORLD,stderr,"Matrix is diagonal.\n");		
		} else {
			PetscFPrintf(PETSC_COMM_WORLD,stderr,"Matrix is not diagonal.\n");		
		}
	}

	deallocate_row_space(&row_vals);
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "dummy_representative_check"
/*! \brief This function is a dummy function which is used when translational invariance is not being used. */
PetscTruth dummy_representative_check(struct parameter_struct *parameters, uPetscInt index) {
	return PETSC_TRUE;
}

#undef __FUNCT__
#define __FUNCT__ "dummy_group_representative"
/*! \brief This function is a dummy function which is used when translational invariance is not being used. */
uPetscInt  dummy_group_representative(struct parameter_struct *parameters, uPetscInt basis, PetscTruth phase_check, PetscScalar *phase,PetscTruth *group_exists) {
	*phase = 1.0;
	*group_exists = PETSC_TRUE;
	return basis;
}

#undef __FUNCT__
#define __FUNCT__ "dummy_group_normal"
/*! \brief This function is a dummy function which is used when translational invariance is not being used. */
PetscReal  dummy_group_normal(struct parameter_struct *parameters, uPetscInt index) {
	PetscReal normal = 1.0;
	return normal;
}

#undef __FUNCT__
#define __FUNCT__ "create_matrix_communicator"
/*! \brief Create a new communicator with the subset of processes that own some elements. Only for small basis sizes will the new communicator be different than PETSC_COMM_WORLD. */
PetscInt create_matrix_communicator(struct parameter_struct *parameters) {
	MPI_Group processes_being_used;
	MPI_Group global_group;
	MPI_Comm  comm_world = PETSC_COMM_WORLD;
	int number_relevant_ranks = 0;
	int *relevant_ranks;
	for ( int proc = 0 ; proc < parameters->size ; proc++ ) {
		if ( get_rank_chunk_size(parameters->current_basis_size,proc,parameters->size) > 0 ) {
			number_relevant_ranks++; 
		}
	}
	
	if ( get_rank_chunk_size(parameters->current_basis_size,parameters->rank,parameters->size) > 0 ) {
		parameters->in_communicator = PETSC_TRUE; 
	} else {
		parameters->in_communicator = PETSC_FALSE; 
	}
	
	
	if ( number_relevant_ranks > 0 ) {
		PetscErrorCode ierr = PetscMalloc(number_relevant_ranks*sizeof(int),&relevant_ranks);CHKERRQ(ierr);
		int rank_idx = 0;
		//PetscFPrintf(PETSC_COMM_SELF,stderr,"Including %d tasks in new communicator.\n",number_relevant_ranks);
		for ( int proc = 0 ; proc < parameters->size ; proc++ ) {
			if ( get_rank_chunk_size(parameters->current_basis_size,proc,parameters->size) > 0 ) {
				//PetscFPrintf(PETSC_COMM_SELF,stderr,"Including rank %d as element %d.\n",proc,rank_idx);
				relevant_ranks[rank_idx++] = proc; 
			}
		}
		MPI_Comm_group(comm_world,&global_group);
		MPI_Group_incl(global_group,number_relevant_ranks,relevant_ranks,&processes_being_used);
		MPI_Comm_create(comm_world,processes_being_used,&(parameters->mat_solver_comm));
		MPI_Group_free(&processes_being_used);
		MPI_Group_free(&global_group);
		PetscFree(relevant_ranks);
	}
	parameters->restricted_size = number_relevant_ranks;
	return 0;
}
