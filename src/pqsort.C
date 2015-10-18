#include "pqsort.h"

/*! \brief Function to perform parallel quicksort to sort an array of unsigned PetscInts.*/
int parallel_qsort_uPetscInt(struct parameter_struct *parameters, MPI_Comm *subcomm, uPetscInt *array, PetscInt llength) {
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
						PetscInt *send_buf;
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
					PetscInt *send_buf;
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
					PetscInt *send_buf;
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
				PetscInt *send_buf;
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
				PetscInt *send_buf;
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
					PetscInt *send_buf;
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
					parallel_qsort_uPetscInt(parameters,&left_subcomm,array,pivot_pos_within_chunk);
				}
			}
			if ( number_relevant_ranks_right ) { 
				if ( pivot_pos_within_chunk != lengths[pivot_process] ) {
				// 	cout << "Pivot process going to right partition of size " << number_relevant_ranks_right << " and local array size " << llength-pivot_pos_within_chunk << ": ";
// 					for ( int i = 0 ; i < llength - pivot_pos_within_chunk ; i++ ) cout << array[pivot_pos_within_chunk + i] << " " ;
// 					cout << endl;
					parallel_qsort_uPetscInt(parameters,&right_subcomm,&(array[pivot_pos_within_chunk]),llength-pivot_pos_within_chunk);
				}
			}
		} else if ( rank < pivot_process) {
			if ( number_relevant_ranks_left ) {
				//cout << "Process left of pivot partition of size " << number_relevant_ranks_left << " and local array size " << llength << endl;
				parallel_qsort_uPetscInt(parameters,&left_subcomm,array,llength);
			}
		} else if ( rank > pivot_process ) {
			if ( number_relevant_ranks_right ) { 
				//cout << "Process right of pivot partition of size " << number_relevant_ranks_right << " and local array size " << llength << ": " ;
				//for ( int i = 0 ; i < llength  ; i++ ) cout << array[i] << " " ;
				//cout << endl;
				parallel_qsort_uPetscInt(parameters,&right_subcomm,array,llength);
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
