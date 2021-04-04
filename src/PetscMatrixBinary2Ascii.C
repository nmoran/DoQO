#include "slepceps.h"
#include "slepcsvd.h"
#include <fstream>
#include <iostream>

#define MAX_STRING 200

using namespace std;

typedef PetscBool PetscTruth;

int main( int argc, char **argv ) {
    PetscErrorCode ierr; 
    string in_mat, out_mat;    
    PetscTruth use_complex = PETSC_FALSE, is_ascii = PETSC_FALSE;
    PetscMPIInt rank, size;
    stringstream ss;           
    PetscTruth flg;
    
    /*------------------------------------------------------------------------------------------------
    Initialise some values.
    ------------------------------------------------------------------------------------------------*/
    char help[] = "Program that converts a binary Petsc vector to Ascii.\n";
    
    SlepcInitialize(&argc,&argv,(char*)0,help);	
    
    //get rank of this process and number of processes 
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr); 
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);        
    
    char buf[MAX_STRING];
    ierr = PetscOptionsGetString(NULL,NULL,"-i",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	PetscPrintf(PETSC_COMM_WORLD,"Must specify input matrix on the command line with switch -i.\n");
	return 1;
      }
    PetscPrintf(PETSC_COMM_WORLD,"Input file: %s\n",buf);
    in_mat.assign(buf);
    
    ierr = PetscOptionsGetString(NULL,NULL,"-o",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	PetscPrintf(PETSC_COMM_WORLD,"Must specify output matrix on the command line with switch -o.\n");
	return 1;
      }
    PetscPrintf(PETSC_COMM_WORLD,"Output file: %s\n",buf);
    out_mat.assign(buf);
                
    ierr = PetscOptionsGetString(NULL,NULL,"-ascii",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	is_ascii = PETSC_FALSE;
      } 
    else 
      {
	if ( strcmp(buf,"true") == 0 ) 
	  {
	    is_ascii = PETSC_TRUE;
	  } 
	else 
	  {	
	    is_ascii = PETSC_FALSE;
	  }
      }
      
    use_complex = PETSC_FALSE;
    ierr = PetscOptionsGetString(NULL,NULL,"-complex",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE) 
      {
	if ( strcmp(buf,"false") == 0 ) 
	  {
	    use_complex = PETSC_FALSE;
	  } 
	else if ( strcmp(buf,"true") == 0 ) 
	  {	
	    use_complex = PETSC_TRUE;	    
	    PetscPrintf(PETSC_COMM_WORLD,"Using complex part of wavefunction also.\n");		
	  }
      }
      
    PetscTruth use_compact = PETSC_FALSE;
    ierr = PetscOptionsGetString(NULL,NULL,"-compact",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_TRUE) 
      {
	if ( strcmp(buf,"false") == 0 ) 
	  {
	    use_complex = PETSC_FALSE;
	  } 
	else if ( strcmp(buf,"true") == 0 ) 
	  {	
	    use_compact = PETSC_TRUE;	    
	    PetscPrintf(PETSC_COMM_WORLD,"Using compact output format.\n");		
	  }
      }    
       
    //load in the vector as local vector on each processor.
    PetscViewer viewer;
    PetscViewerCreate(PETSC_COMM_SELF,&viewer);
    Mat in_mat_mat;
    
    if ( is_ascii ) 
      {
	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, in_mat.c_str() ,&viewer);CHKERRQ(ierr); //ascii not working well for large sizes at the moment.
      } 
    else 
      {
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, in_mat.c_str() ,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      }
    //ierr = MatLoad(viewer,MATSEQDENSE,&in_mat_mat);CHKERRQ(ierr);
    ierr = MatLoad(in_mat_mat, viewer);CHKERRQ(ierr);
    PetscViewerDestroy(&viewer);
        
    PetscInt mat_size_n, mat_size_m;
    MatGetSize(in_mat_mat,&mat_size_n,&mat_size_m);
    
    ofstream out_mat_file(out_mat.c_str());
    if ( out_mat_file.is_open() )
      {
	for ( PetscInt row = 0 ; row < mat_size_n ; row ++ ) 
	  {
	    if ( use_compact == PETSC_FALSE ) 	  
	      {
		out_mat_file.precision(14);
	      }
	    else 
	      {
		out_mat_file.precision(2);
	      }
	    const PetscScalar *vals;
	    MatGetRow(in_mat_mat, row, PETSC_NULL, PETSC_NULL, &vals);	    
	    if ( use_complex  == PETSC_FALSE )
	      {
	      for ( int i = 0 ; i < mat_size_m;  i++ ) 
		{
		  out_mat_file << scientific << PetscRealPart(vals[i])  ; 
		  if ( i != ( mat_size_m - 1 ) ) out_mat_file << " ";
		}
	      }
	    else 
	      {
		for ( int i = 0 ; i < mat_size_m;  i++ ) 
		{
		  if ( PetscAbsScalar(vals[i]) > 1e-8 ) 
		    {
		      if ( PetscAbsScalar(PetscImaginaryPart(vals[i])) > 1e-9 )
			{
			  out_mat_file << scientific << "(" << PetscRealPart(vals[i]) << "," << PetscImaginaryPart(vals[i]) << ")" ; 
			}
		      else
			{
			  out_mat_file << scientific << PetscRealPart(vals[i])  ; 
			}
		    }
		  else
		    {
		      out_mat_file << "0" ;
		    }
		  if ( i != ( mat_size_m - 1 ) ) out_mat_file << " ";
		}
	      }	    
	    out_mat_file << endl ;  
	    MatRestoreRow(in_mat_mat, row, &mat_size_m, PETSC_NULL, &vals);	    
	  }
	out_mat_file.close();
      }
    else
      {
	PetscPrintf(PETSC_COMM_WORLD,"Error opening file for writing %s,\n", out_mat.c_str());
	return 1;
      }       
    
    ierr = SlepcFinalize();CHKERRQ(ierr);
    return 0;
}
