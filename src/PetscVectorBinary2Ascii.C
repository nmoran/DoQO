#include "slepceps.h"
#include "slepcsvd.h"
#include <fstream>
#include <iostream>

#define MAX_STRING 200

using namespace std;

typedef PetscBool PetscTruth;

int main( int argc, char **argv ) {
    PetscErrorCode ierr; 
    string in_vec, out_vec;    
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
	PetscPrintf(PETSC_COMM_WORLD,"Must specify input vector on the command line with switch -i.\n");
	return 1;
      }
    PetscPrintf(PETSC_COMM_WORLD,"Input file: %s\n",buf);
    in_vec.assign(buf);
    
    ierr = PetscOptionsGetString(NULL,NULL,"-o",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	PetscPrintf(PETSC_COMM_WORLD,"Must specify output vector on the command line with switch -o.\n");
	return 1;
      }
    PetscPrintf(PETSC_COMM_WORLD,"Output file: %s\n",buf);
    out_vec.assign(buf);
                
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
       
    //load in the vector as local vector on each processor.
    PetscViewer viewer;
    PetscViewerCreate(PETSC_COMM_SELF,&viewer);
    Vec in_vec_vec;
    
    if ( is_ascii ) 
      {
	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, in_vec.c_str() ,&viewer);CHKERRQ(ierr); //ascii not working well for large sizes at the moment.
      } 
    else 
      {
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, in_vec.c_str() ,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      }
    //ierr = VecLoad(viewer,VECSEQ,&in_vec_vec);CHKERRQ(ierr);
    ierr = VecLoad(in_vec_vec, viewer);CHKERRQ(ierr);
    PetscViewerDestroy(&viewer);
        
    PetscInt vec_size;
    VecGetSize(in_vec_vec,&vec_size);
    
    ofstream out_vec_file(out_vec.c_str());
    if ( out_vec_file.is_open() )
      {
	out_vec_file.precision(14);
	PetscScalar *vals;
	VecGetArray(in_vec_vec, &vals);
	if ( use_complex  == PETSC_FALSE )
	  {
	  for ( int i = 0 ; i < vec_size;  i++ ) 
	    {
	      out_vec_file << scientific << PetscRealPart(vals[i]) << endl ; 
	    }
	  }
	else 
	  {
	    for ( int i = 0 ; i < vec_size;  i++ ) 
	    {
	      out_vec_file << scientific << "(" << PetscRealPart(vals[i]) << "," << PetscImaginaryPart(vals[i]) << ")" << endl ; 
	    }
	  }
	VecRestoreArray(in_vec_vec, &vals);	
	out_vec_file.close();
      }
    else
      {
	PetscPrintf(PETSC_COMM_WORLD,"Error opening file for writing %s,\n", out_vec.c_str());
	return 1;
      }       
    
    ierr = SlepcFinalize();CHKERRQ(ierr);
    return 0;
}
