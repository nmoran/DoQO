#include "slepceps.h"
#include "slepcsvd.h"
#include <fstream>
#include <iostream>

#define MAX_STRING 200

using namespace std;

typedef PetscBool PetscTruth;

int trim(const string str, string *trimmed) {
	size_t start = str.find_first_not_of(" \t\n\r");
	if(start == string::npos) return 1;
	trimmed->assign(str,start,str.find_last_not_of(" \t\n\r") - start + 1);
	return 0; 
}

int main( int argc, char **argv ) {
    PetscErrorCode ierr;     
    PetscTruth use_complex = PETSC_FALSE, is_ascii = PETSC_FALSE;
    PetscMPIInt rank, size;
    stringstream ss;           
    PetscTruth flg;
    
    /*------------------------------------------------------------------------------------------------
    Initialise some values.
    ------------------------------------------------------------------------------------------------*/
    char help[] = "Program that sums vectors given in file according to weight given in file also.";
    
    SlepcInitialize(&argc,&argv,(char*)0,help);	
    
    //get rank of this process and number of processes 
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr); 
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);        
    
    
    char buf[MAX_STRING];
    ierr = PetscOptionsGetString(NULL,"-vecs",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	PetscPrintf(PETSC_COMM_WORLD,"Must specify file which lists vectors and weights with -vecs option.\n");
	return 1;
      }
    PetscPrintf(PETSC_COMM_WORLD,"Input file: %s\n",buf);
    string vecs_file;
    vecs_file.assign(buf);
    
    ierr = PetscOptionsGetString(NULL,"-o",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE) 
      {
	PetscPrintf(PETSC_COMM_WORLD,"Must specify output vector on the command line with switch -o.\n");
	return 1;
      }
    PetscPrintf(PETSC_COMM_WORLD,"Output file: %s\n",buf);
    string out_vec;
    out_vec.assign(buf);
                
    ierr = PetscOptionsGetString(NULL,"-ascii",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    if ( flg == PETSC_FALSE ) 
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
          
    ierr = PetscOptionsGetString(NULL,"-overlap_vec",buf,MAX_STRING,&flg);CHKERRQ(ierr);
    string overlap_vec_filename;
    bool use_overlap = false;
    if ( flg == PETSC_TRUE) 
      {	  
	PetscPrintf(PETSC_COMM_WORLD,"Overlap vec: %s\n",buf);	
	overlap_vec_filename.assign(buf);  	
	use_overlap = true;
      }
    
      
    /*use_complex = PETSC_FALSE;
    ierr = PetscOptionsGetString(NULL,"-complex",buf,MAX_STRING,&flg);CHKERRQ(ierr);
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
      }*/    
       
    char **VecNames; 
    PetscScalar *VecWeights; 
    //read in information from file
    int NumberVectors = 0;
    ifstream in_file( vecs_file.c_str(), ios_base::in);
    if ( in_file.is_open() )
      {
	string line;
	while( getline(in_file, line) )
	  {
	    PetscPrintf(PETSC_COMM_WORLD,"%s\n",line.c_str());
	    NumberVectors++;
	  }		  	  	  
	in_file.close();
	ierr = PetscMalloc(NumberVectors*sizeof(char*),&VecNames);CHKERRQ(ierr);
	ierr = PetscMalloc(NumberVectors*sizeof(PetscScalar),&VecWeights);CHKERRQ(ierr);
	for (int i = 0 ; i < NumberVectors ; i++ ) 
	  {
	    ierr = PetscMalloc(MAX_STRING*sizeof(char),&VecNames[i]);CHKERRQ(ierr);
	  }
	in_file.open( vecs_file.c_str(), ios_base::in);
	int count = 0;
	if ( in_file.is_open() )
	  {
	    while( getline(in_file, line) )
	      {
		line = line.substr(line.find_first_not_of(" \t\n\r"),line.length()-line.find_first_not_of(" \t\n\r")); //trim whitespace at start.
		string leftside = line.substr(0,line.find_first_of(" \t\n\r"));
		string rightside = line.substr(line.find_first_of(" \t\n\r"),line.length()-line.find_first_of(" \t\n\r"));
		rightside = rightside.substr(rightside.find_first_not_of(" \t\n\r"),rightside.length()-rightside.find_first_not_of(" \t\n\r")); //trim whitespace at start.
		PetscPrintf(PETSC_COMM_WORLD,"Line: %d, <%s>,<%s>\n", count, leftside.c_str(), rightside.c_str());
		strcpy(VecNames[count], leftside.c_str());
		VecWeights[count] = atof(rightside.c_str()); 
		PetscPrintf(PETSC_COMM_WORLD,"Line: %d, <%s>,<%.14g>\n", count, VecNames[count], PetscRealPart(VecWeights[count]));
		count++;		
	      }		  	  	  
	  }
      }
    else
      {
	PetscPrintf(PETSC_COMM_SELF,"Problem opening file %s. Quiting\n", vecs_file.c_str());
	return 1;
      }
                  
    Vec *vecs_array;
    ierr = PetscMalloc(NumberVectors*sizeof(Vec), &vecs_array);CHKERRQ(ierr);
    
    for ( int vec_idx = 0 ; vec_idx < NumberVectors ; vec_idx++ ) 
      {
	PetscViewer viewer;
	PetscViewerCreate(PETSC_COMM_SELF,&viewer);	
	
	if ( is_ascii ) 
	  {
	    ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, VecNames[vec_idx] ,&viewer);CHKERRQ(ierr); //ascii not working well for large sizes at the moment.
	  } 
	else 
	  {
	    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, VecNames[vec_idx] ,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
	  }
	//ierr = VecLoad(viewer,VECSEQ, &vecs_array[vec_idx]);CHKERRQ(ierr);
	ierr = VecLoad(vecs_array[vec_idx], viewer);CHKERRQ(ierr);
	PetscViewerDestroy(&viewer);
      }
      
    Vec overlap_vec;   
    if ( use_overlap ) 
      {
	PetscViewer viewer;
	PetscViewerCreate(PETSC_COMM_SELF,&viewer);	
	
	if ( is_ascii ) 
	  {
	    ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, overlap_vec_filename.c_str() ,&viewer);CHKERRQ(ierr); //ascii not working well for large sizes at the moment.
	  } 
	else 
	  {
	    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, overlap_vec_filename.c_str() ,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
	  }
	//ierr = VecLoad(viewer,VECSEQ, &overlap_vec);CHKERRQ(ierr);
	ierr = VecLoad(overlap_vec, viewer);CHKERRQ(ierr);
	PetscViewerDestroy(&viewer);
	
	PetscPrintf(PETSC_COMM_WORLD, "Overlap weights\n");
	for ( int vec_idx = 0 ; vec_idx < NumberVectors ; vec_idx++ ) 
	  {
	    ierr = VecDot(vecs_array[vec_idx], overlap_vec, &VecWeights[vec_idx]);CHKERRQ(ierr);
	    PetscPrintf(PETSC_COMM_WORLD,"Line: %d, <%s>,<%.14g>\n", vec_idx, VecNames[vec_idx], PetscRealPart(VecWeights[vec_idx]));
	  }
	VecDestroy(&overlap_vec);
      }
           
    Vec C;
    ierr = VecDuplicate(vecs_array[0], &C); CHKERRQ(ierr);
    VecZeroEntries(C);
    VecMAXPY(C, NumberVectors, VecWeights, vecs_array);    
    
    PetscReal norm;
    ierr = VecNormalize(C, &norm);CHKERRQ(ierr);
    
    PetscViewer viewer;
    PetscViewerCreate(PETSC_COMM_SELF,&viewer);	
    
    if ( is_ascii ) 
      {
	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, out_vec.c_str() ,&viewer);CHKERRQ(ierr); //ascii not working well for large sizes at the moment.
      } 
    else 
      {
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, out_vec.c_str() ,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      }
    ierr = VecView	(C, viewer);CHKERRQ(ierr);
    PetscViewerDestroy(&viewer);
    
    for ( int vec_idx = 0 ; vec_idx < NumberVectors  ; vec_idx++ ) 
     {
       VecDestroy(&vecs_array[vec_idx]); 
       PetscFree(VecNames[vec_idx]);
     }
    PetscFree(VecNames);
    PetscFree(VecWeights);
    VecDestroy(&C);
        
    
    ierr = SlepcFinalize();CHKERRQ(ierr);
    return 0;
}
