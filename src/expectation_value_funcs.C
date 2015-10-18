#include "expectation_value_funcs.h"

#include "utils.h"

/* 
 * Function to read the expectation values input file and set the relevant parameters in the structure.
 */
int read_expectation_values_file(struct parameter_struct *parameters)
{
    if ( parameters->calculate_expectation_values == PETSC_TRUE )
      {
	//initialise the structure values
	parameters->exp_vals_info.number_vectors =0;
	parameters->exp_vals_info.momentum_sector = -1;
	parameters->exp_vals_info.rotation_sector = -1;
	parameters->exp_vals_info.filling_sector = -1;
	parameters->exp_vals_info.parity_sector = -1;			
	
	string filename(parameters->input_file);
	TiXmlDocument doc(parameters->expectation_values_file);
	bool loadOkay = doc.LoadFile();
	
	if ( !loadOkay )
	  {
	    PetscPrintf(PETSC_COMM_WORLD,"Could not load expectaction values file '%s'. Error='%s'. Exiting.\n", parameters->expectation_values_file , doc.ErrorDesc() );
	    return 1;
	  }
	  TiXmlHandle docHandle( &doc );
	
	TiXmlElement* element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("FILLING").Element();
	if ( element)
	  {
	    if ( element->Attribute("sector") != NULL ) 
	      {
		parameters->exp_vals_info.filling_sector = (PetscInt)atoi(element->Attribute("sector"));		
	      }
	  }
	  
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("ROTATION").Element();  
	if ( element)
	  {
	    if ( element->Attribute("sector") != NULL ) 
	      {
		parameters->exp_vals_info.rotation_sector = (PetscInt)atoi(element->Attribute("sector"));		
	      }
	  }      	    
	
	element = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("EIGEN_PAIRS").Element();  
	if ( element )
	  {
	    if ( element->Attribute("number") != NULL ) 
	      {
		parameters->exp_vals_info.number_vectors = (PetscInt)atoi(element->Attribute("number"));
	      }
	    if ( parameters->exp_vals_info.number_vectors > 0 ) 
	      {				
		if ( element->Attribute("spaces") != NULL ) 
		  {
		    parameters->exp_vals_info.number_spaces = atoi(element->Attribute("spaces"));		     
		  }
		if ( parameters->exp_vals_info.number_spaces > 0 ) 
		  {
		    PetscMalloc(parameters->exp_vals_info.number_spaces*sizeof(PetscInt),&parameters->exp_vals_info.number_vectors_per_space);
		    //PetscMalloc(parameters->exp_vals_info.number_spaces*sizeof(string*),&parameters->exp_vals_info.vectors);
		    parameters->exp_vals_info.vectors = new string*[parameters->exp_vals_info.number_spaces];		    
		    for ( int space = 0 ; space < parameters->exp_vals_info.number_spaces ; space ++ ) 
		      {
			TiXmlElement* element2 = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("EIGEN_PAIRS").Child("EIGENSPACE",space).Element();  
			if ( element2 ) 
			  {
			    parameters->exp_vals_info.number_vectors_per_space[space] = 0;
			    if ( element2->Attribute("degeneracy") != NULL ) 
			      {
				parameters->exp_vals_info.number_vectors_per_space[space] = atoi(element2->Attribute("degeneracy"));				
			      }
			    parameters->exp_vals_info.vectors[space] = new string[parameters->exp_vals_info.number_vectors_per_space[space]];
			    for ( int j = 0 ; j < parameters->exp_vals_info.number_vectors_per_space[space] ; j++ ) 
			      {
				TiXmlElement* element3 = docHandle.FirstChild( "SIMULATION" ).FirstChild( "PARAMETERS").FirstChild("EIGEN_PAIRS").Child("EIGENSPACE",space).Child("EIGENSTATE",j).Element();  
				if ( element3 ) 
				  {
				    if ( element3->Attribute("state_vector") != NULL ) 
				      {
					parameters->exp_vals_info.vectors[space][j].assign(element3->Attribute("state_vector"));
				      }				
				  }
			      }
			  }
		      }
		  }
	      }
	  }
      }
    return 0;
}

int calculate_expectation_values(struct parameter_struct* parameters, Mat *A)
{
    char op_filename[MAX_STRING];
    //stringstream ss();      
    //ss << parameters->output_prefix << "_expectation." << parameters->current_output_additional_prefix << ".xml" ;    
    sprintf(op_filename,"%s_expectation.%s.xml", parameters->output_prefix, parameters->current_output_additional_prefix);
    TiXmlDocument* doc;
    TiXmlElement * element;
    
    if (parameters->rank == 0 ) 
      {
	doc = new TiXmlDocument(op_filename);

	TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "", "" );
	doc->LinkEndChild( decl );

	element = new TiXmlElement( "EXPECTATION_VALUES" );
	element->SetAttribute("number", parameters->exp_vals_info.number_vectors);
	doc->LinkEndChild( element );
  
      }
      
    PetscViewer viewer;
    PetscErrorCode ierr;
    
    bool calc_overlap = false;
    Vec overlap ;
    Vec total_vec ; 
    
    
    if ( strcmp(parameters->exp_vals_info.overlap_vec.c_str(),"") != 0 ) 
      {
	calc_overlap = true;
	ierr = PetscViewerCreate(parameters->mat_solver_comm,&viewer);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(parameters->mat_solver_comm, parameters->exp_vals_info.overlap_vec.c_str() ,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
	//ierr = VecLoad(viewer,VECMPI,&overlap);CHKERRQ(ierr)
	ierr = VecLoad(overlap, viewer); CHKERRQ(ierr);
	//PetscViewerDestroy(viewer);
	PetscViewerDestroy(&viewer);
	ierr = VecCreate(parameters->mat_solver_comm,&total_vec);CHKERRQ(ierr);
	ierr = VecSetSizes(total_vec, get_rank_chunk_size(parameters->current_basis_size, parameters->rank, parameters->size), parameters->current_basis_size); CHKERRQ(ierr);
	ierr = VecSetType(total_vec, VECMPI); CHKERRQ(ierr);
	ierr = VecZeroEntries(total_vec);CHKERRQ(ierr);
      }
      
    
    for ( int vec_space = 0 ; vec_space < parameters->exp_vals_info.number_spaces ; vec_space++ ) 
      {
	TiXmlElement* element2;
	PetscInt dim = parameters->exp_vals_info.number_vectors_per_space[vec_space];
	if (parameters->rank == 0 ) 
	  {	
	    element2 = new TiXmlElement("EIGENSPACE");	    
	    element2->SetAttribute("dimension",dim);
	    element->LinkEndChild(element2);
	  }
	if ( dim == 1) 
	  {	    	    
	    Vec b,c;	
	    ierr = PetscViewerCreate(parameters->mat_solver_comm,&viewer);CHKERRQ(ierr);
	    ierr = PetscViewerBinaryOpen(parameters->mat_solver_comm, parameters->exp_vals_info.vectors[vec_space][0].c_str() ,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
	    //ierr = VecLoad(viewer,VECMPI,&b);CHKERRQ(ierr)
	    ierr = VecLoad(b, viewer); CHKERRQ(ierr);
	    //PetscViewerDestroy(viewer);
	    PetscViewerDestroy(&viewer);
	    PetscInt vec_size;
	    ierr = VecGetSize(b,&vec_size);CHKERRQ(ierr);
	    ierr = VecDuplicate(b, &c);CHKERRQ(ierr);
	    if ( vec_size != parameters->current_basis_size ) 
	      {
		PetscPrintf(parameters->mat_solver_comm, "Vector length does not match operator dimension.\n");
		return 1;
	      }
	      	      
	    ierr = MatMult(*A, b, c);CHKERRQ(ierr);
	    PetscScalar val; 
	    ierr = VecDot(c, b, &val); CHKERRQ(ierr);
	    
	    char state_filename[MAX_STRING];
	    sprintf(state_filename, "%s.%s_space_%d_vec_%d.vec" , parameters->output_prefix, parameters->current_output_additional_prefix, (int)vec_space, 0) ;
	    PetscViewerCreate(parameters->mat_solver_comm, &viewer);
	    PetscViewerBinaryOpen(parameters->mat_solver_comm, state_filename,FILE_MODE_WRITE,&viewer);		
	    VecView(b,viewer);
	    //PetscViewerDestroy(viewer);
	    PetscViewerDestroy(&viewer);
	    
	    sprintf(state_filename, "%s.%s_space_%d_vec_%d_ascii.vec" , parameters->output_prefix, parameters->current_output_additional_prefix, (int)vec_space, 0) ;
	    PetscViewerCreate(parameters->mat_solver_comm, &viewer);
	    PetscViewerASCIIOpen(parameters->mat_solver_comm, state_filename,&viewer);
	    VecView(b,viewer);
	    //PetscViewerDestroy(viewer);
	    PetscViewerDestroy(&viewer);
	    
	    PetscScalar overlap_val;	 
	    if ( calc_overlap ) 
	      {
		ierr = VecDot(overlap, b, &overlap_val);CHKERRQ(ierr);
		if ( PetscAbsScalar(overlap_val) > 1e-12 && PetscAbsScalar(val) < 1e-12 ) 
		  {
		      VecScale(b, overlap_val);
		      PetscScalar * vec_entries1, *vec_entries2;
		      ierr = VecGetArray(b, &vec_entries1);CHKERRQ(ierr);
		      ierr = VecGetArray(total_vec, &vec_entries2);CHKERRQ(ierr);
		      for ( PetscInt k = 0; k < get_rank_chunk_size(parameters->current_basis_size, parameters->rank, parameters->size); k++ ) 
			{
			  vec_entries2[k] += vec_entries1[k];
			}
		      ierr = VecRestoreArray(b, &vec_entries1);CHKERRQ(ierr);
		      ierr = VecRestoreArray(total_vec, &vec_entries2);CHKERRQ(ierr);		      
		  }
	      }
	      
	    if ( parameters->rank == 0 ) 
	      {
		TiXmlElement* element3 = new TiXmlElement("VECTOR");
		element3->SetAttribute("name", parameters->exp_vals_info.vectors[vec_space][0].c_str());
		element3->SetDoubleAttribute("val_real", PetscRealPart(val));
		element3->SetDoubleAttribute("val_imag", PetscImaginaryPart(val)); 
		if ( calc_overlap ) 
		  {		    
		    TiXmlElement* element4 = new TiXmlElement("OVERLAP");
		    stringstream ss;
		    ss.str("");
		    ss.setf(ios::scientific,ios::floatfield);
		    ss.precision(16);
		    ss << PetscRealPart(overlap_val);
		    TiXmlText *text = new TiXmlText(ss.str());		    
		    element4->LinkEndChild(text);element3->LinkEndChild(element4);
		  }
		element2->LinkEndChild(element3);
	      }	    
	    
	      
	    VecDestroy(&b);
	    VecDestroy(&c);
	  } 
	else 
	  {	    
	    //craete a matrix with the expectation values between vectors.
	    Mat E;
	    MatCreate(PETSC_COMM_SELF, &E);
	    MatSetSizes(E, PETSC_DECIDE, PETSC_DECIDE, dim, dim);
	    MatSetType(E, MATSEQDENSE);
	    
	    Vec *b,c;
	    
	    ierr = PetscMalloc(dim * sizeof(Vec), &b);CHKERRQ(ierr);      
	    for ( PetscInt x = 0 ; x < dim ; x++ ) 
	      {
		ierr = PetscViewerCreate(parameters->mat_solver_comm,&viewer);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(parameters->mat_solver_comm, parameters->exp_vals_info.vectors[vec_space][x].c_str() ,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
		//ierr = VecLoad(viewer,VECMPI,&b[x]);CHKERRQ(ierr)
		ierr = VecLoad(b[x], viewer);CHKERRQ(ierr);
		PetscViewerDestroy(&viewer);
	      }
	    ierr = VecDuplicate(b[0], &c);CHKERRQ(ierr);
	    
	    for ( PetscInt x = 0 ; x < dim ; x++ ) 
	      {		
		ierr = MatMult(*A, b[x], c);CHKERRQ(ierr);
		for ( PetscInt y = x ; y < dim ; y++ ) 
		  {		    
		    PetscScalar val; 
		    ierr = VecDot(c, b[y], &val); CHKERRQ(ierr);
		    ierr = MatSetValue(E, x, y, val, INSERT_VALUES);CHKERRQ(ierr);
		    if ( x != y ) 
		      {
			ierr = MatSetValue(E, y, x,PetscConj(val), INSERT_VALUES);CHKERRQ(ierr);
		      }		    		    
		  }		
	      }
	    ierr = MatAssemblyBegin(E, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	    ierr = MatAssemblyEnd(E, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	    //finished creating and population this matrix.
	    
	    //Now diagonalise this matrix.
	    EPS solver;  
	    ierr = EPSCreate(PETSC_COMM_SELF,&solver);CHKERRQ(ierr);
	    ierr = EPSSetOperators(solver,E,PETSC_NULL);CHKERRQ(ierr);   		
	    ierr = EPSSetProblemType(solver, EPS_HEP);CHKERRQ(ierr);
	    ierr = EPSSetType(solver, EPSLAPACK);CHKERRQ(ierr);	    
	    ierr = EPSSolve(solver);CHKERRQ(ierr);
	    
	    //Now take eigenvalues and use the eigenvectors to rotate the original vectors such that they are eigenstates of the overall operator here. 
	    PetscInt converged;
	    ierr = EPSGetConverged(solver,&converged);CHKERRQ(ierr);
	    if ( converged != dim ) PetscPrintf(parameters->mat_solver_comm,"Number of converged eigenvalues (%d) less than matrix dimensions(%d).", (int)converged, (int)dim);    	    
	    for (PetscInt vec_idx = 0 ; vec_idx < converged ; vec_idx++ ) 
	      {	
		ierr = VecZeroEntries(c);CHKERRQ(ierr);		
		PetscScalar valr, vali, *vec_entries, *tmp_vec_entries, *sol_vals;
		Vec solr, soli, tmpvec; 
		ierr = VecCreate(PETSC_COMM_SELF,&solr);CHKERRQ(ierr);
		ierr = VecSetSizes(solr, dim, dim); CHKERRQ(ierr);
		ierr = VecSetType(solr, VECSEQ); CHKERRQ(ierr);
		ierr = VecDuplicate(solr, &soli); CHKERRQ(ierr);		
		
		ierr = EPSGetEigenpair(solver,vec_idx,&valr,&vali,solr,soli);CHKERRQ(ierr);
				 
		PetscPrintf(PETSC_COMM_WORLD,"Eig val: %g \n",PetscRealPart(valr));
		//ierr = VecGetArray(c, &vec_entries); CHKERRQ(ierr);
		ierr = VecGetArray(solr, &sol_vals); CHKERRQ(ierr);		
		//ierr = VecMAXPY(c, dim, sol_vals, b); CHKERRQ(ierr);
		for ( int j = 0 ; j < dim ; j++ ) 
		  {
		    ierr = VecAXPY(c, PetscConj(sol_vals[j]), b[j]); CHKERRQ(ierr);
		    //ierr = VecGetArray(b[j], &tmp_vec_entries); CHKERRQ(ierr);
		    
		    //for ( PetscInt k = 0; k < get_rank_chunk_size(parameters->current_basis_size, parameters->rank, parameters->size); k++ ) 
		    //  {
			//vec_entries[k] += sol_vals[j] * tmp_vec_entries[k];
		      //}		    
		    PetscPrintf(PETSC_COMM_WORLD,"%g + i%g\n",PetscRealPart(sol_vals[j]), PetscImaginaryPart(sol_vals[j]));
		    //ierr = VecRestoreArray(b[j], &tmp_vec_entries); CHKERRQ(ierr);
		  }		
		//ierr = VecRestoreArray(c, &vec_entries); CHKERRQ(ierr);
		ierr = VecRestoreArray(solr, &sol_vals); CHKERRQ(ierr);		
		
		 {
		    PetscScalar *vals;
		    PetscReal l_max = 0.0; 
		    PetscInt max_idx; 
		    VecGetArray(c,&vals);
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
		    VecRestoreArray(c, &vals);
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
		    VecScale(c, 1.0/max_val);
		    ierr = PetscFree(all_vals); CHKERRQ(ierr);
		    ierr = PetscFree(send_buf); CHKERRQ(ierr);			    
		    VecNormalize(c, &l_max);
		  }
		
		char state_filename[MAX_STRING];
		sprintf(state_filename, "%s.%s_space_%d_vec_%d.vec" , parameters->output_prefix, parameters->current_output_additional_prefix, (int)vec_space, (int)vec_idx) ;
		PetscViewerCreate(parameters->mat_solver_comm, &viewer);
		PetscViewerBinaryOpen(parameters->mat_solver_comm, state_filename,FILE_MODE_WRITE,&viewer);		
		VecView(c,viewer);
		PetscViewerDestroy(&viewer);
		
		sprintf(state_filename, "%s.%s_space_%d_vec_%d_ascii.vec" , parameters->output_prefix, parameters->current_output_additional_prefix, (int)vec_space, (int)vec_idx) ;
		PetscViewerCreate(parameters->mat_solver_comm, &viewer);
		PetscViewerASCIIOpen(parameters->mat_solver_comm, state_filename,&viewer);
		VecView(c,viewer);
		PetscViewerDestroy(&viewer);
		
		PetscScalar overlap_val;	 
		if ( calc_overlap ) 
		  {
		    ierr = VecDot(overlap, c, &overlap_val);CHKERRQ(ierr);		    		    
		    if ( PetscAbsScalar(overlap_val) > 1e-10 && PetscAbsScalar(valr) < 1e-10 ) 
		      {
			VecScale(c, overlap_val);
			PetscScalar *vec_entries1, *vec_entries2;
			ierr = VecGetArray(c, &vec_entries1);CHKERRQ(ierr);
			ierr = VecGetArray(total_vec, &vec_entries2);CHKERRQ(ierr);
			for ( PetscInt k = 0; k < get_rank_chunk_size(parameters->current_basis_size, parameters->rank, parameters->size); k++ ) 
			  {
			    vec_entries2[k] += vec_entries1[k];
			  }
			ierr = VecRestoreArray(c, &vec_entries1);CHKERRQ(ierr);
			ierr = VecRestoreArray(total_vec, &vec_entries2);CHKERRQ(ierr);		      
		      }		    		    
		  }				
		
		if ( parameters->rank == 0 ) 
		  {				    		    
		    TiXmlElement* element3 = new TiXmlElement("VECTOR");
		    element3->SetAttribute("name", parameters->exp_vals_info.vectors[vec_space][vec_idx].c_str());
		    element3->SetDoubleAttribute("val_real", PetscRealPart(valr));
		    element3->SetDoubleAttribute("val_imag", PetscImaginaryPart(valr));		    
		    if ( calc_overlap ) 
		      {
			TiXmlElement* element4 = new TiXmlElement("OVERLAP");
			stringstream ss;
			ss.str("");
			ss.setf(ios::scientific,ios::floatfield);
			ss.precision(16);
			ss << PetscRealPart(overlap_val);
			TiXmlText *text = new TiXmlText(ss.str());			
			element4->LinkEndChild(text);element3->LinkEndChild(element4);
		      }
		    element2->LinkEndChild(element3);		
		  }		  
		VecDestroy(&solr);
		VecDestroy(&soli);
	      }
	    EPSDestroy(&solver);	    
	    MatDestroy(&E);
	    for ( int i = 0 ; i < dim ; i++ ) VecDestroy(&b[i]);
	    VecDestroy(&c);			          
	  }	   	   	  	
      }
    
    if ( calc_overlap ) 
      {		
	  {
	    PetscScalar *vals;
	    PetscReal l_max = 0.0; 
	    PetscInt max_idx; 
	    VecGetArray(total_vec,&vals);
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
	    VecRestoreArray(total_vec, &vals);
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
	    VecScale(total_vec, 1.0/max_val);
	    ierr = PetscFree(all_vals); CHKERRQ(ierr);
	    ierr = PetscFree(send_buf); CHKERRQ(ierr);			    
	    VecNormalize(total_vec, &l_max);			    
	  }
	
	char state_filename[MAX_STRING];
	sprintf(state_filename, "%s.%s_total.vec" , parameters->output_prefix, parameters->current_output_additional_prefix) ;
	PetscViewerCreate(parameters->mat_solver_comm, &viewer);
	PetscViewerBinaryOpen(parameters->mat_solver_comm, state_filename,FILE_MODE_WRITE,&viewer);
	VecView(total_vec,viewer);
	PetscViewerDestroy(&viewer);
	
	sprintf(state_filename, "%s.%s_total_ascii.vec" , parameters->output_prefix, parameters->current_output_additional_prefix) ;
	PetscViewerCreate(parameters->mat_solver_comm, &viewer);
	PetscViewerASCIIOpen(parameters->mat_solver_comm, state_filename, &viewer);
	VecView(total_vec,viewer);
	PetscViewerDestroy(&viewer);
	
	VecDestroy(&overlap);
	VecDestroy(&total_vec);
      }
     
    if ( parameters->rank == 0 ) 
      {
	if ( !doc->SaveFile()) return 1 ;
      }
    return 0;
}


int deallocate_expectaction_value_space(struct parameter_struct *parameters)
{
  if ( parameters->calculate_expectation_values == PETSC_TRUE && parameters->exp_vals_info.number_vectors )
    {
      if ( parameters->exp_vals_info.number_spaces > 0 ) 
	{
	  for (int space = 0 ; space < parameters->exp_vals_info.number_spaces ; space++ ) 
	    {
	      if ( parameters->exp_vals_info.number_vectors_per_space[space] > 0 ) 
		{
		  delete[] parameters->exp_vals_info.vectors[space];
		}	      
	    }
	  delete[] parameters->exp_vals_info.vectors;  
	  PetscFree(parameters->exp_vals_info.number_vectors_per_space);
	}      
    }
  return 0;
}
