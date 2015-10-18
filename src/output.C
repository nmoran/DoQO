/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 * 
 * File containing functions called from the main file solver.c  
 * 
 * Version: $Id$ 
 * 
*/

#include "output.h" 
#include "conservation.h"

#undef __FUNCT__
#define __FUNCT__ "write_xml_data"
/*! \brief Funcion that writes the xml file containing the results of an individual diagonalisation. */
int write_xml_data(struct parameter_struct *parameters, struct solver_results *results){
	char op_filename[MAX_STRING];

	sprintf(op_filename,"%s.%s.xml",parameters->output_prefix,parameters->current_output_additional_prefix);
	
	if ( parameters->rank == 0 ) {
		
		
		TiXmlDocument doc(op_filename);
		//fp = fopen(op_filename,"w");
	
		TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "", "" );
		doc.LinkEndChild( decl );
	
		TiXmlElement * element = new TiXmlElement( "SIMULATION" );
		doc.LinkEndChild( element );
		
		TiXmlElement * element2 = new TiXmlElement( "PARAMETERS" );
		element->LinkEndChild( element2 );
	
		TiXmlElement * element3 = new TiXmlElement( "TASK");
		stringstream ss;
		ss << parameters->current_task;
		TiXmlText *text = new TiXmlText(ss.str());
		element3->LinkEndChild(text);element2->LinkEndChild(element3);
		
		element3 = new TiXmlElement( "PROCESSES");
		ss.str("");ss << parameters->size;
		text = new TiXmlText(ss.str());
		element3->LinkEndChild(text);element2->LinkEndChild(element3);
		
		if ( parameters->use_parity_sectors == PETSC_TRUE ) { 
			PetscInt real_parity_sector  = parameters->parity_sector;
			if ( parameters->parity_info.num_relevant_sectors > 0 ) { 
					real_parity_sector = parameters->parity_info.relevant_sectors[parameters->parity_sector];
			}
			element3 = new TiXmlElement( "PARITY");
			element3->SetAttribute("sector",real_parity_sector);
			element3->SetAttribute("number_operators",parameters->parity_info.number_masks);
			element2->LinkEndChild(element3);
			for ( int i = 0 ; i < parameters->parity_info.number_masks ; i++ ) {
				TiXmlElement *element4 = new TiXmlElement("OPERATOR");
				ss.str("");ss << "P" << i ;
				element4->SetAttribute("name",ss.str());
				ss.str("");ss <<  get_mask_parity(real_parity_sector,i);
				text = new TiXmlText(ss.str());
				element4->LinkEndChild(text);element3->LinkEndChild(element4);
			}
		}
		
		if ( parameters->use_conservation_sectors == PETSC_TRUE ) { 
			PetscInt real_conservation_sector  = parameters->conservation_sector;
			if ( parameters->parity_info.num_relevant_sectors > 0 ) { 
					real_conservation_sector = parameters->parity_info.relevant_sectors[parameters->conservation_sector];
			}
			element3 = new TiXmlElement( "FILLING");
			element3->SetAttribute("sector",real_conservation_sector);
			element3->SetAttribute("number_operators",parameters->parity_info.number_masks);
			element2->LinkEndChild(element3);
			for ( int i = 0 ; i < parameters->parity_info.number_masks ; i++ ) {
				TiXmlElement *element4 = new TiXmlElement("OPERATOR");
				ss.str("");ss << "C" << i ;
				element4->SetAttribute("name",ss.str());
				ss.str("");ss << get_mask_filling(parameters,real_conservation_sector,i);
				text = new TiXmlText(ss.str());
				element4->LinkEndChild(text);element3->LinkEndChild(element4);
			}
		}
		
		if ( parameters->use_momentum_sectors == PETSC_TRUE )  {
			PetscInt real_momentum_sector  = parameters->momentum_sector;
			if ( parameters->momentum_info.num_relevant_sectors > 0 ) { 
				real_momentum_sector = parameters->momentum_info.relevant_sectors[parameters->momentum_sector];
			}
			element3 = new TiXmlElement( "MOMENTUM");
			element3->SetAttribute("sector",real_momentum_sector);
			element3->SetAttribute("number_dimensions",2);
			element2->LinkEndChild(element3);
			for ( int i = 0 ; i < 2 ; i++ ) {
				TiXmlElement *element4 = new TiXmlElement("DIRECTION");
				ss.str("");ss << i;
				element4->SetAttribute("name",ss.str());
				PetscReal mom ;
				if ( parameters->use_spectral_flow == PETSC_FALSE ) {				
					if ( i == 0 ) {
						mom = (real_momentum_sector / parameters->momentum_info.num_sectors[1]) * 2.0 * PETSC_PI / parameters->momentum_info.num_sectors[0];
					} else { 
						mom = (real_momentum_sector % parameters->momentum_info.num_sectors[1]) * 2.0 * PETSC_PI / parameters->momentum_info.num_sectors[1];
					}
				} else {
					if ( i == 0 ) {
						mom = ((PetscReal)(real_momentum_sector / parameters->momentum_info.num_sectors[1]) / (PetscReal)parameters->momentum_info.num_sectors[0] + (parameters->spectral_info.alpha[0]*(PetscReal)parameters->real_conservation_sector)/ (PetscReal)parameters->momentum_info.num_sectors[0] )* 2.0 * PETSC_PI ;
					} else { 
						//mom = (real_momentum_sector % parameters->momentum_info.num_sectors[1]) * 2.0 * PETSC_PI / parameters->momentum_info.num_sectors[1];
						mom = ((PetscReal)(real_momentum_sector % parameters->momentum_info.num_sectors[1]) / (PetscReal)parameters->momentum_info.num_sectors[0] + (parameters->spectral_info.alpha[1]*(PetscReal)parameters->real_conservation_sector)/ (PetscReal)parameters->momentum_info.num_sectors[1] )* 2.0 * PETSC_PI ;
					}
				}
				ss.str(""); ss << mom;
				text = new TiXmlText(ss.str());
				element4->LinkEndChild(text);element3->LinkEndChild(element4);
			}
		}
		
		
		//add something to the output file identifying the spectral flow point that is being used.
		if ( parameters->use_spectral_flow  ){
			element3 = new TiXmlElement( "SPECTRAL_FLOW");	
			element3->SetAttribute("sector",parameters->spectral_info.real_current_point[1]*parameters->spectral_info.number_points[0] + parameters->spectral_info.real_current_point[0]);
			element2->LinkEndChild(element3);
		}
		
		if ( parameters->use_rotation_invariance )  {
			element3 = new TiXmlElement( "ROTATION");
			if ( parameters->rotation_number_relevant_sectors > 0 ) 
			{
			    element3->SetAttribute("sector",parameters->rotation_relevant_sectors[parameters->rotation_current_sector]);
			} 
			else 
			{
			    element3->SetAttribute("sector",parameters->rotation_current_sector);
			}
			element2->LinkEndChild(element3);
			/*PetscReal rot ;
			rot = ((PetscReal)parameters->rotation_info.real_sector / (PetscReal)parameters->rotation_info.num_sectors) * 2.0 * PETSC_PI ;
			ss.str(""); ss << rot;
			text = new TiXmlText(ss.str());
			element3->LinkEndChild(text);*/
		}
					
		element3 = new TiXmlElement( "BASIS_SIZE");
		ss.str(""); ss << parameters->current_basis_size;
		text = new TiXmlText(ss.str());
		element3->LinkEndChild(text);element2->LinkEndChild(element3);	
			
		element3 = new TiXmlElement("EIGENVALUES_CONVERGED");
		ss.str(""); ss << results->eigenvalues_converged;
		text = new TiXmlText(ss.str());
		element3->LinkEndChild(text);element2->LinkEndChild(element3);
		
		element3 = new TiXmlElement("ITERATIONS_TAKEN");
		ss.str("");ss << results->iterations_taken;
		text = new TiXmlText(ss.str());
		element3->LinkEndChild(text);element2->LinkEndChild(element3);
		
		element3 = new TiXmlElement("TIME_TAKEN");
		ss.str("");ss << results->time_taken;
		text = new TiXmlText(ss.str());
		element3->LinkEndChild(text);element2->LinkEndChild(element3);
		
		element3 = new TiXmlElement("DEGENERATE_BANDS");
		ss.str("");ss << results->degenerate_bands;
		text = new TiXmlText(ss.str());
		element3->LinkEndChild(text);element2->LinkEndChild(element3);
		
		if ( parameters->calculate_correlations == PETSC_TRUE ) {
			element3 = new TiXmlElement("CORRELATIONS");
			element3->SetAttribute("number",parameters->number_ops);
			element2->LinkEndChild(element3);
			
			//if we have a trial wave function and correlations
			if ( parameters->prod_wf_overlap ) {
				element3 = new TiXmlElement("TRIAL_CORRELATIONS");
				element3->SetAttribute("number",parameters->number_ops);
				element2->LinkEndChild(element3);
				for ( int j = 0 ; j < parameters->number_ops ; j++ ) {
					TiXmlElement *element4 = new TiXmlElement("CORRELATION");
					ss.str("");
					ss << j; 
					element4->SetAttribute("number",ss.str());
					ss.str("");
					ss << PetscRealPart(results->trial_correlations[j]) << " + i" << PetscImaginaryPart(results->trial_correlations[j]);
					text = new TiXmlText(ss.str());
				 	element4->LinkEndChild(text);
 					element3->LinkEndChild(element4);
 				}
			}
		}
		
		//put in values of the model parameters for this run
		for ( int i = 0 ; i < parameters->hamiltonian.number_parameters ; i++ ) { 
			element3 = new TiXmlElement("TASK_PARAMETER");
			element3->SetAttribute("name",parameters->tasks.parameter_names[parameters->hamiltonian.parameter_indices[i]]);
			ss.setf(ios::scientific,ios::floatfield);
			ss.precision(16);
			#if defined(PETSC_USE_COMPLEX) 
			ss.str(""); ss << PetscRealPart(parameters->tasks.task_parameter_values[parameters->current_task][parameters->hamiltonian.parameter_indices[i]]) << " + i" 
			<< PetscImaginaryPart(parameters->tasks.task_parameter_values[parameters->current_task][parameters->hamiltonian.parameter_indices[i]])  ;
			#else
			ss.str(""); ss << parameters->tasks.task_parameter_values[parameters->current_task][parameters->hamiltonian.parameter_indices[i]] ;
			#endif
			text = new TiXmlText(ss.str());
			element3->LinkEndChild(text);element2->LinkEndChild(element3);
		}

		element3 = new TiXmlElement("EIGEN_PAIRS");	
		element3->SetAttribute("number",results->eigenvalues_converged);
		element3->SetAttribute("spaces",results->degenerate_bands);
		element2->LinkEndChild(element3);
			
		for ( int i = 0 ; i < results->degenerate_bands ; i++ ) {
			TiXmlElement *element4 = new TiXmlElement("EIGENSPACE");
			element4->SetAttribute("number",i);
			element4->SetAttribute("degeneracy",results->degeneracies[i+1] - results->degeneracies[i]);
			ss.str("");
			ss.setf(ios::scientific,ios::floatfield);
			ss.precision(16);
			ss << results->eigenvalues[results->degeneracies[i]];
			//sprintf(buf,"%le",results->eigenvalues[results->degeneracies[i]]);
			element4->SetAttribute("energy",ss.str());
			element3->LinkEndChild(element4);
			
		  	//fprintf(fp,"\t\t<EIGENSPACE number=\"%d\" degeneracy=\"%d\" energy=\"%.16e\">\n",i,results->degeneracies[i+1] - results->degeneracies[i],results->eigenvalues[results->degeneracies[i]]);
			for ( int j = 0 ; j < results->degeneracies[i+1] - results->degeneracies[i] ; j++ ) {
				TiXmlElement *element5 = new TiXmlElement("EIGENSTATE");
				element5->SetAttribute("number",results->degeneracies[i]+j);
				//sprintf(buf,"%le",results->errors[results->degeneracies[i]+j]);
				ss.str(""); ss << results->errors[results->degeneracies[i]+j];
				element5->SetAttribute("error",ss.str());
				if ( parameters->save_states == PETSC_TRUE ) {
					ss.str("") ; ss << parameters->output_prefix << "." << parameters->current_output_additional_prefix << "_state_" << results->degeneracies[i]+j  << ".vec";
					//sprintf(buf, "%s.%s_state_%d.vec" , parameters->output_prefix, parameters->current_output_additional_prefix, (int)(results->degeneracies[i]+j) ) ;
					element5->SetAttribute("state_vector",ss.str());
				}
				if ( parameters->prod_wf_overlap ) {
					ss.str("");ss << results->overlaps[results->degeneracies[i]+j];
					element5->SetAttribute("overlap",ss.str());
					//cout << "Overlap : " << ss.str() << endl;
				}
	
				if ( parameters->other_wf_overlap ) {
					ss.str("");ss << results->other_overlaps[results->degeneracies[i]+j];
					element5->SetAttribute("other_overlap",ss.str());
				}
				
				ss.str(""); ss << results->eigenvalues[results->degeneracies[i]+j];
				//sprintf(buf,"%le",results->eigenvalues[results->degeneracies[i]+j]);
				text = new TiXmlText(ss.str());
				element5->LinkEndChild(text);element4->LinkEndChild(element5);
			 }//end loop over states in degenerate subspace.
			 
			 if ( parameters->calculate_correlations == PETSC_TRUE ) {
				 TiXmlElement *element6 = new TiXmlElement("CORRELATIONS");
				 ss.str("");
				 ss << parameters->number_ops ;  
				 element6->SetAttribute("number",ss.str());
				 for ( int k = 0 ; k < parameters->number_ops ; k++ ) {
				 	int dim = results->degeneracies[i+1] - results->degeneracies[i];
				 	TiXmlElement *element7 = new TiXmlElement("CORRELATION");
				 	ss.str("");
				 	ss << k; 
				 	element7->SetAttribute("number",ss.str());
				 	ss.str("");
				 	ss <<  dim;
				 	element7->SetAttribute("dimension",ss.str());
				 	ss.str("");
				 	for ( int l = 0 ; l < dim ; l++ ) {
				 		for ( int m = 0 ; m < dim ; m++ ) {
				 			if ( m > 0 ) {
				 				ss << ", " ; 
				 			}
				 			ss << PetscRealPart(results->correlations[i][k][l*dim+m]) << " + i" << PetscImaginaryPart(results->correlations[i][k][l*dim+m]) ;
				 		}
				 		ss << "; ";
				 	}
				 	text = new TiXmlText(ss.str());
				 	element7->LinkEndChild(text);
				 	element6->LinkEndChild(element7);
				 }
				 element4->LinkEndChild(element6);
			 }
		} //end loop over degenerate bands.
		
		if ( !doc.SaveFile()) return 1 ;
	} //if rank == 0
	

	return 0 ;
}

#undef __FUNCT__
#define __FUNCT__ "write_main_xml_output"
/*! \brief Funcion that writes the main xml output file  . */
int write_main_xml_output(struct parameter_struct *parameters){
	char buf[MAX_STRING];
	sprintf(buf,"%s.output.xml",parameters->output_prefix); 
	
	TiXmlDocument doc(buf);
	
	TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "", "" );
	doc.LinkEndChild( decl );
	
	TiXmlElement * element = new TiXmlElement( "SIMULATION" );
	doc.LinkEndChild( element );
	
	TiXmlElement * element2 = new TiXmlElement( "PARAMETERS" );
	element->LinkEndChild( element2 );
	
	TiXmlElement * element3 = new TiXmlElement( "MODEL_TYPE" );
	TiXmlText *text;
	if ( parameters->model_type == SPIN_HALF ) { 
		 text = new TiXmlText("SPIN_HALF");
	} else if ( parameters->model_type == FERMIONIC ) { 
		 text = new TiXmlText("FERMIONIC");
	} else if ( parameters->model_type == SPIN_ONE ) { 
		 text = new TiXmlText("SPIN_ONE");
	} else {
		text = new TiXmlText("UNKNOWN");
	} 
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );
	
	element3 = new TiXmlElement( "CODE_REVISION" );
	text = new TiXmlText(CODE_REV); 
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );
	
	element3 = new TiXmlElement( "SCALAR_TYPE" );
#if defined(PETSC_USE_COMPLEX) 
			text = new TiXmlText("complex"); 
#else
			text = new TiXmlText("real"); 
#endif
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );

	element3 = new TiXmlElement( "INTEGER_SIZE" );
#if defined(PETSC_USE_64BIT_INDICES)
		text = new TiXmlText("64bits"); 
#else
		text = new TiXmlText("32bits"); 
#endif
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );
	
	
	element3 = new TiXmlElement( "MODEL_FILE" );
	text = new TiXmlText(parameters->hamiltonian.filename); 
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );
	
	element3 = new TiXmlElement( "TASK_LIST" );
	text = new TiXmlText(parameters->tasks.filename); 
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );
	
	element3 = new TiXmlElement( "SOLVER_TYPE" );
	text = new TiXmlText(parameters->method); 
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );
	
	element3 = new TiXmlElement( "MAX_ITERATIONS" );
	sprintf(buf,"%d",(int)parameters->max_its);
	text = new TiXmlText(buf); 
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );
	
	element3 = new TiXmlElement( "SOLVER_TOLERANCE" );
	sprintf(buf,"%le",parameters->tol);
	text = new TiXmlText(buf); 
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );
	
	element3 = new TiXmlElement( "DEGENERACY_TOLERANCE" );
	sprintf(buf,"%le",parameters->deg_tol);
	text = new TiXmlText(buf); 
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );
	
	element3 = new TiXmlElement( "PHASE_TOLERANCE" );
	sprintf(buf,"%le",parameters->phase_tol);
	text = new TiXmlText(buf); 
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );

	element3 = new TiXmlElement( "MOMENTUM_SECTORS" );
	if ( parameters->use_momentum_sectors == PETSC_TRUE ) {
		text = new TiXmlText("enabled"); 
		sprintf(buf,"%d",(int)(parameters->momentum_info.num_sectors[0]*parameters->momentum_info.num_sectors[1]));
		element3->SetAttribute("file",parameters->momentum_info_file);
		element3->SetAttribute("number",buf);
	} else {
		text = new TiXmlText("disabled"); 
	}
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );
	
	element3 = new TiXmlElement( "PARITY_SECTORS" );
	if ( parameters->use_parity_sectors == PETSC_TRUE ) {
		text = new TiXmlText("enabled"); 
		sprintf(buf,"%d",parameters->parity_info.number_parity_sectors);
		element3->SetAttribute("file",parameters->parity_masks);
		element3->SetAttribute("number",buf);
	} else {
		text = new TiXmlText("disabled"); 
	}
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );
	
	element3 = new TiXmlElement( "FILLING_SECTORS" );
	if ( parameters->use_conservation_sectors == PETSC_TRUE ) {
		text = new TiXmlText("enabled"); 
		sprintf(buf,"%d",parameters->parity_info.number_parity_sectors);
		element3->SetAttribute("file",parameters->conservation_masks);
		element3->SetAttribute("number",buf);
	} else {
		text = new TiXmlText("disabled"); 
	}
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );	
	
	if ( parameters->use_spectral_flow == PETSC_TRUE ) {
		element3 = new TiXmlElement( "SPECTRAL_FLOW" );
		text = new TiXmlText("enabled"); 
		sprintf(buf,"%d",(int)parameters->spectral_info.number_points[0]);
		element3->SetAttribute("number_points_0",buf);
		sprintf(buf,"%d",(int)parameters->spectral_info.number_points[1]);
		element3->SetAttribute("number_points_1",buf);
		element3->LinkEndChild(text); element2->LinkEndChild( element3 );
	}

	element3 = new TiXmlElement( "ROTATION_SECTORS" );
	if ( parameters->use_rotation_invariance == PETSC_TRUE ) {	
		text = new TiXmlText("enabled"); 
		if ( parameters->rotation_number_relevant_sectors != -1 ) 
		{
			sprintf(buf,"%d",(int)parameters->rotation_number_relevant_sectors);
		} else sprintf(buf,"%d",(int)parameters->number_rotation_sectors);
		element3->SetAttribute("number_sectors",buf);				
	}  else {
		text = new TiXmlText("disabled"); 
	}
	element3->LinkEndChild(text); element2->LinkEndChild( element3 );
		
	element2 = new TiXmlElement("RESULTS");
	element->LinkEndChild( element2 );
	
	int num_parity_sectors,num_conservation_sectors,num_momentum_sectors,num_rotation_sectors;
	//work out how many sectors there are taking into consideration parity, conservation and momentum
	if ( parameters->use_parity_sectors == PETSC_TRUE ) {
		num_parity_sectors = parameters->parity_info.number_parity_sectors;
		if ( parameters->parity_info.num_relevant_sectors > 0 ) num_parity_sectors = parameters->parity_info.num_relevant_sectors;
	} else { 
		num_parity_sectors = 1;
	}
	
	if ( parameters->use_conservation_sectors == PETSC_TRUE ) {
		num_conservation_sectors = parameters->parity_info.number_parity_sectors;
		if ( parameters->parity_info.num_relevant_sectors > 0 ) num_conservation_sectors = parameters->parity_info.num_relevant_sectors;
	} else { 
		num_conservation_sectors = 1;
	}
	
	if ( parameters->use_momentum_sectors == PETSC_TRUE ) { 
		num_momentum_sectors = parameters->momentum_info.num_sectors[0] * parameters->momentum_info.num_sectors[1];
		if ( parameters->momentum_info.num_relevant_sectors > 0 ) num_momentum_sectors = parameters->momentum_info.num_relevant_sectors ; 
	} else {
		num_momentum_sectors = 1;
	}

	if ( parameters->use_rotation_invariance == PETSC_TRUE ) { 
		if ( parameters->rotation_number_relevant_sectors != -1 ) 
		{
			num_rotation_sectors = parameters->rotation_number_relevant_sectors ;
		} else num_rotation_sectors = parameters->number_rotation_sectors ;
	} else {
		num_rotation_sectors = 1;
	}
	
	
	//loop over tasks
	for ( int task = 0 ; task < parameters->tasks.number_tasks ; task++ ) {
		for ( int parity_sector = 0 ; parity_sector < num_parity_sectors ; parity_sector++ ) { 
			for ( int conservation_sector = 0 ; conservation_sector < num_conservation_sectors ; conservation_sector++ ) { 
				for ( int current_point_1 = 0 ; current_point_1 < parameters->spectral_info.number_relevant_points[1]; current_point_1++) { 
					for ( int current_point_0 = 0 ; current_point_0 < parameters->spectral_info.number_relevant_points[0]; current_point_0++) { 
						for ( int momentum_sector = 0 ; momentum_sector < num_momentum_sectors ; momentum_sector++ ) { 
							for ( int rotation_sector = 0 ; rotation_sector < num_rotation_sectors ; rotation_sector++ ) { 
		
								stringstream ss;
								ss << parameters->output_prefix << ".task_" << task ;
								element3 = new TiXmlElement("OUTPUT");
								element3->SetAttribute("task",task);
								
								if ( parameters->use_parity_sectors ) {
									int real_parity_sector = ( parameters->parity_info.num_relevant_sectors > 0 ) ? parameters->parity_info.relevant_sectors[parity_sector]: parity_sector;
									ss << "_parity_" << real_parity_sector;
									element3->SetAttribute("parity_sector",real_parity_sector);
								} else if ( parameters->use_conservation_sectors ) {
									int real_conservation_sector = ( parameters->parity_info.num_relevant_sectors > 0 ) ? parameters->parity_info.relevant_sectors[conservation_sector]: conservation_sector;
									ss << "_filling_" << real_conservation_sector;
									element3->SetAttribute("filling_sector",real_conservation_sector);
								}
								
								if ( parameters->use_spectral_flow ) {
									int real_current_0 = ( parameters->spectral_info.number_relevant_points[0] < parameters->spectral_info.number_points[0] ) ? parameters->spectral_info.relevant_points[0][current_point_0]: current_point_0;
									int real_current_1 = ( parameters->spectral_info.number_relevant_points[1] < parameters->spectral_info.number_points[1] ) ? parameters->spectral_info.relevant_points[1][current_point_1]: current_point_1;
									ss << "_spectral_" << real_current_1 * parameters->spectral_info.number_points[0] + real_current_0;
									element3->SetAttribute("spectral_flow_point",real_current_1 * parameters->spectral_info.number_points[0] + real_current_0);
								}
								
								if ( parameters->use_momentum_sectors ) {
									int real_momentum_sector = ( parameters->momentum_info.num_relevant_sectors > 0 ) ? parameters->momentum_info.relevant_sectors[momentum_sector]: momentum_sector;
									ss << "_momentum_" <<  real_momentum_sector;
									element3->SetAttribute("momentum_sector",real_momentum_sector);
								}
								
								if ( parameters->use_rotation_invariance ) {
									if ( parameters->rotation_number_relevant_sectors > 0 ) 
									{
										ss << "_rotation_" <<  parameters->rotation_relevant_sectors[rotation_sector];
										element3->SetAttribute("rotation_sector",parameters->rotation_relevant_sectors[rotation_sector]);
									} 
									else 
									{
										ss << "_rotation_" <<  rotation_sector;
										element3->SetAttribute("rotation_sector",rotation_sector);
									}									
								}
								ss << ".xml" ;
								text = new TiXmlText(ss.str().c_str());
								element3->LinkEndChild(text); element2->LinkEndChild(element3);
							}//end loop over rotation sectors							
						}//end loop over momentum sectors
					} //end loop over spectral flow sectors in direction 0
				} //end loop over spectral flow sectors in direction 1
			} // end loop over conservation sectors
		} //end loop over parity sectors
	} //end loop over tasks

	

	if ( !doc.SaveFile() ) {
		return 1;
	}
	
	return 0;
}
