/* 
 * Copyright (C) 2010    Niall Moran, Graham Kells, Jiri Vala
 * 
 * 
*/

#include "tasks.h"


#undef __FUNCT__
#define __FUNCT__ "parse_task_file"
/*! \brief Function that reads task values and names.*/ 
int parse_task_file(struct task_parameters *tasks, struct parameter_struct *parameters, struct matrix_meta_data *md ) {
	int i ,j  ;
	PetscErrorCode ierr;
	string line;
	vector<string> lines; 

	//count the number of tasks.
	ifstream taskfile;
	taskfile.open(tasks->filename,ifstream::in);
	
	if (!taskfile.is_open()) {
		return 1;
	}
	
	tasks->number_tasks = 0;
	while ( getline(taskfile,line) ) {
		tasks->number_tasks++;
		lines.push_back(line);
	}
	taskfile.close();

	//allocate space for parameter values
	ierr = PetscMalloc(tasks->number_tasks * sizeof(PetscScalar*),&tasks->task_parameter_values);CHKERRQ(ierr);
	
	for ( i = 0 ; i < tasks->number_tasks ; i++ ) { 
		ierr = PetscMalloc(md->number_parameters * sizeof(PetscScalar),&tasks->task_parameter_values[i]);CHKERRQ(ierr);
		//initialise all parameters to 0
		for ( j = 0 ; j < md->number_parameters ; j++ ) {
			tasks->task_parameter_values[i][j] = 0.0; 
		}
	}
	//done allocating space
	
	//now read in each line and process it.
	int task = 0;
	for ( vector<string>::iterator line_iter = lines.begin() ; line_iter != lines.end(); line_iter++ ) {
		vector<string> parts ; 
		string delim(",");
		tokenize_string(*line_iter,delim,&parts); //split line at each comma into a vector of substrings.
		//for ( vector<string>::iterator iter = parts.begin() ; iter != parts.end() ; iter++ ) { //iterate over each substring and set appropriate parameter value. 
		//	cout << *iter << endl ;
		//}
		for ( vector<string>::iterator iter = parts.begin() ; iter != parts.end() ; iter++ ) { //iterate over each substring and set appropriate parameter value. 
			vector<string> parts2; 
			string delim2("=");
			tokenize_string(*iter,delim2,&parts2); //split the part at the equals sign. Parameter name will be on the left and value on the right. 
			vector<string>::iterator iter2 = parts2.begin();
			string param_name;
			trim(*iter2,&param_name); //trim any excess space
			//now find the index of the parameter and leave as -1 if not found.
			int param_num = -1;
			for ( i = 0 ; i < tasks->number_parameters  ; i++ ) {
				if ( param_name.compare(tasks->parameter_names[i]) == 0 ) {
					param_num = i ;
				}
			}
			
			//now read the value
			*iter2++;
			PetscScalar param_value = str2scalar(*iter2);
			if ( param_num != -1 ) {
				tasks->task_parameter_values[task][param_num] = param_value;
			}
		} //end loop over separate parameter assignment values.
		task += 1; 
	} //end while loop over lines 

	return 0;
}




#undef __FUNCT__
#define __FUNCT__ "print_task_info"
/*! /brief Function that prints out information on the tasks.*/
void print_task_info(struct task_parameters *tasks , struct matrix_meta_data *md) {
	int i, j ; 
	
	PetscFPrintf(PETSC_COMM_WORLD, stderr,"%d Tasks\n", tasks->number_tasks ) ;
	for ( i = 0 ; i < tasks->number_tasks; i++ ) {
		PetscFPrintf(PETSC_COMM_WORLD, stderr,"Task %d\n--------\n", i ) ;
		for ( j = 0 ; j < md->number_parameters ; j++ ) {
			//PetscFPrintf(PETSC_COMM_WORLD,stderr, "%s: %lf + i%lf\n", md->parameter_names[j], tasks->task_parameter_values[i][j],tasks->imag_parameter_values[i][j]) ;
		}
	}
	PetscFPrintf(PETSC_COMM_WORLD, stderr,"\n" ) ;
}

#undef __FUNCT__
#define __FUNCT__ "apply_task_parameters"
/*! \brief Use the task parameters for task with index 'task'. */
int apply_task_parameters(struct matrix_meta_data *md, struct task_parameters *tasks , int task ) {
	int j, k ;
	
	for ( j = 0 ; j < md->number_terms ; j++ ) {
		if ( md->number_term_params[j] > 0 ) {
			md->term_multipliers[j] = tasks->task_parameter_values[task][md->term_parameters[j][0]] ;	
		}
		
		for ( k = 1 ; k < md->number_term_params[j] ; k++ ) {
			md->term_multipliers[j] *=  tasks->task_parameter_values[task][md->term_parameters[j][k]] ;
		}
	}

	return 0 ;
}

#undef __FUNCT__
#define __FUNCT__ "deallocate_task_space"
/*! \brief Deallocate the space used to store task values. */
int deallocate_task_space(struct task_parameters *tasks,struct matrix_meta_data *md ) {
	int task;
	PetscErrorCode ierr;
	
	for ( task = 0 ; task < tasks->number_tasks ; task++ ) {
		if (md->number_parameters > 0 ) {
			ierr = PetscFree(tasks->task_parameter_values[task]);CHKERRQ(ierr);
		}
	}

	if ( tasks->number_tasks > 0 ) {
		ierr = PetscFree(tasks->task_parameter_values);CHKERRQ(ierr);
	}
	
	return 0 ;
}
