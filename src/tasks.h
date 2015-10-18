/* 
 * Copyright (C) 2010    Niall Moran, Graham Kells, Jiri Vala
 * 
 * 
*/

#ifndef _TASKS_H_
#define _TASKS_H_


#include "doqo.h" 
#include "hamiltonian.h"
#include "tasks.h"

int parse_task_file(struct task_parameters *tasks, struct parameter_struct *parameters, struct matrix_meta_data *md ) ;
int apply_task_parameters(struct matrix_meta_data *md, struct task_parameters *tasks , int task ) ;
int deallocate_task_space(struct task_parameters *tasks, struct matrix_meta_data *md);
void print_task_info(struct task_parameters *tasks , struct matrix_meta_data *md);

#endif
