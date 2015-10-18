/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 *
 * Header file for utility functions  
 * 
 * Version: $Id$ 
 * 
*/



#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include "doqo.h"
#include "symmetry.h"


//void write_parameter_element(FILE *fp, char *name, char *type, char *value);


int write_xml_data(struct parameter_struct *parameters,struct solver_results *results);
int write_main_xml_output(struct parameter_struct *parameters);

#endif

