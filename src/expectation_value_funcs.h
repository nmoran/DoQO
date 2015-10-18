#ifndef _EXPECTATION_VALUE_FUNCS_H_
#define _EXPECTATION_VALUE_FUNCS_H_

#include "doqo.h"

int read_expectation_values_file(struct parameter_struct *parameters);
int calculate_expectation_values(struct parameter_struct* parameters, Mat *A);
int deallocate_expectaction_value_space(struct parameter_struct *parameters);

#endif