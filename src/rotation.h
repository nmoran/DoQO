/*
 * Copyright (C) 2008    Niall Moran, Graham Kells, Jiri Vala
 *
 * Header file for rotation functions  
 * 
 * Version: $Id$ 
 * 
*/



#ifndef _ROTATION_H_
#define _ROTATION_H_

#include "doqo.h"
#include "hamiltonian.h"
#include "conservation.h"
#include "momentum.h"
#include "tasks.h"
#include "utils.h"

int parse_rotation_info(struct parameter_struct *parameters,TiXmlHandle *docHandle);
int read_rotation_information(struct parameter_struct *parameters);
uPetscInt apply_rotation(struct parameter_struct *parameters, uPetscInt basis,int *number, int rot_op_idx);
uPetscInt apply_rotation(struct parameter_struct *parameters, uPetscInt basis,int *number, int rot_op_idx, PetscScalar *phase);
PetscScalar get_rotation_phase(struct parameter_struct *parameters, PetscInt sector);
void calculate_rotation_phases(struct parameter_struct *parameters, int rot_op_idx);
uPetscInt  rotation_group_representative(struct parameter_struct *parameters, uPetscInt basis, PetscTruth phase_check, PetscScalar *phase, PetscTruth *group_exists);
PetscReal  rotation_group_normal(struct parameter_struct *parameters, uPetscInt index);
PetscTruth rotation_representative_check(struct parameter_struct *parameters, uPetscInt index);
void set_rotation_sector(struct parameter_struct *parameters, PetscInt sector);
void deallocate_rotation_memory(struct parameter_struct *parameters);
#endif

