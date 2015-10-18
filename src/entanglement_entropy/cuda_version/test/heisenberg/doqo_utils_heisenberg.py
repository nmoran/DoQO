#!/usr/bin/python

import sys
sys.path.append('/home/niall/local/Solver/branches/fermionic/models/')

from doqo_utils import * 
from numpy import *

if sys.argv.__len__() >= 2 : 
	sites = int(sys.argv[1])
else :
	print '%d args but 1 expected\n'%(sys.argv.__len__())
	quit()


new_sim = DoQO_simulation(sites)
new_sim.conf['model_type'] = "SPIN_HALF"
new_sim.conf['output_prefix'] = 'Heisenberg_N_%d'%(sites)
new_sim.conf['operator_file'] = '%s.ham'%(new_sim.conf['output_prefix'])
new_sim.conf['num_eigenvalues'] = 3
new_sim.conf['tasks_file'] = '%s_tasks'%(new_sim.conf['output_prefix'])
new_sim.conf['task_parameters'] = ['J']
new_sim.conf['tasks_values'] = [[0.5]]
new_sim.write_task_file()
new_sim.conf['filling_symmetry'] = True
new_sim.conf['filling_file'] = '%s.cons'%(new_sim.conf['output_prefix'])
new_sim.conf['filling_details'] = dict()
new_sim.conf['filling_details']['number'] = 1
new_sim.conf['filling_details']['sectors'] = [(sites/2)]
new_sim.write_generic_filling_file()
new_sim.conf['use_disk'] = False
new_sim.conf['save_matrix'] = False
new_sim.conf['save_matrix_format'] = 'ascii'
new_sim.conf['save_states'] = True
new_sim.conf['save_states_format'] = 'ascii'
new_sim.conf['save_states_real'] = True
new_sim.conf['save_basis'] = True
new_sim.conf['benchmark'] = False 
new_sim.conf['verbosity'] = 5
new_sim.conf['rotation_symmetry'] = True
new_sim.conf['rotation_file'] = '%s.rot'%(new_sim.conf['output_prefix'])
new_sim.conf['rotation_details'] = dict()
new_sim.conf['rotation_details']['number'] = 1
new_sim.conf['rotation_details']['sectors'] = [0]
write_chain_translation(sites, '%s.rot'%(new_sim.conf['output_prefix']))
new_sim.write_input_file()

my_unit_cell = Unit_cell(2,1)
my_lattice = Heisenberg_Lattice(sites,my_unit_cell,pbc=True)
my_lattice.add_interactions_and_params()
my_lattice.write_model_file(new_sim.conf['operator_file'])
	
