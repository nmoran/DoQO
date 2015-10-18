#!/usr/bin/python

import numpy as np
import sys
import pdb
import matplotlib.pyplot as p
import pickle

def read_matlab_matrix(filename):
	try:
		f = open(filename,"r")
	except IOError:
		print "Error opening file %s."%(filename)
		exit()
	
	line = f.readline()
	rows = 0 
	cols = 0
	while line != "" : 
		pos = line.find("zeros")
		if pos != -1 :
			rows =  int(line.split("(")[1].split(")")[0].split(",")[0])
			cols =  int(line.split("(")[1].split(")")[0].split(",")[1])
			break
		line = f.readline()
	
	#create empty matrix 
	mat = np.matrix(np.zeros([rows,cols]))	
	
	#now read the matrix entries 
	#first line is just assignment so ignore
	line = f.readline() 
	
	row = 0
	while row < rows : 
		line = f.readline() # this reads in a row of values
		vals = line.split(" ") #splits by space
		col = 0
		while col < cols :
			mat[row,col] = float(vals[col])
			col += 1
		row += 1 
	
	return mat
	

def calculate_ee_from_rho_a(mat):
	#where rho a is the reduced density matrix
	#first we diagonalise
	(d,v) = np.linalg.eig(mat)
	
	v = np.matrix(v,dtype=np.complex128)

	#compare v^* * d * v to mat
	diff = mat - (v * np.matrix(np.diag(d)) * v.transpose().conjugate())
	
	#print "Max is %lf\n"%(diff.max()) 
	#print "Min is %lf\n"%(diff.min())

	d_new = np.ndarray([d.__len__()],dtype=np.complex128)
	
	#then replace and zeros in the diagonal with minimum float
	idx = 0
	while idx < d.__len__() :
		d_new[idx] = d[idx]
		idx += 1 
		
	idx = 0
	while idx < d_new.__len__() :
		if abs(d_new[idx]) == 0.0:
			d_new[idx] = complex(1e-200,0.0)
		idx += 1 
		
	log_d  = 	np.log(d_new)
	log_mat = v *  np.matrix(np.diag(log_d),dtype=np.complex128) * v.transpose().conjugate()   
	
	entropy = np.trace(-mat*log_mat)
		
	return log_d, log_mat , d_new , v ,entropy , diff
	
	
def calculate_es_from_W(mat):
	u,s,v = np.linalg.svd(mat)

	#normalise
	tot = s.sum
	
	es = np.zeros([s.__len__()])
	idx = 0 
	while idx < s.__len__()  :
		es[idx] = -np.log(s[idx]) * 2.0
		idx += 1
	
	es_normed = np.zeros([s.__len__()])
	tot = 0.0
	for val in es :
		tot += np.exp(-val)
		
	log_tot = np.log(tot)
	#print log_tot
	
	es_normed = np.zeros([s.__len__()])
	idx = 0 
	while idx < es_normed.__len__()  :
		es_normed[idx] = es[idx] + log_tot
		idx += 1
		
	tot = 0
	for val in es_normed :
		tot += np.exp(-val)
		
	#print tot 
		
	entropy = 0.0
	for val in es :
		entropy += val * np.exp(-val)
		
	entropy_normed = 0.0
	for val in es_normed :
		entropy_normed += val * np.exp(-val)
		
	return es, es_normed, entropy, entropy_normed


use_rho = False
use_W = False
if sys.argv.__len__() >= 2 :
	rho_a_file = sys.argv[1]
	use_rho = True
	if sys.argv.__len__() >= 3 :
		W_file = sys.argv[2]
		use_W = True
else :
	print "Two files must be specifiied.\n"
	quit()

if use_rho : 
	rho_a = read_matlab_matrix(rho_a_file)
	logd ,log_mat, d , v , entropy, diff = calculate_ee_from_rho_a(rho_a)
	idx = 0 
	while idx < logd.__len__() :
		print '%.4e		%.4f'%(d[idx],-logd[idx])
		idx+=1

	print "Entropy using regular method: %f\n"%(entropy)

#W = np.zeros([2,2])
#W[0,0] = 1/np.sqrt(2)
#W[1,1] = 1/np.sqrt(2)

if use_W : 
	W = read_matlab_matrix(W_file)
	
	es, es_normed, en2, en2_normed = calculate_es_from_W(W)
	
	print "Entropy using newer method: %f (normed %f)\n"%(en2,en2_normed)
	es_file = "%s"%(W_file + '.es')
	print "Writing entanglement spectrum to file: %s\n"%(es_file)
	
	f = open(es_file,'wb')
	pickle.dump(es,f)
	f.close()

	esn_file = "%s"%(W_file + '.esn')
	print "Writing entanglement spectrum (normed) to file: %s\n"%(esn_file)

	f = open(esn_file,'wb')
	pickle.dump(es_normed,f)
	f.close()
