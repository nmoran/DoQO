#! /usr/bin/env python

from xml.sax import make_parser
from xml.sax import saxutils
from xml.sax import ContentHandler 
from xml.sax.handler import feature_namespaces
import sys, getopt, string
from numpy import *
import os.path
from optparse import OptionParser
import pdb

def normalize_whitespace(text):
	"Remove redundant whitespace from a string"
	return ' '.join(text.split())


class ostream:
	def __init__(self, file):
		self.file = file
		
	def __lshift__(self, obj):
		self.file.write(str(obj));
		return self

class ExtractOutputInfo(ContentHandler):
	def __init__(self):
		self.inOutput = 0
		self.files = list()
		self.tasks = list()

	def startElement(self, name, attrs):
		if name == 'OUTPUT':
			#print "in output"
			self.inOutput = 1;
			self.task = int(attrs.get('task', None))
			self.tasks.append(self.task)
			self.momentum_sector = attrs.get('momentum_sector', None)
			self.buffer = '';
	
	def characters(self, ch):
		if self.inOutput:
			self.buffer = self.buffer + ch
			

	def endElement(self, name):
		if name == 'OUTPUT' and self.inOutput :
			#print "end output"
			self.inOutput = 0
			self.buffer = normalize_whitespace(self.buffer)
			self.files.append(self.buffer);


class ExtractSpectrumInfo(ContentHandler):
	def __init__(self):
		self.inMomentum = 0
		self.inDirection = 0
		self.inEigenstate = 0
		self.inEigenspace = 0
		self.files = list()
		self.momentum =  zeros(2)
		self.rotation = 0 
		self.error = float(0)
		self.states = list()
		self.correlations = list()
		self.correlation_matrices = list()
		self.correlations_present = False
		self.filling = 0
		self.parity = 0 
		self.spectral = 0

	def startElement(self, name, attrs):
		
		if name == 'MOMENTUM':
			#print "in output"
			self.inMomentum = 1;
		if name == 'ROTATION':
			#print "in output"
			self.rotation = int(attrs.get('sector', None))  
		elif name == 'DIRECTION':
			#print "in output"
			self.inMomentum = 1;
			self.inDirection = int(attrs.get('name', None)) + 1 
			self.buffer = '';
		elif name == 'EIGENSTATE':
			self.inEigenstate = 1
			self.buffer = 0
			self.error = float(attrs.get('error',None))
			self.buffer = '';
			if  attrs.get('overlap',None) :
				self.overlap = float(attrs.get('overlap',None))
			else: self.overlap = 0.0
			if  attrs.get('other_overlap',None) :
				self.other_overlap = float(attrs.get('other_overlap',None))
			else: self.other_overlap = 0.0
		elif name == 'EIGENSPACE' :
			self.inEigenspace = 1;
			self.space_correlations = list()
			self.cmatrices = list()
		elif name == 'CORRELATION' :
			self.correlations_present = True
			if self.inEigenspace == 1 :
				self.in_correlation = 1
				self.buffer = ''
				self.dim = int(attrs.get('dimension',None))
				self.corr_number = int(attrs.get('number',None))
		elif name == 'FILLING':
			self.filling = int(attrs.get('sector', None)) 
		elif name == 'PARITY':
			self.parity = int(attrs.get('sector', None)) 
		elif name == 'SPECTRAL_FLOW':
			self.spectral = int(attrs.get('sector', None)) 	
				
	
	def characters(self, ch):
		if self.inDirection > 0 or self.inEigenstate :
			self.buffer = self.buffer + ch
			

	def endElement(self, name):
		if name == 'DIRECTION'  :
			#print "end output"
			self.buffer = normalize_whitespace(self.buffer)
			self.momentum[self.inDirection-1] = float(self.buffer) 
			self.inDirection =0
		elif name == 'EIGENSTATE':
			self.buffer = normalize_whitespace(self.buffer)
			tmp = EigenState()
			tmp.energy = float(self.buffer)
			tmp.error = self.error
			tmp.momentum[0] = self.momentum[0];tmp.momentum[1] = self.momentum[1]
			tmp.rotation = self.rotation
			tmp.occupation = self.filling
			tmp.parity = self.parity
			tmp.spectral = self.spectral
			tmp.overlap = self.overlap
			tmp.other_overlap = self.other_overlap
			self.states.append(tmp)
		elif name == 'EIGENSPACE':
			self.inEigenspace = 0 
			self.correlations.append(self.space_correlations)
			self.correlation_matrices.append(self.cmatrices)
			#print 'All correlations:' ,self.correlations
		elif name == 'CORRELATION' :
			if self.inEigenspace == 1 :
				self.buffer = normalize_whitespace(self.buffer)
				mat = zeros([self.dim,self.dim])
				rows = self.buffer.split(';')
				row = 0
				while row < self.dim :
					cols = rows[row].split(',')
					col = 0 
					while col < self.dim : 
						mat[row,col] = float(cols[col].split()[0])
						col +=1
					row += 1
				(d,v) = linalg.eig(mat)
				#print mat
				corr = list()
				d.sort()
				for e in d: 
					corr.append(e)
				#print 'corr list:', corr
				self.space_correlations.append(corr)
				self.cmatrices.append(mat)
				#print self.state_correlations
				#print self.correlations
			
	def finish(self) :
		eigenstate_idx = 0
		#print self.correlations
		if self.correlations_present :
			for eigspace in self.correlations :
				dim = eigspace[0].__len__()  #the dimension of the eigenspace
				#print dim
				eig_idx = 0 
				while eig_idx < dim : 
					corrs = list()  #create a new list for the eigenstate
					for correlations in eigspace :
						corrs.append(correlations[eig_idx])
					eig_idx += 1  
					self.states[eigenstate_idx].correlations.extend(corrs)
					self.states[eigenstate_idx].correlation_matrices.append(self.correlation_matrices)
					eigenstate_idx += 1 
			
			
			
def eigenstate_compare(x,y):
	if  (x.energy - y.energy) < -1e-10 : 
		return -1 
	elif (x.energy - y.energy ) > 1e-10 :
		return 1 
	else :
		return 0 


def eigenstate_spectral_compare(x,y):
	if  (x.spectral - y.spectral) < -1e-10 : 
		return -1 
	elif (x.spectral - y.spectral ) > 1e-10 :
		return 1 
	else :
		return 0 
	
	
class EigenState : 
	def __init__(self):
		self.energy = 0.0
		self.error = 0.0
		self.momentum = zeros(2)
		self.rotation = 0
		self.occupation = 0
		self.correlations = list()
		self.correlation_matrices = list()
		self.parity = 0 
		self.spectral = 0 	
		self.overlap = 0
		self.other_overlap = 0
		
	def print_state(self):
		print '%.20lf , %lf , %lf , %lf, %d, %d, %d, %d, %le, %le'%(self.energy,self.error,self.momentum[0],self.momentum[1],self.occupation,self.parity, self.spectral,self.rotation, self.overlap, self.other_overlap)
		#print self.correlations
	

def contains(theString, theQueryValue):
	return theString.find(theQueryValue) > -1


if __name__ == '__main__':
	# Create a parser

	optparser = OptionParser()
	optparser.add_option("-p", "--pattern", dest="pattern",
                  help="Pattern to check for in the output filenames.", default="")
	optparser.add_option("-f", "--filling", dest="filling",
                  help="Restrict to specified filling.", default=-1)

	(optlist, args) = optparser.parse_args()

	
	#print args
	#print optlist
	#print optlist.pattern
	
	if args.__len__() >= 1 :
		in_file = args[0] # this is the output file which lists other output files for each task and sector.
		if os.path.isfile(in_file) == False :
			print "%s does not exists."%(in_file)
			quit()
	else : 
		print "Usage: \n extract_xml_results.py <output filename>  [<gap/dispersion/correlations/spectral>]"
		quit()
	
	
	style = 'plain'
	if args.__len__() >= 2 :
		if args[1] == 'gap' :
			style = 'gap'
		elif args[1] == 'dispersion' :
			style = 'dispersion'
		elif args[1] == 'correlations' :
			style = 'correlations'
		elif args[1] == 'spectral' :
			style = 'spectral'
		elif args[1] == 'plain' :
			style = 'plain'
		elif args[1] == 'zero_degs' :
			style = 'zero_degs'
		else :
			print "Output style '%s' not available."%(args[1])
			quit()
	
	chosen_task = -1
	if args.__len__() >= 3 :
		chosen_task = int(args[2])


	parser = make_parser()
	parser.setFeature(feature_namespaces, 0)
	dh = ExtractOutputInfo()
	parser.setContentHandler(dh)
	parser.parse(in_file) 
	#print dh.tasks
	
	all_eigenstates = dict()
	fidx = 0 
	eig_task_idx = 0 
	task = dh.tasks[0]
	all_eigenstates[task] = list()
	while fidx < dh.files.__len__() :
		file = dh.files[fidx]
		if contains(file,optlist.pattern):
			if dh.tasks[fidx] != task : 
				task = dh.tasks[fidx]
				all_eigenstates[task] = list()
			
			indiv_parser = make_parser()
			dh2 = ExtractSpectrumInfo()
			indiv_parser.setContentHandler(dh2)
			indiv_parser.parse(file)
			dh2.finish()
		
			for state in dh2.states : 
				if int(optlist.filling) == -1 or state.occupation == int(optlist.filling) : 
					all_eigenstates[task].append(state)
		fidx += 1 
		
	#now we will output depending on the style specified
	if style == 'plain' : 		
		task_idx = 0 
		while task_idx <= task :
			if chosen_task == -1 or task_idx == chosen_task : 
				all_eigenstates[task_idx].sort(eigenstate_compare)
				for eig_state in all_eigenstates[task_idx] :
					eig_state.print_state()
				print ""
			task_idx += 1 
	elif style == 'gap' : 
		task_idx = 0 
		while task_idx <= task :
			if chosen_task == -1 or task_idx == chosen_task : 
				if all_eigenstates[task_idx].__len__() > 0 :
					all_eigenstates[task_idx].sort(eigenstate_compare)
					gs_energy = all_eigenstates[task_idx][0].energy
					first_energy = gs_energy
					eidx = 1
					while eidx < all_eigenstates[task_idx].__len__() :
						if abs(all_eigenstates[task_idx][eidx].energy - gs_energy) > 1.0e-5 :
							first_energy = all_eigenstates[task_idx][eidx].energy
							eidx = all_eigenstates[task_idx].__len__()
						eidx += 1 
					print "%d   %lf"%(task_idx,first_energy-gs_energy) 
			task_idx += 1
	elif style == 'dispersion' :
		task_idx = 0 
		while task_idx <= task :
			all_eigenstates[task_idx].sort(eigenstate_compare)
			gs_energy = all_eigenstates[task_idx][0].energy
			
			eidx = 0
			while eidx < all_eigenstates[task_idx].__len__() :
				print "%lf   %lf"%(all_eigenstates[task_idx][eidx].momentum[0],all_eigenstates[task_idx][eidx].energy-gs_energy) 	
				eidx += 1 
			print ""
			task_idx += 1
	elif style == 'correlations' :
		task_idx = 0 
		for task_idx in all_eigenstates.keys() :		
			if chosen_task == -1 or task_idx == chosen_task :
				#print 'Tasks %d\n-------'%(task_idx)
				eidx = 0
				while eidx < all_eigenstates[task_idx].__len__() :
					#print 'State %d\n-------'%(eidx)
					corr_idx = 1
					for corr in all_eigenstates[task_idx][eidx].correlations : 
						print '%d %d %d %.9e'%(task_idx,eidx,corr_idx,corr) 
						#print all_eigenstates[task_idx][eidx].correlation_matrices
						corr_idx += 1 
					eidx += 1 
			task_idx += 1
			
	elif style == 'spectral' :
		print "@default linestyle 1"
		print "@default symbol size 0.41"
		task_idx = 0 
		while task_idx <= task :
			all_eigenstates[task_idx].sort(eigenstate_spectral_compare)
			
			eidx = 0
			spectral = all_eigenstates[task_idx][eidx].spectral
			print '@	s%d symbol 1'%(all_eigenstates[task_idx][eidx].spectral)
			print '@	s%d line type 0'%(all_eigenstates[task_idx][eidx].spectral)
			print '@    s%d symbol linewidth 1.5'%(all_eigenstates[task_idx][eidx].spectral)
			while eidx < all_eigenstates[task_idx].__len__() :
				if spectral != all_eigenstates[task_idx][eidx].spectral :  
					print ""
					print '@	s%d symbol 1'%(all_eigenstates[task_idx][eidx].spectral)
					print '@	s%d line type 0'%(all_eigenstates[task_idx][eidx].spectral)
					print '@    s%d symbol linewidth 1.5'%(all_eigenstates[task_idx][eidx].spectral)
					spectral = all_eigenstates[task_idx][eidx].spectral
				print "%lf   %lf"%(all_eigenstates[task_idx][eidx].momentum[0],all_eigenstates[task_idx][eidx].energy) 	
				eidx += 1 
			print ""
			task_idx += 1
	elif style == 'zero_degs' :
		task_idx = 0
		while task_idx <= task :
			all_eigenstates[task_idx].sort(eigenstate_spectral_compare)
			#find max filling
			max_filling = 0 
			for eigenstate in all_eigenstates[task_idx] :
				if eigenstate.occupation > max_filling:
					max_filling = eigenstate.occupation
			#create array of zeros of this size
			degs = zeros(max_filling+1)
			for eigenstate in all_eigenstates[task_idx] :
				if abs(eigenstate.energy) < 1.0e-10 :
					degs[eigenstate.occupation] += 1
			#print out the fillings
			for deg in degs :
				print "%d, "%(deg),
			task_idx += 1

						
			
			
 

	
		
	
