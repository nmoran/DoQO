#! /usr/bin/env python

from xml.sax import make_parser
from xml.sax import saxutils
from xml.sax import ContentHandler 
from xml.sax.handler import feature_namespaces
import sys, getopt, string
from numpy import *


class ostream:
    def __init__(self, file):
        self.file = file
        
    def __lshift__(self, obj):
        self.file.write(str(obj));
        return self

def normalize_whitespace(text):
    "Remove redundant whitespace from a string"
    return ' '.join(text.split())

 
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
        self.files = list()
        self.momentum =  zeros(2)
        self.error = float(0)
        self.energies = list()
        self.errors = list()

    def startElement(self, name, attrs):
        
        if name == 'MOMENTUM':
            #print "in output"
            self.inMomentum = 1;
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
            self.energies.append(float(self.buffer))
            self.errors.append(self.error);
            
     

if __name__ == '__main__':
    # Create a parser
    parser = make_parser()
    in_file = sys.argv[1]
    task_to_process = -1
    try:
        task_to_process = int(sys.argv[2])
    except IndexError:
        pass
#    print 'File to parse is ', in_file 
    parser.setFeature(feature_namespaces, 0)
    dh = ExtractOutputInfo()
    parser.setContentHandler(dh)
    parser.parse(in_file)
 #   print dh.tasks
    all_eigenvalues = list()
    all_errors = list()
    all_momentums = list()
    fidx = 0 
    task = dh.tasks[0]
    while fidx < dh.files.__len__() :
        file = dh.files[fidx]
        
        if dh.tasks[fidx] != task : 
            #simple bubble sort. good for pretty small numbers of eigenstates
            idx = 0 
            while idx < all_eigenvalues.__len__() :
                idx2 = idx + 1
                while idx2 < all_eigenvalues.__len__():
                    if all_eigenvalues[idx2] < all_eigenvalues[idx]:
                        #swap them
                        val = all_eigenvalues[idx2]
                        all_eigenvalues[idx2] = all_eigenvalues[idx]
                        all_eigenvalues[idx] = val
                        #swap errors and momentum also
                        val = all_errors[idx2]
                        all_errors[idx2] = all_errors[idx]
                        all_errors[idx] = val
                        mom = all_momentums[idx2]
                        all_momentums[idx2] = all_momentums[idx]
                        all_momentums[idx] = mom
                    idx2 += 1
                idx+=1    
                
            first_gap = all_eigenvalues[1]-all_eigenvalues[0]
            idx = 1
            while idx < all_eigenvalues.__len__():
                #all_eigenvalues[idx] =  (all_eigenvalues[idx]-all_eigenvalues[0])/first_gap
                all_eigenvalues[idx] =  (all_eigenvalues[idx]-all_eigenvalues[0])
                idx +=1
            all_eigenvalues[0] = 0
                    
            idx = 0 
#            print 'task = %d'%(task)
            if task_to_process == -1 or task_to_process == task : 
                while idx < all_eigenvalues.__len__():
                   print '%lf , %lf , %lf , %lf'%(all_momentums[idx][0],all_momentums[idx][1],all_eigenvalues[idx],all_errors[idx])	
                   #print '%lf , %lf , %lf , %lf'%(all_eigenvalues[idx],all_errors[idx],all_momentums[idx][0],all_momentums[idx][1])
                   idx+=1
                print ''
                print ''
                
            all_eigenvalues = list()
            all_errors = list()
            all_momentums = list()
            task = dh.tasks[fidx]
        
        #print file
        indiv_parser = make_parser()
        dh2 = ExtractSpectrumInfo()
        indiv_parser.setContentHandler(dh2)
        indiv_parser.parse(file)
        #print dh2.momentum
        #print dh2.energies
        #print dh2.errors
        idx = 0
        while idx < dh2.energies.__len__() :
            all_eigenvalues.append(dh2.energies[idx])
            all_errors.append(dh2.errors[idx])
            all_momentums.append(dh2.momentum)
            idx += 1

        fidx +=1
        
    idx = 0 
    while idx < all_eigenvalues.__len__() :
        idx2 = idx + 1
        while idx2 < all_eigenvalues.__len__():
            if all_eigenvalues[idx2] < all_eigenvalues[idx]:
                #swap them
                val = all_eigenvalues[idx2]
                all_eigenvalues[idx2] = all_eigenvalues[idx]
                all_eigenvalues[idx] = val
                #swap errors and momentum also
                val = all_errors[idx2]
                all_errors[idx2] = all_errors[idx]
                all_errors[idx] = val
                mom = all_momentums[idx2]
                all_momentums[idx2] = all_momentums[idx]
                all_momentums[idx] = mom
            idx2 += 1
        idx+=1    
        
    first_gap = all_eigenvalues[1]-all_eigenvalues[0]
    idx = 1
    while idx < all_eigenvalues.__len__():
        #all_eigenvalues[idx] =  (all_eigenvalues[idx]-all_eigenvalues[0])/first_gap
        all_eigenvalues[idx] =  (all_eigenvalues[idx]-all_eigenvalues[0])
        idx +=1
    all_eigenvalues[0] = 0
    
 #   print 'task = %d'%(task)
    if task_to_process == -1 or task_to_process == task : 
        idx = 0 
        while idx < all_eigenvalues.__len__():
            print '%lf , %lf , %lf , %lf'%(all_momentums[idx][0],all_momentums[idx][1],all_eigenvalues[idx],all_errors[idx])
            idx+=1
