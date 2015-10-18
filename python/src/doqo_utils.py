from numpy import *
import itertools
import re


# Class to manage a DoQO simulation.
class DoQO_simulation:
    def __init__(self,sites):
        self.conf = dict()
        self.conf['output_prefix'] = "doqo_output"
        self.conf['operator_file'] = ""
        self.conf['num_eigenvalues'] = 4
        self.conf['verbosity'] = 5
        self.sites = sites
        
    def write_input_file(self):
        filename = "%s.input.xml" % (self.conf['output_prefix'])
        fp = open(filename,"w")
        fp.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fp.write('<SIMULATION>\n')
        fp.write('\t<PARAMETERS>\n')
        if self.conf.has_key('model_type') : fp.write('\t\t<MODEL_TYPE>%s</MODEL_TYPE>\n'%(self.conf['model_type']))
        fp.write('\t\t<MODEL_FILE>%s</MODEL_FILE>\n'%(self.conf['operator_file']))
        fp.write('\t\t<TASK_LIST>%s</TASK_LIST>\n'%(self.conf['tasks_file']))
        fp.write('\t\t<OUTPUT_PREFIX>%s</OUTPUT_PREFIX>\n'%(self.conf['output_prefix']))
        fp.write('\t\t<EIGENVALUES>%d</EIGENVALUES>\n'%(self.conf['num_eigenvalues']))
        if self.conf.has_key('save_states') and self.conf['save_states'] : 
            if self.conf.has_key('save_states_real') and self.conf['save_states_real'] : 
                fp.write('\t\t<SAVE_STATES format="%s" real="true">true</SAVE_STATES>\n'%(self.conf['save_states_format']))
            else:
                fp.write('\t\t<SAVE_STATES format="%s">true</SAVE_STATES>\n'%(self.conf['save_states_format']))
        if self.conf.has_key('save_matrix') and self.conf['save_matrix'] : 
            fp.write('\t\t<SAVE_MATRIX format="%s" >true</SAVE_MATRIX>\n'%(self.conf['save_matrix_format']))
        if self.conf.has_key('use_disk') and self.conf['use_disk'] : fp.write('\t\t<USE_DISK>true</USE_DISK>\n')
        if self.conf.has_key('save_basis') and self.conf['save_basis'] : fp.write('\t\t<SAVE_BASIS>true</SAVE_BASIS>\n')
        if self.conf.has_key('use_bst') and self.conf['use_bst'] : fp.write('\t\t<USE_BST>true</USE_BST>\n')
        if self.conf.has_key('nn_exclusion') and self.conf['nn_exclusion'] : 
            fp.write('\t\t<NN_EXCLUSION recursive="%s" adjacency_file="%s">true</NN_EXCLUSION>\n'%(self.conf['nn_exclusion_recursive'],self.conf['adj_filename']))
        if self.conf.has_key('prod_wf') and self.conf['prod_wf'] and self.conf.has_key('prod_wf_file') :
            fp.write('\t\t<PROD_WF_OVERLAP file="%s" />\n'%(self.conf['prod_wf_file']))
        fp.write('\t\t<SOLVER_TYPE>krylovschur</SOLVER_TYPE>\n')
        fp.write('\t\t<MAX_ITERATIONS>10000</MAX_ITERATIONS>\n')
        fp.write('\t\t<VERBOSITY>%d</VERBOSITY>\n'%(self.conf['verbosity']))
        if self.conf.has_key('expectation_values_file') :
            fp.write('\t\t<EXPECTATION_VALUES file="%s" '%(self.conf['expectation_values_file']))
            if self.conf.has_key('expectation_values_overlap_vec') :
                fp.write(' overlap_vec="%s" '%(self.conf['expectation_values_overlap_vec']))
            fp.write('/>\n')
        if self.conf.has_key('benchmark') and self.conf['benchmark'] :
            fp.write('\t\t<BENCHMARK>true</BENCHMARK>\n')
        fp.write('\t\t<SYMMETRIES>\n')
        if self.conf.has_key('filling_symmetry') and self.conf['filling_symmetry'] : 
            fp.write('\t\t\t<FILLING file="%s">\n'%(self.conf['filling_file']))
            if self.conf.has_key('filling_details') :    
                fp.write('\t\t\t\t<RELEVANT_SECTORS number="%d">\n'%(self.conf['filling_details']['number']))
                sector = 0
                while sector < self.conf['filling_details']['number'] :
                    fp.write('\t\t\t\t\t<SECTOR>%d</SECTOR>\n'%(self.conf['filling_details']['sectors'][sector]))
                    sector +=1
                fp.write('\t\t\t\t</RELEVANT_SECTORS>\n')
            fp.write('\t\t\t</FILLING>\n')
        if self.conf.has_key('parity_symmetry') and self.conf['parity_symmetry'] : 
            fp.write('\t\t\t<PARITY file="%s">\n'%(self.conf['parity_file']))
            if self.conf.has_key('parity_details') :    
                fp.write('\t\t\t\t<RELEVANT_SECTORS number="%d">\n'%(self.conf['parity_details']['number']))
                sector = 0
                while sector < self.conf['parity_details']['number'] :
                    fp.write('\t\t\t\t\t<SECTOR>%d</SECTOR>\n'%(self.conf['parity_details']['sectors'][sector]))
                    sector +=1
                fp.write('\t\t\t\t</RELEVANT_SECTORS>\n')
            fp.write('\t\t\t</PARITY>\n')
        if self.conf.has_key('translation_symmetry') and self.conf['translation_symmetry'] : 
            fp.write('\t\t\t<MOMENTUM file="%s">\n'%(self.conf['translation_file']))
            fp.write('\t\t\t</MOMENTUM>\n')
        if self.conf.has_key('rotation_symmetry') and self.conf['rotation_symmetry'] : 
            fp.write('\t\t\t<ROTATION file="%s">\n'%(self.conf['rotation_file']))
            if self.conf.has_key('rotation_details') :    
                fp.write('\t\t\t\t<RELEVANT_SECTORS number="%d">\n'%(self.conf['rotation_details']['number']))
                sector = 0
                while sector < self.conf['rotation_details']['number'] :
                    fp.write('\t\t\t\t\t<SECTOR>%d</SECTOR>\n'%(self.conf['rotation_details']['sectors'][sector]))
                    sector +=1
                fp.write('\t\t\t\t</RELEVANT_SECTORS>\n')
            fp.write('\t\t\t</ROTATION>\n')
        fp.write('\t\t</SYMMETRIES>\n')
        if self.conf.has_key('calculate_correlations') and self.conf['calculate_correlations'] : 
            correlations = self.conf['correlations']
            fp.write('\t\t<CALCULATE_CORRELATIONS number="%d">\n'%(correlations.__len__()))
            corr_idx = 0
            while corr_idx < correlations.__len__() : 
                fp.write('\t\t\t<CORRELATOR number="%d">\n'%(correlations[corr_idx].__len__()))
                corr_site = 0
                while corr_site < correlations[corr_idx].__len__() :
                    fp.write('\t\t\t\t<SITE>%d</SITE>\n'%(correlations[corr_idx][corr_site]))
                    corr_site += 1 
                fp.write('\t\t\t</CORRELATOR>\n')
                corr_idx += 1
            fp.write('\t\t</CALCULATE_CORRELATIONS>\n')
        fp.write('\t</PARAMETERS>\n')
        fp.write('</SIMULATION>\n')
        fp.close()                                        
        
    def write_generic_filling_file(self):
        fp = open(self.conf['filling_file'],"w")
        i = 0
        while i < self.sites : 
            fp.write('1')
            i+= 1
        fp.write('\n')
        fp.close()

    def write_generic_parity_file(self):
        fp = open(self.conf['parity_file'],"w")
        i = 0
        while i < self.sites : 
            fp.write('1')
            i+= 1
        fp.write('\n')
        fp.close()
        
    def write_task_file(self):
        fp = open(self.conf['tasks_file'],"w")
        task_params = self.conf['task_parameters']
        tasks_values = self.conf['tasks_values']
        for vals in tasks_values :
            idx = 0 
            while idx < task_params.__len__() : 
                fp.write('%s = %.20lf'%(task_params[idx],vals[idx]))
                if idx < (task_params.__len__() - 1) :
                    fp.write(',')
                idx+=1
            fp.write('\n')
        fp.close()

    #vec is the name of the state file we are writing the output file for. 
    #filling is the filling
    #rotation is the rotation sector
    def write_expectation_value_file(self, vec, filling, rotation):          
        filename = re.sub('.vec$','.xml', vec)
        fp = open(filename,"w")
        fp.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fp.write('<SIMULATION>\n')
        fp.write('\t<PARAMETERS>\n')
        fp.write('\t\t<FILLING sector="%d" />\n'%(filling))
        fp.write('\t\t<ROTATION sector="%d" />\n'%(rotation))
        fp.write('\t\t<EIGEN_PAIRS number="1" spaces="1">\n')
        fp.write('\t\t\t<EIGENSPACE number="0" degeneracy="1" >\n')
        fp.write('\t\t\t\t<EIGENSTATE number="0" state_vector="%s" />\n'%(vec))
        fp.write('\t\t\t</EIGENSPACE>\n')
        fp.write('\t\t</EIGEN_PAIRS>\n')
        fp.write('\t</PARAMETERS>\n')
        fp.write('</SIMULATION>\n')
        fp.close()        

class Lattice: 
    def __init__(self,sites,unit_cell,pbc=True):
        self.sites = sites
        self.adjacency = zeros([sites,sites])
        self.edges = 0
        self.unit_cell = unit_cell
        self.interactions = list()
        self.pbc = pbc
        self.parameters = list()
        
    def set_adjacency(self,adjacency_mat):
        self.adjacency = adjacency_mat
        i = 0
        self.edges = 0
        while i < self.sites: 
            j=i+1
            while j < self.sites:
                if self.adjacency[i][j] == 'I' :
                    self.edges +=1 
                j+=1
            i+=1
        
    def get_neighbours(self,site):
        neighbours = list()
        i = 0
        while i < self.sites : 
            if self.adjacency[i][site] != 'o' :
                neighbours.append(i)
            i+=1
        return neighbours
        
    def write_adjacency_file(self,filename):
        fp = open(filename,"w")
        fp.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fp.write('<EDGES number="%d">\n'%(self.edges))
        i = 0
        while i < self.sites: 
            j=i+1
            while j < self.sites:
                if self.adjacency[i][j] == 'I' :
                    from_site = i
                    to_site = j
                    fp.write('<EDGE from="%d" to="%d" ></EDGE>\n'%(from_site+1,to_site+1))
                j+=1
            i+=1        
        fp.write('</EDGES>\n');
        fp.close()
        
    def update_parameters(self):
        for inter in self.interactions : 
            for param in inter.params: 
                try:
                    idx = self.parameters.index(param)
                except ValueError:
                    self.parameters.append(param)
                    
    def write_model_file(self,filename):
        fp = open(filename,"w")
        fp.write('SITES %d\n'%(self.sites))
        fp.write('PARAMETERS\n')
        for param in self.parameters: 
            fp.write('%s\n'%(param))
        fp.write('TERMS\n')
        for inter in self.interactions: 
            fp.write('%s\n'%(inter.interaction_string()))
        fp.write('\n')
        fp.close()            
    
class SUSY_Lattice(Lattice):
    def __init__(self,sites,unit_cell,pbc=True,weights=None):
        Lattice.__init__(self,sites,unit_cell,pbc=pbc)
        if weights != None :
            self.weights = weights
        else:  
            self.weights = list()
            for i in arange(0, sites):
                  self.weights.append('I')
  
    def add_adjacency_for_square_octagon_lattice(self, rows, cols, rowbc, colbc):
        amat = list()
        i = 0 
        while i < self.sites :
            j = 0
            foo = list()
            while j < self.sites : 
                foo.append('o')
                j+=1
            amat.append(foo)
            i+=1
        
        i = 0 
        while i < rows :
            j = 0 
            while j < cols :
                base = i*cols*4 + j*4
                #these are the four links on each square
                amat[base][base+1] = 'I';amat[base+1][base] = 'I' ;
                amat[base][base+2] = 'I';amat[base+2][base] = 'I' ;
                amat[base+1][base+3] = 'I';amat[base+3][base+1] = 'I' ;
                amat[base+2][base+3] = 'I';amat[base+3][base+2] = 'I' ;
                if i < (rows - 1):
                    up_base = (i+1)*cols*4 + j*4
                    amat[base+1][up_base+2] = 'I'; amat[up_base+2][base+1] = 'I';  
                elif colbc == 'periodic':
                    up_base = ((i+1)%rows)*cols*4 + j*4
                    amat[base+1][up_base+2] = 'I'; amat[up_base+2][base+1] = 'I';
                elif colbc == 'moebius':
                    up_base = ((i+1)%rows)*cols*4 + (cols-j-1)*4
                    amat[base+1][up_base+2] = 'I'; amat[up_base+2][base+1] = 'I';
                if j < (cols - 1):
                    right_base = i*cols*4 + (j+1)*4
                    amat[base+3][right_base] = 'I'; amat[right_base][base+3] = 'I';
                elif rowbc == 'periodic':
                    right_base = i*cols*4 + ((j+1)%cols)*4
                    amat[base+3][right_base] = 'I'; amat[right_base][base+3] = 'I';
                elif rowbc == 'moebius':
                    right_base = (rows-i-1)*cols*4 + ((j+1)%cols)*4
                    amat[base+3][right_base] = 'I'; amat[right_base][base+3] = 'I';
                j+=1
            i+=1
        self.set_adjacency(amat)
        
        
    def add_adjacency_for_triangular_lattice(self,rows,cols):
        if self.sites != rows * cols : 
            print "Sites does not match product of rows and cols : %d * %d != %d\n"%(rows,cols,self.sites)
        amat = list()
        i = 0 
        while i < self.sites :
            j = 0
            foo = list()
            while j < self.sites : 
                foo.append('o')
                j+=1
            amat.append(foo)
            i+=1
        
        i = 0 
        while i < rows :
            j = 0 
            while j < cols :
                lower_left = i * cols + j 
                lower_right = i * cols + ((j+1)%cols)
                upper_left =  ((i+1)%rows) * cols + j
                upper_right = ((i+1)%rows) * cols + ((j+1)%cols)
            
                #the lower left of each square has three links
                amat[lower_left][lower_right] = 'I';amat[lower_right][lower_left] = 'I' ;
                amat[lower_left][upper_right] = 'I';amat[upper_right][lower_left] = 'I' ;
                amat[lower_left][upper_left] = 'I' ;amat[upper_left][lower_left] = 'I'  ;
                j+=1
            i+=1
        self.set_adjacency(amat)    
        
        
    def add_adjacency_for_open_chain_lattice(self):
        length = self.sites
        amat = list()
        i = 0 
        while i < self.sites :
            j = 0
            foo = list()
            while j < self.sites : 
                foo.append('o')
                j+=1
            amat.append(foo)
            i+=1
        
        i = 0 
        while i < self.sites :
            if i < (self.sites - 1) :
                left = i
                right = i+1
                amat[left][right] = 'I'; amat[right][left] = 'I';            
            i+=1
        self.set_adjacency(amat)    
        
    def add_defect(self,source,dest,symbol):
        amat = self.adjacency
        amat[source-1][dest-1] = symbol; amat[dest-1][source-1] = symbol
        self.set_adjacency(amat)
    
    
    def add_interactions(self):
        site = 1 
        while site <= self.sites :
            neighbours = self.get_neighbours(site-1) #get list of neighbours of 'site'
            
            #first do the kinetic hopping terms from the given site to each of its neighbours
            for neighbour in neighbours : 
                #get neighbouring sites of each neighbour
                neighbours2 = self.get_neighbours(neighbour)
                
                local_iters = list() 
                int_sites = [neighbour+1,site]
                int_types = ['C','A']
                int_params = ['kin_const']
                if self.weights[site-1] != 'I' or self.weights[neighbour] != 'I':
                      if self.weights[site-1] == self.weights[neighbour] :
                          tmpstr = '%s2'%(self.weights[site-1])
                          int_params.append(tmpstr)
                      else:
                          if self.weights[site-1] != 'I':
                              int_params.append(self.weights[site-1])
                          if self.weights[neighbour] != 'I':
                              int_params.append(self.weights[neighbour])
                              
                inter = Interaction(int_sites,int_types,int_params)
                local_iters.append(inter)
                
                link_strength = self.adjacency[site-1][neighbour]
                if link_strength != 'I' : 
                    # for iter in local_iters:
#                         iter.params.append(link_strength)
                    new_iters = list()
                    for iter in local_iters :
                        new_iter = Interaction(iter)
                        new_iters.append(new_iter)
                    for iter in local_iters :                        
                        iter.params.append(link_strength)
                        iter.params.append('TWO')
                    for iter in new_iters :
                        iter.params.append('%s2'%(link_strength))
                        #iter.params.append('%s'%(link_strength))
                        iter.params.append('minus')
                    local_iters.extend(new_iters)
                    
                for a in neighbours :
                    if a != neighbour :
                        link_strength = self.adjacency[site-1][a]
                        if link_strength != 'I' :
                            new_iters = list()
                            for iter in local_iters :
                                new_iter = Interaction(iter)
                                new_iter.params.append(link_strength)
                                new_iter.params.append('minus')
                                new_iter.sites.extend([a+1,a+1])
                                new_iter.types.extend(['C','A'])
                                new_iters.append(new_iter)    
                            local_iters.extend(new_iters)
                            
                for a in neighbours2 :
                    if a != site-1 :
                        link_strength = self.adjacency[neighbour][a]
                        if link_strength != 'I' :
                            new_iters = list()
                            for iter in local_iters :
                                new_iter = Interaction(iter)
                                new_iter.params.append(link_strength)
                                new_iter.params.append('minus')
                                new_iter.sites.extend([a+1,a+1])
                                new_iter.types.extend(['C','A'])
                                new_iters.append(new_iter)                
                            local_iters.extend(new_iters)
                self.interactions.extend(local_iters)
                    
            #next do the potential terms for the given site. 
            num = 0 
            while num <= neighbours.__len__() : 
                combs = itertools.combinations(neighbours,num) #find all combinations of neighbours of length num
                for comb in combs : 
                    int_sites = list()
                    int_types = list()
                    int_params = list()
                    if self.weights[site-1] != 'I':
                        tmpstr = '%s2'%(self.weights[site-1])
                        int_params.append(tmpstr)
                    sign = 1
                    for s in comb :
                        sign *= -1
                        int_sites.append(s+1);int_sites.append(s+1)
                        int_types.append('C');int_types.append('A');
                        if self.adjacency[site-1][s] != 'I' : 
                            int_params.append(self.adjacency[site-1][s])
                    if sign == -1 :
                        int_params.append('minus')
                    int_params.append('pot_const')
                    inter = Interaction(int_sites,int_types,int_params)
                    self.interactions.append(inter)        
                num += 1         
            site += 1
            
#    def add_staggered_interactions(self, StaggeredSites):
#        """
#        Method to add interactions for the staggered case.
#        
#        Parameters
#        -----------
#        StaggeredSites: list
#            List of the site indices which should be staggered.
#        """        
#            
            
    def write_axial_square_octagon_rotation(self,rows,cols,filename):
        sites = rows * cols * 4
        fp = open(filename,"w")
        fp.write('<ROTATION_OPS number="1">\n')
        fp.write('\t<ROTATION_OP number="%d">\n'%(sites))
        site = 0
        while site < sites : 
          from_site = site + 1
          from_plaquette = site / 4 
          from_site_within_plaquette = site % 4
          to_site_within_plaquette = from_site_within_plaquette
          if from_site_within_plaquette == 1 :
              to_site_within_plaquette = 2
          if from_site_within_plaquette == 2 :
              to_site_within_plaquette = 1
          from_plaquette_row = from_plaquette / cols
          from_plaquette_col = from_plaquette - from_plaquette_row * cols
          to_plaquette_row = rows - from_plaquette_row - 1
          to_plaquette_col = from_plaquette_col
          to_plaquette = to_plaquette_row * cols + to_plaquette_col
          to_site = to_plaquette*4 + to_site_within_plaquette + 1 
          fp.write('\t\t<MAPPING from="%d" to="%d" />\n'%(from_site, to_site))
          site +=1    
        fp.write('\t</ROTATION_OP>\n')
        fp.write('</ROTATION_OPS>\n')    
        fp.close()    
        
class Unit_cell:
    def __init__(self,sites,dimension):
        self.dimension = dimension
        self.sites = sites
        
class Interaction:
    def __init__(self,sites,types=[],params=[]):
        if isinstance(sites,Interaction):
            self.sites = list(sites.sites)
            self.types = list(sites.types)
            self.params = list(sites.params)
        else:
            self.sites = sites #list of vectors
            self.types = types  #list of characters representing operator on term
            self.params = params
            
    def interaction_string(self):
        int_str = str()
        i = 0 
        while i < self.sites.__len__() :
            if i > 0 : int_str += ','
            int_str += '%d %c'%(self.sites[i],self.types[i])
            i+=1
        int_str += '* '
        i = 0 
        while i < self.params.__len__() :
            if i > 0 : int_str += ','
            int_str += '%s'%(self.params[i])
            i+=1 
        return int_str            
        
class Operator_matrix:
    def __init__(self,filename):
        #first work out the dimension of the matrix
        self.dimension = 0 
        fp = open(filename,"r")
        line = fp.readline()
        while ( line.__len__() > 0 ) : 
            nums = extract_numbers(line)
            if self.dimension < nums[0] : 
                self.dimension = nums[0]
            idx = 1 
            while idx < nums.__len__() : 
                if self.dimension < nums[idx] : 
                    self.dimension = nums[idx]
                idx += 2 
            line = fp.readline()
        fp.close()
        
        self.dimension += 1
        #create matrix
        self.matrix = zeros([self.dimension,self.dimension],dtype=complex64)
        
        #populate the matrix 
        fp = open(filename,"r")
        line = fp.readline()
        while ( line.__len__() > 0 ) : 
            nums = extract_numbers(line)
            idx = 1
            while idx < nums.__len__() : 
                self.matrix[nums[0],nums[idx]] = nums[idx+1] 
                idx += 2 
            line = fp.readline()
        fp.close()
        
        
class ThinTorus_Lattice(Lattice): 
    def add_interactions_and_params(self):
        self.parameters.extend(['Jx','Jy','Jz'])
        self.length = self.sites/2  #the length is half the number of sites.
        rung = 0
        while rung < self.length :
            iter = Interaction([2*rung+1,2*rung+2],['Z','Z'],['Jz'])
            self.interactions.append(iter)
            if rung != (self.length-1) : 
                iter = Interaction([2*rung+1,2*(rung+1)+2],['X','X'],['Jx'])
                self.interactions.append(iter)
                iter = Interaction([2*rung+2,2*(rung+1)+1],['Y','Y'],['Jy'])
                self.interactions.append(iter)
            elif self.pbc :
                iter = Interaction([2*rung+1,2],['X','X'],['Jx'])
                self.interactions.append(iter)
                iter = Interaction([2*rung+2,1],['Y','Y'],['Jy'])
                self.interactions.append(iter)
            rung+=1
            
            
class Heisenberg_Lattice(Lattice): 
    def add_interactions_and_params(self):
        self.parameters.extend(['Jx','Jy','Jz'])        
        site = 1
        while site <= self.sites :
            next_site = (site % self.sites) + 1 
            iter = Interaction([site,next_site],['X','X'],['Jx'])
            self.interactions.append(iter)
            iter = Interaction([site,next_site],['Y','Y'],['Jy'])
            self.interactions.append(iter)
            iter = Interaction([site,next_site],['Z','Z'],['Jz'])
            self.interactions.append(iter)
            site+=1
            
class Heisenberg_Spin_One_Lattice(Lattice): 
    def add_interactions_and_params(self):
        self.parameters.extend(['J'])
        self.parameters.extend(['Jz'])
        site = 1
        while site <= self.sites :
            if self.sites != 2 or site != 2 :
                next_site = (site % self.sites) + 1 
                iter = Interaction([site,next_site],['R','L'],['J'])
                self.interactions.append(iter)
                iter = Interaction([site,next_site],['L','R'],['J'])
                self.interactions.append(iter)
                iter = Interaction([site,next_site],['S','S'],['Jz'])
                self.interactions.append(iter)
            site+=1            
            
    def task_values(self):        
        parameter_vals = dict()        
        parameter_vals['J'] = 1.0
        parameter_vals['Jz'] = 1.0
        return parameter_vals
        
        
class AKLT_Spin_One_Lattice(Lattice): 
    def add_interactions_and_params(self):
        self.parameters.extend(['J1'])
        self.parameters.extend(['J2'])
        self.parameters.extend(['J3'])
        iter = Interaction([],[],['J3'])
        self.interactions.append(iter)
        site = 1
        while site <= self.sites :
            if self.sites != 2 or site != 2 :
                next_site = (site % self.sites) + 1 
                iter = Interaction([site,next_site],['R','L'],['J1'])
                self.interactions.append(iter)
                iter = Interaction([site,next_site],['L','R'],['J1'])
                self.interactions.append(iter)
                iter = Interaction([site,next_site],['S','S'],['J1'])
                self.interactions.append(iter)
                
                iter = Interaction([site,next_site,site,next_site],['R','L','R','L'],['J2'])
                self.interactions.append(iter)
                iter = Interaction([site,next_site,site,next_site],['R','L','L','R'],['J2'])
                self.interactions.append(iter)
                iter = Interaction([site,next_site,site,next_site],['R','L','S','S'],['J2'])
                self.interactions.append(iter)
                iter = Interaction([site,next_site,site,next_site],['L','R','R','L'],['J2'])
                self.interactions.append(iter)
                iter = Interaction([site,next_site,site,next_site],['L','R','L','R'],['J2'])
                self.interactions.append(iter)
                iter = Interaction([site,next_site,site,next_site],['L','R','S','S'],['J2'])
                self.interactions.append(iter)
                iter = Interaction([site,next_site,site,next_site],['S','S','R','L'],['J2'])
                self.interactions.append(iter)
                iter = Interaction([site,next_site,site,next_site],['S','S','L','R'],['J2'])
                self.interactions.append(iter)
                iter = Interaction([site,next_site,site,next_site],['S','S','S','S'],['J2'])
                self.interactions.append(iter)
            site+=1            
            
    def task_values(self):        
        parameter_vals = dict()        
        parameter_vals['J1'] = 0.5
        parameter_vals['J2'] = 1.0/6.0
        parameter_vals['J3'] = float(self.sites)/3.0
        return parameter_vals
        
            
                                
class Haldane_Shastry_Chain(Lattice): 
    def __init__(self,sites,unit_cell,pbc=True, Delta=1.0):
        Lattice.__init__(self, sites, unit_cell, pbc)
        self.Delta = Delta

    def add_interactions_and_params(self):        
        half = self.sites / 2
        l = 1
        while l <= half:          
          self.parameters.extend(['J%d'%(l)])        
          l+=1
        self.parameters.extend(['Delta'])
        int_array = zeros([self.sites,self.sites])        
        site = 0
        while site < self.sites :
            separation = 1
            while separation <= half :
                next_site = (site + separation) % self.sites 
                #if abs(next_site - site) == separation  or (abs(next_site - site)%half) == separation:
                int_array[site,next_site] = separation
                int_array[next_site,site] = separation                  
                #else : 
                #      print '|%d - %d| != %d'%(site, next_site , separation)
                separation +=1
            site+=1          
        site = 0
        print int_array
        while site < self.sites :
            next_site = site
            while next_site < self.sites :
                separation = int_array[site,next_site]
                if separation > 0 :
                      interaction = 'J%d'%(separation)
                      iter = Interaction([site+1,next_site+1],['X','X'],[interaction])
                      self.interactions.append(iter)
                      iter = Interaction([site+1,next_site+1],['Y','Y'],[interaction])
                      self.interactions.append(iter)
                      iter = Interaction([site+1,next_site+1],['Z','Z'],[interaction, 'Delta'])
                      self.interactions.append(iter)
                next_site += 1
            site+=1
                                    
    def task_values(self):
        half = self.sites / 2
        N = float(self.sites)
        parameter_vals = dict()
        parm = 1
        parameter_vals['Delta'] = self.Delta
        while parm <= half:
            parameter_vals['J%d'%(parm)] = ((N/pi)*sin(pi*float(parm)/N))**(-2)
            #parameter_vals['J%d'%(parm)] = (sin(pi*float(parm)/N))**(-2)
            parm+=1
        return parameter_vals
        
        
class Haldane_Shastry_Chain_Spin_One(Lattice): 
    def add_interactions_and_params(self):        
        half = self.sites / 2
        l = 1
        while l <= half:          
          self.parameters.extend(['J%d'%(l)])        
          l+=1
        int_array = zeros([self.sites,self.sites])        
        site = 0
        while site < self.sites :
            separation = 1
            while separation <= half :
                next_site = (site + separation) % self.sites 
                #if abs(next_site - site) == separation  or (abs(next_site - site)%half) == separation:
                int_array[site,next_site] = separation
                int_array[next_site,site] = separation                  
                #else : 
                #      print '|%d - %d| != %d'%(site, next_site , separation)
                separation +=1
            site+=1          
        site = 0
        print int_array
        while site < self.sites :
            next_site = site
            while next_site < self.sites :
                separation = int_array[site,next_site]
                if separation > 0 :
                      interaction = 'J%d'%(separation)
                      iter = Interaction([site+1,next_site+1],['R','L'],[interaction])
                      self.interactions.append(iter)
                      iter = Interaction([site+1,next_site+1],['L','R'],[interaction])
                      self.interactions.append(iter)
                      iter = Interaction([site+1,next_site+1],['S','S'],[interaction])
                      self.interactions.append(iter)
                next_site += 1
            site+=1
                                    
    def task_values(self):
        half = self.sites / 2
        N = float(self.sites)
        parameter_vals = dict()
        parm = 1
        while parm <= half:
            parameter_vals['J%d'%(parm)] = ((N/pi)*sin(pi*float(parm)/N))**(-2)
            #parameter_vals['J%d'%(parm)] = (sin(pi*float(parm)/N))**(-2)
            parm+=1
        return parameter_vals        

class SpinHalfSSquared_Lattice(Lattice): 
    def add_interactions_and_params(self):
        self.parameters.extend(['N','J'])        
        iter = Interaction([],[],['N'])
        self.interactions.append(iter)
        site = 1
        while site <= self.sites :
            other_site = site + 1 
            while other_site <= self.sites :                                                  
                  iter = Interaction([site,other_site],['X','X'],['J'])
                  self.interactions.append(iter)
                  iter = Interaction([site,other_site],['Y','Y'],['J'])
                  self.interactions.append(iter)
                  iter = Interaction([site,other_site],['Z','Z'],['J'])
                  self.interactions.append(iter)
                  other_site+=1
            site+=1
            
    def task_values(self):        
        N = float(self.sites)
        parameter_vals = dict()        
        parameter_vals['N'] = (3.0*N)/4.0
        parameter_vals['J'] = 1.0/2.0        
        return parameter_vals


class SpinOneSSquared_Lattice(Lattice): 
    def add_interactions_and_params(self):
        self.parameters.extend(['N','J'])        
        #iter = Interaction([],[],['N'])
        #self.interactions.append(iter)
        site = 1
        while site <= self.sites :
            other_site = 1
            while other_site <= self.sites :                                                  
                  iter = Interaction([site,other_site],['R','L'],['J'])
                  self.interactions.append(iter)
                  iter = Interaction([site,other_site],['L','R'],['J'])
                  self.interactions.append(iter)
                  iter = Interaction([site,other_site],['S','S'],['J'])
                  self.interactions.append(iter)
                  other_site+=1
            site+=1
            
    def task_values(self):        
        N = float(self.sites)
        parameter_vals = dict()        
        parameter_vals['N'] = (3.0*N)/4.0
        parameter_vals['J'] = 1.0#/2.0        
        return parameter_vals

            
def extract_numbers(str):
    idx = 0
    nums = list()
    #first split at the colon
    parts = str.split(':')
    
    #divide first part which should be of format 'row %d' where %d is the row number
    subparts = parts[0].split(' ') 
    nums.append(int(subparts[1])) #append the row number
    
    parts[1]
    subparts = parts[1].lstrip().rstrip().split('(')
    
    for part in subparts :
        try :
            idx = part.index(')') #this just ensures it is an actual pair
            part = part.rstrip().rstrip(')')
            subsubparts = part.split(',')
            col = int(subsubparts[0])
            val = to_complex(str)
            nums.append(col)
            nums.append(val)
        except     ValueError:
            pass
    return nums
    
def to_complex(str):
    num = complex(0,0) #create a zero valued complex number
    str = str.lstrip().rstrip() #strip whitespace from sites
    idx = 0 
    buf = ''
    imag = False
    in_num = False
    while idx < str.__len__() : 
        if (str[idx] == '+' or str[idx] == '-')  and in_num == False : 
            buf += str[idx]
            in_num = True
        elif (str[idx] >= '0' and str[idx] <= '9') or str[idx] == 'e' or str[idx] == '.' : 
            buf += str[idx]
            in_num = True
        elif (str[idx] == 'j' or str[idx] == 'J' or str[idx] == 'i' or str[idx] == 'I' ) : 
            imag = True
        elif (str[idx] == '+' or str[idx] == '-' or str[idx] == ' ')  and in_num == True : 
            num = float(buf)
            if imag :
                num += complex(0,num)
            else  :
                num += complex(num,0)
            
            if str[idx] == '-' :
                buf = '-'
                in_num = True
            else :
                buf = ''
                in_num = False
        idx += 1 
        if in_num : 
            num = float(buf)
            if imag :
                num += complex(0,num)
            else  :
                num += complex(num,0)

    return num
    
    
def write_momentum_file(vecs,norm_vecs,dims,filename):
    fp = open(filename,"w")
    vec_idx = 0
    while vec_idx < 2 :
        fp.write('LATTICE_VECTOR%d = %d,%d\n'%(vec_idx+1,vecs[vec_idx][0],vecs[vec_idx][1]))
        vec_idx+=1
    vec_idx = 0
    while vec_idx < 2 :
        fp.write('NORM_VECTOR%d = %f,%f\n'%(vec_idx+1,norm_vecs[vec_idx][0],vecs[vec_idx][1]))
        vec_idx+=1
    fp.write('LATTICE_DIMENSIONS = %d, %d\n'%(dims[0],dims[1]))
    fp.close()            


def write_chain_translation(sites,filename):
    fp = open(filename,"w")
    fp.write('<ROTATION_OPS number="1">\n')
    fp.write('\t<ROTATION_OP number="%d">\n'%(sites))
    site = 0
    while site < sites : 
      from_site = site + 1
      to_site = ((site + 1) % sites) + 1 
      fp.write('\t\t<MAPPING from="%d" to="%d" />\n'%(from_site, to_site))
      site +=1    
    fp.write('\t</ROTATION_OP>\n')
    fp.write('</ROTATION_OPS>\n')    
    fp.close()            
