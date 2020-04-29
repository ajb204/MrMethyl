#!/usr/bin/python

import sys,os,re
from commAux import CommAux
import numpy,scipy,copy,random
from datetime import datetime
import matplotlib.pyplot as plt
from scipy.linalg import expm
from timeit import default_timer
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy import mat,zeros
from scipy.optimize import leastsq
from scipy.optimize import least_squares


############################################################
# START FROM HERE!
#dIS is hbar gammaIgammaS / 4pi e0 r^6
#csa is hbar ((parallel-permendicular) * gamma * B0)/3
#J(teff,w) is teff/(1+wc^2teff^2) - Kay type - numerical constants are as stated
#P20=1/2(3cos^2b-1) #polynomial
#P21=(3sinbcosb) #polynomial
#P22=(3sin^2b)    #polynomial
#
#To use:
# first make a basis.
# then construct the desired Hamiltonian from dipolar and CSA terms
# then evaluate the spectral density functions. Options here include:
#    -kay True -approximate the functions with (3cos^2t-1) angles
#    -kay False -use explicit geometry
#    -macro True - take macromolecular limit - tauM terms-> taum, high frequencies  to zero.
#    -macro False - don't take macrmolecular limit
#    -taum 'taum' - allow for methyl tumbling
#    -taum 0     -don't allow for methyl tumbling
# After spectral density funcitons have been constructed, calculate desired relaxation rates.
# Have fun!
############################################################

class mr_noesy():
    def __init__(self,peak_file_dir,peak_list_dir,pdb_file_dir,macro=False, ax=True, xx=True, cx=True, ca=True, sparse=True, sym=True, irreducible=False, num=False):
        #THIS IS THE INIT CLASS FOR SET-UP:
        #CREATE A CommAux class instance to deal with cross-relaxation calculations
        self.cross=CommAux('CH3CH3',latex_file='output/eps/cross_latex.tex',sparse=sparse,sym=sym)
        #OPEN LOG FILE FOR WRITING TO:
        print_log=open('output/log/log_Comm','w')
        
        #INITIALIZE CommAux class self.cross with spins&interactions:
        if(ax):
            self.cross.GetDip('H2','C1','d1',tag='single')
            self.cross.GetDip('H3','C1','d1',tag='single')
            self.cross.GetDip('H4','C1','d1',tag='single')
            self.cross.GetDip('H6','C5','d2',tag='single')
            self.cross.GetDip('H7','C5','d2',tag='single')
            self.cross.GetDip('H8','C5','d2',tag='single')
        if(xx):
            self.cross.GetDip('H2','H3','e1',tag='single')
            self.cross.GetDip('H2','H4','e1',tag='single')
            self.cross.GetDip('H3','H4','e1',tag='single')
            self.cross.GetDip('H6','H7','e1',tag='single')
            self.cross.GetDip('H6','H8','e1',tag='single')
            self.cross.GetDip('H7','H8','e1',tag='single')
        if(cx):
            self.cross.GetCSA('H2','cX1',tag='single')
            self.cross.GetCSA('H3','cX1',tag='single')
            self.cross.GetCSA('H4','cX1',tag='single')
            self.cross.GetCSA('H6','cX2',tag='single')
            self.cross.GetCSA('H7','cX2',tag='single')
            self.cross.GetCSA('H8','cX2',tag='single')
        if(ca):
            self.cross.GetCSA('C1','cA1',tag='static')
            self.cross.GetCSA('C5','cA2',tag='static')
        for i in (2,3,4):
            i_str='H'+str(i)
            for j in (6,7,8):
                j_str='H'+str(j)
                self.cross.GetDip(i_str,j_str,'m',tag='double')
        
        #Set self.cross.macro and call CrossCorrHam() to calculate cross-relaxation terms:
        self.cross.macro=macro
        self.cross.CrossCorrHam()
        #Load up the NMR peak list file:
        cross_matrix,reflected,d1,d2,c1,h1,c2,h2,map_key_nmr,map_index_nmr,pairs,mask_matrix=self.peak_file_parser(peak_file_dir,peak_list_dir)
        #Generate the "wobble mask" which screens out guys who are dynamic:
        wobble_mask=self.wobble_masker(cross_matrix,map_index_nmr,map_key_nmr,c1,h1,c2,h2)
        #Generate data mask:
        self.mask=numpy.logical_and(wobble_mask,mask_matrix)
        #Mask out the diagonals: (I think we want to do this?)
        diag_mask=numpy.ones(self.mask.shape,dtype=bool)
        numpy.fill_diagonal(diag_mask,False) 
        self.mask=numpy.logical_and(self.mask,diag_mask)
        print 'NUMBER OF PEAKS INCLUDED IN OPTIMIZATION: ',numpy.sum(self.mask)
        print 'NUMBER OF RESIDUES INCLUDED IN WOBBLE-MASK:',numpy.sqrt(numpy.sum(wobble_mask))
        print 'NUMBER OF PEAKS INCLUDED IN DATA SET MASK:',numpy.sum(mask_matrix)
        #Load up the PDB file:
        location,orientation,sites,self.map_key_xrd,self.map_index_xrd,self.residues=self.get_pdb_coord_methyl(pdb_file_dir)
        #Set-up other parameters:
        self.cross.SetFreq(600)
        cA=20.
        cX=1.
        self.tau=400*1E-3
        self.num_samples=2000
        #Assemble the basis set:
        self.cross.Assemble(('H2z','H3z','H4z'),'M1z')
        self.cross.Assemble(('H6z','H7z','H8z'),'M2z')
        #SetPars section: (Initializes values in CommAux module cross for order parameter calculations)
        num=location.shape[0]
        i,j=numpy.arange(num),numpy.arange(num)
        ii,jj=numpy.meshgrid(i,j,indexing='ij')
        self.cross.pars['v_m1'],self.cross.pars['v_m2']=orientation[i,:],orientation[j,:]
        self.cross.pars['c1_site'],self.cross.pars['c2_site']=location[i,:],location[j,:]
        self.cross.pars['sites1'],self.cross.pars['sites2']=sites[i,:,:],sites[j,:,:]
        self.cross.pars['cA'],self.cross.pars['cX']=cA,cX
        v_mm=location[jj,:]-location[ii,:]
        self.cross.pars['v_mm']=v_mm
        self.cross.SetPars('CH3CH3')
        #Set-up indexing matrices to convert between XRD and NMR indexing:
        xrd_indices,nmr_indices=[],[]
        for key in map_key_nmr:
            xrd_indices.append(self.map_key_xrd[key])
            nmr_indices.append(map_key_nmr[key])
        xrd_indices,nmr_indices=numpy.array(xrd_indices),numpy.array(nmr_indices)
        self.ii_xrd,self.jj_xrd=numpy.meshgrid(xrd_indices,xrd_indices,indexing='ij')
        self.ii_nmr,self.jj_nmr=numpy.meshgrid(nmr_indices,nmr_indices,indexing='ij')
        #Create and calculate distance matrix:
        self.dist=numpy.zeros(self.ii_nmr.shape)
        dist_matrix=numpy.sqrt(numpy.sum(v_mm*v_mm,axis=2))
        self.dist[self.ii_nmr,self.jj_nmr]=dist_matrix[self.ii_xrd,self.jj_xrd]
        #Create and fill the NMR peak matrix:
        self.nmr_peaks=numpy.zeros(self.ii_nmr.shape)
        self.nmr_peaks[self.ii_nmr,self.jj_nmr]=cross_matrix[self.ii_nmr,self.jj_nmr]
        #Create the XRD peak matrix: (This will be populated in the optimization loop, when we do a function call)
        self.xrd_peaks=numpy.zeros(self.ii_nmr.shape)
        #Initialize the initial intensities as the diagonal NMR peaks (this should do as a halfway decent proxy)
        d=numpy.arange(len(nmr_indices))
        self.init_intensities=numpy.zeros(self.nmr_peaks.shape)
        self.init_intensities[d,d]=self.nmr_peaks[d,d]
        #Initialize t_m,t_c:
        self.cross.tm=25.E-12
        self.cross.tc=25.E-9
        
        
        #Run the optimization:
        x_init=self.pack()
        print x_init
        print x_init.shape
        x_opt=leastsq(self.chi_func,x_init)
        self.unpack(x_opt[0])
        
    
    

    
    def unpack(self,x):
#         self.cross.tm=x[0,1]
#         self.cross.tc=x[1,1]
#         self.init_intensities=x[:,0]
        self.cross.tm=x[0]
        self.cross.tc=x[1]
        d=numpy.arange(self.nmr_peaks.shape[0])
        self.init_intensities[d,d][self.mask[d,d]]=x[2:]
    
    def pack(self):
        d=numpy.arange(self.nmr_peaks.shape[0])
        x=numpy.array([self.cross.tm,self.cross.tc])
        x=numpy.append(x,self.init_intensities[d,d][self.mask[d,d]])
#         x=numpy.zeros((len(self.init_intensities),2))
#         x[:,0]=self.init_intensities
#         x[0,1],x[1,1]=self.cross.tm,self.cross.tc
        return x
        
    
    def chi_func(self,x):
        self.unpack(x)
        self.CalcModel()
        print 'chi2',numpy.sum((self.nmr_peaks[self.mask]-self.xrd_peaks[self.mask])**2.)
        print 'difference matrix = \n',self.nmr_peaks[self.mask]-self.xrd_peaks[self.mask]
        return (self.nmr_peaks[self.mask]-self.xrd_peaks[self.mask]).flatten()
        
    
    def CalcModel(self):
        #Calculate the rates:
        r_cross=self.cross.CalcRate('M1z', 'M2z')
        r_auto=self.cross.CalcRate('M1z', 'M1z')
        #Create matrices to store results:
        num_methyls=r_cross.shape[0]
        result_matrix=numpy.zeros((num_methyls,num_methyls,self.num_samples))
        #Calculate intensities by repeating for 2000 combinations:
        for rep in numpy.arange(self.num_samples):
            list_methyls=self.methyl_list_generator(self.residues)
            index_list=numpy.array([self.map_key_xrd[x] for x in list_methyls])
            rate_matrix=numpy.zeros((len(index_list),len(index_list)))
            exp_matrix=numpy.zeros((len(index_list),len(index_list)))
            #Have 2 indexing regimes. x,y are for the compacted rate_matrix. This is for the calculation here
            x,y=numpy.arange(len(index_list)),numpy.arange(len(index_list))
            xx,yy=numpy.meshgrid(x,y,indexing='ij')
            #i,j are for the expanded matrix. Maps the compacted matrix onto the expanded matrix. 
            i,j=index_list,index_list
            ii,jj=numpy.meshgrid(i,j,indexing='ij')
            #Load up and take matrix exponential:
            rate_matrix[xx,yy]=r_cross[ii,jj]
            rate_matrix[x,y]=numpy.sum((r_auto[ii,jj]-r_auto[ii,ii]),axis=1)+r_auto[i,j]
            exp_matrix[:,:]=expm(-rate_matrix[:,:]*self.tau)
            result_matrix[ii,jj,rep]=exp_matrix[xx,yy]

        #AVERAGE OVER 2000 samples:
        sim_matrix=numpy.average(result_matrix[:,:,:], axis=2)
        #MATCH SIMULATED XRD DATA ONTO NMR INDICES:
        self.xrd_peaks[self.ii_nmr,self.jj_nmr]=sim_matrix[self.ii_xrd,self.jj_xrd]*self.init_intensities[self.ii_nmr,self.jj_nmr]
    
    

        
        
        
        
        
        
        
        
        
        
    #NOW HAVE ALL THE OTHER RANDOM FUNCTIONS:    
    def wobble_masker(self,cross_matrix,map_index_nmr,map_key_nmr,c1,h1,c2,h2,lim=0.80):
        indices=[]
        for i in numpy.arange(cross_matrix.shape[0]):
            res=map_index_nmr[i]
            if(res[0:3]=='ALA'):
                i=map_key_nmr[res]
                if(c1[i,i]!=0.):
                    indices.append(i)
            elif(res[0:3]=='ILE'):
                p_g=(14.8-c1[i,i])/5.5
                if(c1[i,i]!=0. and (p_g<1.-lim or p_g>lim)):
                    indices.append(i)
            elif(res[0:3]=='LEU'):
                res_a,res_b=res[0:-1]+'A',res[0:-1]+'B'
                i_a,i_b=map_key_nmr[res_a],map_key_nmr[res_b]
                c_a,c_b=c2[i_a,i_a],c2[i_b,i_b]
                p_t=0.5+(c_a-c_b)/10.
                if(c_a!=0 and c_b!=0 and (p_t<1.-lim or p_t>lim)):
                    indices.append(i)
            elif(res[0:3]=='VAL'):
                res_a,res_b=res[0:-1]+'A',res[0:-1]+'B'
                i_a,i_b=map_key_nmr[res_a],map_key_nmr[res_b]
                c_a,c_b=c2[i_a,i_a],c2[i_b,i_b]
                if(c_a!=0. and c_b!=0. and ((c_a>21.0 and c_b>21.0) or (c_a<20. and c_b<20.))):
                    indices.append(i)
        i,j=numpy.array(indices),numpy.array(indices)
        ii,jj=numpy.meshgrid(i,j,indexing='ij')
        wobble_mask=numpy.zeros(cross_matrix.shape,dtype=bool)
        wobble_mask[ii,jj]=True #This will then be *AND* with the data mask
        return wobble_mask    
        

    def string_ab_swap(self,string):
        if string[-1]=='A':
            return string[:-1]+'B'
        elif string[-1]=='B':
            return string[:-1]+'A'    

    def get_pdb_coord_methyl(self,pdb_file):
        #Initialize stuff - Set which file to open; initilaize lists for all atoms in file, and then result list - methyl_list
        print pdb_file
        file = open(pdb_file, 'r')
        atom_list = []
        methyl_list = {}

        #Iterate over lines in file, and fill atom_list with position, type & residue (type & number) for each atom:
        for line in file:
            if(line[0:4]=='ATOM'):
                atom_number = line[6:12]
                atom_name = line[13:16]
                residue_name = line[17:20]
                residue_number = int(line[22:26].strip())
                x_coord = float(line[30:38])
                y_coord = float(line[38:46])
                z_coord = float(line[46:54])

                atom={}
                atom['num']=atom_number
                atom['ele']=atom_name
                atom['residue']=residue_name+' '+str(residue_number)
                atom['x'], atom['y'], atom['z']=x_coord,y_coord,z_coord
                atom_list.append(atom)

        #Now iterate over all atoms in atom_list and group in residue (type & number):
        residue_list={}
        for atom in atom_list:
            if(atom['residue'] in residue_list):
                residue_list[atom['residue']][atom['ele']]=(atom['x'],atom['y'],atom['z'])
            else:
                residue_list[atom['residue']]={}
                residue_list[atom['residue']][atom['ele']]=(atom['x'],atom['y'],atom['z'])

        #Now iterate over residues and generate methyls from ILE, LEU, VAL:
        #General strategy - Each methyl consits of 3 vectors - c_m, c_a, c_b. These totally specify a methyl groups position & orientation, *including* psi angle
        residues=[]
        count_ile,count_leu,count_val,count_ala = 0,0,0,0
        for key in residue_list:
            residue = residue_list[key]
            if(key[0:3]=='ILE'):
                residues.append(key)
                methyl=self.CalcMetGeom(residue['CD1'],residue['CG1'],residue['CB '])
                methyl_list[key+' A']=methyl
                count_ile+=1
            elif(key[0:3]=='LEU'):
                residues.append(key)
                methyl=self.CalcMetGeom(residue['CD1'],residue['CG '],residue['CB '])
                methyl_list[key+' A']=methyl
                methyl=self.CalcMetGeom(residue['CD2'],residue['CG '],residue['CB '])
                methyl_list[key + ' B']=methyl
                count_leu+=1
            elif(key[0:3]=='VAL'):
                residues.append(key)
                methyl=self.CalcMetGeom(residue['CG1'],residue['CB '],residue['CA '])
                methyl_list[key+' A']=methyl
                methyl=self.CalcMetGeom(residue['CG2'],residue['CB '],residue['CA '])
                methyl_list[key + ' B']=methyl
                count_val+=1
            elif(key[0:3]=='ALA'):
                residues.append(key)
                methyl=self.CalcMetGeom(residue['CB '],residue['CA '],residue['C  '])
                methyl_list[key+' A']=methyl
                count_ala+=1

        #Want to return: mapping dictionaries from methyl residue name & number (EG. LEU 113 A) to numpy index, and index to name.
        #3 numpy arrays - N*3 array containing coordinates of methyl C (c_m), N*3 for orientation (v_m), and N*3*3 positions of methyl H (sites)
        num=len(methyl_list)
        location=numpy.zeros((num, 3))
        orientation=numpy.zeros((num, 3))
        sites=numpy.zeros((num, 3, 3))
        map_key={}
        map_index={}
        i=0 #Index which incremented to slot methyls into numpy arrays
        print 'count_leu = ',count_leu
        print 'count_val = ',count_val
        print 'count_ile = ',count_ile
        print 'count_ala = ',count_ala
        for key in methyl_list:
            location[i, :]=methyl_list[key]['loc']
            orientation[i, :]=methyl_list[key]['orientation']
            sites[i, :, :]=methyl_list[key]['sites']
            map_key[key]=i
            map_index[i]=key
            i+=1
        return location, orientation, sites, map_key, map_index, list(set(residues))

    
    def peak_list_parser(self,peak_list_dir):
        file=open(peak_list_dir,'r')
        peak_list=[]
        for line in file:
            res1,res2=line.split()[0],line.split()[1]
            res1,res2=self.peak_list_str_parser(res1),self.peak_list_str_parser(res2)
            peak_list.append((res1,res2))
            peak_list.append((res1,res1))
            peak_list.append((res2,res2))
            if(res1[0:3]=='LEU' or res1[0:3]=='VAL'):
                peak_list.append((res1[0:-1]+'B',res2))
                peak_list.append((res1[0:-1]+'B',res1[0:-1]+'B'))
            if(res2[0:3]=='LEU' or res2[0:3]=='VAL'):
                peak_list.append((res1,res2[0:-1]+'B'))
                peak_list.append((res2[0:-1]+'B',res2[0:-1]+'B'))
        return peak_list

    
    def peak_list_str_parser(self,res_str):    
        res_letter=res_str[0]
        integers=[int(s) for s in re.findall(r'-?\d+\.?\d*', res_str)]
        if(res_letter=='I'):
            name='ILE'+' '+str(integers[0])+' A'
        elif(res_letter=='A'):
            name='ALA'+' '+str(integers[0])+' A'
        elif(res_letter=='L'):
            name='LEU'+' '+str(integers[0])+' A'
        elif(res_letter=='V'):
            name='VAL'+' '+str(integers[0])+' A'
        return name

    
    def peak_file_parser(self,peak_file_dir,peak_list_dir):
        peak_list=self.peak_list_parser(peak_list_dir)
        file=open(peak_file_dir, 'r')
        cross,reflected={},{}
        diag_1,diag_2={},{}
        c1_shift,c2_shift={},{}
        h1_shift,h2_shift={},{}
        for line in file:
            data=line.split()
            res1=self.residue_parser(data[0])
            res2=self.residue_parser(data[1])
            c1,h1=data[2],data[3]
            c2,h2=data[4],data[5]
            v1,v2=data[6],data[7]
            a1,a2=data[14],data[15]
            if not(res1 in cross):
                cross[res1]={}
                reflected[res1]={}
                diag_1[res1]={}
                diag_2[res1]={}
                c1_shift[res1],c2_shift[res1]={},{}
                h1_shift[res1],h2_shift[res1]={},{}
            cross[res1][res2]=v1
            reflected[res1][res2]=v2
            diag_1[res1][res2]=a1
            diag_2[res1][res2]=a2
            c1_shift[res1][res2],c2_shift[res1][res2]=c1,c2
            h1_shift[res1][res2],h2_shift[res1][res2]=h1,h2
        map_key,map_index={},{}
        i=0
        for res1 in cross:
            map_key[res1]=i
            map_index[i]=res1
            i+=1

        num=len(cross)
        cross_matrix,reflected_matrix,d1_matrix,d2_matrix=numpy.zeros((num,num)),numpy.zeros((num,num)),numpy.zeros((num,num)),numpy.zeros((num,num))
        c1_matrix,h1_matrix,c2_matrix,h2_matrix=numpy.zeros((num,num)),numpy.zeros((num,num)),numpy.zeros((num,num)),numpy.zeros((num,num))
        mask_matrix=numpy.zeros((num,num),dtype=bool)

        peak_list_return=[]
        for peak in peak_list:
            res1,res2=peak
            try:
                i1,i2=map_key[res1],map_key[res2]
                cross_matrix[i1,i2]=cross[res1][res2]
                reflected_matrix[i1,i2]=reflected[res1][res2]
                d1_matrix[i1,i2]=diag_1[res1][res2]
                d2_matrix[i1,i2]=diag_2[res1][res2]
                c1_matrix[i1,i2],c2_matrix[i1,i2]=c1_shift[res1][res2],c2_shift[res1][res2]
                h1_matrix[i1,i2],h2_matrix[i1,i2]=h1_shift[res1][res2],h2_shift[res1][res2]
                if(cross_matrix[i1,i2]!=0 and reflected_matrix[i1,i2]!=0. and d1_matrix[i1,i2]!=0. and d2_matrix[i1,i2]!=0. and c1_matrix[i1,i2]!=0. and c2_matrix[i1,i2]!=0. and h1_matrix[i1,i2]!=0. and h2_matrix[i1,i2]!=0.):
                    mask_matrix[i1,i2]=True
                    peak_list_return.append((res1,res2))
            except KeyError:
                pass
        return cross_matrix,reflected_matrix,d1_matrix,d2_matrix,c1_matrix,h1_matrix,c2_matrix,h2_matrix,map_key,map_index,peak_list_return,mask_matrix

    
    def nmr_masker(self,h1,h2,c1,c2,map_key,pairs,threshold=0.01):
        #SETUP MATRICES TO STORE DATA:
        num=h1.shape[0]
        mask=numpy.zeros((num,num),dtype=bool)
        i_m,j_m=[],[]
        for peak1 in pairs:
            for peak2 in pairs:
                i,j,k,l=map_key[peak1[0]],map_key[peak1[1]],map_key[peak2[0]],map_key[peak2[1]]
                diff=2.*numpy.sqrt(((h1[i,j]-h1[k,l])/(h1[i,j]+h1[k,l]))**2 + ((h2[i,j]-h2[k,l])/(h2[i,j]+h2[k,l]))**2 + ((c1[i,j]-c1[k,l])/(c1[i,j]+c1[k,l]))**2 + ((c2[i,j]-c2[k,l])/(c2[i,j]+c2[k,l]))**2)
                if (diff<threshold) and (peak1[0]!=peak2[0]) and (peak1[1]!=peak2[1]):
                    i_m.append(i)
                    i_m.append(k)
                    j_m.append(j)
                    j_m.append(l)
    #                 mask[i,j],mask[k,l]=True,True
    #                 print 'RESIDUES WHICH ARE MASKED ARE = '
    #                 print peak1,': ',h1[i,j],c1[i,j],h2[i,j],c2[i,j],'DIST = ',dist[i,j]
    #                 print peak2,': ',h1[k,l],c1[k,l],h2[k,l],c2[k,l],'DIST = ',dist[k,l]
    #                 print ' DIFF = ',diff,'\n'
        i_m,j_m=numpy.array(i_m),numpy.array(j_m)
        ii,jj=numpy.meshgrid(i_m,j_m,indexing='ij')
        mask[ii,jj]=True
        return mask

    
    def residue_parser(self,residue_str):
        res_letter=residue_str[0]
        integers=[int(s) for s in re.findall(r'-?\d+\.?\d*', residue_str)]
        out_str=''
        if(res_letter=='I'):
            out_str='ILE'+' '+str(integers[0])+' A'
        elif(res_letter=='A'):
            out_str='ALA'+' '+str(integers[0])+' A'
        elif(res_letter=='L'):
            if('D1' in residue_str):
                out_str='LEU'+' '+str(integers[0])+' A'
            elif('D2' in residue_str):
                out_str='LEU'+' '+str(integers[0])+' B'
        elif(res_letter=='V'):
            if('G1' in residue_str):
                out_str='VAL'+' '+str(integers[0])+' A'
            elif('G2' in residue_str):
                out_str='VAL'+' '+str(integers[0])+' B'
        if(out_str==''):
            print 'BLANK AMINO ACID CODE. AMINO ACID STRING TO BE PARSED = ',residue_str
            sys.exit(100)
        return out_str


    def methyl_list_generator(self,residue_list):
        #EXAMPLES OF RESIDUE NAMES:
        #'ILE   57', 'LEU  123', 'VAL  176'
        methyl_list=[]
        for residue in residue_list:
            if(residue[0:3]=='ILE' or residue[0:3]=='ALA'):
                methyl_list.append(residue+' A')
            if(residue[0:3]=='VAL' or residue[0:3]=='LEU'):
                ab = random.choice([' A', ' B'])
                methyl_list.append(residue+ab)
        return methyl_list


    def CalcMetGeom(self,c_m, c_a, c_b, r_m=1.09, theta_m=numpy.arccos(1./3.)):
        #This function takes the position of three carbon atoms - cm_pos is position of methyl carbon; ca_pos, cb_pos are position of the other 2 carbon atoms
        #This function returns the position of the three hydrogen atoms.
        c_m=numpy.array(c_m)
        c_a=numpy.array(c_a)
        c_b=numpy.array(c_b)

        v_m=numpy.array(c_m-c_a)
        v_m=v_m/numpy.sqrt(numpy.dot(v_m,v_m))
        #Calculate phi_rot, theta_rot terms. These terms are the terms which do the rotation of the methyl group axes, IE to re-orientate the methyl groups
        r=numpy.sqrt(numpy.dot(v_m,v_m))
        theta_rot=numpy.arccos(v_m[2]/r)
        phi_rot=numpy.arctan2(v_m[0],v_m[1])

        #This calculate psi angle. First calculate the rotation matrix to re-orientate methyl group:
        rot=numpy.dot(self.DoRot((0,0,1),phi_rot),self.DoRot((1,0,0),theta_rot))
        inv_rot=numpy.linalg.inv(rot)
        v_a=c_b-c_a

        v_a=numpy.dot(inv_rot,v_a)
        psi_i=numpy.arctan2(v_a[0,0],v_a[0,1])+numpy.pi/3.

        sites=3 #Keep below code generalized so we can change stuff up for diffusive model/n-site model if we want later
        delta_psi=2.*numpy.pi/sites
        psi=psi_i+(delta_psi*numpy.linspace(0,sites-1,sites))

        #Calculate these factors as are re-used (hang-over from when these were used many, many times in random walk simulations)
        sin_psi,cos_psi=r_m*numpy.sin(psi),r_m*numpy.cos(psi)
        sin_theta_m,cos_theta_m=numpy.sin(theta_m),numpy.cos(theta_m)

        #Calculate the Cartesian coordinates of the methyl hydrogens. Remember - methyl hydrogens are in a ring, but using C as ccoordinate centre, *not* centre of mass
        sites=numpy.matrix([numpy.array(sin_psi)*sin_theta_m, numpy.array(cos_psi)*sin_theta_m, numpy.full(len(sin_psi), r_m*cos_theta_m)])
        sites=numpy.dot(rot, sites)        

        #This displaces the methyl group by offset vector.
        sites[0, :]+=c_m[0]
        sites[1, :]+=c_m[1]
        sites[2, :]+=c_m[2]

        methyl = {'loc':c_m, 'orientation':v_m, 'sites':numpy.transpose(numpy.array(sites))}
        return methyl

   
    def DoRot(self,k,phi):
        k=k/numpy.sqrt(numpy.dot(k, k))
        a=numpy.cos(phi/2.)
        b=-k[0]*numpy.sin(phi/2.)
        c=-k[1]*numpy.sin(phi/2.)
        d=-k[2]*numpy.sin(phi/2.)
        r=mat(zeros((3,3)))
        r[0,0]=a**2.+b**2.-c**2.-d**2.
        r[0,1]=2*(b*c-a*d)
        r[0,2]=2*(b*d+a*c)
        r[1,0]=2*(b*c+a*d)
        r[1,1]=a**2.+c**2.-b**2.-d**2.
        r[1,2]=2*(c*d-a*b)
        r[2,0]=2*(b*d-a*c)
        r[2,1]=2*(c*d+a*b)
        r[2,2]=a**2.+d**2.-b**2.-c**2.
        return r

    
    def get_pdb_coord_methyl(self,pdb_file):
        #Initialize stuff - Set which file to open; initilaize lists for all atoms in file, and then result list - methyl_list
        print pdb_file
        file = open(pdb_file, 'r')
        atom_list = []
        methyl_list = {}

        #Iterate over lines in file, and fill atom_list with position, type & residue (type & number) for each atom:
        for line in file:
            if(line[0:4]=='ATOM'):
                atom_number=line[6:12]
                atom_name=line[13:16]
                residue_name=line[17:20]
                residue_number=int(line[22:26].strip())
                x_coord=float(line[30:38])
                y_coord=float(line[38:46])
                z_coord=float(line[46:54])
                atom={}
                atom['num']=atom_number
                atom['ele']=atom_name
                atom['residue']=residue_name+' '+str(residue_number)
                atom['x'], atom['y'], atom['z']=x_coord, y_coord, z_coord
                atom_list.append(atom)

        #Now iterate over all atoms in atom_list and group in residue (type & number):
        residue_list={}
        for atom in atom_list:
            if(atom['residue'] in residue_list):
                residue_list[atom['residue']][atom['ele']]=(atom['x'],atom['y'],atom['z'])
            else:
                residue_list[atom['residue']]={}
                residue_list[atom['residue']][atom['ele']]=(atom['x'],atom['y'],atom['z'])

        #Now iterate over residues and generate methyls from ILE, LEU, VAL:
        #General strategy - Each methyl consits of 3 vectors - c_m, c_a, c_b. These totally specify a methyl groups position & orientation, *including* psi angle
        #REF:
    #     def CalcMetGeom(c_m, c_a, c_b, r_m=1.09, theta_m=numpy.arccos(1./3.)):
        residues=[]
        count_ile, count_leu, count_val,count_ala=0,0,0,0
        for key in residue_list:
            residue=residue_list[key]
            if(key[0:3]=='ILE'):
                residues.append(key)
                methyl=self.CalcMetGeom(residue['CD1'], residue['CG1'], residue['CB '])
                methyl_list[key+' A'] = methyl
                count_ile+=1
            elif(key[0:3]=='LEU'):
                residues.append(key)
                methyl=self.CalcMetGeom(residue['CD1'], residue['CG '], residue['CB '])
                methyl_list[key+' A'] = methyl
                methyl=self.CalcMetGeom(residue['CD2'], residue['CG '], residue['CB '])
                methyl_list[key+' B'] = methyl
                count_leu+=1
            elif(key[0:3]=='VAL'):
                residues.append(key)
                methyl=self.CalcMetGeom(residue['CG1'], residue['CB '], residue['CA '])
                methyl_list[key+' A'] = methyl
                methyl=self.CalcMetGeom(residue['CG2'], residue['CB '], residue['CA '])
                methyl_list[key+' B'] = methyl
                count_val+=1
            elif(key[0:3]=='ALA'):
                residues.append(key)
                methyl=self.CalcMetGeom(residue['CB '],residue['CA '],residue['C  '])
                methyl_list[key+' A']=methyl
                count_ala+=1

        #Want to return: mapping dictionaries from methyl residue name & number (EG. LEU 113 A) to numpy index, and index to name.
        #3 numpy arrays - N*3 array containing coordinates of methyl C (c_m), N*3 for orientation (v_m), and N*3*3 positions of methyl H (sites)
    #     methyl = {'loc':c_m, 'orientation':v_m, 'sites':numpy.array(sites)}
        num=len(methyl_list)
        location=numpy.zeros((num, 3))
        orientation=numpy.zeros((num, 3))
        sites=numpy.zeros((num, 3, 3))
        map_key={}
        map_index={}
        i=0 #Index which incremented to slot methyls into numpy arrays
        print 'count_leu = ', count_leu
        print 'count_val = ', count_val
        print 'count_ile = ', count_ile
        print 'count_ala = ', count_ala
        for key in methyl_list:
            location[i, :]=methyl_list[key]['loc']
            orientation[i, :]=methyl_list[key]['orientation']
            sites[i, :, :]=methyl_list[key]['sites']
            map_key[key]=i
            map_index[i]=key
            i+=1
        return location, orientation, sites, map_key, map_index, list(set(residues))

#     self,peak_file_dir,peak_list_dir,pdb_file_dir,macro=False, ax=True, xx=True, cx=True, ca=True, sparse=True, sym=True, irreducible=False, num=False):
mr_noesy('input/correlate.3','input/ein_noes','input/1EZA.pdb')
# def DoTestSystem(macro=True, ax=True, xx=True, cx=True, ca=True, sparse=True, sym=True, irreducible=False, num=False):
#     #First set up the spin system and then cross-correlate:
#     print 'Calculating methyl-methyl rates'
#     cross = CommAux('CH3CH3', latex_file='output/eps/cross_latex.tex', sparse=sparse, sym=sym)
#     print_log=open('output/log/log_Comm','w')
#     if(ax):
#         cross.GetDip('H2', 'C1', 'd1', tag ='single')
#         cross.GetDip('H3', 'C1', 'd1', tag ='single')
#         cross.GetDip('H4', 'C1', 'd1', tag ='single')
#         cross.GetDip('H6', 'C5', 'd2', tag ='single')
#         cross.GetDip('H7', 'C5', 'd2', tag ='single')
#         cross.GetDip('H8', 'C5', 'd2', tag ='single')
#     if(xx):
#         cross.GetDip('H2', 'H3', 'e1', tag ='single')
#         cross.GetDip('H2', 'H4', 'e1', tag ='single')
#         cross.GetDip('H3', 'H4', 'e1', tag ='single')
#         cross.GetDip('H6', 'H7', 'e2', tag ='single')
#         cross.GetDip('H6', 'H8', 'e2', tag ='single')
#         cross.GetDip('H7', 'H8', 'e2', tag ='single')
#     if(cx):
#         cross.GetCSA('H2', 'cX1', tag='single')
#         cross.GetCSA('H3', 'cX1', tag='single')
#         cross.GetCSA('H4', 'cX1', tag='single')
#         cross.GetCSA('H6', 'cX2', tag='single')
#         cross.GetCSA('H7', 'cX2', tag='single')
#         cross.GetCSA('H8', 'cX2', tag='single')
#     if(ca):
#         cross.GetCSA('C1', 'cA1', tag='static')
#         cross.GetCSA('C5', 'cA2', tag='static')

#     index = numpy.arange(3)
#     for i in (2,3,4):
#         i_str = 'H'+str(i)
#         for j in (6,7,8):
#             j_str = 'H'+str(j)
#             cross.GetDip(i_str, j_str, 'm', tag='double')
    
#     cross.macro=macro
#     cross.CrossCorrHam()
    
#     cross.tm=10E-12
#     cross.tc=10E-9
#     t_point = 400*1E-3
#     cross.SetFreq(600)  #MHz. Set spectrometer proton frequency
#     cross.Assemble(('H2z','H3z','H4z'),'M1z')
#     cross.Assemble(('H6z','H7z','H8z'),'M2z')
    
#     cA=20.
#     cX=1.
    
#     #Iterate over all of the pairings:
#     #So each methyl in "methyls" will consists of 2 3-vectors - one giving position and one giving orientation
#     #These three variables hold - r_single holds auto-relaxation of methyl i on it's own, r_auto holds auto-relaxation of i with j present; r_cross holds cross-relaxation from j to i.
    
    
#     num=2
#     location = numpy.zeros((num, 3))
#     orientation = numpy.zeros((num, 3))
#     sites = numpy.zeros((num, 3, 3))
    
# #      methyl = {'loc':c_m, 'orientation':v_m, 'sites':numpy.transpose(numpy.array(sites))}
# #     def CalcMetGeom(c_m, c_a, c_b, r_m=1.09, theta_m=numpy.arccos(1./3.)):
# #     m1,m2=CalcMetGeom([0.,0.,1.],[0.,0.,0.],[1.,0.,0.]),CalcMetGeom([0.,5.,1.],[0.,5.,0.],[1.,5.,0.])
# #     m1,m2=CalcMetGeom([0.,0.,1.],[0.,0.,0.],[1.,0.,0.]),CalcMetGeom([0.,0.,6.],[0.,0.,5.],[1.,0.,5.])
#     m1,m2=CalcMetGeom([0.,0.,1.],[0.,0.,0.],[1.,0.,0.]),CalcMetGeom([0.,0.,4.],[0.,0.,5.],[1.,0.,5.])
#     location[0,:],location[1,:]=m1['loc'],m2['loc']
#     orientation[0,:],orientation[1,:]=m1['orientation'],m2['orientation']
#     sites[0,:,:],sites[1,:,:]=m1['sites'],m2['sites']
    
#     i = numpy.arange(num)
#     j = numpy.arange(num)
#     ii, jj = numpy.meshgrid(i, j,indexing='ij')
#     c_m1 = location[i, :]
#     c_m2 = location[j, :]
#     v_m1 = orientation[i, :]
#     v_m2 = orientation[j, :]
#     v_mm = location[jj, :] - location[ii, :]
#     sites1 = sites[i, :, :]
#     sites2 = sites[j, :, :]
    
#     #SetPars and parameter set-up:
#     cross.pars['c2_site']=c_m2
#     cross.pars['c1_site']=c_m1
#     cross.pars['sites1']=sites1
#     cross.pars['v_m1']=v_m1
#     cross.pars['v_mm']=v_mm
#     cross.pars['cA']=cA
#     cross.pars['cX']=cX
#     cross.pars['v_m2']=v_m2
#     cross.pars['sites2']=sites2
#     cross.SetPars('CH3CH3')
    
#     r_cross = cross.CalcRate('M1z', 'M2z')
#     r_auto = cross.CalcRate('M1z', 'M1z')
#     rate_matrix=numpy.zeros((num,num))
#     i,j=numpy.arange(num),numpy.arange(num)
#     ii,jj=numpy.meshgrid(i,j,indexing='ij')
#     rate_matrix[ii,jj]=r_cross[ii,jj]
#     rate_matrix[i,j]=numpy.sum((r_auto[ii,jj]-r_auto[ii,ii]),axis=1)+r_auto[i,j]
#     exp_matrix=expm(-rate_matrix[:,:]*t_point)
    
#     #DEBUGGING:
#     numpy.savetxt('output/log/test_cross',r_cross)
#     numpy.savetxt('output/log/test_auto',r_auto)
#     numpy.savetxt('output/log/test_exp',exp_matrix)