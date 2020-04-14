#!/usr/bin/python

import sys,os,re
from commAux import CommAux
from datetime import datetime
import numpy
import matplotlib.pyplot as plt
import random
from scipy.linalg import expm
from timeit import default_timer
from scipy.stats import pearsonr
from scipy.stats import spearmanr

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




def DoCH3_ExtM(macro=True, ax=True, xx=True, cx=True, ca=True, sparse=True, sym=True, irreducible=False, num=False):
    #First set up the spin system and then cross-correlate:
    print 'Calculating methyl-methyl rates'
    auto = CommAux('CH3', latex_file='eps/auto_latex.tex', sparse=sparse, sym=sym)
    cross = CommAux('CH3CH3', latex_file='eps/cross_latex.tex', sparse=sparse, sym=sym)
    
    if(ax):
        auto.GetDip('H2', 'C1', 'd1', tag ='single')
        auto.GetDip('H3', 'C1', 'd1', tag ='single')
        auto.GetDip('H4', 'C1', 'd1', tag ='single')
        
        cross.GetDip('H2', 'C1', 'd1', tag ='single')
        cross.GetDip('H3', 'C1', 'd1', tag ='single')
        cross.GetDip('H4', 'C1', 'd1', tag ='single')
        cross.GetDip('H6', 'C5', 'd2', tag ='single')
        cross.GetDip('H7', 'C5', 'd2', tag ='single')
        cross.GetDip('H8', 'C5', 'd2', tag ='single')
    if(xx):
        auto.GetDip('H2', 'H3', 'e1', tag ='single')
        auto.GetDip('H2', 'H4', 'e1', tag ='single')
        auto.GetDip('H3', 'H4', 'e1', tag ='single')

        cross.GetDip('H2', 'H3', 'e1', tag ='single')
        cross.GetDip('H2', 'H4', 'e1', tag ='single')
        cross.GetDip('H3', 'H4', 'e1', tag ='single')
        cross.GetDip('H6', 'H7', 'e2', tag ='single')
        cross.GetDip('H6', 'H8', 'e2', tag ='single')
        cross.GetDip('H7', 'H8', 'e2', tag ='single')
    if(cx):
        auto.GetCSA('H2', 'cX1', tag='single')
        auto.GetCSA('H3', 'cX1', tag='single')
        auto.GetCSA('H4', 'cX1', tag='single')
        
        cross.GetCSA('H2', 'cX1', tag='single')
        cross.GetCSA('H3', 'cX1', tag='single')
        cross.GetCSA('H4', 'cX1', tag='single')
        cross.GetCSA('H6', 'cX2', tag='single')
        cross.GetCSA('H7', 'cX2', tag='single')
        cross.GetCSA('H8', 'cX2', tag='single')
    if(ca):
        auto.GetCSA('C1', 'cA1', tag='static')
        cross.GetCSA('C1', 'cA1', tag='static')
        
        cross.GetCSA('C5', 'cA2', tag='static')

    index = numpy.arange(3)
    for i in (2,3,4):
        i_str = 'H'+str(i)
        for j in (6,7,8):
            j_str = 'H'+str(j)
            cross.GetDip(i_str, j_str, 'm', tag='double')
    
    auto.macro=macro
    auto.CrossCorrHam()
    
    cross.macro=macro
    cross.CrossCorrHam()
    
    #LOAD UP THE METHYLS FROM THE PDB FILE:
    location,orientation,sites,map_key_xrd,map_index_xrd,residues=get_pdb_coord_methyl('1EZA.pdb')
    
    #LOAD UP THE PEAK LIST FILE:
    cross_matrix,reflected_matrix,d1_matrix,d2_matrix,map_key_nmr,map_index_nmr,mask_matrix=peak_file_parser('correlate.3')
    
    #DEBUGGING:
    print cross_matrix.shape
    print location.shape
    for key in map_key_xrd:
        if key not in map_key_nmr:
            print 'KEY = ',key,' IS NOT IN NMR DATA'
            
    for key in map_key_nmr:
        if key not in map_key_xrd:
            print 'KEY = ',key,' IS NOT IN XRD DATA'
    
    start=True
    for key in map_key_xrd:
        if start:
            start=False
            test=key
        if(int(re.search('\d+',key).group())>int(re.search('\d+',test).group())):
            test=key
    print 'MAXIMUM AMINO ACID VALUE FOR XRD DATA IS = ',test
           
    test=''
    start=True
    for key in map_key_nmr:
        if start:
            start=False
            test=key
        try:
            if(int(re.search('\d+',key).group())>int(re.search('\d+',test).group())):
                test=key
        except AttributeError:
            print 'ATTRIBUTE ERROR, key = ',key,' test = ',test
            
    print 'MAXIMUM AMINO ACID VALUE FOR NMR DATA IS = ',test
    
    #Now we have the spin system set-up, put into it the numerical values:
    cross.tm=0.5*10E-12
    cross.tc=2.5*10E-9
    t_point = 400*1E-3
    cross.SetFreq(600)  #MHz. Set spectrometer proton frequency
    cross.Assemble(('H2z','H3z','H4z'),'M1z')
    cross.Assemble(('H6z','H7z','H8z'),'M2z')
    
    cA=20.
    cX=1.
    
    #Iterate over all of the pairings:
    #So each methyl in "methyls" will consists of 2 3-vectors - one giving position and one giving orientation
    #These three variables hold - r_single holds auto-relaxation of methyl i on it's own, r_auto holds auto-relaxation of i with j present; r_cross holds cross-relaxation from j to i.
    
    #rate_matrix is an n*n rate matrix
    num = location.shape[0]
    rate_matrix = numpy.zeros((num, num))
    dist_matrix = numpy.zeros((num, num))
    
    i = numpy.arange(num)
    j = numpy.arange(num)
    ii, jj = numpy.meshgrid(i, j,indexing='ij')
    c_m1 = location[i, :]
    c_m2 = location[j, :]
    v_m1 = orientation[i, :]
    v_m2 = orientation[j, :]
    v_mm = location[jj, :] - location[ii, :]
    sites1 = sites[i, :, :]
    sites2 = sites[j, :, :]
    
    
    #SetPars and parameter set-up:
    cross.pars['c2_site']=c_m2
    cross.pars['c1_site']=c_m1
    cross.pars['sites1']=sites1
    cross.pars['v_m1']=v_m1
    cross.pars['v_mm']=v_mm
    cross.pars['cA']=cA
    cross.pars['cX']=cX
    cross.pars['v_m2']=v_m2
    cross.pars['sites2']=sites2
    cross.SetPars('CH3CH3')
    
    r_cross = cross.CalcRate('M1z', 'M2z')
    r_auto = cross.CalcRate('M1z', 'M1z')
    
    #DEBUGGING:
    print 'r_cross = \n',r_cross.shape,'\n',r_cross
    print 'r_auto = \n',r_auto.shape,'\n',r_auto
    
    num_samples=800
    num_methyls = location.shape[0]
    num_residues = len(residues)
    result_matrix=numpy.zeros((num_methyls,num_methyls,num_samples))
    
    
    for rep in numpy.arange(num_samples):
        list_methyls = methyl_list_generator(residues)
        index_list = numpy.array([map_key_xrd[x] for x in list_methyls])
        rate_matrix = numpy.zeros((len(index_list),len(index_list)))
        exp_matrix = numpy.zeros((len(index_list),len(index_list)))
        
        #Have 2 indexing regimes. x,y are for the compacted rate_matrix. This is for the calculation here
        x,y=numpy.arange(len(index_list)),numpy.arange(len(index_list))
        xx,yy=numpy.meshgrid(x,y,indexing='ij')
        
        #i,j are for the expanded matrix. Maps the compacted matrix onto the expanded matrix. 
        i,j=index_list,index_list
        ii,jj=numpy.meshgrid(i,j,indexing='ij')
        
        rate_matrix[xx,yy]=r_cross[ii,jj]
        rate_matrix[x,y]=numpy.sum((r_auto[ii,jj]-r_auto[ii,ii]),axis=1)+r_auto[i,j]
        exp_matrix[:,:] = expm(-rate_matrix[:,:] * t_point)
        result_matrix[ii,jj,rep] = exp_matrix[xx,yy]
    
    sim_matrix = numpy.average(result_matrix[:,:,:], axis=2)
    dist_matrix = numpy.sqrt(numpy.sum(v_mm * v_mm, axis=2))
    
    print 'sim_matrix = \n',sim_matrix
    print 'dist_matrix = \n',dist_matrix
    
    #CAST DATA SETS INTO THE SAME INDEXING SCHEME. Choose NMR data indexing scheme as the smaller of the two.
    xrd,xrd_ab=[],[]
    nmr,nmr_ab=[],[] 
    #So iterate over the indices for the NMR peak list:
    for key in map_key_nmr:
        if(key[0:3]=='ALA' or key[0:3]=='ILE'):
            xrd.append(map_key_xrd[key])
            nmr.append(map_key_nmr[key])
        if(key[0:3]=='LEU' or key[0:3]=='VAL'):
            xrd_ab.append((map_key_xrd[key],map_key_xrd[string_ab_swap(key)]))
            nmr_ab.append((map_key_nmr[key],map_key_nmr[string_ab_swap(key)]))
    
    #NOW CONVERT INDICES FROM XRD TO NMR:
    xrd_indices,nmr_indices=[],[]
    for key in map_key_nmr:
        xrd_indices.append(map_key_xrd[key])
        nmr_indices.append(map_key_nmr[key])
    
    #Create square matrices which map residues we've picked above onto either the NMR or XRD data:
    xrd,nmr=numpy.array(xrd),numpy.array(nmr)
    xrd_ab,nmr_ab=numpy.array(xrd_ab),numpy.array(nmr_ab)
    
    
    xrd_indices,nmr_indices=numpy.array(xrd_indices),numpy.array(nmr_indices)
    print 'xrd_indices = ',xrd_indices
    print 'nmr_indices = ',nmr_indices
    ii_xrd,jj_xrd=numpy.meshgrid(xrd_indices,xrd_indices,indexing='ij')
    ii_nmr,jj_nmr=numpy.meshgrid(nmr_indices,nmr_indices,indexing='ij')
    
    #MAP THE MASK FOR MISSING DATA DOWN ONTO THE NMR DATA
    mask=numpy.zeros(ii_nmr.shape)
    mask[ii_nmr,jj_nmr]=mask_matrix
    mask=numpy.logical_not(mask)
    print 'THERE ARE = ',numpy.sum(mask),' NMR PEAKS'
    
    #MAP XRD DATA ONTO THE NMR DATA:
    xrd_peaks=numpy.zeros(ii_nmr.shape)
    nmr_peaks=numpy.zeros(ii_nmr.shape)
    dist=numpy.zeros(ii_nmr.shape)
    
    nmr_peaks[ii_nmr,jj_nmr]=cross_matrix
    xrd_peaks[ii_nmr,jj_nmr]=sim_matrix[ii_xrd,jj_xrd]
    dist[ii_nmr,jj_nmr]=dist_matrix[ii_xrd,jj_xrd]
    
    #CALCULATE CORRELATION BETWEEN DATA ROWS, AND CHECK A/B ASSIGNMENTS:
    #ITERATE OVER THE INDICES WHICH APPLY TO RESIDUES WITH A/B POSSIBILITIES:
    i,j=numpy.arange(nmr_peaks.shape[0]),numpy.arange(nmr_peaks.shape[0])
    
    print 'nmr and xrd data: \n'
    print nmr_peaks[0:16,0:16][mask[0:16,0:16]]
    print xrd_peaks[0:16,0:16][mask[0:16,0:16]]
    
    for i_n in numpy.arange(nmr_ab.shape[0]):
        i_a,i_b=nmr_ab[i_n,0],nmr_ab[i_n,1]
        corr_1=pearsonr(nmr_peaks[i_a,:][mask[i_a,:]],xrd_peaks[i_a,:][mask[i_a,:]])[0]+pearsonr(nmr_peaks[i_b,:][mask[i_b,:]],xrd_peaks[i_b,:][mask[i_b,:]])[0]
        corr_2=pearsonr(nmr_peaks[i_a,:][mask[i_a,:]],xrd_peaks[i_b,:][mask[i_a,:]])[0]+pearsonr(nmr_peaks[i_b,:][mask[i_b,:]],xrd_peaks[i_a,:][mask[i_b,:]])[0]
#         print i_a,i_b,corr_1,corr_2
#         print nmr_peaks[i_a,:][mask[i_a,:]]
#         print xrd_peaks[i_a,:][mask[i_a,:]]
#         print nmr_peaks[i_b,:][mask[i_b,:]]
#         print xrd_peaks[i_b,:][mask[i_b,:]]
#         print '\n'
        if corr_2>corr_1:
            print 'NEED TO SWAP ',i_a,' AND ',i_b
            i[i_a],i[i_b]=i_b,i_a
            j[i_a],j[i_b]=i_b,i_a
    ii,jj=numpy.meshgrid(i,j,indexing='ij')
    nmr_peaks=nmr_peaks[ii,jj]
    xrd_peaks=xrd_peaks[ii,jj]
    dist=dist[ii,jj]
    print 'ii = \n',ii
    print 'jj = \n',jj
    print 'TOTAL CORRELATION = ',pearsonr(nmr_peaks[mask].flatten(),xrd_peaks[mask].flatten())
    plt.scatter(nmr_peaks[mask].flatten(),xrd_peaks[mask].flatten())
    plt.show()
    
    print 'REPEAT - SEE IF WE DO ANY MORE SWAPS: '
    i,j=numpy.arange(nmr_peaks.shape[0]),numpy.arange(nmr_peaks.shape[0])
    for i_n in numpy.arange(nmr_ab.shape[0]):
        i_a,i_b=nmr_ab[i_n,0],nmr_ab[i_n,1]
        corr_1=pearsonr(nmr_peaks[i_a,:][mask[i_a,:]],xrd_peaks[i_a,:][mask[i_a,:]])[0]+pearsonr(nmr_peaks[i_b,:][mask[i_b,:]],xrd_peaks[i_b,:][mask[i_b,:]])[0]
        corr_2=pearsonr(nmr_peaks[i_a,:][mask[i_a,:]],xrd_peaks[i_b,:][mask[i_a,:]])[0]+pearsonr(nmr_peaks[i_b,:][mask[i_b,:]],xrd_peaks[i_a,:][mask[i_b,:]])[0]
#         print i_a,i_b,corr_1,corr_2
#         print nmr_peaks[i_a,:][mask[i_a,:]]
#         print xrd_peaks[i_a,:][mask[i_a,:]]
#         print nmr_peaks[i_b,:][mask[i_b,:]]
#         print xrd_peaks[i_b,:][mask[i_b,:]]
#         print '\n'
        if corr_2>corr_1:
            print 'NEED TO SWAP ',i_a,' AND ',i_b
            i[i_a],i[i_b]=i_b,i_a
            j[i_a],j[i_b]=i_b,i_a
    ii,jj=numpy.meshgrid(i,j,indexing='ij')
    nmr_peaks=nmr_peaks[ii,jj]
    xrd_peaks=xrd_peaks[ii,jj]
    dist=dist[ii,jj]
    print 'ii = \n',ii
    print 'jj = \n',jj
    
    

    
    
    
def string_ab_swap(string):
    if string[-1]=='A':
        return string[:-1]+'B'
    elif string[-1]=='B':
        return string[:-1]+'A'    
    
def peak_file_parser(peak_file_dir):
    file = open(peak_file_dir, 'r')
    cross,reflected={},{}
    diag_1,diag_2={},{}
    for line in file:
        data = line.split()
        res1 = residue_parser(data[0])
        res2 = residue_parser(data[1])
        v1,v2=data[6],data[7]
        a1,a2=data[14],data[15]
        if not(res1 in cross):
            cross[res1]={}
            reflected[res1]={}
            diag_2[res1]={}
            diag_1[res1]={}
        cross[res1][res2]=v1
        reflected[res1][res2]=v2
        diag_1[res1][res2]=a1
        diag_2[res1][res2]=a2
    map_key,map_index={},{}
    i=0
    for res1 in cross:
        map_key[res1]=i
        map_index[i]=res1
        i+=1
    
    num = len(cross)
    cross_matrix, reflected_matrix, d1_matrix, d2_matrix = numpy.zeros((num,num)),numpy.zeros((num,num)),numpy.zeros((num,num)),numpy.zeros((num,num))
    mask_matrix=numpy.ones((num,num),dtype=bool)
    for res1 in cross:
        for res2 in cross[res1]:
            i1,i2=map_key[res1],map_key[res2]
            cross_matrix[i1,i2]=cross[res1][res2]
            reflected_matrix[i1,i2]=reflected[res1][res2]
            d1_matrix[i1,i2]=diag_1[res1][res2]
            d2_matrix[i1,i2]=diag_2[res1][res2]
            mask_matrix[i1,i2]=False
    
    return cross_matrix,reflected_matrix,d1_matrix,d2_matrix,map_key,map_index,mask_matrix
    
        
def residue_parser(residue_str):
    res_letter = residue_str[0]
    integers = [int(s) for s in re.findall(r'-?\d+\.?\d*', residue_str)]
    out_str = ''
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
    
def methyl_list_generator(residue_list):
    #EXAMPLES OF RESIDUE NAMES:
    #'ILE   57', 'LEU  123', 'VAL  176'
    methyl_list = []
    for residue in residue_list:
        if(residue[0:3]=='ILE' or residue[0:3]=='ALA'):
            methyl_list.append(residue+' A')
        if(residue[0:3]=='VAL' or residue[0:3]=='LEU'):
            ab = random.choice([' A', ' B'])
            methyl_list.append(residue+ab)
    return methyl_list


def CalcMetGeom(c_m, c_a, c_b, r_m=1.09, theta_m=numpy.arccos(1./3.)):
    #This function takes the position of three carbon atoms - cm_pos is position of methyl carbon; ca_pos, cb_pos are position of the other 2 carbon atoms
    #This function returns the position of the three hydrogen atoms.
    c_m = numpy.array(c_m)
    c_a = numpy.array(c_a)
    c_b = numpy.array(c_b)
    
    v_m = numpy.array(c_m - c_a)
    v_m = v_m/numpy.sqrt(numpy.dot(v_m,v_m))
    #Calculate phi_rot, theta_rot terms. These terms are the terms which do the rotation of the methyl group axes, IE to re-orientate the methyl groups
    r = numpy.sqrt(numpy.dot(v_m, v_m))
    theta_rot = numpy.arccos(v_m[2]/r)
    phi_rot = numpy.arctan2(v_m[0], v_m[1])
    
    #This calculate psi angle. First calculate the rotation matrix to re-orientate methyl group:
    rot = numpy.dot(DoRot((0,0,1), phi_rot), DoRot((1,0,0), theta_rot))
    inv_rot = numpy.linalg.inv(rot)
    v_a = c_b - c_a
    
    v_a = numpy.dot(inv_rot, v_a)
    psi_i = numpy.arctan2(v_a[0,0],v_a[0,1]) + numpy.pi/3.
    
    sites=3 #Keep below code generalized so we can change stuff up for diffusive model/n-site model if we want later
    delta_psi = 2. * numpy.pi/sites
    psi = psi_i + (delta_psi * numpy.linspace(0, sites-1, sites))

    #Calculate these factors as are re-used (hang-over from when these were used many, many times in random walk simulations)
    sin_psi, cos_psi = r_m * numpy.sin(psi), r_m * numpy.cos(psi)
    sin_theta_m, cos_theta_m = numpy.sin(theta_m), numpy.cos(theta_m)



    #Calculate the Cartesian coordinates of the methyl hydrogens. Remember - methyl hydrogens are in a ring, but using C as ccoordinate centre, *not* centre of mass
    sites = numpy.matrix([numpy.array(sin_psi) * sin_theta_m, numpy.array(cos_psi) * sin_theta_m, numpy.full(len(sin_psi), r_m*cos_theta_m)])
    sites = numpy.dot(rot, sites)        
    
    #This displaces the methyl group by offset vector.
    sites[0, :]+=c_m[0]
    sites[1, :]+=c_m[1]
    sites[2, :]+=c_m[2]
    
    methyl = {'loc':c_m, 'orientation':v_m, 'sites':numpy.transpose(numpy.array(sites))}
    
    return methyl

def DoRot(k,phi):
    from scipy import mat,zeros
    k = k / numpy.sqrt(numpy.dot(k, k))
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
    
def get_pdb_coord_methyl(pdb_file):
    #Initialize stuff - Set which file to open; initilaize lists for all atoms in file, and then result list - methyl_list
    print pdb_file
    file = open(pdb_file, 'r')
    atom_list = []
    methyl_list = {}
    
    #Iterate over lines in file, and fill atom_list with position, type & residue (type & number) for each atom:
    for line in file:
        if(line[0:4]=='ATOM'):
#             print line
            atom_number = line[6:12]
            atom_name = line[13:16]
            residue_name = line[17:20]
            residue_number = int(line[22:26].strip())
            x_coord = float(line[30:38])
            y_coord = float(line[38:46])
            z_coord = float(line[46:54])
            
            atom = {}
            atom['num']  = atom_number
            atom['ele'] = atom_name
            atom['residue'] = residue_name + ' ' + str(residue_number)
            atom['x'], atom['y'], atom['z'] = x_coord, y_coord, z_coord
            atom_list.append(atom)

#                 print 'ATOM = ', atom_name
#                 print 'RESIDUE = ', residue_name
#                 print 'X = ', x_coord
#                 print 'Y = ', y_coord
#                 print 'Z = ', z_coord
#                 print line
    
    #Now iterate over all atoms in atom_list and group in residue (type & number):
    residue_list = {}
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
    residues = []
    count_ile, count_leu, count_val,count_ala = 0,0,0,0
    for key in residue_list:
        residue = residue_list[key]
        if(key[0:3]=='ILE'):
            residues.append(key)
            methyl = CalcMetGeom(residue['CD1'], residue['CG1'], residue['CB '])
            methyl_list[key+' A'] = methyl
            count_ile+=1
        elif(key[0:3]=='LEU'):
            residues.append(key)
            methyl = CalcMetGeom(residue['CD1'], residue['CG '], residue['CB '])
            methyl_list[key + ' A'] = methyl
            methyl = CalcMetGeom(residue['CD2'], residue['CG '], residue['CB '])
            methyl_list[key + ' B'] = methyl
            count_leu+=1
        elif(key[0:3]=='VAL'):
            residues.append(key)
            methyl = CalcMetGeom(residue['CG1'], residue['CB '], residue['CA '])
            methyl_list[key + ' A'] = methyl
            methyl = CalcMetGeom(residue['CG2'], residue['CB '], residue['CA '])
            methyl_list[key + ' B'] = methyl
            count_val+=1
        elif(key[0:3]=='ALA'):
            residues.append(key)
            methyl=CalcMetGeom(residue['CB '],residue['CA '],residue['C  '])
            methyl_list[key+' A']=methyl
            count_ala+=1

    #Want to return: mapping dictionaries from methyl residue name & number (EG. LEU 113 A) to numpy index, and index to name.
    #3 numpy arrays - N*3 array containing coordinates of methyl C (c_m), N*3 for orientation (v_m), and N*3*3 positions of methyl H (sites)
#     methyl = {'loc':c_m, 'orientation':v_m, 'sites':numpy.array(sites)}
    num = len(methyl_list)
    location = numpy.zeros((num, 3))
    orientation = numpy.zeros((num, 3))
    sites = numpy.zeros((num, 3, 3))
    map_key = {}
    map_index = {}
    i = 0 #Index which incremented to slot methyls into numpy arrays
    print 'count_leu = ', count_leu
    print 'count_val = ', count_val
    print 'count_ile = ', count_ile
    print 'count_ala = ', count_ala
    for key in methyl_list:
        location[i, :] = methyl_list[key]['loc']
        orientation[i, :] = methyl_list[key]['orientation']
        sites[i, :, :] = methyl_list[key]['sites']
        map_key[key] = i
        map_index[i] = key
        i+=1
        
    return location, orientation, sites, map_key, map_index, list(set(residues))
    

start=datetime.now()
DoCH3_ExtM(macro=False, ax=True, xx=True, cx=True, ca=True, sparse=True, sym=True, irreducible=False, num=False)
print 'TIME TAKEN = ', datetime.now()-start
