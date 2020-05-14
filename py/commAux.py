#!/usr/bin/python

import sys,os,re
import string,copy,numpy
from math import sqrt,fabs
import math
import numpy
from string import digits
from fractions import Fraction,gcd

from scipy import mat,zeros
from basis import comm,HermConj,compare,OuterProd,OuterProdSparse,InitBasis,Trace,TraceSparse,IsZero,IsZeroSparse,GetBasis
from scipy.linalg import expm
from scipy import linalg
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist,pdist,squareform

from numpy import tan as Tan
from numpy import cos as Cos
from numpy import sin as Sin
from numpy import sqrt as Sqrt
from numpy import pi as Pi
from numpy import arccos as ACos
from numpy import arctan2 as ATan2
from numpy import arcsin as ASin
from numpy import average as Ave
from numpy import exp as Exp

from datetime import datetime

import matplotlib.pyplot as plt


#from array import array
#import operator

class CommAux():
    def __init__(self,typ,latex_file='output/eps/latex.tex',sparse=True,sym=True,latex=True):
        if(os.path.exists('redfield')==0):
            print 'Warning: No libraries found'
        if(os.path.exists('output')==0):
            os.system('!mkdir output')
        if(os.path.exists('output/pdf')==0):
            os.system('!mkdir output/pdf')
            os.makedirs('output/pdf')
        if(os.path.exists('output/eps')==0):
            os.system('!mkdir output/eps')
            os.makedirs('output/eps')
        if(os.path.exists('output/log')==0):
            os.system('!mkdir output/log')
            os.makedirs('output/log')
        if(os.path.exists('output/log/log_Aux')):
            os.remove('output/log/log_Aux')
        self.ncpus=2
        
        print_log=open('output/log/log_Aux','w')
        print_log.write('THIS IS A LOG FILE WHICH SAYS LOTS OF COOL STUFF: \n')
        #DS hello
        self.GNUFILE='\'C:\Program Files\gnuplot\bin\gnuplot.exe\''
        
        self.SPARSE=sparse  #use scipy sparse routines
        self.SYM=sym        #use symbolic representation of matricies
        self.LATEX=latex    #make latex outputs?
        self.POINT=False #pointless. and wrong.

        self.EXT=False #is there an external proton in the model? will get set automatically.
        self.MET=False #Is there an external methyl group to worry about ??
        self.CYCLIC=False #when assembling symmetric basis, make one spin the focus
        self.SQUARE=True #complete the square on TROSY results
        
        self.baseHalf,self.baseOne=InitBasis(sparse=sparse)  #get basis
        self.baseTag=typ
        self.GetBaseType()  #get the basetype (CHHH..)
        #sys.stdout=open('log/'+self.baseTag+'.log','w') #redirect log.

        print_log.write('Running RelCalc')
        print_log.write('---------------')
        print_log.write('Type:     '+typ)
        print_log.write('baseType: '+self.baseType)
        print_log.close()
        self.baseSize=len(self.baseType)
        self.en=int(re.findall(r'\d+',typ)[0]) #number of protons

        #all sorted by coherence 
        self.base={}   #will store matrix representation
        self.baseOp={} #will store expanded coherence representation

        self.Aint=[] #list of terms from interacting Hamiltonian
        self.Atyp=[] #list types of interacting Hamiltonian and whether or not they are rotating.
        self.resSet={} #will store final symbolic results
        self.resCalc={}
        self.consts={} #will store constants
        self.pars={}   #This is basically a transfer room for variables from commTest
        self.Cint={}

        self.edict={}  #relates dipolar constants to AX distance

        self.macro=False #flag for macromolecular limit
        self.kay=False   #flag for Kay and Bull representation
        self.taum=0      #flag for methyl rotation

        #factors to map spin operator to irreducible basis
        #used when setting up Hamiltonians
        self.o={}
        self.o[0]={}
        self.o[0][0]= (1./2.)*numpy.sqrt(8./3.) #for zz terms#
        self.o[0][1]=-(1./8.)*numpy.sqrt(8./3.) #for +- and -+
        self.o[1]= (1./2.)                      #for + and -
        self.o[2]= (1./2.)                      #for ++ and --

        #constants
        self.mu0d4pi= 1E-7    #Tm-1A-1
        self.hbar=  1.05457173E-34 #m2 kg s-1
        self.gammaH= 2.675E8   #rad s-1 T-1
        self.gammaC= 6.72828E7   #rad s-1 T-1
        self.gammaF= 2.51662E8 #rad s-1 T-1
        self.gammaD= 4.1065E7 #rad s-1 T-1

        self.sfrq=600.   #spectrometer frequency in MHz
        self.SetFreq(self.sfrq) #set spectrometer frequencies
        self.tc=1E-8  #initialise rotational correlation time
        self.tm=1E-12 #initialise symmetry axis correlation time
        self.Woessner=True   #Woessner 1/0
        self.n=3          #symmetry axis number
        #self.Woessner=1
        
        #don't forget to set frequency
        #set calculate dipolar and csa constants
        #Plib['P2']='3sin^2b' 
        #Plib['P1']='3sinbcosb' 
        #Plib['P0']='1/2(3cos^2b-1)'


        #paths:
        self.LATEXFILE=latex_file

        #choose which functions to use if sparse
        if(self.SPARSE):   
            self.IsZero=IsZeroSparse
            self.Trace=TraceSparse
            self.OuterProd=OuterProdSparse
        else:  #don't use sparse routines
            self.IsZero=IsZero
            self.Trace=Trace
            self.OuterProd=OuterProd
            
        #if symbolic, then read the libraries
        if(self.SYM):
            self.ReadCommDictSingle() #read the two spin commutator index
            self.ReadTraceDict()      #read the trace dictionary

        #if making latex outputs, initialise output
        if(self.LATEX):
            self.InitLatex() #initialise latex output

    #expand tag to get full basis
    #creates self.baseType and self.baseList
    def GetBaseType(self):
        self.baseType=''
        for i in range(len(self.baseTag)): #loop over specification
            if(self.baseTag[i] in '1234567890'): #is current entry a digit?
                if(self.baseTag[i-1] in '1234567890'): #is the last element a digit?
                    for j in range(int(self.baseTag[i-1]+self.baseTag[i])-1):
                        self.baseType+=self.baseTag[i-2]
                else:
                    for j in range(int(self.baseTag[i])-1):
                        self.baseType+=self.baseTag[i-1]
            else:
                self.baseType+=self.baseTag[i]
        self.baseList=numpy.array(list(self.baseType))

    #set spectrometer frequency and angular frequencies of various nuclei
    def SetFreq(self,sfrq):
        self.sfrq=sfrq
        self.omegaH=2*numpy.pi*self.sfrq*1E6  #proton spectrometer frequency in rad s-1
        self.B0=self.omegaH/self.gammaH       #B field T
        self.omegaC=2*numpy.pi*self.sfrq*1E6/self.gammaH*self.gammaC  #carbon spectrometer frequency in rad s-1
        self.omegaF=2*numpy.pi*self.sfrq*1E6/self.gammaH*self.gammaF  #carbon spectrometer frequency in rad s-1
        self.omegaD=2*numpy.pi*self.sfrq*1E6/self.gammaH*self.gammaD  #carbon spectrometer frequency in rad s-1
    
    
    
    def DoRot(self,k,phi):
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
    
    
    
    #get relevant gyromagnetic ratios
    def GetNucGamma(self,n):
        if(n=='H'):
            return self.gammaH
        if(n=='C'):
            return self.gammaC
        if(n=='D'):
            return self.gammaD
        if(n=='F'):
            return self.gammaF
        print 'no idea which gamma. add this to the lookup table:'
        print n
        sys.exit(100)

    def CalcConstsDip(self, tag, n1, n2, sites1, sites2, typ='single', transpose=False):
        num = sites1.shape[0]
        g1=self.GetNucGamma(n1)
        g2=self.GetNucGamma(n2)
        self.consts[tag]=self.mu0d4pi*self.hbar*g1*g2/((1.E-10)**3.)
        sites = numpy.zeros((num, num, 3, 3))
        xx,yy,zz = numpy.meshgrid(numpy.arange(num), numpy.arange(num), numpy.arange(3))
        sites[xx,yy,zz,:] = sites1[xx,zz,:] - sites2[xx,zz,:] #So take difference between vectors within *same* methyl group (indexed by xx). zz indixes 3 possible sites
        if(transpose):
            self.Cint[tag]=typ,numpy.swapaxes(sites,0,1)
        else:
            self.Cint[tag]=typ,sites

    def CalcConstsCSA(self, tag, n1, csa, sites1, typ='single', transpose=False):
        self.consts[tag]=self.GetNucGamma(n1)*self.sfrq/self.gammaH*2*numpy.pi*csa/3.
        if(typ=='single'):
            num = sites1.shape[0]
            sites = numpy.zeros((num, num, 3, 3))
            xx,yy,zz = numpy.meshgrid(numpy.arange(num), numpy.arange(num), numpy.arange(3), indexing='ij')
            norm = numpy.sqrt(numpy.sum(sites1 * sites1, axis=2))
            sites=numpy.divide(sites1[xx,zz,:], norm[xx,zz, numpy.newaxis])
        elif(typ=='static'):
            num = sites1.shape[0]
            sites = numpy.zeros((num, num, 3))
            xx,yy = numpy.meshgrid(numpy.arange(num), numpy.arange(num))
            norm = numpy.sqrt(numpy.sum(sites1 * sites1, axis=1))
            sites[xx,yy,:]=numpy.divide(sites1[xx,:], norm[xx, numpy.newaxis])
        if(transpose):
            self.Cint[tag]=typ, numpy.swapaxes(sites,0,1)
        else:
            self.Cint[tag]=typ, sites
    
    def CalcConstsMet(self, tag, n1, n2, sites1, sites2, typ='double'):
        g1=self.GetNucGamma(n1)
        g2=self.GetNucGamma(n2)
        self.consts[tag]=self.mu0d4pi*self.hbar*g1*g2/((1.E-10)**3.)
        num = sites1.shape[0]
        sites = numpy.zeros((num,num,3,3,3))
        xx,yy,uu,vv = numpy.meshgrid(numpy.arange(num),numpy.arange(num),numpy.arange(3),numpy.arange(3))
        sites[xx,yy,uu,vv,:] = sites1[xx,uu,:] - sites2[yy,vv,:] #Take difference between vectors within *different* methyl groups (indexed by xx,yy). uu,vv index different configuration of each methyl group
        self.Cint[tag]=typ, sites
                
    def CalcConsts(self, tag_k, tag_l):
        #This is a new function which will calculate ensemble-averaged correlation function values.
        #Will need different cases for different types of k,l interactions - whether "static", "single", "double"
        #Will also need the interaction values which the interactions can take
        #Then we can ensemble average over the start values. And then the end values depending on the population factors.
        #Which would give us a correlation function
        #We want S^2 factors - so go do the maths.
        typ_k, site_k = self.Cint[tag_k]
        typ_l, site_l = self.Cint[tag_l]
        tag = self.GetCont(tag_k, tag_l)
        
        num = site_k.shape[0]
        C0 = numpy.zeros((num, num))
        S2 = numpy.zeros((num, num))
        x,y,u,v,w,z=numpy.arange(num),numpy.arange(num),numpy.arange(3),numpy.arange(3),numpy.arange(3),numpy.arange(3)
        if(typ_l=="static"):
            if(typ_k=="static"):
                xx, yy = numpy.meshgrid(x,y, indexing='ij')
                C0 = self.CalcZeta(site_k[xx,yy,:], site_l[xx,yy,:])
                self.consts['C('+tag+')']= numpy.nan_to_num(C0)
            elif(typ_k=="single"):
                xx, yy, uu = numpy.meshgrid(x,y,u, indexing='ij')
                chi = numpy.zeros((num, num, 3))
                chi[xx,yy,uu]=self.CalcZeta(site_k[xx,yy,uu,:], site_l[xx,yy,:])
                xx,yy=numpy.meshgrid(x,y,indexing='ij')
                C0 = 1./3. * numpy.sum(chi[xx,yy,:], axis=2)
                self.consts['C('+tag+')'] = numpy.nan_to_num(C0)
            elif(typ_k=="double"):
                xx, yy, uu, vv = numpy.meshgrid(x,y,u,v, indexing='ij')
                chi = numpy.zeros((num, num, 3, 3))
                chi[xx,yy,uu,vv] = self.CalcZeta(site_k[xx,yy,uu,vv,:], site_l[xx,yy,:])
                xx,yy=numpy.meshgrid(x,y,indexing='ij')
                C0 = 1./9. * numpy.sum(chi[xx,yy,:,:], axis=(2,3))
                self.consts['C('+tag+')'] = numpy.nan_to_num(C0)
        elif(typ_l=="single"):
            if(typ_k=="static"):
                xx,yy,uu = numpy.meshgrid(x,y,u, indexing='ij')
                chi = numpy.zeros((num, num, 3))
                chi[xx,yy,uu] = self.CalcZeta(site_k[xx,yy,:], site_l[xx,yy,uu,:])
                C0 = 1./3. * numpy.sum(chi, axis=2)
                self.consts['C('+tag+')'] = numpy.nan_to_num(C0)
            elif(typ_k=="single"):
                xx,yy,uu,vv = numpy.meshgrid(x,y,u,v, indexing='ij')
                chi = numpy.zeros((num, num, 3, 3))
                chi[xx,yy,uu,vv] = self.CalcZeta(site_k[xx,yy,uu,:], site_l[xx,yy,vv,:])
                xx,yy,uu = numpy.meshgrid(x,y,u, indexing='ij')
                C0 = 1./3. * numpy.sum(chi[xx,yy,uu,uu], axis=2) #Sum over the diagonal elements of chi
                S2 = 1./9. * numpy.sum(chi, axis=(2,3)) #Sum over *all* elements of chi
                self.consts['C('+tag+')'] = numpy.nan_to_num(C0)
                self.consts['S('+tag+')^2'] = numpy.nan_to_num(S2/C0)
                self.consts['(1-S('+tag+')^2)'] = numpy.nan_to_num(1. - S2/C0)
            elif(typ_k=="double"):
                chi = numpy.zeros((num,num,3,3,3))
                xx,yy,uu,vv,ww = numpy.meshgrid(x,y,u,v,w, indexing='ij')
                chi[xx,yy,uu,vv,ww] = self.CalcZeta(site_k[xx,yy,uu,vv,:], site_l[xx,yy,ww,:])
                xx,yy,uu,vv = numpy.meshgrid(x,y,u,v, indexing='ij')
                C0 = 1./9. * numpy.sum(chi[xx,yy,uu,vv,uu], axis =(2,3)) #*Should* be a sum over the diagonal elements of chi
                S2 = 1./27. * numpy.sum(chi, axis = (2,3,4)) #Sum over *all* the elements of chi
                self.consts['C('+tag+')'] = numpy.nan_to_num(C0)
                self.consts['S('+tag+')^2'] = numpy.nan_to_num(S2/C0)
                self.consts['(1-S('+tag+')^2)'] = numpy.nan_to_num(1. - S2/C0)
        elif(typ_l=="double"):
            if(typ_k=="static"):
                xx,yy,uu,vv = numpy.meshgrid(x,y,u,v, indexing='ij')
                chi = numpy.zeros((num, num, 3, 3))
                chi[xx, yy, uu, vv] = self.CalcZeta(site_k[xx,yy,:], site_l[xx,yy,uu,vv,:])
                C0 = 1./9. * numpy.sum(chi, axis=(2,3)) #Sum over *all* elements of chi
                self.consts['C('+tag+')'] = numpy.nan_to_num(C0)
            elif(typ_k=="single"):
                chi = numpy.zeros((num,num,3,3,3))
                xx,yy,uu,vv,ww = numpy.meshgrid(x,y,u,v,w, indexing='ij')
                chi[xx,yy,uu,vv,ww] = self.CalcZeta(site_k[xx,yy,uu,:], site_l[xx,yy,vv,ww,:])
                xx,yy,uu,vv = numpy.meshgrid(x,y,u,v, indexing='ij')
                C0 = 1./9. * numpy.sum(chi[xx,yy,uu,uu,vv], axis=(2,3)) #Sum over diagonal elements of chi
                S2 = 1./27. * numpy.sum(chi, axis=(2,3,4))
                self.consts['C('+tag+')']= numpy.nan_to_num(C0)
                self.consts['S('+tag+')^2']= numpy.nan_to_num(S2/C0)
                self.consts['(1-S('+tag+')^2)'] = numpy.nan_to_num(1. - S2/C0)
            elif(typ_k=="double"):
                xx,yy,uu,vv,ww,zz = numpy.meshgrid(x,y,u,v,w,z, indexing='ij')
                chi = numpy.zeros((num, num,3,3,3,3))
                chi[xx,yy,uu,vv,ww,zz]=self.CalcZeta(site_k[xx,yy,uu,vv,:], site_l[xx,yy,ww,zz,:])

                chi_0 = numpy.zeros((num, num))
                chi_1 = numpy.zeros((num, num))
                chi_2 = numpy.zeros((num, num))
                
                uu, vv= numpy.meshgrid(u,v, indexing='ij')
                chi_0 = numpy.sum(chi[:,:,uu,vv,uu,vv], axis=(2,3))
                #chi_1 is sum_i sum_j sum_k chi_i,j,i,k + chi_i,j,k,j
                xx,yy,uu,vv,ww = numpy.meshgrid(x,y,u,v,w, indexing='ij')
                chi_1 = numpy.sum(chi[xx,yy,uu,vv,uu,ww], axis=(2,3,4)) + numpy.sum(chi[xx,yy,uu,vv,ww,vv], axis=(2,3,4)) - 2*chi_0
                chi_2 = numpy.sum(chi, axis=(2,3,4,5)) - chi_1 - chi_0
                
                C0 = 1/9. * chi_0
                S_inf  = 1./81. * (chi_0 + chi_1 + chi_2)
                S_tm   = 1./81. * (4.*chi_0 + chi_1 - 2.*chi_2)
                S_tm_2 = 1./81. * (4.*chi_0 - 2.*chi_1 + chi_2)
                
                self.consts['C('+tag+')'] = numpy.nan_to_num(C0)
                self.consts['S_inf('+tag+')^2'] = numpy.nan_to_num(S_inf/C0)
                self.consts['S_tm('+tag+')^2'] = numpy.nan_to_num(S_tm/C0)
                self.consts['S_tm_2('+tag+')^2'] = numpy.nan_to_num(S_tm_2/C0)
        else:
            print_log=open('output/log/log_Aux','a')
            print_log.write("MUPPET. Need type of interaction as static, single, double. Check code & commTest set-up \n")
            print_log.close()
            
    def PolarToCartN(self,v):
        vnew=v.copy()
        vnew[:,0]=v[:,0]*Sin(v[:,1])*Cos(v[:,2])
        vnew[:,1]=v[:,0]*Sin(v[:,1])*Sin(v[:,2])
        vnew[:,2]=v[:,0]*Cos(v[:,1])
        return vnew
    def GetSphPars(self,v1,v2):
        #AA1,BB1,CC1,dzz1=self.GetSphPars(thetaE1,thetaD1)
        AA= v1[0]**2. + v2[0]**2. - 2. *v1[0]* v2[0]*Cos(v1[1]) *Cos(v2[1])
        BBB=2.*v1[0]* v2[0]*Sin(v2[1]) *Sin(v1[1]) 
        BB =   BBB* Cos(v2[2] - v1[2]) 
        CC = - BBB* Sin(v2[2] - v1[2])
        dzz = v1[0]*Cos(v1[1])-v2[0]*Cos(v2[1])  #1-2
        return AA,BB,CC,dzz
    
    def GetSph2(self,cosT,RED,phiT):
        sinT=Sqrt(1-cosT**2.)
        P=self.GetP(cosT,sinT)
        #spherical Harmonics:
        i=complex(0,1.)
        yr={}
        yr[0] =Sqrt(6)*P[0] /RED**3. *complex(1.,0)
        yr[1] =        P[1] /RED**3. * Exp(phiT*i)
        yr[2] =1./2.*  P[2] /RED**3. * Exp(2.*phiT*i) 
        yr[-1]=        P[1] /RED**3. * Exp(-phiT*i)
        yr[-2]=1./2.*  P[2] /RED**3. * Exp(-2.*phiT*i)
        return yr
    
    def GetP(self,cosT,sinT):
        P={}
        P[0]=3./2.*cosT**2.-1/2.
        P[1]=3.*cosT*sinT
        P[2]=3.*sinT**2.
        return P
    #get dipolar spherical harmonic from
    #E-D consants, polar vectors in the ring, and xy co-ordinates of external
    
    def Ydip(self,AA,BB,CC,dzz,phiM,v2,xyE):
        RED = Sqrt(AA - BB*Cos(phiM) - CC*Sin(phiM)) #new value of RED
        cosT=dzz/RED
        xyD=self.GetXY((v2[0],v2[1],v2[2]+phiM)) #get set of cartesian co-ordinates for methyl ring
        phiT=ATan2(xyE[1]-xyD[1],xyE[0]-xyD[0])  #get Phi E1->D1
        return self.GetSph2(cosT,RED,phiT) #get spherical harmonics interaction 1
    
    def GetXY(self,v):
        x=v[0]*Sin(v[1])*Cos(v[2])
        y=v[0]*Sin(v[1])*Sin(v[2])
        #z=v[0]*Cos(v[1])
        return x,y
    
    def OrderCalc(self,typ,yr1,yr2):
        S={}
        if(typ==0):
            S[0]=Ave(yr1[0]*yr2[0]) #<y1(0)y2(0)>
            for key in 1,2: # <y1(q) y2(-q)> +  <y1(-q) y2(q)> 
                S[key]=Ave(yr1[key]*yr2[-key])+Ave(yr1[-key]*yr2[key])
        else:
            S[0]=Ave(yr1[0])*Ave(yr2[0]) #<y1(0)><y2(0)>
            for key in 1,2: # <y1(q)><y2(-q)> +  <y1(-q)><y2(q)> 
                S[key]=Ave(yr1[key])*Ave(yr2[-key])+Ave(yr1[key])*Ave(yr2[-key])            
        #print Ave(yr1[1]*yr2[-1]),Ave(yr1[-1]*yr2[1])
        return (S[0]+S[1]+S[2]).real,S            

    
    def CalcZeta(self, v1, v2):
        #v1, v2 are arrays which are of a variety of dimensions! Woo!
        #The *last* dimension is the one we care about - IE if we have (x1,x2,...xn), we care about xn as the vector dimension. Rest is *which* vector!
        v1 = numpy.array(v1)
        v2 = numpy.array(v2)
        #This function takes 2 vectors which orientate interactions in space. 
        #We then take the dot-product...
        #...and calculate 1/r**3 terms & product of these
        #...and then calculate P_0(cos(zeta))=1/2(3x**2 - 1)
        #...and then return that value
        #We need to compare this against OrderCalc and spherical terms
        #But if this is the case we can do away with spherical harmonics altogether!
        r1 = numpy.sqrt(numpy.sum(v1 * v1, axis=-1))
        r2 = numpy.sqrt(numpy.sum(v2 * v2, axis=-1))
        zeta = numpy.sum(v1 * v2, axis=-1)/(r1 * r2) #This gives cos(zeta) term, which is all we need
        return 6./r1**3. * 1./r2**3. * 0.5 * (3.*zeta**2. - 1.)
        

    #calculate q=0 spectral density function
    #def Jomega(self,om): #om=self.omega0*fac
    #    return (self.tc/(1+(om)**2.*self.tc**2.))
    #calculate q!=0 spectral density function
    def Jomega(self,k,om):
        #print 'komega'
        #k=Fraction(k)
        #print k,om
        if(self.Woessner):  #if using jumping
            if(k!=0):
                k=(self.n) #(Woessner)
        taueff=1./(1./self.tc + k*1./(self.tm))
        #print 'kay'
        #print taueff
        #print self.tm
        #print 'a' ,taueff/(1+(om)**2.*(taueff)**2.),self.tm*k
        return taueff/(1.+(om)**2.*(taueff)**2.)


    #angular terms for rotating Hamiltonian into molecular frame
    #Yq=sqrt(6)D(q,0)=sqrt(6)sqrt((2-m)!/(2+m)!) P|q|(cos(theta)) e^{-iq alpha}
    #where P|q|(cos(theta)) is associated Legendre polynomial
    def GetSph(self,beta,alpha):
        Y={}
        #try to evaluate the 'pi's
        try:
            if(alpha[:2]=='pi'):
                alpha=alpha.replace('pi','1pi')
        except:
            pass
        try:
            if(beta[:2]=='pi'):
                beta=beta.replace('pi','1pi')
        except:
            pass
        try:
            alpha=eval(alpha.replace('pi','*'+str(numpy.pi)))
        except:
            pass
        try:
            beta=eval(beta.replace('pi','*'+str(numpy.pi)))
        except:
            pass

        #leading constant is the factor that converts associated
        #legendre polynomials to spherical harmonics.
        #def Getccc(l,m):
        #    from math import factorial as ff
        #    return sqrt((2*l+1)*ff(l-m)/(4*numpy.pi*ff(l+m)))
        #print (Getccc(2,2)/Getccc(2,2))
        #print (Getccc(2,1)/Getccc(2,2))
        #print (Getccc(2,0)/Getccc(2,2))

        #constants= sqrt(6)D(q,0)=sqrt(6 (2-q)!/(2+q)!)Pq e(-iaq)
        #spherical harmonics: numerical factor  * associated legendre polynomial * complex alpha rotation
        if(type(beta)==str): #deal with the betas
            Y[-2]=[0.5,'P2']
            Y[-1]=[1.,'P1']
            Y[0]=[sqrt(6),'P0']
            Y[1]=[1.,'P1']
            Y[2]=[0.5,'P2']
        else:
            Y2=3*(Sin(beta)**2. ) 
            Y1=3*Sin(beta)*Cos(beta)
            Y0=(3.*Cos(beta)**2.-1)*0.5
            Y[-2]=[0.5*Y2,'']
            Y[-1]=[Y1,''] #ignore stupid sign problem
            Y[0]=[sqrt(6)*Y0,'']
            Y[1]=[Y1,'']            
            Y[2]=[0.5*Y2,'']

        if(type(alpha)==str): #now add the alphas
            Y[-2].append('exp(-2ialpha)')
            Y[-1].append('exp(-1ialpha)')
            Y[ 0].append('')
            Y[ 1].append('exp(1ialpha)')
            Y[ 2].append('exp(2ialpha)')
        else:
            Y[-2].append(Exp(-2.*complex(0,alpha)))
            Y[-1].append(Exp(-1.*complex(0,alpha)))
            Y[ 0].append(Exp( 0.*complex(0,alpha)))
            Y[ 1].append(Exp( 1.*complex(0,alpha)))
            Y[ 2].append(Exp( 2.*complex(0,alpha)))

        #Plib={}
        #Plib['P2']='3sin^2b'
        #Plib['P1']='3sinbcosb'
        #Plib['P0']='1/2(3cos^2b-1)'
        return Y


    #split string into numerical segments, store them in a dictionary
    #and assemble a full set of operators (eg E,x,x,E,E,E...)
    def GetOp(self,coh):
        n=''
        v=''
        cc={}
        for c in coh: #for each character...
            if(c in '1234567890'): #if it's a number...
                n+=c  #increment the number counter
            elif(len(n)!=0): #we've finished a character.
                cc[int(n)]=c #this is the operator for this number
                n='';v=''  #reset counters
            else:   #if we have no number and c is not a character
                v+=c  #this is the name of the nucleus.
        op=[] #reset list.
        for i in range(self.baseSize): #for each basis member...
            if(i+1 in cc.keys()): #add operator we've found in dictionary
                op.append(cc[i+1])  
            else: #otherwise this element must have the identity
                op.append('E')
        return numpy.array(op) #return list of operators spanning whole basis.


    #create matrix form of operator with direct products
    #or if symbolic, create full operator representation
    def MakeOp(self,coh,save=True):
        #print 'making',coh,self.SYM
        if(coh in self.baseOp.keys()):
            #print coh,'already in basis'
            if(save==False and self.SYM==False):
                return self.base[coh]
            else:
                return

        op=self.GetOp(coh)

        mask=op!='E'
        active=op[mask]
        if(len(active)==1): #if we generate the 1 spin commutators this is unneccesasry
            if(mask[0]==False):
                mask[0]=True
            else:
                mask[1]=True
        #active="".join(op[mask])
        active=op[mask].tostring()
        #if(len(active)!=2):
        #    print 'shit'
        #    sys.exit(100)
        self.baseOp[coh]=(1,op,mask,active,op.tostring()),  #save operator representation
        if(self.SYM):
            return


        #print 'adding:',coh,op,self.baseSize
        #print baseHalf.keys()

        for i in range(self.baseSize): #for each member in the basis...
            if(self.baseType[i]=='D'): #add relevant basis (spin1)for nucleus i
                if(op[i]=='E'):
                    targ=self.baseOne['E']
                else:
                    targ=self.baseOne['C1'+op[i]]
            else: #add relevant basis (spin1/2) for nucleus i
                if(op[i]=='E'):
                    targ=self.baseHalf['E']
                else:
                    targ=self.baseHalf['C1'+op[i]]

            if(i==0):#expand current entry by the new one by outer products
                #test=copy.deepcopy(targ)
                test=targ.copy()
            else:
                test=self.OuterProd(targ,test)
        #print coh,test,test.shape

        if(self.SPARSE): #make sure matrix is sparse if needed
            test=csr_matrix(test) 
            test.eliminate_zeros()
        if(save):
            #self.base[coh]=copy.deepcopy(test) #csr_matrix(test) #save matrix representation
            self.base[coh]=test.copy()
        else:
            return test

    #show current memory and cpu usage statistics
    def ShowMem(self):
        import psutil
        #print psutil.cpu_percent()
        cok= dict(psutil.virtual_memory()._asdict())
        for key,vals in cok.items():
            print key,vals/1024./1024./1024.,'Gb'
                    

    #Dipolar Hamiltonian
    #d is gammaIgammaS / 4pi e0 r^6
    #requires o to map to irreducible basis, operators and spherical harmonics
    #H(mol)= sum_i   o O Yq .Carry full set of Y(q_i) with each entry
    #H(mol)=sum q YqTq
    #H(lab)=sum q sum q2 Yq2 D_q,q2 T_q  #so need to carry around full set of Y terms
    #H(lab)=sum i sum q2 Yq2 D(qi,q2) Oi oi  #so need to carry around full set of Y terms

#     inst.GetDip('H5','H2','f',beta=0,alpha=0,rot='ext')

    def GetDip(self,a,b,ii,tag='static'):
        #DOCUMENTATION:
        #Arguments:
            #a,b - the spins which the dipolar interaction is between
            #ii - the relevent interaction constant - 'd', 'm' etc.
            #tag = 'static', 'single', 'double' -> This indicates the number of timescales under which the interaction evolves
            #As an example of each:
                #'static' - applies to the case where there is a fixed dipolar interaction between 2 spins, no *local* tumbling.
                    #EG. Between 2 fixed backbone amide protons.
                #'single' - applies to the case where one of the spins is on a methyl ring, and the other is fixed.
                    #EG. Between a methyl hydrogen and a fixed backbone amide proton
                    #EG. Between two methyl hydrogens on the *same* ring
                #'double' - applies to case where both of the spins are on methyl rings.
                    #EG. Between two methyl hydrogens on *different* methyl rings
        
        #Interaction terms look like:
            #Aint term - these is a list of different 2-spin operators, with a bunch of terms. Each 2-spin operators looks like:
                #[0.5, 'H5mH2m', 'f', -2, {0: [2.449489742783178, '', (1+0j)], 1: [0.0, '', (1+0j)], 2: [0.0, '', (1+0j)], -1: [0.0, '', (1+0j)], -2: [0.0, '', (1+0j)]}]
                #Form is:
                    #0.5 - o_i factor (5 irrep T_q terms -> 2-spin operators)
                    #'H5mH2m' - String specifies operators
                    #'f' - Relevent interaction constant tag
                    #-2 - Cohberence order
                    #Tuple - Spherical harmonics term. Gives the spatial oreintation of dipolar interaction within molecular frame. In our updated case this will be the *initial*  value of the sph_harm term without having moved. This is ideal
            #Atyp term - Takes the form: 
                #['dipolar', tag, initial, sites]:
                    #'dipolar' - Gives the type of the interaction - dipolar, CSA, quadrupolar
                    #tag = 'static', 'single', 'double'
        
        print 'Including dipolar',a,b
        Anew=[]
        #o,operator,const,coherence order,rotations
        if(a[0]==b[0]): #is this a homonuclear term?
            homo=True
        else:
            homo=False
        Anew.append([ self.o[0][0],a+'z'+b+'z',ii, 0]) #self.o[0][0] = 'sqrt(6) 1/2(1 - 3cos^2t)'
        #only include if not in macromolecular mode. If term is homonuclear,
        #then include anyway as these will be J(0)
        if(homo==True or self.macro==False): #skip only if macromolecular mode and heteronuclear
            Anew.append([ self.o[0][1],a+'p'+b+'m',ii, 0]) #'sqrt(6) 1/2(1 - 3cos^2t)'
            Anew.append([ self.o[0][1],a+'m'+b+'p',ii, 0]) #'sqrt(6) 1/2(1 - 3cos^2t)'
        Anew.append([ self.o[1]   ,a+'z'+b+'p',ii, 1])# '    3 sint cost exp(+p)'
        Anew.append([ self.o[1]   ,a+'z'+b+'m',ii,-1])# '    3 sint cost exp(+p)'
        Anew.append([ self.o[1]   ,a+'p'+b+'z',ii, 1])# '    3 sint cost exp(+p)'
        Anew.append([ self.o[1]   ,a+'m'+b+'z',ii,-1])# '    3 sint cost exp(+p)'
        if(self.macro==False): #include only if not in macromolecular mode
            Anew.append([ self.o[2]   ,a+'p'+b+'p',ii, 2]) #  '1/2 3 sin^2t exp(+2p)'
            Anew.append([ self.o[2]   ,a+'m'+b+'m',ii,-2])  # '1/2 3 sin^2t exp(+2p)'
        self.Aint.append(Anew)
        self.Atyp.append(('dipolar',tag)) #increment type of motion and interaction
        for i in range(len(Anew)): #create the required operators
            self.MakeOp(Anew[i][1])
#         if(self.LATEX):#make pretty latex output
#             alphaVal,betaVal=self.FormatAlphaBeta(alpha,beta)
#             if(ii not in self.edict.keys()):
#                 eVal='R_'+ii
#             else:
#                 eVal='%.2fR_d' % self.edict[ii]
#             outy=open(self.LATEXFILE,'a')
#             outy.write('Dipolar between %s and %s $\\beta$= $%s$ $\\alpha$= $%s$ distance $%s$ constant $%s=\\frac{\\hbar\\gamma_%s\\gamma_%s}{4\\pi\\epsilon_0 (%s)^3}$\n\n' % (a,b,betaVal,alphaVal,eVal,ii,a[0],b[0],eVal))
#             outy.close()

    #take strings and make nice latex output
    def FormatAlphaBeta(self,alpha,beta):
        betaVal=beta
        alphaVal=alpha
        try:
            betaVal=str('%.2f\\pi' % (beta/Pi))
        except:
            betaVal=betaVal.replace('beta','\\beta').replace('pi','\\pi')
        try:
            alphaVal=str('%.2f\\pi' % (alpha/Pi))
        except:
            alphaVal=alphaVal.replace('alpha','\\alpha').replace('pi','\\pi')
        return alphaVal,betaVal

    #introduce a quadruoplar interaction
    def GetQuad(self,a,ii,beta=0,alpha=0,tag='single'):
        print 'Including quad',a
        Y=self.GetSph(beta,alpha) #get spherical harmonics
        Anew=[]
        #o,operator,const,coherence order,rotations
        Anew.append([ self.o[0][0],a+'z'+a+'z',ii, 0,Y]) #'sqrt(6) 1/2(1 - 3cos^2t)'
        Anew.append([ self.o[0][1],a+'p'+a+'m',ii, 0,Y]) #'sqrt(6) 1/2(1 - 3cos^2t)'
        Anew.append([ self.o[0][1],a+'m'+a+'p',ii, 0,Y]) #'sqrt(6) 1/2(1 - 3cos^2t)'
        Anew.append([ self.o[1]   ,a+'z'+a+'p',ii, 1,Y])# '    3 sint cost exp(+p)'
        Anew.append([ self.o[1]   ,a+'z'+a+'m',ii,-1,Y])# '    3 sint cost exp(+p)'
        Anew.append([ self.o[1]   ,a+'p'+a+'z',ii, 1,Y])# '    3 sint cost exp(+p)'
        Anew.append([ self.o[1]   ,a+'m'+a+'z',ii,-1,Y])# '    3 sint cost exp(+p)'
        if(self.macro==False):
            Anew.append([ self.o[2]   ,a+'p'+a+'p',ii, 2,Y]) #  '1/2 3 sin^2t exp(+2p)'
            Anew.append([ self.o[2]   ,a+'m'+a+'m',ii,-2,Y])  # '1/2 3 sin^2t exp(+2p)'
        self.Aint.append(Anew)
        self.Atyp.append(('quad',rot))

        for i in range(len(Anew)): #create the required operators
            self.MakeOp(Anew[i][1][:3]) #1st term  #NOTE! will be tripped by number>9
            self.MakeOp(Anew[i][1][3:]) #2nd term  #NOTE! will be tripped by number>9
        #products of these two...
        for o1,o2 in (('z','z'),('m','p'),('p','m'),('z','p'),('p','z'),('m','z'),('z','m'),('p','p'),('m','m')):
            self.base[a+o1+a+o2]=numpy.dot(self.base[a+o1],self.base[a+o2])
        if(self.LATEX):#make pretty latex output
            alphaVal,betaVal=self.FormatAlphaBeta(alpha,beta)
            outy=open(self.LATEXFILE,'a')
            outy.write('Quadrupolar for %s $\\beta$= $%s$ $\\alpha$= $%s$ constant $%s=\\hbar\\gamma_%sB_0 \\Delta\\delta/3$\n\n' % (a,betaVal,alphaVal,ii,a[0]))
            outy.close()
        

    #csa constant is (parallel-perpendicular) * gamma * B0/3)
    def GetCSA(self,a,ii,tag='single'):
        #TERMS LIKE:
#         [0.5, 'H4p', 'cX', 1, {0: [2.449489742783178, 'P0', (1+0j)], 1: [1.0, 'P1', (1+0j)], 2: [0.5, 'P2', (1+0j)], -1: [1.0, 'P1', (1+0j)], -2: [0.5, 'P2', (1+0j)]}]
        #Again:
            #1st term - o_i factor. 5 T_q terms -> 9 Zeeman 2-spin operators
            #2nd term - Associated Zeeman single spin operator?
            #3rd term - Tag which describes this as a H CSA (CSA for H4)
            #4th term - Coherence order
            #5th term - Spherical harmonic which orientates interaction within the molecular frame
        #Again - we append this interaction Hamiltonian to a list Aint
        #The description of this Hamiltonian - "CSA", rot - is appended to list Atyp
        #And we make the single spin operators which we need
        print 'Including CSA',a
        
        Anew = []
        #o,operator,const,coherence order,rotations
        Anew.append([ self.o[0][0],a+'z',ii, 0])#'sqrt(6) 1/2(1 - 3cos^2t)'
        Anew.append([ self.o[1]   ,a+'p',ii, 1])# '    3 sint cost exp(+p)'
        Anew.append([ self.o[1]   ,a+'m',ii,-1])# '    3 sint cost exp(+p)'
        self.Aint.append(Anew)
        self.Atyp.append(('CSA',tag)) #increment type
        for i in range(len(Anew)): #create the required operators
            self.MakeOp(Anew[i][1])

            
#         #make pretty latex output
#         if(self.LATEX):#make pretty latex output
#             alphaVal,betaVal=self.FormatAlphaBeta(alpha,beta)
#             outy=open(self.LATEXFILE,'a')
#             outy.write('CSA for %s $\\beta$= $%s$ $\\alpha$= $%s$ constant $%s=\\hbar\\gamma_%sB_0 \\Delta\\delta/3$\n\n' % (a,betaVal,alphaVal,ii.replace('c','c_'),a[0]))
#             outy.close()



    #work out frequency from coherence order
    #text string created from numbers
    def ParseSum(self,freqstr,ssum,freq,sgn):
        ssum=ssum*sgn
        if(ssum!=0):
            if(numpy.fabs(ssum)==1):
                if(ssum>0):
                    freqstr+='+'+freq
                else:
                    freqstr+='-'+freq
            else:
                if(ssum>0):
                    freqstr+='+'+str(ssum)+freq
                else:
                    freqstr+='-'+str(ssum)+freq
        return freqstr

    #work out frequency of operator from its name
    def GetFreq(self,coh):
        #TESTING:
        mask=self.baseList=='H'
        Hsum=self.baseOp[coh][0][1][mask].tostring().count('p')
        Hsum-=self.baseOp[coh][0][1][mask].tostring().count('m')

        mask=self.baseList=='D'
        Dsum=self.baseOp[coh][0][1][mask].tostring().count('p')
        Dsum-=self.baseOp[coh][0][1][mask].tostring().count('m')

        mask=self.baseList=='C'
        Csum=self.baseOp[coh][0][1][mask].tostring().count('p')
        Csum-=self.baseOp[coh][0][1][mask].tostring().count('m')

        q=numpy.array((Csum,Hsum,Dsum)) #save frequency (C,H,D)

        #1. make D positive.
        #2. make C positive if D is zero.
        #3. make H positive if D and C are zero.
        freqstr=''
        if(Dsum==0 and Hsum==0 and Csum==0):
            freqstr+=str(0)
        if(Dsum!=0):
            if(Dsum>0):
                sgn=1
            else:
                sgn=-1
        elif(Csum!=0):
            if(Csum>0):
                sgn=1
            else:
                sgn=-1
        elif(Hsum!=0):
            if(Hsum>0):
                sgn=1
            else:
                sgn=-1
        else:
            sgn=0

        freqstr=self.ParseSum(freqstr,Dsum,'wD',sgn)
        freqstr=self.ParseSum(freqstr,Csum,'wC',sgn)
        freqstr=self.ParseSum(freqstr,Hsum,'wH',sgn)
    
        if(freqstr!='0'):
            final=freqstr[1:]
        else:
            final=freqstr
        if(freqstr==''):
            print 'Cannot determine frequency!',Csum,Hsum,Dsum
            sys.exit(100)
        #print coh,final,q
        return final,q


    #typeset for printing to screen
    def FormatSpecSh(self,face,j1):
        if(face>0):
            stry=' +'
        else:
            stry=' -'
        if(face<0):
            face=face*-1
        if(face.denominator==1):
            face=face.numerator
            stry+= '('+str(face.numerator)+') '+j1
        else:
            stry+= '('+str(face.numerator)+'/'+str(face.denominator)+') '+j1
        return stry

    #typeset for printing to screen
    def FormatSpecLatex(self,face,j1):
        j1=j1.replace('tc,','') #no need for this.
        j1=j1.replace('w','\\omega_')
        j1=j1.replace('(','\\left(')
        j1=j1.replace(')','\\right)')
        #j1=j1.replace('tc','\\tau_c')
        
        j1=j1.replace('tm/1','\\tau_m')
        j1=j1.replace('tm/4','\\frac{\\tau_m}{4}')
        j1=j1.replace('tm','\\tau_m')

        stry=self.FormatLatexFrac(face)
        return stry+' '+j1
        """
        if(face>0):
            stry=' +'
        else:
            stry=' -'
        if(face<0):
            face=face*-1
        if(face.denominator==1):
            face=face.numerator
            stry+= ' '+str(face.numerator)+' '+j1
        else:

            if(len(str(face.numerator))>3):
                facenew=(face*face).limit_denominator()
                if(len(str(facenew.numerator))<len(str(face.numerator))):
                    if((facenew.denominator**0.5).is_integer()):
                        stry+= '\\frac{\sqrt{'+str(facenew.numerator)+'}}{'+str(int(facenew.denominator**0.5))+'} '+j1
                    else:
                        stry+= '\\left(\\frac{'+str(facenew.numerator)+'}{'+str(facenew.denominator)+'}\\right)^{1/2} '+j1
                else:
                    stry+= '\\frac{'+str(face.numerator)+'}{'+str(face.denominator)+'} '+j1
            else:

                stry+= '\\frac{'+str(face.numerator)+'}{'+str(face.denominator)+'} '+j1
        return stry
        """
    def FormatLatexFrac(self,face,plus=True):
        face=Fraction(face).limit_denominator()
        stry=''
        if(face>0):
            if(plus):
                stry=' +'
        else:
            stry=' -'
        if(face<0):
            face=face*-1
        if(face.denominator==1):
            face=face.numerator
            stry+= ' '+str(face.numerator)
        else:

            if(len(str(face.numerator))>3):
                facenew=(face*face).limit_denominator()
                if(len(str(facenew.numerator))<len(str(face.numerator))):
                    

                    if((facenew.numerator**0.5).is_integer()):
                        stry+= '\\frac{'+str(int(facenew.numerator**0.5))+'}'

                        #{'+str(int(facenew.denominator**0.5))+'}'
                        #nom=facenew.numerator**0.5
                    else:
                        stry+= '\\frac{\\sqrt{'+str(facenew.numerator)+'}}'
                        #{'+str(facenew.denominator)+'}}'
                        #nom=facenew.numerator**0.5

                    if((facenew.denominator**0.5).is_integer()):
                        #stry+= '\\frac{\sqrt{'+str(facenew.numerator)+'}}
                        stry+='{'+str(int(facenew.denominator**0.5))+'}'
                        #den=facenew.denominator**0.5
                    else:
                        #stry+= '\\sqrt{\\frac{'+str(facenew.numerator)+'}
                        stry+='{\\sqrt{'+str(facenew.denominator)+'}}'
                else:
                    stry+= '\\frac{'+str(face.numerator)+'}{'+str(face.denominator)+'}'
            else:

                stry+= '\\frac{'+str(face.numerator)+'}{'+str(face.denominator)+'}'
        return stry

    #take conjugate of operator (easy in shift basis)
    def GetConjOp(self,ref,p):
        n=''
        for j in range(len(ref)):
            if(ref[j]=='m'):
                n+='p'
            elif(ref[j]=='p'):
                n+='m'
            else:
                n+=ref[j]
        return n,-p

    #evaluate double commutators.
    #sorted by O1, then O2. So if first commutator is zero for one of a set
    #no point calculating any more. Same for O2.
    #Once double commutator is done, append all constant/frequency terms that apply.
    def EvalEntry(self,Pre,Posts):
        result={}
        pre=self.base[Pre]
        norm=self.Norms[Pre] #set normalisation constant to current prefactor 
        for O1,d1 in self.Aind.items():
#             print 'O1: ',O1
#             print "dl: ", d1
            #prefac d1 d2 J(j1) <post | [O2,[O1,pre]]> / <post|pre>
            a=comm(self.base[O1],pre)
            if(self.IsZero(a)): #if comm is zero...
                continue #no need to try any other commutators
            for O2,d2 in d1.items():
#                 print 'O2:', O2
                b=comm(self.base[O2],a)
#                 print "b: ", b
                if(self.IsZero(b)): #if double comm is zero...
                    continue #no need to try any other commutators.
                #doing this now saves a ton of memory (no need to save double comm matrix)
                for Post in Posts: #for each operator, close the result.
                    post=self.base[Post] #post operator
                    tr=self.Trace(post,b)  #finish the expression
                    if(fabs(tr)<1E-6):#if trace is zero or effectively zero, abort.
                        continue      #skip rest
                    #keep!
                    normPost=self.Norms[Post] #post normalisation constant
                    trD=normPost*norm #normalisation
                    for dd,d3 in d2.items(): #if nonzero, save.
                        for j1,prefac in d3.items(): #if nonzero, save.
                            try: #index by post...
                                ccc=result[Post]
                            except:
                                result[Post]={}
                            try: #index by dd...
                                ccc=result[Post][dd]
                            except:
                                result[Post][dd]={}
                            try: #and j1...
                                ccc=result[Post][dd][j1]
                            except: 
                                result[Post][dd][j1]=0
                            result[Post][dd][j1]+= (tr*prefac/trD) #increment final number
                            
                            if(j1=='J(0,wH)' and dd=='e^2' and Post=='CzHzHz' and Pre=='CzHzHz'):
                                print '<',Post,'[',O2,'[',O1,',',Pre,']]>=',Fraction(tr*prefac/trD).limit_denominator(),j1,dd,tr,prefac,trD
                            #num,com,op,ma=compare(a,self.base)
                            #print num,op
                            #num,com,op,ma=compare(b,self.base)
                            #print num,op
                            #print tr,prefac,trD
                            #if(tr!=0.25 and tr!=-0.25):
                            #    sys.exit(100)
                            """
                            if(dd=='f^2'):
                                if(Pre=='aaa' and Post=='aab'):
                                    print '1',Pre,Post,O1,O2,dd,j1,tr,norm,normPost,trD,prefac
                                if(Pre=='bbb' and Post=='abb'):
                                    print '2',Pre,Post,O1,O2,dd,j1,tr,norm,normPost,trD,prefac
                                if(j1=='J(tc,0)'):
                                    if(Pre=='aaa'):
                                        ii=0
                                    elif(Pre=='aab'):
                                        ii=1
                                    elif(Pre=='abb'):
                                        ii=2
                                    elif(Pre=='bbb'):
                                        ii=3
                                    if(Post=='aaa'):
                                        jj=0
                                    elif(Post=='aab'):
                                        jj=1
                                    elif(Post=='abb'):
                                        jj=2
                                    elif(Post=='bbb'):
                                        jj=3
                                    
                                    self.tempy[ii,jj]+=tr/(len(self.baseOp[Pre]))
                           """
        try:
            ccc=self.resSet[Pre]
        except:
            self.resSet[Pre]={}
        for Post in Posts: #for each operator, close the result.
            try:
                self.resSet[Pre][Post]=self.CleanUp(result[Post]) #save result
            except:
                self.resSet[Pre][Post]={}

    #work out symmetry operators that relate 
    #positions of O1 and O2
    def GetSym(self,O1,O2):
        o1=self.baseOp[O1][0][2][1:]*1.
        o2=self.baseOp[O2][0][2][1:]*1.
        foggle1=[]
        for key,val in self.maps.items(): #for each symmetry op
            o2t=numpy.dot(val,o2) #apply symmetry operator
            if((o2t==o1).all()):
                newOp=key.split('_')[0]
                foggle1.append(newOp)
        return sorted(foggle1)


    #symbolically calculate the pre operators (double commutators)
    #do all in one function (slow because not consolidating trace operations)
    def EvalEntrySym(self,Pre,Posts):
        self.norm=self.Norms[Pre] #set normalisation constant to current prefactor (used in evalpost)
        baseops=self.baseOp[Pre]
        result={}

        #parallelise over this loop.
        for O1,d1 in self.Aind.items():
            #print 'prefac d1 d2 J(j1) <post | [O2,[O1,pre]]> / <post|pre'
            #print 'O1:',O1            
            a=self.CommSym(self.baseOp[O1],baseops) #[O1,pre]
            if(a==0): #if commutator is zero...
                continue #no need to try any other commutators
            for O2,d2 in d1.items():
                b=self.CommSym(self.baseOp[O2],a) #[O2,[O1,pre]]
                if(b==0): #if double comm is zero...
                    continue

                #refmat=self.CommVal(b)
                #comm2=comm(self.base[O2],comm1)
                #if(self.IsZero(refmat-comm2)):
                #    print 'good'
                #else:
                #    print 'shit'
                print "b: ", b
                for bvals in b: #go over each item in double commutator...
                    for Post in Posts: #for each operator, close the result.
                        print "POST: ", post
                        postOp=self.baseOp[Post] #unpack
                        tr=self.TraceSym(postOp,((1.,bvals[1]),))  #finish the expression. list of operators in bval[0]
                        if(fabs(tr)<1E-6):#if trace is zero or effectively zero, abort.
                            continue      #skip rest

                        normPost=self.Norms[Post] #normalisation constant for post
                        trD=normPost*self.norm #normalisation

                        for dd,d3 in d2.items(): #if nonzero, save.
                            for j1,prefac in d3.items():
                                try: #index by post...
                                    ccc=result[Post]
                                except:
                                    result[Post]={}
                                try: #index by dd...
                                    ccc=result[Post][dd]
                                except:
                                    result[Post][dd]={}
                                try: #and j1...
                                    ccc=result[Post][dd][j1]
                                except: 
                                    result[Post][dd][j1]=0
                                result[Post][dd][j1]+= (tr*prefac*bvals[0]/trD) #increment final number

        try:
            ccc=self.resSet[Pre]
        except:
            self.resSet[Pre]={}
        for Post in Posts: #for each operator, close the result.
            try:
                self.resSet[Pre][Post]=self.CleanUp(result[Post]) #save result
            except:
                self.resSet[Pre][Post]={}






    #symbolically calculate the pre operators (double commutators)
    def EvalPreSym(self,Pre,Aind):
        self.norm=self.Norms[Pre] #set normalisation constant to current prefactor (used in evalpost)
        self.bvals={} #store non-zero double commutators
        self.blist={}
        baseops=self.baseOp[Pre]

        for O1,d1 in Aind.items():
#             print 'prefac d1 d2 J(j1) <post | [O2,[O1,pre]]> / <post|pre'
#             print 'O1:',O1            
            a=self.CommSym(self.baseOp[O1],baseops) #[O1,pre]
            if(a==0): #if commutator is zero...
                continue #no need to try any other commutators
            for O2,d2 in d1.items():
                b=self.CommSym(self.baseOp[O2],a) #[O2,[O1,pre]]
                if(b==0): #if double comm is zero...
                    continue

                #refmat=self.CommVal(b)
                #comm2=comm(self.base[O2],comm1)
                #if(self.IsZero(refmat-comm2)):
                #    print 'good'
                #else:
                #    print 'shit'
                for bvals in b:
                    test=bvals[4] #compressed representation
                    try:#surprisingly fast way to do this
                        ccc=self.blist[test]
                    except:
                        self.blist[test]={}
                        self.blist[test]['val']=(1.,bvals[1]),
                        self.blist[test]['arg']=[]
                    #this is very slow #if(test not in self.blist.keys()): 
                    for dd,d3 in d2.items(): #if nonzero, save.
                        for j1,prefac in d3.items():
                            self.blist[test]['arg'].append((j1,bvals[0]*prefac,dd))

    #symbolically calculate the pre operators (double commutators)
    def EvalPreSymParallel(self,Pre):
        #print 'bef',len(self.blist.keys())

        self.norm=self.Norms[Pre] #set normalisation constant to current prefactor (used in evalpost)
        self.blist={}

        baseops=self.baseOp[Pre]


        print 'items:',len(self.Aind.keys())

        #split Aind over multiple jobs.
        jabs=[]
        self.jobs={}
        cnt=0
        run=[]
        for O1,d1 in self.Aind.items():
            jabs.append((O1,d1))
            if(len(jabs)>len(self.Aind.keys())/self.ncpus):
                self.jobs[cnt]=jabs
                jabs=[]
                run.append((cnt,baseops))
                cnt+=1
        if(len(jabs)>0):
            self.jobs[cnt]=jabs
            run.append((cnt,baseops))
            cnt+=1
            jabs=[]
            

        import pp
        ppservers = ()
        job_server = pp.Server(self.ncpus, ppservers=ppservers)
        #print jobs
        #print len(jobs)
        #jobs = [(cnt, job_server.submit(self.DoEvalPreSymParallel,(cnt,baseops,), globals=globals())) for cnt,baseops in run]
        jobs = [(cnt, job_server.submit(self.DoEvalPreSymParallel,(cnt,baseops,),globals=globals())) for cnt,baseops in run]
        for input,job in jobs:
            self.blist.update(job())
            job()
            
        job_server.print_stats()
        job_server.destroy()
        #job1=job_server.submit(self.DoEvalPreSymParallel,run, globals=globals())

    def DoEvalPreSymParallel(self,val,baseops):
        #n3=datetime.now()
        blist={}
        print 'processor:',val,len(self.jobs[val])
        for O1,d1 in self.jobs[val]:
            #d1=self.Aind[O1]
            #print 'prefac d1 d2 J(j1) <post | [O2,[O1,pre]]> / <post|pre'
            #print 'O1:',O1            
            a=self.CommSym(self.baseOp[O1],baseops) #[O1,pre]
            if(a==0): #if commutator is zero...
                continue #no need to try any other commutators
            for O2,d2 in d1.items():
                b=self.CommSym(self.baseOp[O2],a) #[O2,[O1,pre]]
                if(b==0): #if double comm is zero...
                    continue
                for bvals in b:
                    test=bvals[4] #compressed representation
                    try:#surprisingly fast way to do this
                        ccc=blist[test]
                    except:
                        blist[test]={}
                        blist[test]['val']=(1.,bvals[1]),
                        blist[test]['arg']=[]
                    #this is very slow #if(test not in self.blist.keys()): 
                    for dd,d3 in d2.items(): #if nonzero, save.
                        for j1,prefac in d3.items():
                            blist[test]['arg'].append((j1,bvals[0]*prefac,dd))
        #print 'tt',val,datetime.now()-n3
        return blist





    #find non-zero elements and work out corresponding operator
    #experimental. should loop over matrix and find non-zero elements.
    def WhoAmI(self,a,verb='n'):

        #print a.shape
        rows,cols=a.nonzero()
        if(verb=='y'):
            print
            print 'analysing'
            print a.shape
            print a
        reps=[]
        for row,col in zip(rows,cols):
            aa=row
            bb=col
            #print 'non-zero element:',aa,bb,a[aa,bb]
            sz=a.shape[0]
            bis=[]
            while(sz>=2):
                #print '   ',aa,bb,bis,sz
                if(aa>sz/2-1):
                    aa-=(sz/2)
                    if(bb>sz/2-1):
                        bb-=(sz/2)
                        bis.append('b')
                        #print 'b'
                    else:
                        bis.append('m')
                        #print 'm'
                else:
                    if(bb>sz/2-1):
                        bb-=(sz/2)
                        bis.append('p')
                        #print 'p'
                    else:
                        bis.append('a')
                        #print 'a'
                sz=sz/2
            bis=numpy.array(bis)
            if(numpy.fabs(a[row,col])>1E-6):
                if(verb=='y'):
                    print 'element:',bis.tostring(),'values:',row,col,a[row,col]
                #labby=self.GetLabelPure(bis[1:],self.baseSize-1)            
                reps.append((a[row,col],bis,1,1,bis.tostring()))
        return reps

    #close double commutator with required operators
    #symbolically
    def EvalPostSym(self,Pre,Post): 
        try: #save result
            ccc=self.resSet[Pre]
        except:
            self.resSet[Pre]={}
        try: #if there is already a pre/post entry, skip
            cc=self.resSet[Pre][Post]
            return
        except:
            pass



        result={}       #store the final result in dictionary, indexed by prefactor and frequency
        if(self.POINT): #get set of operators that constitute 'post'
            postOp=self.baseOp[Post+'Post']
        else:
            postOp=self.baseOp[Post]

        #print 'Number of double commutators to close:',len(self.blist.keys())
        #for i,pist in enumerate(postOp):
        #print len(postOp),len(self.figgle[Post])

        result={}
        normPost=self.Norms[Post] #normalisation constant for post
        trD=normPost*self.norm #normalisation
        for b in self.blist.keys():
            tr=self.TraceSym(postOp,self.blist[b]['val'])  #finish the expression. list of operators in bval[0]
            #tr=self.Trace(self.CommVal(postOp),self.CommVal(self.blist[b]['val']))
            #print tr,self.Trace(self.base[Post],self.CommVal(bval[0])),self.Trace(self.base[Post],bval[4])
            if(fabs(tr)<1E-6):#if trace is zero or effectively zero, abort.
                continue      #skip rest
            #result={}
            #otherwise, maybe keep!
            for j1,prefac,dd in self.blist[b]['arg']:#unpack frequencies and constants
                try: #index by dd...(symbolic)
                    ccc=result[dd]
                except:
                    result[dd]={}
                try: #and j1... (symbolic)
                    ccc=result[dd][j1]
                except:
                    result[dd][j1]=0
                result[dd][j1]+=(tr*prefac/trD) #increment (value)
        self.resSet[Pre][Post]=self.CleanUp(result) #save result


    #close double commutator with required operators
    #symbolically
    def EvalPostSymParallel(self,Pre,Post): 
        
        try: #save result
            ccc=self.resSet[Pre]
        except:
            self.resSet[Pre]={}
        """
        try: #if there is already a pre/post entry, skip
            cc=self.resSet[Pre][Post]
            return
        except:
            pass
        """



        result={}       #store the final result in dictionary, indexed by prefactor and frequency
        if(self.POINT): #get set of operators that constitute 'post'
            self.postOp=self.baseOp[Post+'Post']
        else:
            self.postOp=self.baseOp[Post]

        print 'Number of double commutators to close:',len(self.blist.keys())
        #for i,pist in enumerate(postOp):
        #print len(postOp),len(self.figgle[Post])

        result={}
        normPost=self.Norms[Post] #normalisation constant for post
        self.trD=normPost*self.norm #normalisation

        #split Aind over multiple jobs.
        jabs=[]
        self.jobs={}
        cnt=0
        run=[]
        for b in self.blist.keys():
            jabs.append(b)
            if(len(jabs)>len(self.blist.keys())/self.ncpus):
                self.jobs[cnt]=jabs
                jabs=[]
                run.append((cnt,))
                cnt+=1
        if(len(jabs)>0):
            self.jobs[cnt]=jabs
            run.append((cnt,))
            cnt+=1
            jabs=[]
            

        import pp
        ppservers = ()
        job_server = pp.Server(self.ncpus, ppservers=ppservers)
        #print jobs
        #print len(jobs)


        jobs = [(cnt, job_server.submit(self.DoEvalPostSymParallel,(cnt,),globals=globals())) for cnt in run]
        for input,job in jobs:
            print 'fockles',input
            result.update(job())


        self.resSet[Pre][Post]=self.CleanUp(result) #save result

        job_server.print_stats()
        job_server.destroy()


    def DoEvalPostSymParallel(self,val):
        #print 'fcukes'
        #print val
        #print 'processor:',val[0],len(self.jobs[val[0]])
        result={}
        for b in self.jobs[val[0]]:
            #print b
            #print self.blist[b]['val']
            tr=self.TraceSym(self.postOp,self.blist[b]['val'])  #finish the expression. list of operators in bval[0]
            #tr=self.Trace(self.CommVal(postOp),self.CommVal(self.blist[b]['val']))
            #print tr,self.Trace(self.base[Post],self.CommVal(bval[0])),self.Trace(self.base[Post],bval[4])
            if(math.fabs(tr)<1E-6):#if trace is zero or effectively zero, abort.
                continue      #skip rest
            #result={}
            #otherwise, maybe keep!
            for j1,prefac,dd in self.blist[b]['arg']:#unpack frequencies and constants
                try: #index by dd...(symbolic)
                    ccc=result[dd]
                except:
                    result[dd]={}
                try: #and j1... (symbolic)
                    ccc=result[dd][j1]
                except:
                    result[dd][j1]=0
                result[dd][j1]+=(tr*prefac/self.trD) #increment (value)
        return result

    #group terms and remove zeros
    def CleanUp(self,result):
        for dd,d1 in result.items(): #for constants...
            for j1,val in d1.items(): #for frequencies...
                if(numpy.fabs(val.imag)>1E-6):
                    print 'value is imaginary...'
                    print dd,j1,val
                    sys.exit(100)
                if(numpy.fabs(val.real)<1E-6):#clean up residual zero values
                    del result[dd][j1] #kill the entry
                else:
                    result[dd][j1]=val.real
            if(len(result[dd].keys())==0): #make sure something is left
                del result[dd]
        return result



    #make sure pres and posts are both lists.
    #even if just one entry is added.
    def CheckLists(self,Pres,Posts):
        a='aaa'
        if(type(Pres)==type(a)):
            Pres=Pres,
        if(type(Posts)==type(a)):
            Posts=Posts,
        return Pres,Posts

    #calculate normalisation constants
    #and matrix operators required for a calculation
    def GetNorms(self,Pres,Posts):
        #print 'Calculating matricies and normalisation constants...'
        self.Norms={} #calculate normalisation constants for post operators
        for Post in Posts: #for each operator, close the result.
            self.MakeOp(Post) #create pre operator, if not already there
            post=self.base[Post]
            if(Post not in self.Norms.keys()):
                self.Norms[Post]=Sqrt(self.Trace(post,post))
        for Pre in Pres: #for each operator, close the result.
            self.MakeOp(Pre) #create pre operator, if not already there
            pre=self.base[Pre]
            if(Pre not in self.Norms.keys()):
                self.Norms[Pre]=Sqrt(self.Trace(pre,pre))
        #for key,vals in self.Norms.items():
        #    print key,vals
        #    print self.baseOp[key]
        #sys.exit(100)

    #calculate normalisation constants
    #and matrix operators required for a calculation
    def GetNormsSym(self,Pres,Posts):
        #print 'Calculating matricies and normalisation constants...'
        self.Norms={} #calculate normalisation constants for post operators
        for Post in Posts: #for each operator, close the result.
            self.MakeOp(Post) #create pre operator, if not already there
            if(Post not in self.Norms.keys()):
                self.Norms[Post]=Sqrt(self.TraceSym(self.baseOp[Post],self.baseOp[Post]))
        for Pre in Pres: #for each operator, close the result.
            self.MakeOp(Pre) #create pre operator, if not already there
            if(Pre not in self.Norms.keys()):
                self.Norms[Pre]=Sqrt(self.TraceSym(self.baseOp[Pre],self.baseOp[Pre]))




    #evaulte rate for pre and post on correlated Hamiltonians
    #saves in self.resSet, indexed by pre and post.
    #results indexed by symbolic constant, and symbolic frequency.
    def EvalRate(self,Pres,Posts,verb='y',calc='n'):
        #Pres & posts are lists of operators (in string representation), to evaluate the relaxation rate between. So evaluate rate between every Pre & every Post
        #self.SYM indicates whether this is to be done symbolically or using matrices
        #calc indicates whether we should calculate numerical results
        Pres, Posts=self.CheckLists(Pres, Posts) #make sure pre and post are lists

        #when considering Pres and Posts combinations:
        #if coherence order of Post does not match pre, then will have rate=0.
        #need to implement this in evalpostsym (symbolic) and evalentry (matrix)
        #code will work this out, but there's no reason to because answer is known.
        if(self.SYM):
            self.GetNormsSym(Pres, Posts) #calculate normalisation constants and matricies for pre/post
        else:
            self.GetNorms(Pres, Posts) #calculate normalisation constants and matricies for pre/post
        from datetime import datetime

        for Pre in Pres: #loop over Pre vals
            if(self.SYM):


                n=datetime.now()
                #if(len(Pre.split('H'))==1 and self.CYCLIC):
                #    self.EvalPreSym(Pre,self.Aind2) #work out the double commutators.
                #else:
                self.EvalPreSym(Pre, self.Aind) #work out the double commutators.
                #print 'A',datetime.now()-n


                #n2=datetime.now()
                #self.EvalPreSymParallel(Pre) #work out the double commutators.
                #print 'B',datetime.now()-n2


                #n3=datetime.now()
                for Post in Posts: #for each operator, close the result.
                    #n=datetime.now()
                    self.EvalPostSym(Pre, Post) #evaluate and save final rate expression
                #print '33',datetime.now()-n3

                #n4=datetime.now()
                #for Post in Posts: #for each operator, close the result.
                #    #n=datetime.now()
                #    self.EvalPostSymParallel(Pre,Post) #evaluate and save final rate expression
                #print '44',datetime.now()-n4

                #print 'Time',datetime.now()-n
            else:
                self.EvalEntry(Pre,Posts) #work out the double commutators.

            #this is basically a bad idea...
            #for symbolic much better to consolidate elements for quicker traces.
            #n2=datetime.now()
            #self.EvalEntrySym(Pre,Posts)
            #print 'new',datetime.now()-n2

            for Post in Posts: #for each operator, close the result.
                if(verb=='y'): #print symbolic result to screen
                    self.ShowRate(Pre,Post)
                    if(self.LATEX): #print symbolic result to screen
                        #self.ShowRate(Pre,Post,latex=True)
                        self.ShowRate(Pre,Post,latex=True,sq=self.SQUARE)
                if(calc=='y'): #calculate rate using current constants
                    self.CalcRate(Pre,Post,verb='y')


    #format constant term so that it looks nice in latex
    def FormatDLatex(self,dd):
        #if(len(dd.split('P0'))>2):
        #    print dd
        #    #print dd.replace('P0P0'
        #    sys.exit(10)
        dd=dd.replace('c','c_')
        dd=dd.replace('frac_','frac')
        dd=dd.replace('P1P1','P_1^2')
        dd=dd.replace('P0P0','P_0^2')
        dd=dd.replace('P2P2','P_2^2')
        dd=dd.replace('P1','P_1')
        dd=dd.replace('P0','P_0')
        dd=dd.replace('P2','P_2')
        dd=dd.replace('^(1/2)','^{1/2}')
        return dd

    #try to simplify fractions
    #will try squaring and factorising
    #value to get something that looks nice.
    def AdjustFraction(self,frac,deno,op=True):
        #print 'frac',frac
        #print 'deno',deno
        #print frac*deno
        #print (frac*deno)**0.5
        if(op==False):
            test=Fraction( (frac*deno)).limit_denominator()
            #print 'testA',test
            if(len(str(test.numerator))>4):
                return '('+str((frac*deno))+')',test
            else:
                #print frac*deno,Fraction(frac*deno).limit_denominator()
                #print 'returning:', str(((frac*deno)))
                if(test.denominator==1):
                    if(test==1):
                        return '',test
                    else:
                        return str(int(round(frac*deno))),test
                else:
                    return str(test.numerator)+'/'+str(test.denominator),test
        #print 'test',test,test.numerator*1./test.denominator
        test=Fraction( (frac*deno)**0.5).limit_denominator()
        
        if(len(str(test.numerator))>4):
            return '('+str((frac*deno))+')^(1/2)',test
        else:
            #print 'returning:', str(int((frac*deno)**0.5))
            return str(int((frac*deno)**0.5)),test


    #adjust the value in the results dict
    #remove the entry if its effectively zero
    def AdjustVal(self,dic,a,b,num):
        dic[a][b]+=num
        if(numpy.fabs(dic[a][b])<1E-6):
            del dic[a][b]

    #format term to get the non-square version
    #for completing the square
    def FormatText(self,term):
        test=term.split("^")
        ret=test[0]
        if(len(test)>1):
            tast=test[1].split("P")
            if(len(tast)>1):
                ret+="P"+tast[1]
        return ret

    #take a mask and try to complete the square
    #to reveal TROSY effects
    def FactorTROSY(self,res,terms,freq):
        da =Fraction(res[terms[0]][freq]).limit_denominator() #dipolar term ^2
        ca =Fraction(res[terms[1]][freq]).limit_denominator() #csa term^2
        dca=Fraction(res[terms[2]][freq]).limit_denominator() #mixed term
        if(ca.numerator<0):
            bigFac=-1.
            da*=-1
            ca*=-1
            dca*=-1
        else:
            bigFac=1.
            #return False,0,0,0,0,0,0

        if(dca<0):
            facA=1.
            tA='-'
        else:
            facA=-1.
            tA='+'

        #dcanew=danewnum*canewnum*2.-dca*denoa*facA
        #if(numpy.fabs(dcanew)<1E-6): 
        #    print 'square complete!'

        #dipolar coefficient  a1*d^2
        #mix coefficient      a2*dc
        #csa coefficient      a3*c^2
        #( a2/(2(sqrt(a3)) +  sqrt(a3) )^2
        # (d1 d+d2 c)**2.
        # d1^2 d^2 + 2 d1 d2 dc + d2^2 c^2
        #d2=sqrt(a3)
        #a2=2 d1 d2
        #d1=a2/(2sqrt(a3))
            
        #then factor out what you need for J until the coefficients are nice.
        #a2=dca
        #a3=ca
        #d1=a2/2(sqrt(a3))
        
        #subtract from dipolar:
        d1=numpy.fabs(dca.numerator)/(dca.denominator)/(2.*ca**0.5) #needs a positive ca
        d2=sqrt(ca)

        #should be zero
        #print 'a3',d2**2.,ca
        #print 'a2',-1*facA*d1*2*d2,dca
        #print 'a1',d1**2.,da

        #subtract the following factors:
        self.AdjustVal(res,terms[0],freq,-d1**2.*bigFac)
        self.AdjustVal(res,terms[1],freq,-ca*bigFac )
        self.AdjustVal(res,terms[2],freq,-1*dca*bigFac )  

        #try and factor out a constant factor to make a nicer fraction.
        to=1.
        while(1==1):
            #print to
            d4,d4num=self.AdjustFraction(Fraction(d1*sqrt(to*1.)).limit_denominator(),1,op=False) #factor this from chemical shifts...
            d5,d5num=self.AdjustFraction(Fraction(ca**0.5*sqrt(to*1.)).limit_denominator(),1,op=False) #factor this from chemical shifts...
            if(len(str(d4num.numerator))<4 and len(str(d5num.numerator))<4):
                break
                #print to,d4num,d5num,d4,d5
            if(to==60.):
                to=1.
                d4,d4num=self.AdjustFraction(Fraction(d1*sqrt(to*1.)).limit_denominator(),1,op=False) #factor this from chemical shifts...
                d5,d5num=self.AdjustFraction(Fraction(ca**0.5*sqrt(to*1.)).limit_denominator(),1,op=False) #factor this from chemical shifts...
                break
            to+=1.

        denoa=d4num.denominator
        if(len(str(denoa))>4):
            denoa=1.

        if(denoa==1 and to==1): #amounts to failure
            d4num=Fraction(d1).limit_denominator()
            d4f=Fraction(d1**2.).limit_denominator()
            d4g=Fraction(d4f*d4f).limit_denominator()
            if(len(str(d4g.numerator))>len(str(d4f.numerator))):
                if((d4f.denominator**0.5).is_integer()):
                    #d4='('+str(d4f.numerator)+'/'+str(d4f.denominator)+')^(1/2)'
                    d4='\\frac{\\sqrt{'+str(d4f.numerator)+'}}{'+str(d4f.denominator)+'}'
                else:
                    d4='\\left(\\frac{'+str(d4f.numerator)+'}{'+str(d4f.denominator)+'}\\right)^(1/2)'
            else:
                #d4=str(d4g.numerator)+'/'+str(d4g.denominator)
                d4='\\frac{'+str(d4g.numerator)+'}{'+str(d4g.denominator)+'}'

            d5num=Fraction(ca**0.5).limit_denominator()
            d5f=Fraction(ca).limit_denominator()
            d5g=Fraction(d5f*d5f).limit_denominator()
            if(len(str(d5g.numerator))>len(str(d5f.numerator))):
                if((d5f.denominator**0.5).is_integer()):
                    #stry+= '\\frac{\sqrt{'+str(facenew.numerator)+'}}{'+str(int(facenew.denominator**0.5))+'} '+j1
                    d5='\\frac{\sqrt{'+str(d5f.numerator)+'}}{'+str(int(d5f.denominator**0.5))+'}'
                else:
                    d5='(\\frac{'+str(d5f.numerator)+'}{'+str(d5f.denominator)+'})^(1/2)'
                #d5='('+str(d5f.numerator)+'/'+str(d5f.denominator)+')^(1/2)'
            else:
                d5='\\frac{'+str(d5g.numerator)+'}{'+str(d5g.denominator)+'}'
                #d5=str(d5g.numerator)+'/'+str(d5g.denominator)

            #d5='('+str(d5f.numerator)+'/'+str(d5f.denominator)+')^(1/2)'
            #d4='('+str(d4f.numerator)+'/'+str(d4f.denominator)+')^(1/2)'
            #print d4,d4num,d4f,d4g
            #print d5,d5num,d5f,d5g
        else:
            d4,d4num=self.AdjustFraction(Fraction(d1*sqrt(to*1.)*(denoa)).limit_denominator(),1,op=False) #factor this from chemical shifts...
            d5,d5num=self.AdjustFraction(Fraction(ca**0.5*sqrt(to*1.)*(denoa)).limit_denominator(),1,op=False) #factor this from chemical shifts...


        jay=1./(1.*to*denoa**2.)*bigFac #set J factor
        #if(jay==-1):
        #    sys.exit(100)
        #print to,denoa,da,dca,ca,d1,Fraction(d1**2.).limit_denominator()
        #print 1./(1.*to*denoa**2.)

        #factor out common factors
        if(d4num.denominator==d5num.denominator):
            n=gcd(d4num.numerator,d5num.numerator) #find greatest common divisor
            if(n>1):
                #ay*=n**2.
                #print d4num,d5num,n
                d4,d4num=self.AdjustFraction(d4num/n,1,op=False) #factor this from chemical shifts...
                d5,d5num=self.AdjustFraction(d5num/n,1,op=False) #factor this from chemical shifts...
                jay*=n**2.
                #sys.exit(100)



        if(numpy.fabs(d4num**2*jay-d1**2.*bigFac)>1E-6 or numpy.fabs(-1*facA*d4num*d5num*2.*jay-dca*bigFac)>1E-6 or numpy.fabs(d5num**2.*jay-ca*bigFac)>1E-6):
            print 'shit'
            print 'a1',d4num**2./to/denoa**2.,d1**2.
            print 'a2',-1*facA*d4num*d5num*2./(to*1.)/denoa**2.,dca
            print 'a3',d5num**2./to/denoa**2.,ca
            sys.exit(100)

        return d4,d4num,tA,d5,d5num,jay





    #aim to complete the square
    #between dipolar and CSA terms
    #to emphasise TROSY effects
    def CompleteSquare(self,result,terms,freqs):

        #make sure all terms are present
        for term in terms:
            if(term not in result.keys()):
                return result
            for freq in freqs:
                if(freq not in result[term].keys()):
                    return result
                    
        res=copy.deepcopy(result) #lets get to work.

        dipVal,dipNum,tA,csaVal,csaNum,jfac=self.FactorTROSY(res,terms,freqs[0])


        csa=self.FormatText(terms[1])
        dip=self.FormatText(terms[0])

        if(len(freqs)==1):
            newkeyA=str('\\left('+dipVal+dip+tA+csaVal+csa+'\\right)^2') #J(wc)
            res[newkeyA]={}
            res[newkeyA][freqs[0]]=jfac #/denoa #danewnum**2./denoa
            #print newkeyA,'=',jfac
            return res

        dipVal2,dipNum2,tA2,csaVal2,csaNum2,jfac2=self.FactorTROSY(res,terms,freqs[1])


        #newkeyA=str('('+dipVal+dip+tA+csaVal+csa+')^2') #J(wc)
        #res[newkeyA]={}
        #res[newkeyA][freqs[0]]=jfac #/denoa #danewnum**2./denoa


        d6,d6num=self.AdjustFraction(csaNum/dipNum*dipNum2,1,op=False) #factor this from chemical shifts...

        newkeyA=str('\\left('+dipVal2+dip+tA+d6+csa+'\\right)^2') #J(wc)
        res[newkeyA]={}
        res[newkeyA][freqs[0]]=jfac*dipNum**2./dipNum2**2 #/denoa #danewnum**2./denoa
        #print newkeyA,'=',res[newkeyA][freqs[0]],freqs[0]

        newkeyB=str('\\left('+dipVal2+dip+tA+csaVal2+csa+'\\right)^2') #J(wc)
        if(newkeyB not in res.keys()):
            res[newkeyB]={}
        res[newkeyB][freqs[1]]=jfac2 #/denoa #danewnum**2./denoa
        #print newkeyB,'=',res[newkeyB][freqs[1]],freqs[1]
        #sys.exit(10)

        return res



    #print symbolic rate expression
    def ShowRate(self,Pre,Post,latex=False,sq=False):
        result=self.resSet[Pre][Post]
        lines=[]
        if(latex):

            lon='\nRate for $'+Post.replace('a','\\alpha ').replace('b','\\beta ')+'$ on $'+Pre.replace('a','\\alpha ').replace('b','\\beta ')+'$ :'
            lines.append(lon)
            lon='\\begin{equation}\\begin{array}{rl}'
            lines.append(lon)
            lon='R_{'+Pre.replace('a','\\alpha ').replace('b','\\beta ')+','+Post.replace('a','\\alpha ').replace('b','\\beta ')+'}='
            lines.append(lon)
        else:
            lines.append('Rate for '+Post+' on '+Pre+':')
            #print 'Rate for',Post,'on',Pre,':'

        if(sq):
            #group specified gropus of terms and  frequencies.
            #complex the square by taking the coefficients
            #of the 2nd and 3rd terms.

            #res=self.CompleteSquare(res,('d^2P0P0','cA^2','dcAP0'),('J(0,wC)','J(0,0)'))

            terms=[]
            terms.append('d^2P0P0')
            terms.append('cA^2')
            terms.append('dcAP0')
            freqs=[]
            freqs.append('J(0,wC)')
            freqs.append('J(0,0)')
            res=self.CompleteSquare(result,terms,freqs)


            terms=[]
            terms.append('e^2')
            terms.append('cX^2P2P2')
            terms.append('ecXP2')
            freqs=[]
            freqs.append('J(4,wH)')
            res=self.CompleteSquare(res,terms,freqs)

            terms=[]
            terms.append('e^2')
            terms.append('cX^2P2P2')
            terms.append('ecXP2')
            freqs=[]
            freqs.append('J(4,0)')
            res=self.CompleteSquare(res,terms,freqs)


            terms=[]
            terms.append('e^2')
            terms.append('cX^2P1P1')
            terms.append('ecXP1')
            freqs=[]
            freqs.append('J(1,wH)')
            res=self.CompleteSquare(res,terms,freqs)


            terms=[]
            terms.append('e^2')
            terms.append('cX^2P0P0')
            terms.append('ecXP0')
            freqs=[]
            freqs.append('J(0,wH)')
            res=self.CompleteSquare(res,terms,freqs)


            terms=[]
            terms.append('d^2')
            terms.append('cA^2')
            terms.append('dcA')
            freqs=[]
            freqs.append('J(0,0)')
            freqs.append('J(0,wC)')
            res=self.CompleteSquare(res,terms,freqs)

            terms=[]
            terms.append('e^2')
            terms.append('cX^2')
            terms.append('ecX')
            freqs=[]
            freqs.append('J(0,wH)')
            freqs.append('J(4,wH)')
            res=self.CompleteSquare(res,terms,freqs)


            terms=[]
            terms.append('d^2P0P0')
            terms.append('cX^2P0P0')
            terms.append('dcXP0P0')
            freqs=[]
            freqs.append('J(0,wH)')
            freqs.append('J(0,0)')
            res=self.CompleteSquare(res,terms,freqs)
            terms=[]
            terms.append('d^2P2P2')
            terms.append('cX^2P2P2')
            terms.append('dcXP2P2')
            freqs=[]
            freqs.append('J(4,wH)')
            freqs.append('J(4,0)')
            res=self.CompleteSquare(res,terms,freqs)
            terms=[]
            terms.append('d^2P1P1')
            terms.append('cX^2P1P1')
            terms.append('dcXP1P1')
            freqs=[]
            freqs.append('J(1,wH)')
            freqs.append('J(1,0)')
            res=self.CompleteSquare(res,terms,freqs)


            terms=[]
            terms.append('e^2')
            terms.append('f^2')
            terms.append('fe')
            freqs=[]
            freqs.append('J(0,wH)')
            res=self.CompleteSquare(res,terms,freqs)

            terms=[]
            terms.append('f^2')
            terms.append('g^2')
            terms.append('gf')
            freqs=[]
            freqs.append('J(0,wH)')
            res=self.CompleteSquare(res,terms,freqs)


            terms=[]
            terms.append('e^2')
            terms.append('cX^2')
            terms.append('ecX')
            freqs=[]
            freqs.append('J(0,wH)')
            res=self.CompleteSquare(res,terms,freqs)

            terms=[]
            terms.append('e^2')
            terms.append('cX^2')
            terms.append('ecX')
            freqs=[]
            freqs.append('J(0,0)')
            res=self.CompleteSquare(res,terms,freqs)


            terms=[]
            terms.append('f^2')
            terms.append('cX^2')
            terms.append('fcX')
            freqs=[]
            freqs.append('J(0,wH)')
            res=self.CompleteSquare(res,terms,freqs)

            terms=[]
            terms.append('f^2')
            terms.append('cX^2')
            terms.append('fcX')
            freqs=[]
            freqs.append('J(0,0)')
            res=self.CompleteSquare(res,terms,freqs)




        else: #use the existing dictionary.
            res=result

        cnt=0 #count how many active lines
        for dd in res.keys(): #for constants...
            fin=''
            for j1 in res[dd].keys(): #for frequencies...
                val=res[dd][j1] #get numerical values...
                if(latex):
                    fin+=self.FormatSpecLatex(Fraction(val.real).limit_denominator(),j1) #format forprinting
                else:
                    fin+=self.FormatSpecSh(Fraction(val.real).limit_denominator(),j1) #format forprinting
            if(len(fin)!=0):
                cnt+=1
                if(latex):
                    dl=self.FormatDLatex(dd)
                    lon='& + '+dl+' \\left['+fin+'\\right]\\\\'
                    #lon=lon.replace('J(0,0)','tc')
                    lines.append(lon)
                else:
                    lines.append('+ '+dd+' [ '+fin+' ]')
                    #print '+ '+dd,'[',
                    #print fin,' ]'


        if(cnt==0):
            if(latex):
                lon='0\\\\'
                lines.append(lon)
            else:
                lines.append('0')
                #print 0


        if(latex):
            lon='\end{array}\end{equation}'
            lines.append(lon)
            #print lines
            outy=open(self.LATEXFILE,'a')

            #write the first three then sort the rest.
            outy.write('%s\n' % lines[0])
            lines.pop(0)
            outy.write('%s\n' % lines[0])
            lines.pop(0)
            outy.write('%s\n' % lines[0])
            lines.pop(0)

            lines=sorted(lines)

            for i,line in enumerate(lines):

                line=line.replace('J\\left(0,0\\right)','\\tau_c')
                
                #DS NEW:
                line=line.replace('S_inf', 'S_{\\infty}')
                line=line.replace('S_tm_2', 'S_{\\frac{\\tau_{m}}{2}}')
                line=line.replace('S_tm', 'S_{\\tau_{m}}')
                
                
                #print line
                outy.write('%s\n' % line)
            outy.close()

        else:
            print lines[0]
            lines.pop(0)
            lines=sorted(lines)
            for line in lines:
                print line.replace('J(0,0)','tc')



    def AddCalc(self,Pre,Post,jj,dd,v):
        if(dd not in self.resCalc[Pre][Post].keys()):
            self.resCalc[Pre][Post][dd]={}
        if(jj not in self.resCalc[Pre][Post][dd].keys()):
            self.resCalc[Pre][Post][dd][jj]={}
        self.resCalc[Pre][Post][dd][jj]=v


    def CalcCarbonRate(self,Pre,Post,typ):

        n=Pre.count('a')
        m=Pre.count('b')
        N=n+m

        if(Pre not in self.resCalc.keys()):
            self.resCalc[Pre]={}
        if(Post not in self.resCalc[Pre].keys()):
            self.resCalc[Pre][Post]={}

        if(typ=='rot'):

            phiM=2.*Pi/N
            #constant for all XX vectors for rotating AXN
            psX={}
            psX[0]=(3*Cos(Pi/2.)**2.-1)/2.
            psX[1]=(3*Cos(Pi/2.)*Sin(Pi/2.))
            psX[2]=(3*Sin(Pi/2.)*Sin(Pi/2.))
            for q in 0,1,2:
                y=2./(7.*q**2.-3.*q+2.)  #1,1/3,1/12
                ps='P'+str(q)+'P'+str(q) 

                if(q==0):
                    fac=1.
                elif(q==2 and N==2):
                    fac=1.
                else:
                    #fac=Sin(Pi*q)/Tan(Pi*q/2.) * (-1)**q-1
                    fac=-1.

                #for e->cX cross correlation. #NEED TO CHECK THIS.
                #phiM2=30.
                phiM2=2.*Pi/N/(N-1.) #works for N=2 and N=3

                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',0)'    ,'d^2'+ps,y*1/5.*(N+((n-m)**2.-N)*Cos(q*phiM)) )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC)'   ,'d^2'+ps,y*3/20.*(N+((n-m)**2.-N)*Cos(q*phiM)))
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC-wH)','d^2'+ps,y*1/20.*(N))
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'   ,'d^2'+ps,y*3/20.*(N))
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC+wH)','d^2'+ps,y*6/20.*(N))
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',0)'    ,'cA^2',y*4/5.)
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC)'   ,'cA^2',y*3/5.)
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',0)'    ,'dcA'+ps[:-2],y*4/5.*( ((n-m))*Cos(q*phiM)))
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC)'   ,'dcA'+ps[:-2],y*3/5.*( ((n-m))*Cos(q*phiM)))
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',2wH)'  ,'e^2',psX[q]**2.*y*3/10.*(n*(n-1)+m*(m-1)))
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'  ,'cX^2'+ps,y*3/5.*(N))
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'  ,'e^2',psX[q]**2.*y*3/10.*(N*(N-1)*0.5+ 0.5*Cos(q*phiM)*(N-2)*((n-m)**2.-N)))


                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'   ,'ecX'+ps[:-2],psX[q]*y*4.*3./20.*(((n*(n-1)-m*(m-1))*fac*(Cos(q*(phiM2))))))

        else: #expressions for static geometry, using co-ordinates.

            P0ax=-1/(N-1.)  #(109.47o)
            #P0xx60=-1/8.  #(60o)
            #P0xx45= 1/4.  #(60o)
            #P0xx180= 1.  #(180o)
            #P0xx90= -1/2.  #(90o)
            #P0xc=(3*Cos(((180-109.47122063449068)/2 )/180.*Pi)**2.-1)*0.5 # (35.3o)
            #P0sqrt3=0.     # (54.7)
            #P0sq2sq3=0.5   #  (35.2o)

            self.AddCalc(Pre,Post,'J(0,0)'    ,'d^2',1/5. *(N+((n-m)**2.-N)*P0ax) )
            self.AddCalc(Pre,Post,'J(0,wC)'   ,'d^2',3/20.*(N+((n-m)**2.-N)*P0ax) )
            self.AddCalc(Pre,Post,'J(0,wC-wH)','d^2',1/20.*N)
            self.AddCalc(Pre,Post,'J(0,wH)'   ,'d^2',3/20.*N)
            self.AddCalc(Pre,Post,'J(0,wC+wH)','d^2',6/20.*N)
            self.AddCalc(Pre,Post,'J(0,0)'    ,'cA^2',4/5.)
            self.AddCalc(Pre,Post,'J(0,wC)'   ,'cA^2',3/5.)
            self.AddCalc(Pre,Post,'J(0,0)'    ,'dcA',4/5.*( (n-m)*P0ax))
            self.AddCalc(Pre,Post,'J(0,wC)'   ,'dcA',3/5.*( (n-m)*P0ax))
            self.AddCalc(Pre,Post,'J(0,wH)'   ,'cX^2',3/5.*(N))

            #next, terms involving proton-proton dipolar
            #these are calculaed from the geometry.
            #we need to calculate propertie of one corner.
            #we need to know how many XX vectors of which types converge at a corner (total=(N-1)*N/2 )
            #we need to know the angles between the different XX vector types at the corner. (xx-xx) (total=(N-1)*(N-2))
            #we need to know the angle between the XX vector and the AX vector (xx-xCSA) (total=(N-1)*N/2)

            #print 'doing 2wH'
            for key,val in self.CROSSdiag.items(): #there are N-1 of these
                #print key,val,(n*(n-1)+m*(m-1))/(N-1)
                self.AddCalc(Pre,Post,'J(0,2wH)'  ,key,3/10.*(val*(n*(n-1)+m*(m-1))/(N-1.) ) )
            #print 'doing wH'
            for i,v1 in enumerate(self.edict.keys()):
                for j,v2 in enumerate(self.edict.keys()):
                    if(j>=i):
                        if(i==j):
                            lab=v1+'^2'
                        else:
                            if(v1>v2):
                                lab=v1+v2
                            else:
                                lab=v2+v1
                        num=0
                        #print lab,
                        if(lab in self.CROSSdiag.keys()): #there are N-1 of these
                            num+=self.CROSSdiag[lab]*N*(N-1.)/2./(N-1.)
                        if(lab in self.CROSSxx.keys()): #there are (N-1)(N-2) of these
                            num+=self.CROSSxx[lab]*((n-m)**2.-N)*(N-2)/((N-1.)*(N-2.))
                        #print num
                        self.AddCalc(Pre,Post,'J(0,wH)'   ,lab,3/10.*(num))
            for key,val in self.CROSSax.items(): #val/N-1 is average value.
                self.AddCalc(Pre,Post,'J(0,wH)'   ,key+'cX',3/5.*(((n*(n-1)-m*(m-1))*(val)/(N-1.))))
            #print self.resCalc[Pre][Post]
            
            """
            if(self.baseTag=='CH4'): #tetrahedron
                self.AddCalc(Pre,Post,'J(0,2wH)'  ,'e^2',3/10.*(n*(n-1)+m*(m-1)))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'e^2',3/20.*(N*(N-1)+ P0xx60*(N-2)*((n-m)**2.-N)))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'ecX',3/5.*(((n*(n-1)-m*(m-1))*P0xc)))

            if(self.baseTag=='CH6'): #octahedron
                #cube e is long, f is short
                self.AddCalc(Pre,Post,'J(0,2wH)'  ,'e^2',3/10.*(1.*(n*(n-1)+m*(m-1))/(N-1.)))
                self.AddCalc(Pre,Post,'J(0,2wH)'  ,'f^2',3/10.*(4.*(n*(n-1)+m*(m-1))/(N-1.)))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'f^2',3/10.*(4*(N/2.)+ (4*P0xx60+2*P0xx90)*((n-m)**2.-N)/(N-1.) ))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'fe' ,3/10.*(          (4*P0xx45         )*((n-m)**2.-N)/(N-1.) ))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'e^2',3/10.*(1*(N/2.)  ))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'ecX',3/5.*(((n*(n-1)-m*(m-1))*(2*P0sq2sq3)/(N-1))))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'fcX',3/5.*(((n*(n-1)-m*(m-1))*( P0xx180  )/(N-1))))
                

            if(self.baseTag=='CH8'): #cube
                #cube e is short, f is middle, g is long
                self.AddCalc(Pre,Post,'J(0,2wH)'  ,'e^2',3/10.*(3. *(n*(n-1)+m*(m-1))/(N-1)))
                self.AddCalc(Pre,Post,'J(0,2wH)'  ,'f^2',3/10.*(3. *(n*(n-1)+m*(m-1))/(N-1)))
                self.AddCalc(Pre,Post,'J(0,2wH)'  ,'g^2',3/10.*(1. *(n*(n-1)+m*(m-1))/(N-1)))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'e^2',3/10.*(3.*(N/2.)+ (3.*P0xx90  )*((n-m)**2.-N)/(N-1.) ))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'f^2',3/10.*(3.*(N/2.)+ (3.*P0xx60  )*((n-m)**2.-N)/(N-1.) ))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'g^2',3/10.*(1.*(N/2.)    ))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'fe',3/10.*(            (6*P0xx45+3*P0xx90)*((n-m)**2.-N)/(N-1.) ))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'gf',3/10.*(            (3*P0sq2sq3) *((n-m)**2.-N)/(N-1.) ))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'eg',3/10.*(            (3*P0sqrt3)*((n-m)**2.-N)/(N-1.) ))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'ecX',3/5.*((3*P0sqrt3)*((n*(n-1)-m*(m-1))/(N-1))))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'fcX',3/5.*((3*P0sq2sq3)*((n*(n-1)-m*(m-1))/(N-1))))
                self.AddCalc(Pre,Post,'J(0,wH)'   ,'gcX',3/5.*(( P0xx180  )*((n*(n-1)-m*(m-1))/(N-1))))

            """


    def CalcProtonRate(self,Pre,Post,typ):
        if(Pre not in self.resCalc.keys()):
            self.resCalc[Pre]={}
        if(Post not in self.resCalc[Pre].keys()):
            self.resCalc[Pre][Post]={}

        n=0;m=0
        if(Pre.count('a')>0):
            n=1
        if(Pre.count('b')>0):
            m=1
        N=self.baseType.count('H')

        if(typ=='rot'):
            #constant for all XX vectors for rotating AXN
            psX={}
            psX[0]=(3*Cos(Pi/2.)**2.-1)/2.
            psX[1]=(3*Cos(Pi/2.)*Sin(Pi/2.))
            psX[2]=(-3*Sin(Pi/2.)*Sin(Pi/2.))
            phiM=2.*Pi/N
            for q in 0,1,2:
                y=2./(7.*q**2.-3.*q+2.)
                ps='P'+str(q)+'P'+str(q)
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',0)'    ,'d^2'+ps,y*1/5.  )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC)'   ,'d^2'+ps,y*3/20.*N )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC-wH)','d^2'+ps,y*1./20.*(N +(N-1)*Cos(q*phiM))) 
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'   ,'d^2'+ps,y*3/20. )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC+wH)','d^2'+ps,y*6/20.*(N+(N-1)*Cos(q*phiM)) )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',0)'   ,'cX^2'+ps,y*4/5.)
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'   ,'cX^2'+ps,y*3/5.)
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',0)'   ,'dcX'+ps,y*4/5.*(n-m))  #only cross correlation
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'   ,'dcX'+ps,y*3/5.*(n-m)) #only cross correlation
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC)'   ,'cA^2',y*3/5.)
                #deal with XX dipolar and cross correlation.
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',2wH)'  ,'e^2',psX[q]**2.*y*3/10.*(N-1. )) 
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',0)'    ,'e^2',psX[q]**2.*y*9/20.*(N-1. ))
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'   ,'e^2',psX[q]**2.*y*3/10.*10./4.*(N-1.))            

        else: #expressions for static geometry, using co-ordinates.
            self.AddCalc(Pre,Post,'J(0,0)'    ,'d^2',1/5.  )
            self.AddCalc(Pre,Post,'J(0,wC)'   ,'d^2',3/20.*N )
            self.AddCalc(Pre,Post,'J(0,wC-wH)','d^2',1/20.*(N-1) )
            self.AddCalc(Pre,Post,'J(0,wH)'   ,'d^2',3/20. )
            self.AddCalc(Pre,Post,'J(0,wC+wH)','d^2',6/20.*(N-1) )
            self.AddCalc(Pre,Post,'J(0,0)'   ,'cX^2',4/5.)
            self.AddCalc(Pre,Post,'J(0,wH)'   ,'cX^2',3/5.)
            self.AddCalc(Pre,Post,'J(0,0)'   ,'dcX',4/5.*(n-m))  #only cross correlation
            self.AddCalc(Pre,Post,'J(0,wH)'   ,'dcX',3/5.*(n-m)) #only cross correlation
            self.AddCalc(Pre,Post,'J(0,wC)'   ,'cA^2',3/5.)
            #deal with XX dipolar and cross correlation.
            for key,val in self.CROSSdiag.items():
                self.AddCalc(Pre,Post,'J(0,2wH)'  ,key,3/10.*(val ) )
                self.AddCalc(Pre,Post,'J(0,0)'    ,key,9/20.*(val ) )
                self.AddCalc(Pre,Post,'J(0,wH)'   ,key,3/10.*10./4.*(val))            

    def CalcCzHzRate(self,Pre,Post,typ):
        print 'Calc CzHzRate'
        if(Pre not in self.resCalc.keys()):
            self.resCalc[Pre]={}
        if(Post not in self.resCalc[Pre].keys()):
            self.resCalc[Pre][Post]={}

        n=0;m=0
        if(Pre.count('a')>0):
            n=1
        if(Pre.count('b')>0):
            m=1
        N=self.baseType.count('H')

        if(typ=='rot'):
            #constant for all XX vectors for rotating AXN
            psX={}
            psX[0]=(3*Cos(Pi/2.)**2.-1)/2.
            psX[1]=(3*Cos(Pi/2.)*Sin(Pi/2.))
            psX[2]=(3*Sin(Pi/2.)*Sin(Pi/2.))
            phiM=2.*Pi/N
            for q in 0,1,2:
                y=2./(7.*q**2.-3.*q+2.)
                ps='P'+str(q)+'P'+str(q)
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC)'   ,'d^2'+ps,y*3/10.*(N+ 2.*(N-1.)*Cos(q*phiM)) )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC-wH)','d^2'+ps,y*1/10.*(N-1) )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'   ,'d^2'+ps,y*3/10. )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC+wH)','d^2'+ps,y*6/10.*(N-1) )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'   ,'cX^2'+ps,y*6/5.)
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC)'   ,'cA^2',y*6/5.)
                #deal with XX dipolar and cross correlation.
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'   ,'e^2',psX[q]**2.*y*3/10.*(N-1) )           
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',2wH)'  ,'e^2',psX[q]**2.*y*3/10.*4*(N-1)  )

        else: #expressions for static geometry, using co-ordinates.
            self.AddCalc(Pre,Post,'J(0,wC)'   ,'d^2',3/20.*(N-2)*2 )
            self.AddCalc(Pre,Post,'J(0,wC-wH)','d^2',1/10.*(N-1) )
            self.AddCalc(Pre,Post,'J(0,wH)'   ,'d^2',3/10. )
            self.AddCalc(Pre,Post,'J(0,wC+wH)','d^2',6/10.*(N-1) )
            self.AddCalc(Pre,Post,'J(0,wH)'   ,'cX^2',6/5.)
            self.AddCalc(Pre,Post,'J(0,wC)'   ,'cA^2',6/5.)
            #deal with XX dipolar and cross correlation.
            for key,val in self.CROSSdiag.items():
                self.AddCalc(Pre,Post,'J(0,2wH)'  ,key,3/10.*4*(val ) )
                self.AddCalc(Pre,Post,'J(0,wH)'   ,key,3/10.*(val))            


    def CalcCzHzHzRate(self,Pre,Post,typ):
        print 'Calc CzHzHzRate'
        if(Pre not in self.resCalc.keys()):
            self.resCalc[Pre]={}
        if(Post not in self.resCalc[Pre].keys()):
            self.resCalc[Pre][Post]={}

        n=0;m=0
        if(Pre.count('a')>0):
            n=1
        if(Pre.count('b')>0):
            m=1
        N=self.baseType.count('H')

        if(typ=='rot'):
            #constant for all XX vectors for rotating AXN
            psX={}
            psX[0]=(3*Cos(Pi/2.)**2.-1)/2.
            psX[1]=(3*Cos(Pi/2.)*Sin(Pi/2.))
            psX[2]=(-3*Sin(Pi/2.)*Sin(Pi/2.))
            phiM=2.*Pi/N
            for q in 0,1,2:
                y=2./(7.*q**2.-3.*q+2.)
                ps='P'+str(q)+'P'+str(q)
                #P0ax=-1/(N-1.)
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC)'   ,'d^2'+ps,y*3/10.*(N+(N-2)*4.*Cos(q*phiM)) )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'   ,'d^2'+ps,y*3/10.*2. )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC-wH)','d^2'+ps,y*1/10.*(N-2) )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC+wH)','d^2'+ps,y*6/10.*(N-2) )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'   ,'cX^2'+ps,y*12/5.)
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wC)'   ,'cA^2',y*6/5.)
                #deal with XX dipolar and cross correlation.
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',2wH)'  ,'e^2',psX[q]**2.*y*3/10.*8.*(N-2.) )
                self.AddCalc(Pre,Post,'J('+str(int(q**2.))+',wH)'   ,'e^2',psX[q]**2.*y*3/10.*(2*(N-1)+4.*(N-2)*Cos(q*phiM)  ))

        else: #expressions for static geometry, using co-ordinates.
            #P0ax=-1/(N-1.)
            self.AddCalc(Pre,Post,'J(0,wC)'   ,'d^2',3/10.*(N -(N-2)*4./(N-1.)  ))
            self.AddCalc(Pre,Post,'J(0,wH)'   ,'d^2',3/10.*2. )
            self.AddCalc(Pre,Post,'J(0,wC-wH)','d^2',1/10.*(N-2) )
            self.AddCalc(Pre,Post,'J(0,wC+wH)','d^2',6/10.*(N-2) )
            self.AddCalc(Pre,Post,'J(0,wH)'   ,'cX^2',12/5.)
            self.AddCalc(Pre,Post,'J(0,wC)'   ,'cA^2',6/5.)
            #deal with XX dipolar and cross correlation.
            for key,val in self.CROSSdiag.items():
                self.AddCalc(Pre,Post,'J(0,2wH)'  ,key,3/10.*8.*(N-2.)*val/(N-1.) )
            #print 'doing wH'
            for i,v1 in enumerate(self.edict.keys()):
                for j,v2 in enumerate(self.edict.keys()):
                    if(j>=i):
                        if(i==j):
                            lab=v1+'^2'
                        else:
                            if(v1>v2):
                                lab=v1+v2
                            else:
                                lab=v2+v1
                        num=0
                        if(lab in self.CROSSdiag.keys()): #add diagonals..
                            num+=self.CROSSdiag[lab]*2.
                        if(lab in self.CROSSxx.keys()): #add
                            num+=self.CROSSxx[lab]*8./(N-1.)
                        if(num!=0):
                            self.AddCalc(Pre,Post,'J(0,wH)'   ,lab,3/10.*(num))
            

    def CalcCrossCorrVec(self,vv):
        #calculating 
        veco=vv[1:]-vv[0] #all XX vectors
        norma=numpy.linalg.norm(veco,axis=1) #normalise XX vectors
        typp=[] #map back to distances
        for n in norma:
            for key,vals in self.edict.items():
                if(numpy.fabs(vals-n)<1E-6):
                    typp.append(key)
        diag={}
        for v in typp:
            try:
                diag[v+'^2']+=1
            except:
                diag[v+'^2']=1

        print 'Diagonal vectors around an average vertex:'
        for key,vals in diag.items():
            print key,vals


        xx={}
        xx2={}
        for i,v1 in enumerate(veco):
            for j,v2 in enumerate(veco):
                if(i>j):
                    #print i,j,sorted(typp[i]+typp[j]),ACos(numpy.dot(v1,v2)/norma[i]/norma[j])/Pi*180.
                    
                    if(typp[i]==typp[j]):
                        lab=typp[i]+'^2'
                    else:
                        if(typp[i]>typp[j]):
                            lab=typp[i]+typp[j]
                        else:
                            lab=typp[j]+typp[i]

                    arg=numpy.dot(v1,v2)/norma[i]/norma[j]
                    if(arg>1):
                        arg=1
                    if(arg<-1):
                        arg=-1
                    ang=ACos(arg)/Pi*180.

                    if(lab not in xx.keys()):
                        xx[lab]={}
                    try: 
                        xx[lab][ang]+=1
                    except:
                        xx[lab][ang]=1
                    try:
                        xx2[lab]+=(3*Cos(ang/180.*Pi)**2.-1)*0.5
                    except:
                        xx2[lab]=(3*Cos(ang/180.*Pi)**2.-1)*0.5

        print 'Angles between XX dipolar vectors:'
        for key,vals in xx.items():
            if(numpy.fabs(xx2[key])<1E-9):
                del xx2[key]
                del xx[key]
            else:
                print key,vals,'P0CosZeta:',xx2[key]
        
        ax={}
        ax2={}
        for i,v1 in enumerate(veco):
            arg=numpy.dot(v1,vv[0])/norma[i]
            if(arg>1):
                arg=1
            if(arg<-1):
                arg=-1
            ang=ACos(arg)/Pi*180.

            if(typp[i] not in ax.keys()):
                ax[typp[i]]={}
            try:
                ax[typp[i]][ang]+=1
            except:
                ax[typp[i]][ang]=1
            try:
                ax2[typp[i]]+=(3*Cos(ang/180.*Pi)**2-1)*0.5
            except:
                ax2[typp[i]]=(3*Cos(ang/180.*Pi)**2-1)*0.5

        print 'Angles between AX(CSA) and XX vectors:'
        for key,vals in ax.items():
            print key,vals,'P0cosZeta:',ax2[key]
            if(numpy.fabs(ax2[key])<1E-9):
                del ax2[key]

        #save dictionaryies
        self.CROSSdiag=diag
        self.CROSSxx=xx2
        self.CROSSax=ax2



    def TestForm(self,Pre,Post,typ):

        if(Pre.count('H')>0):
            if(Pre.count('z')==2):
                self.CalcCzHzRate(Pre,Post,typ)
            elif(Pre.count('z')==3):
                self.CalcCzHzHzRate(Pre,Post,typ)
            else:
                self.CalcProtonRate(Pre,Post,typ)
        else:
            self.CalcCarbonRate(Pre,Post,typ)
        
        print
        print 'Testing formula...'
        print 'Pre:',Pre,' Post:',Post
        for key,vals in self.resSet[Pre][Post].items():
            for koi,vols in vals.items():
                try:
                    a=Fraction(self.resCalc[Pre][Post][key][koi]).limit_denominator()
                    b=Fraction(vols).limit_denominator()
                    if(a!=b):
                        print 'shit!!!'
                    print key,koi,Fraction(vols).limit_denominator(),Fraction(self.resCalc[Pre][Post][key][koi]).limit_denominator()
                except:
                    print key,koi,Fraction(vols).limit_denominator(),'error'
                    print self.resCalc[Pre][Post].keys()
                    print

    ##JACS Motional parameters paper
    def CalcFlemming(self):
        cc=self.consts['cA']
        dch=self.consts['d']
        dhh=self.consts['e']
        P2=self.consts['P2']
        P1=self.consts['P1']
        P0=self.consts['P0']
        J0=self.Jomega(0,0)
        Jw=self.Jomega(self.omegaC)
        tauMe=self.tm
        
        Raaa=1/20.*(2*cc-3*dch*P0)**2.*(4*J0 + 3*Jw)+ tauMe*(0.5 * dch**2.*P1**2.    + 1/8.*dch**2.*P2**2.   +27./16.*dhh**2)
        Raab=1/20.*(2*cc-  dch*P0)**2.*(4*J0 + 3*Jw)+ tauMe*(29./30. * dch**2.*P1**2.+ 29/120.*dch**2.*P2**2.+99./80.*dhh**2)
        Rabb=1/20.*(2*cc+  dch*P0)**2.*(4*J0 + 3*Jw)+ tauMe*(29./30. * dch**2.*P1**2.+ 29/120.*dch**2.*P2**2.+99./80.*dhh**2)
        Rbbb=1/20.*(2*cc+3*dch*P0)**2.*(4*J0 + 3*Jw)+ tauMe*(0.5 * dch**2.*P1**2.    + 1/8.*dch**2.*P2**2.   +27./16.*dhh**2)
        return Raaa,Raab,Rabb,Rbbb



    #MethylTROSY 
    def CalcTugarinov(self):
        dch=self.consts['d']
        P0=self.consts['P0']
        J0=self.Jomega(0,0)
        Raaa=9/5.*(dch*P0)**2.*J0 
        Raab=1/5.*(dch*P0)**2.*J0 
        Rabb=1/5.*(dch*P0)**2.*J0 
        Rbbb=9/5.*(dch*P0)**2.*J0 
        return Raaa,Raab,Rabb,Rbbb
  
    #Calculate numerical rate for a given coherences start and finish
    def CalcRate(self,start,finish,verb='n'):
        ####check to see if rate is in main dictionary. If not, calculate##
        go=1
        if(start not in self.resSet.keys()):
            print start,'not in',self.resSet.keys()
            go=0
        elif(finish not in self.resSet[start].keys()):
            go=0
            
        if(go==0):
            if(verb=='y'):
                print
                print 'Entry unavailable:', start, finish
                print 'Trying to calculate...'
            self.EvalRate(start, finish)
            return self.CalcRate(start, finish, verb=verb) #lets go around again... should have fixed it...
            #could get stuck in an infinite loop...

        #####Calculate rate######
        result=self.resSet[start][finish]
        tast=0
        for dd in result.keys():
            keyval=self.GetKeyVal(dd)
            tmp=0
            for j1 in result[dd].keys():
                Jval=self.GetJval(j1)
                #print 'stage:',dd,keyval,Jval,result[dd][j1] #*Jval*keyval
                tmp+=float(result[dd][j1])*Jval*keyval
            tast+=tmp
        if(verb=='y'):
            print 'Rate constant:',start,'on',finish,'=',tast
        #print
        #print self.tc,self.tm
        #print tast

        return tast


    #parse frequency string to get value for spectral density
    def GetJval(self,test):
        val=1
        if(test=='tm'):
            return self.tm
        if(test=='tm/2'):
            return self.tm/2.
        if(test=='tc'):
            return self.tc

        if(test[0]=='J'):
            tast=test.split('(')[1].split(')')[0].split(',') #get arguments of J term

            #if(len(tast)==2): #if two entries, tc,omega
            #    omStr=tast[1]
            #    if(omStr=='0'): #J(0). 
            #        return self.Jomega(0)
            #
            #    omStr=omStr.replace('wC',str(self.omegaC))
            #    omStr=omStr.replace('wH',str(self.omegaH))
            #    omStr=omStr.replace('wD',str(self.omegaD))
            #    omStr=omStr.replace('wF',str(self.omegaF))
            #    return self.Jomega(float(eval(omStr)))  #evaluate J(omega)
            #else: #if three entries, tm,tc and omega
            omStr=tast[1]
            tmStr=tast[0]

            #tmStr=tmStr.replace('tm','1.')
            if(omStr=='0'): #J(0)
                return self.Jomega(float(eval(tmStr)),0)

            omStr=omStr.replace('wC',str(self.omegaC))
            omStr=omStr.replace('wH',str(self.omegaH))
            omStr=omStr.replace('wD',str(self.omegaD))
            omStr=omStr.replace('wF',str(self.omegaF))
            
            return self.Jomega(float(eval(tmStr)),float(eval(omStr)))

        print 'Could not parse',test
        sys.exit(100)


    #get numerical value for interaction constant
    def GetKeyVal(self,key):
        keyInit=key
        val=1
        keys=self.consts.keys()
        keys.sort(reverse=True,key=len)
        #print keys
        for test in keys:
            vals=self.consts[test]
            #for test,vals in self.consts.items():
            #replace with 2x repeat, ^2 and single.
            keyB=key.replace(test+test,'')
            if(len(keyB)!=len(key)):
                val*=(vals*vals)
                key=keyB
            keyB=key.replace(test+'^2','')
            if(len(keyB)!=len(key)):
                val*=(vals*vals)
                key=keyB
            keyB=key.replace(test,'')
            if(len(keyB)!=len(key)):
                val*=vals
                key=keyB
            if(len(key)==0): #as soon as no terms left in key, return value
                return val

        #if(len(key)!=0):#there are terms left. bum.
        print 'Could not fully parse constants:'
        print 'remaining:',key,'from'
        print keyInit
        print 'current constants:'
        print self.consts
        sys.exit(100)


    #loop over interactions, and loop over spin operator terms of each.
    #used to calculate library (independent of start and finish coherences)
    def CrossCorrHam(self):
        if(self.LATEX):
            outy=open(self.LATEXFILE,'a');
            outy.write('Macromolecular terms only:  %s\n\n' % (self.macro));
            if(self.taum!=0):
                if(self.Woessner):
                    outy.write('Local rotations: Woessner jumping model, axis order = %i\n\n' % (self.n))
                else:
                    outy.write('Local rotations: Diffusive model \n\n' )
            outy.write('Spectral density functions: ')
            outy.write('$$J\\left(\omega\\right)=\\frac{\\tau_c}{1+\\omega^2\\tau_c^2}  $$    ')
            if(self.taum!=0):
                outy.write('$$J\\left(\\tau_1,\omega\\right)=\\frac{\\frac{1}{\\tau_c}+\\frac{1}{\\tau_1}}{1+\\omega^2\\left(\\frac{1}{\\tau_c}+\\frac{1}{\\tau_1}\\right)^2}  $$')
            outy.close()
            
        #self.CYCLIC2=True
        self.CYCLIC2=False
        self.Aind={} #will hold non-zero terms with spectral density
        self.terms=0 
        for i in range(len(self.Aint)): #for each block of terms..
            for j in range(len(self.Aint)): #for each other block of terms...
                for ii in range(len(self.Aint[i])): #for each individal operator...
                    for jj in range(len(self.Aint[j])): #for each individual operator...
                        #First 2 terms are operators - 2 spin Zeeman-operators
                        #Second 2 terms specify what kind of interaction we're looking at
                        #This takes form of "dd, ext" or similar - dd specifies dipolar-dipolar, which is what I'm interested in
                        #"ext" specifies a spin external to the methyl group is invovled, so there is a rotating interaction
                        #Can have combinations of "rot, rot" OR "static, rot", OR "static, static"
                        #Aint is a list of Hamiltonians. Each Hamiltonian corresponds to an interaction - dipolar, CSA, quadrupolar etc. Indexed by i, j
                        #Each of these Hamiltonian terms consist themselves of a list of Zeeman 2-spin operators, which we loop over. Indexed by ii, jj
                        #We also need to specify how these interactions change over time - this is specified by Atyp. Indexed by i, j, as Hamiltonians as opposed to operators correspond to interactions
                        self.CalcSpectralDensity(self.Aint[i][ii],self.Aint[j][jj],self.Atyp[i],self.Atyp[j])
        #som=0
        if(self.CYCLIC2):
            self.Aind3=copy.deepcopy(self.Aind)
            self.term3=copy.deepcopy(self.terms)
            for O1,O1s in self.Aind.items():
                baseop1=self.baseOp[O1]          #true/false array for O1
                typ1=self.baseList[baseop1[0][2]]#atoms for O1 
                opp1=baseop1[0][3]               #operators for O1
                argy1=numpy.argwhere(baseop1[0][2]==True) #numbers for 1
                for O2,O2s in O1s.items():
                    baseop2=self.baseOp[O2]          #true/false array for O1
                    typ2=self.baseList[baseop2[0][2]]#atoms for O1 
                    opp2=baseop2[0][3]               #operators for O1
                    argy2=numpy.argwhere(baseop2[0][2]==True) #numbers for 1
                    
                    if((argy1==argy2).all() ):

                        for O1f,O1sf in self.Aind.items():
                            baseop1f=self.baseOp[O1f]          #true/false array for O1
                            typ1f=self.baseList[baseop1f[0][2]]#atoms for O1 
                            opp1f=baseop1f[0][3]               #operators for O1
                            argy1f=numpy.argwhere(baseop1f[0][2]==True) #numbers for 1
                            
                            if(opp1==opp1f):
                                
                                for O2f,O2sf in O1sf.items():
                                    baseop2f=self.baseOp[O2f]          #true/false array for O1
                                    typ2f=self.baseList[baseop2f[0][2]]#atoms for O1 
                                    opp2f=baseop2f[0][3]               #operators for O1
                                    argy2f=numpy.argwhere(baseop2f[0][2]==True) #numbers for 1
                                    
                                    if((argy1f==argy2f).all() ):
                            
                                        if(opp2f==opp2 and (O1!=O1f or O2!=O2f) ):
                                            for d,ds in O2s.items():
                                                for j,js in ds.items():
                                                    
                                                    #if( (argy2f==argy1f).all() and (argy1==argy2).all() and opp1==opp1f and (O1!=O1f or O2!=O2f) and opp2f==opp2 ):
                                                    #    print O1,O2,O1f,O2f
                                                
                                                    try:  #combines diagonals. Appears to work.
                                            #if( opp1==opp1f and (O1!=O1f or O2!=O2f) and opp2f==opp2 ):
                                                
                                                        del self.Aind3[O1f][O2f][d][j]
                                                        self.Aind3[O1][O2][d][j]+=self.Aind[O1f][O2f][d][j]
                                                        print 'combining',O1,O2,O1f,O2f
                                                        self.term3-=1
                                                        pass
                                                    except:
                                                        pass

            self.Aind=copy.deepcopy(self.Aind3)
                                        


                    #print O1,O2,self.Aind[O1][O2][dd][cont]
        #                #som+=self.Aind[O1][O2][dd][cont]
        #                #if(cont=='J(0,0)'):
        #                #    print O1,O2,dd,cont,self.Aind[O1][O2][dd][cont]
        #print som

        print 'Library entries to evaluate:',self.terms #,self.term3
        #sys.exit(100)

        if(self.CYCLIC): #make a second dictionary if using cyclic mode - need to remove J(0) for HH dipolar.
            self.Aind2=copy.deepcopy(self.Aind)
            for O1,O2s in self.Aind.items():
                baseop1=self.baseOp[O1]          #true/false array for O1
                typ1=self.baseList[baseop1[0][2]]#atoms for O1 
                if(typ1[0]=='H' and typ1[1]=='H'):
                    for O2,ds in O2s.items():
                        baseop2=self.baseOp[O2]          #true/false array for O1
                        typ2=self.baseList[baseop2[0][2]]#atoms for O1 
                        if(typ2[0]=='H' and typ2[1]=='H'):
                            for d,js in ds.items():
                                for j,val in js.items():
                                    if(j=='J(0,0)'):
                                        del self.Aind2[O1][O2][d][j]


    #go through Hamiltonians and work out who can cross correlate
    #evaluate the spectral density integrals
    def CalcSpectralDensity(self,one,two,typ1,typ2):
# Anew.append([ self.o[0][0],a+'z'+b+'z',ii, 0,Y]) #self.o[0][0] = 'sqrt(6) 1/2(1 - 3cos^2t)'
        prefac1,O1,ii1,p1=one #unpack one
        prefac2,O2,ii2,p2=two #unpack two
        if(p1-p2!=0): #if p1!=-p2 (total frequency is zero - NOTE: have not yet conjugated p2!)
            return

        O2,p2=self.GetConjOp(O2,p2) #take complex conjugate of O2 operator
        j1,qq1=self.GetFreq(O1) #from operator, work out evolution frequency and coherence order
        j2,qq2=self.GetFreq(O2) #from operator, work out evolution frequency and coherence order

        #print O1,O2,j1,qq1,j2,qq2

        #impose secular hypothesis (frequencies need to be the same)
        if(j1!=j2): #secular approximation - frequency evolution of two terms must be equal.
            return

        #Get the relevent interaction constant tag. Takes the form ab. EG m^2, md, dm etc.
        #Order *matters* here - md and dm are different! Will correspond to different situations!! BUT we re-order for reasons
        #And this isn't a problem - while md and dm *are* different, this is due to the different number of timescales
        #and hence we will calculate different S_inf, S^2 etc. terms
        cont = self.GetCont(ii1, ii2) 
        contN = cont+'C('+cont+')' #EG: C(ab)
        
        if(typ2[1]=='static'):
            freq='J('+str(0)+','+j1+')'
            cint=contN
#             print 'Adding term: ', O1,O2,cint,freq,numpy.real(prefac1*prefac2/5.)
            self.AddCross(O1,O2,cint,freq,numpy.real(prefac1*prefac2/5.)) #Adds the cross-correlation term to the dictionary
        elif(typ2[1]=='single'):
            if(typ1=='static'):
                freq='J('+str(0)+','+j1+')'
                cint=contN
#                 print 'Adding term: ', O1,O2,cint,freq,numpy.real(prefac1*prefac2/5.)
                self.AddCross(O1,O2,cint,freq,numpy.real(prefac1*prefac2/5.)) #Adds the cross-correlation term to the dictionary
            else:
                for y in 1,0: #Loop over all timescales - Have 2: infnty, tm:
                    aa=1.
                    freq='J('+str(y)+','+j1+')'
                    if(y==0):
                        cint=contN+'S('+cont+')^2'
                    if(y==1):
                        cint=contN+'(1-S('+cont+')^2)'
                    #Can apply approximations for macromolecular limit:
                    if(self.macro):
                        if(y==1):
                            freq='tm'
                        if(j1 not in('wC', 'wD', '0')):
                            aa=0
                    if(aa!=0):
                        self.AddCross(O1,O2,cint,freq,numpy.real(prefac1*prefac2/5.)) #Adds the cross-correlation term to the dictionary
#                         print 'Adding term: ', O1,O2,cint,freq,numpy.real(prefac1*prefac2/5.)
        elif(typ2[1]=='double'):
            if(typ1[1]=='static'):
                freq='J('+str(0)+','+j1+')'
                cint=contN
#                 print 'Adding term: ', O1,O2,cint,freq,numpy.real(prefac1*prefac2/5.)
                self.AddCross(O1,O2,cint,freq,numpy.real(prefac1*prefac2/5.)) #Adds the cross-correlation term to the dictionary
            elif(typ1[1]=='single'):
                for y in 1,0: #Loop over all timescales - Have 2: infnty, tm:
                    aa=1.
                    freq='J('+str(y)+','+j1+')'
                    if(y==0):
                        cint=contN+'S('+cont+')^2'
                    if(y==1):
                        cint=contN+'(1-S('+cont+')^2)'
                    #Can apply approximations for macromolecular limit:
                    if(self.macro):
                        if(y==1):
                            freq='tm'
                        if(j1 not in('wC', 'wD', '0')):
                            aa=0
                    if(aa!=0):
                        self.AddCross(O1,O2,cint,freq,numpy.real(prefac1*prefac2/5.)) #Adds the cross-correlation term to the dictionary
#                         print 'Adding term: ', O1,O2,cint,freq,numpy.real(prefac1*prefac2/5.)
            elif(typ1[1]=='double'):
                for y in 2,1,0: #Loop over all different timescales. We have three - infinity, tm, tm/2
                    aa=1.
                    freq='J('+str(y)+','+j1+')'
                    if(y==0):
                        cint=contN+'S_inf('+cont+')^2' #Add in non-decaying residual
                    if(y==1):
                        cint=contN+'S_tm('+cont+')^2' #Add in decaying residual with timescale tm
                    if(y==2):
                        cint=contN+'S_tm_2('+cont+')^2' #Add in decaying residual with timescale tm/2
                    #Can apply approximations for macromolecular limit:
                    if(self.macro):
                        if(y==1):
                            freq='tm'
                        if(y==2):
                            freq='tm/2'
                        if(j1 not in('wC', 'wD', '0')): #If not one of these low frequency terms, then kill
                            aa=0
                    if(aa!=0):
#                         print 'Adding term: ', O1,O2,cint,freq,numpy.real(prefac1*prefac2/5.)
                        self.AddCross(O1,O2,cint,freq,numpy.real(prefac1*prefac2/5.)) #Adds the cross-correlation term to the dictionary
        else:
            print 'GONE WRONG! Tag is not one of: static, single, double. Check code'
            print typ1[1], typ2[1]
        
    #take groups of symmetric operators and condense into a single term.
    #will be effective if dealing with symmetric operators.
    #rho start and end will need to also be symmetrically distributed.
    def AdjOp(self,O1,O2):
        n=self.en  #symmetric number taken from tag.

        baseop1=self.baseOp[O1]          #true/false array for O1
        typ1=self.baseList[baseop1[0][2]]#atoms for O1 
        opp1=baseop1[0][3]               #operators for O1

        #if O1 has an H in it... set it to 2, and adjust any in O2 accordingly.
        if('H' not in typ1): #cannot condense.
            return 

        pl=n+1

        argy1=numpy.argwhere(baseop1[0][2]==True) #numbers for 1
        if(typ1[0]=='H'):
            ref=argy1[0]
            O1a=typ1[0],str(int(pl)),opp1[1]
            if(typ1[0]=='H'):
                val2=((argy1[1]-ref+pl)  -2)%n+2  #new position of H in 2
            else:
                val2=argy1[1]+1
            O1b=typ1[1],str(int(val2)),opp1[0]
            

        elif(typ1[1]=='H'):
            O1a=typ1[0],str(int(argy1[0]+1)),opp1[0]
            ref=argy1[1]
            O1b=typ1[1],str(int(pl)),opp1[1]
        else:
            return O1,O2
        
        O1n=''
        if(O1a[2]!='E'):
            O1n+=''.join(O1a)
        if(O1b[2]!='E'):
            O1n+=''.join(O1b)
        if(O1n not in self.baseOp.keys()):
            O1n=''
            if(O1b[2]!='E'):
                O1n+=''.join(O1b)
            if(O1a[2]!='E'):
                O1n+=''.join(O1a)

        baseop2=self.baseOp[O2]          #true/false array for O2
        typ2=self.baseList[baseop2[0][2]]#atoms for O2 
        opp2=baseop2[0][3]               #operators for O2            
        argy2=numpy.argwhere(baseop2[0][2]==True) #numbers for 1

        if(typ2[0]=='H'):
            val1=((argy2[0]-ref+pl)  -2)%n+2  #new position of H in 2
        else:
            val1=argy2[0]+1
        O2a=typ2[0],str(int(val1)),opp2[0]

        if(typ2[1]=='H'):
            val2=((argy2[1]-ref+pl)  -2)%n+2  #new position of H in 2
        else:
            val2=argy2[1]+1
        O2b=typ2[1],str(int(val2)),opp2[1]


        O2n=''
        if(O2a[2]!='E'):
            O2n+=''.join(O2a)
        if(O2b[2]!='E'):
            O2n+=''.join(O2b)
        if(O2n not in self.baseOp.keys()):
            O2n=''
            if(O2b[2]!='E'):
                O2n+=''.join(O2b)
            if(O2a[2]!='E'):
                O2n+=''.join(O2a)



        return O1n,O2n

        """
        if(typ1[0]=='H'):
            val2=((argy[1]-ref+2.)  -2)%n+2  #new position of H in 2
        else:
            val2=argy[1]+1
            Ob1=typ1[1],val2,opp1[0]
            #elif(typ1[1]=='H'):
            O1a=typ1[0],argy[1]+1,opp1[0]
            ref=argy[1]
            O1a=typ1[1],'2',opp1[1]

        #else:
        #    return O1,O2



            
        if(typ1[0]=='C' and typ1[1]=='H'):  #if AX
            if(typ2[0]=='C' and typ2[1]=='H'): #if also AX.
                O1n='' #set proton in O1 to position 2
                if(opp1[1]!='E'):
                    O1n+=typ1[1]+'2'+opp1[1]
                if(opp1[0]!='E'):
                    O1n+=typ1[0]+'1'+opp1[0]
                argy2=numpy.argwhere(baseop2[0][2]==True)[1] #numbers for 2
                argy1=numpy.argwhere(baseop1[0][2]==True)[1] #numbers for 1
                ref2=((argy2-argy1+2.)  -2)%n+2  #new position of H in 2
                O2n='' #adjust positions in O2 to cyclically match O1.
                if(opp2[1]!='E'):
                    O2n+=typ2[1]+str(int(ref2))+opp2[1]
                if(opp2[0]!='E'):
                    O2n+=typ2[0]+'1'+opp2[0]
                O1=O1n
                O2=O2n
        if(typ1[0]=='H' and typ1[1]=='H'):
            if(typ2[0]=='H' and typ2[1]=='H'):

                argy2=numpy.argwhere(baseop2[0][2]==True) #numbers for 2
                argy1=numpy.argwhere(baseop1[0][2]==True) #numbers for 1
                #print argy1,baseop1[0][2],baseop1[0]
                ref=max(argy1) #set O1 to 2, subtract O1 from O2 and cycle.
                rof=min(argy1)
                tt=((rof-ref+2.)-2.)%n+2
                O1e=O1[:1]+'2'+O1[2:4]+str(int(tt))+O1[5:]
                O1n=typ1[1]+'2'+opp1[1]+typ1[0]+str(int(tt))+opp1[0]
                if(O1n not in self.baseOp.keys()):
                    O1n=typ1[1]+str(int(tt))+opp1[0]+typ1[0]+'2'+opp1[1]
                    if(O1n not in self.baseOp.keys()):
                        print 'problem with this term:'
                        print O1,O1n
                        sys.exit(100)

                ref2=max(argy2) #set O1 to 2, subtract O1 from O2 and cycle.
                rof2=min(argy2)

                ref2a=((ref2-ref+2.)  -2)%n+2
                ref2b=((rof2-ref+2.)  -2)%n+2

                if(ref2a>ref2b):
                    #O2n=O2[:1]+str(int(ref2a))+O2[2:4]+str(int(ref2b))+O2[5:]
                    O2n=typ2[1]+str(int(ref2a))+opp2[1]+typ2[0]+str(int(ref2b))+opp2[0]
                else:
                    #O2n=O2[:1]+str(int(ref2b))+O2[2:4]+str(int(ref2a))+O2[5:]
                    O2n=typ2[1]+str(int(ref2b))+opp2[1]+typ2[0]+str(int(ref2a))+opp2[0]
                O1=O1n
                O2=O2n                
        """

    #add entry to L library
    def AddCross(self,O1,O2,cont,freq,number):
        #add entry to dictionary
#         self.AddCross(O1,O2,cint,freq,numpy.real(prefac1*prefac2*aa)) #Adds the cross-correlation term to the dictionary

        if(self.CYCLIC2):
            O1,O2=self.AdjOp(O1,O2) #adjust operators if condensing using symmetry.
            #prevent C-H H-H dipolar cross correlation for single quantum.

            baseop1=self.baseOp[O1]          #true/false array for O1
            typ1=self.baseList[baseop1[0][2]]#atoms for O1 
            baseop2=self.baseOp[O2]          #true/false array for O1
            typ2=self.baseList[baseop2[0][2]]#atoms for O1 
            #print typ1,typ2
            #if(typ1[0]=='C' and typ1[1]=='H'):

            if(typ1[0]=='C' and typ1[1]=='H'):  #if AX
                if(typ2[0]=='H' and typ2[1]=='H'): #if also AX.            
                    print 'kill',typ1,typ2
                    return 
            if(typ2[0]=='C' and typ2[1]=='H'):  #if AX
                if(typ1[0]=='H' and typ1[1]=='H'): #if also AX.            
                    print 'kill',typ1,typ2
                    return 

            #if(typ2[0]=='H' and typ2[1]=='H'):  #if AX
            #    if(typ1[0]=='H' and typ1[1]=='H'): #if also AX.            


        if(O1 not in self.Aind.keys()):
            self.Aind[O1]={}
        if(O2 not in self.Aind[O1].keys()):
            self.Aind[O1][O2]={}
        if(cont not in self.Aind[O1][O2].keys()):
            self.Aind[O1][O2][cont]={}
        if(freq not in self.Aind[O1][O2][cont].keys()):
            self.Aind[O1][O2][cont][freq]=0 #create new entry
            self.terms+=1  #increment number of live terms
        self.Aind[O1][O2][cont][freq]+=number #increment Llibrary

    #kill any terms that have a value of zero
    #doesn't seem to need to be used now.
    def IsCrossZero(self):
        term=0
        for O1 in self.Aind.keys():
            for O2 in self.Aind[O1].keys():
                for cont in self.Aind[O1][O2].keys():
                    for freq in self.Aind[O1][O2][cont].keys():
                        if(self.Aind[O1][O2][cont][freq]==0):
                            del self.Aind[O1][O2][freq]
                        else:
                            term+=1
        print 'termy!',term

    #format constant
    def GetCont(self,ii1,ii2):
        if(ii1==ii2): #work out relevant constant - add legendre polynomials
            cont=ii1+'^2'
        else:
            if(ii1>ii2): #sort for easy grouping later on
                cont=ii1+ii2
            else:
                cont=ii2+ii1
        return cont

    #turn frequency into spectral density function - add tauM if required.
    def GetJfreq(self,j1,y):
        #if(y!=0): #adjust if we need tauM
        Jstr='J('+str(y**2)+','+j1+')'
        if(Jstr=='J(0,0)'):
            return 'tc'
        else:
            return Jstr
        #return 'J(0,'+j1+')'

    #transform matrix into new basis
    def Transform(self,P,mat):
        return numpy.dot(numpy.linalg.inv(P),numpy.dot(mat,P))

    #take manually specified basis and calculate self relaxation rates
    def IreducibleBasis(self,old='n'):

        self.GetSymmetry([],self.baseTag) 


        #print self.maps

        #rep=self.GetReducible(test)  #get reduible representation
        #irreds=self.GetIrreducible(rep) #reduce it...
        #reppy=GetRepresentations(irreds) #get representations.
        #Pold=Pinv
        #print Pinv
        try:
            Hval=int(re.findall('\d+', self.baseTag)[0]) #get the number out of the tag
        except:
            Hval=int(self.baseTag.count('H'))

        if(Hval==2):
            self.GetCharTable('C2v')
        elif(Hval==3):
            self.GetCharTable('C3v')    
        elif(Hval==4):
            self.GetCharTable('Td')    
        elif(Hval==6):
            self.GetCharTable('Oh')    
        elif(Hval==8):
            self.GetCharTable('Oh')    
        elif(Hval==12):
            self.GetCharTable('Ih')    
        elif(Hval==20):
            self.GetCharTable('Ih')    


        #right. this is hell of clunky.
        #1. work out number of permuatations of and b
        newbase=[]
        #assemble all permutations of as and bs to represent the Zeeman basis
        for i in range(Hval+1): #for each distinct chemical shift
            bas=numpy.concatenate((numpy.repeat('a',Hval-i),numpy.repeat('b',i)))
            test=list(more_itertools.distinct_permutations(bas)) #get the permutations
            for te in test:
                newbase.append("".join(te))
        seek=len(newbase) #number of permutations.
        #print newbase

        #CH2
        #P=[[1,0,0,0],[0,1,1,0],[0,1,-1,0],[0,0,0,1]]
        #BUT: need to reorder them to fit in with the basis we're using.
        #so diagonalisiation works correctly.
        #needs to line up with matrix entries, eg CH3:
        #  aaa, baa, aba, bba, aab, bab, abb, bbb
        #1                                    bbb
        #2                +bba     +bab  +abb
        #3                -bba     +bab
        #4                +bba     +bab -2abb
        #5     +baa +aba      +aab
        #6          -aba      +aab
        #7    -2baa +aba      +aab
        #8 aaa
        #
        newbase=[]
        newbase.append('a')
        newbase.append('b')
        while(len(newbase)!=seek):
            newnewbase=[]
            for base in newbase:
                newnewbase.append('a'+base)
                newnewbase.append('b'+base)
            newbase=newnewbase
        #print newbase


        #now reduce representation and save result in matrix sorted by newbase
        sims=[] #symmetries
        Pinv=[] #mapping function
        nums={}
        for i in range(Hval+1): #for each distinct chemical shift
            print 
            bas=numpy.concatenate((numpy.repeat('b',Hval-i),numpy.repeat('a',i)))
            test=list(more_itertools.distinct_permutations(bas)) #get the permutations
            print 'reducing:',bas,len(test)
            rep=self.GetReducibleFunc(test)
            irreds=self.GetIrreducible(rep) #reduce it...
            self.PrintIrreducible(irreds,test)
            reppy=self.GetRepresentationsFunc(irreds,test)
            self.PrintRepresentation(reppy,test)
            #print reppy
            
            for key,vals in reppy.items():
                for val in vals:
                    sims.append(key)
                    po=[]
                    for tast in newbase:
                        #print 'tat',tast,'val',vals[0]
                        if(tast in val.keys()):
                            po.append(val[tast][0])
                        else:
                            po.append(0)
                    po=numpy.array(po)
                    Pinv.append(po)

                    if( (len(Pinv)) not in nums.keys()):
                        nums[len(Pinv)]=[]
                    nums[len(Pinv)].append(i)

                    #if(len(sims)==13):
                    #    print len(sims),po
                    #    print len(Pinv)
                    #    print (po==Pinv[-1]).all()
                    #    print key
                    #    for val in vals:
                    #        print val
                    #    sys.exit(10)

                    #for ii,pe in enumerate(Pinv):
                    #    if (po==pe).all():
                    #        print 'shize. matechy!'
                    #        print key,vals
                    #        print po
                    #        print sims
                    #        print len(sims)
                    #        print ii
                    #        sys.exit(100)





        #print Pold
        #print numpy.array(Pinv).shape
        #print numpy.array(Pold).shape
        #sys.exit(100)
        #print Pee
        #incoming is technically inverse
        #Pinv=Pold

        print 'Irreducible representation:'
        #print Pinv
        Pinv=numpy.array(Pinv)
        #print Pinv.shape
        """
        for i,po in enumerate(Pinv):
            for j,pe in enumerate(Pinv):
                if(j>i):
                    if((po==pe).all()):
                        print 'shitty'
                        print i,j
                        print po
                        print pe
                        sys.exit(100)
                    if((po==pe*-1).all()):
                        print 'shotty'
                        print i,j
                        print po
                        print pe
                        sys.exit(100)
        """
        P=numpy.linalg.inv(Pinv)

        #normalise:
        for i in range(len(P)):
            norm=numpy.sum(P[:,i]**2.)
            P[:,i]=P[:,i]/numpy.sqrt(norm)



        EE=numpy.identity(2)  #expand by S nucleus
        Pbig=numpy.kron(EE,P)  #expand 
        PbigInv=numpy.linalg.pinv(Pbig) #invert

        if(self.LATEX):
            outy=open(self.LATEXFILE,'a')
            outy.write('Irreducible representation in terms of operators and symmetry:\\\\')
            outy.close()

        #print labels
        for i in range(len(P)):
            #print 'PPP',i+1,P[i]
            zz=numpy.zeros(len(P))
            zz[i]=1.
            tast=numpy.diag(zz) #*2*numpy.sqrt(8)
            tast=numpy.kron([[0,1],[0,0]],tast)  #make with 'p' on S
            newy=self.Transform(PbigInv,tast) #convert tast back into Eigenbasis

            if(self.SPARSE):
                newy=csr_matrix(newy)
                newy.eliminate_zeros()

            reps=self.WhoAmI(newy,verb='n')#break representation into single spin operators
            if(self.LATEX):
                stry='%s : %s\\left>\\right<%s =& ' % (sims[i],str(i+1),str(i+1))
                for k,rep in enumerate(reps):
                    if(rep[0]>0):
                        stry+='+'
                    else:
                        stry+='-'
                    if(numpy.fabs(rep[0])==1):
                        num=''
                    else:
                        f1=Fraction(rep[0]).limit_denominator()
                        if(len(str(f1.numerator))>4): #square it and hope the fraction looks nicer
                            f1=Fraction(rep[0]**2.).limit_denominator()
                            num='(\\frac{%s}{%s})^{1/2}' % (numpy.abs(f1.numerator),f1.denominator)
                        else:
                            num='\\frac{%s}{%s}' % (numpy.abs(f1.numerator),f1.denominator)                        

                    stry+='%s %s' % (num,self.ExpandLatex(rep[1]))
                    lim=max(1,int(9-Hval))
                    if((k+1)%lim==0):
                        stry+='\\\\ &'

                        #self.Expand(b)
                
                outy=open(self.LATEXFILE,'a')
                outy.write('\\begin{equation}')
                outy.write('\\begin{array}{ll}')
                outy.write(stry)
                outy.write('\\end{array}')
                outy.write('\\end{equation}')
                outy.close()

            self.baseOp[str(i+1)]=reps     #save string of operators (symbolic calcs)
            if(self.SYM==False): #save matrix if not doing this symbolically
                self.base[str(i+1)]=newy
        #sys.exit(100)
        #run calculations on this.


        for i in range(Hval+1): #for each distinct chemical shift
            testbasis=[]
            for j in range(len(P)):
                if(i in nums[j+1]):
                    testbasis.append(str(j+1))
                    for k in range(len(P)):
                        if( i not in nums[k+1]):
                            if(str(j+1) not in self.resSet.keys()):
                                self.resSet[str(j+1)]={}
                            self.resSet[str(j+1)][str(k+1)]={}
                            #print 'killing:',j+1,k+1
            self.EvalRate(testbasis,testbasis,calc='y')  #correct

        testbasis=[]
        for i in range(len(P)):
            testbasis.append(str(i+1))
        #self.EvalRate(testbasis,testbasis,calc='y')  #correct

        #will store final symbolic results


        #    testbasis.append(str(i+1))
        #self.EvalRate(testbasis,testbasis,calc='y')  #correct

        #self.EvalRate('2','2',calc='y')  #correct
        #sys.exit(100)
        #for te in testbasis:
        #    self.EvalRate(te,te,calc='y')  #correct
        return testbasis


        """
        #sys.exit(100)
        #print self.base[str(2)].toarray()
        #print self.Transform(Pbig,self.base[str(2)].toarray())

        reps='H2pC1z'
        nums=[]
        print 'assembling'

        #self.base['H2zC1z'].eliminate_zeros()
        #print self.base['H2zC1z']

        nums=[1,1]
        self.Assemble(('H2zC1z','H3zC1z'),'test',nums=nums) #this is A

        nums=[1,-1]
        self.Assemble(('H2zC1z','H3zC1z'),'tast',nums=nums) #this is B
        #print csr_matrix(self.base['test']).eliminate_zeros().conjugate().eliminate_zeros().conjugate()
        #print self.base['test'] #.eliminate_zeros()

        #print (self.base['test']*self.base[str(2)])
        print 'a'
        print comm(self.base['test'],self.base[str(1)]) #1 is A
        print 'b'
        print comm(self.base['test'],self.base[str(2)]) #2 is probably A?
        print 'c'
        print comm(self.base['test'],self.base[str(3)]) #3 is probably B?
        print 'd'
        print comm(self.base['test'],self.base[str(4)]) #4 is A
        #so: A representation of zz gives zero comm with A rep of rho
        #    B representation of zz gives zero comm with A and B rep of rho
        sys.exit(100)
        print
        print 'c'
        print comm(self.base['tast'],self.base[str(1)])
        print
        print 'd'
        print comm(self.base['tast'],self.base[str(2)])
        sys.exit(100)
        print
        print
        print 'mat',csr_matrix(self.Transform(Pbig,comm(self.base['test'].toarray(),self.base[str(3)].toarray())))
        print
        print 'mat',csr_matrix(self.Transform(Pbig,comm(self.base['tast'].toarray(),self.base[str(3)].toarray())))
        print

        #print self.WhoAmI(comm(self.base['test'],self.base[str(2)]))
        sys.exit(10)
        """



    #create a new basis by specified entries (un-normalised)
    def Assemble(self,reps,lab,nums=[]):
        #print reps,lab
        #Take in a list of representations, "reps" (strings which represent operators) which are then packed together, labelled with "lab"
        vals=[]
        for i,rep in enumerate(reps):
            op=self.GetOp(rep) #expanded representation
            if(len(nums)==0):
                vals.append((1.,op,1,1,op.tostring()))
            else:
                vals.append((nums[i],op,1,1,op.tostring()))
        self.baseOp[lab]=vals
        
        if(self.SYM): #no need to calculate the matrix
            return

        if(len(nums)==0):
            for i,rep in enumerate(reps):
                if(i==0):
                    if(rep in self.base):
                        #newy=copy.deepcopy(self.base[rep])
                        newy=self.base[rep].copy()
                    else:
                        newy=self.MakeOp(rep,save=False)
                else:
                    if(rep in self.base):
                        #newy+=copy.deepcopy(self.base[rep])
                        newy+=self.base[rep].copy()
                    else:
                        newy+=self.MakeOp(rep,save=False)
        else:
            for i,rep in enumerate(reps):
                if(i==0):
                    if(rep in self.base):
                        #newy=copy.deepcopy(self.base[rep])
                        newy=self.base[rep].copy()*nums[i]
                    else:
                        newy=self.MakeOp(rep,save=False)*nums[i]
                else:
                    if(rep in self.base):
                        #newy+=copy.deepcopy(self.base[rep])
                        newy+=self.base[rep].copy()*nums[i]
                    else:
                        newy+=self.MakeOp(rep,save=False)*nums[i]

        if(self.SPARSE):
            newy=csr_matrix(newy)
        #self.base[lab]=copy.deepcopy(newy)
        self.base[lab]=newy.copy()



    #allow operator to evolve for time 'angle'
    #is an operator in the chemical shift eigenbase?
    def CheckFreePrecession(self,H0base,targ,angle):
        H0=base[H0base]
        H1=expm(-complex(0,1)*H0*angle)
        H2=expm(complex(0,1)*H0*angle)
        #print commAux.base[H0targ]
        baso=numpy.dot(H1,numpy.dot(base[targ],H2))
        #if(resy!=0):
        #    print 'Final form of the matrix:',resy,mat
        if(targ[2]=='p'):
            calc=base[targ]*Exp(-complex(0,1)*angle)
        else:
            calc=base[targ]*Exp(complex(0,1)*angle)
        if(numpy.sum(numpy.abs(calc-baso))<1E-9):
            print 'GREAT: operator',targ,'is in the chemical shift Eigenbasis'
        else:
            print 'PROBLEM: operator',targ,'is not an the Eigenbasis'
        #    print baso
        #    print calc



    #tests if [b,c]=fac*a
    def CheckComm(self,fac,a,b,c):
        comm1=comm(self.base[b],self.base[c])
        frac,resy,baso=compare(comm1,self.base)
        print '[',b,',',c,']=',frac,resy
        if(self.IsZero(baso)):
            print 'done'
            print
            return
    
        comm2=comm(self.base[a],comm1)
        frac2,resy2,baso2=compare(comm2,self.base)
        print '[',a,',[',b,',',c,']=',frac2,resy2
        if(self.IsZero(baso2)):
            print 'done'
            print
            return
            
        tr=self.Trace(self.base[c],comm2)   
        
        print '<',c,'| [',a,',[',b,',',c,'] > = ',tr
        if(tr==0):
            print 'done'
            print 
            return

        trD=Sqrt(self.Trace(self.base[c],self.base[c])*self.Trace(self.base[c],self.base[c]))    
        print '(<',c,'|',c,'><',c,'|',c,'>)^1/2 = ',trD
        print '<',c,'| [',a,',[',b,',',c,'] > / (<',c,'|',c,'><',c,'|',c,'>)^1/2= ',tr/trD
        j1,qq1=self.GetFreq(b) #from operator, work out evolution frequency
        face=Fraction((tr*fac/trD*6/5.)).limit_denominator()  #turn numerical factor into fraction
        print face,'J(',j1,')'
        print

    def CondenseSym(self,listy,basis,fuggle):
        if(len(basis)==1 or len(fuggle)<=1):
            return listy,basis,[],fuggle
        print 'Condensing basis of size:',len(basis)


        fuggle=numpy.array(fuggle)
        basis=numpy.array(basis)
        listy=numpy.array(listy)
        #print fuggle.shape
        uni,argy,inv,cnt=numpy.unique(fuggle,return_index=True,return_inverse=True,return_counts=True)

        #print fuggle
        #print len(listy),len(basis),len(fuggle)
        #print uni
        #print argy
        #print inv
        #print cnt

        print 'start:'
        for i,bas in enumerate(basis):
            print bas,listy[i],fuggle[i]
        #print basis
        #print fuggle
        basisNew=basis[argy]
        numsNew=cnt
        fuggleNew=fuggle[argy]
        listyNew=listy[argy]

        print 'end:'
        for i,bas in enumerate(basisNew):
            print bas,listyNew[i],numsNew[i],fuggleNew[i]

        #print 'fdasfdsafd', basisNew
        return listyNew,basisNew,numsNew,fuggleNew

        #for i,bas in enumerate(basis):
        #    #condense based on transform
        #    if(fuggle[i] not in oppy):
        #        oppy.append(fuggle[i])
        #        basisNew.append(bas)
        #        nums.append(1)

    def GetLabel(self,listy,n,maj):
        bis=''
        for i in range(n-1):
            bis+='H'+str(n+1-i)+str(listy[n-2-i])
        bis+='H'+str(2)+maj
        bis+='C1p'
        return bis

    def GetLabelPure(self,listy,n):
        bis=''
        for i in range(n):
            bis+='H'+str(n+1-i)+str(listy[i])
        bis+='C1p'
        return bis

    def GetTransform(self,listy,maj):
        ref=0 #position of 'reference' atom
        losty=numpy.insert(listy,ref,maj)
        foggle1=[]
        for key,val in self.maps.items(): #for each symmetry op
            test=((losty==maj)*1.)   #which dudes are alphas?
            tran=numpy.dot(val,test) #apply symmetry operator
            #if( (losty==('a','a','b','b','b','b')).all()):
            #    print
            #    print 'test:',losty
            #    print key
            #    print 'from:',test
            #    print 'to   :',tran
            #    print 'stay?:',val[ref,ref]==1.
            #    print val
            if( (test==tran).all()): #if the final result is the same as start...
                #if(tran[ref]==1): #have we transformed onto another alpha?
                newOp=key.split('_')[0]
                if(val[ref,ref]==1): #did we move?
                    newOp+='A'
                else:
                    newOp+='B'

                #if newOp not in foggle1:
                foggle1.append(newOp)

                
                #if( (numpy.dot(val.transpose(),test)==test).all()):
                #    foggle1.append(key)

        #test=((losty=='b')*1.) #do the same to the betas.
        #if( (numpy.dot(val,test)==test).all()):
        #foggle2.append(key) 
        return sorted(foggle1)
        
        #foggle1=sorted(foggle1)
        #foggle2=sorted(foggle2)
        #if(numpy.sum((losty=='a')*1.)==2):
        #if(losty[0]=='a' and losty[1]=='a'):
        #    sys.exit(100)
        #if(foggle1!=foggle2):
        #    print 'problem. check this'
        #    sys.exit(100)
        #self.fuggle.append((foggle1,foggle2))
        #return foggle1

    #create bass for AXn where X is spin half:
    #e.g. abb== abb+bab+bba
    def GetBasisAXn(self,n):
        print 'Creating A+ Eigenbasis'
        basStr=[]
        self.figgle={}
        for i in range(n+1): #for each distinct chemical shift

            if(self.CYCLIC):
                if(i==0 or i==n): #at the ends, no simplification
                    en=n
                    imax=i
                else:  #in the middle, we can simplify by aligning one proton
                    en=n-1
                    imax=i-1
            else:
                en=n
                imax=i

            bas=[]  #assemble string with each combination of as and bs


            if(en==n-1 and self.CYCLIC):#remove one minority letter
                if(i>(n+1.)/2.): #work out the minority letter
                    maj='b'
                else:
                    maj='a'
                print 'minority letter:',maj
                if(maj=='a'): #remove an 'a' from the start
                    for j in range(en):
                        if(j>=imax):
                            bas.append('b')
                        else:
                            bas.append('a')
                    lab='a'+"".join(bas)
                else: #remove a 'b' from the end
                    for j in range(en):
                        if(j>=imax+1):
                            bas.append('b')
                        else:
                            bas.append('a')
                    lab="".join(bas)+'b'
            else:
                for j in range(en):
                    if(j>=imax):
                        bas.append('b')
                    else:
                        bas.append('a')
                lab="".join(bas)
                
            print 'bas',bas
            basStr.append(lab)
            print 'lab',lab

            listyo=list(more_itertools.distinct_permutations(bas)) #get the permutations

            fuggle=[] #store operators for each listy
            basis=[]  #store label for each listy
            for listy in listyo: #these are in the forward direction
                if(self.CYCLIC):
                    if(en==n-1): #dealing with inner line
                        if(self.POINT and self.SYM): #work out transforms if needed
                            fuggle.append(self.GetTransform(listy,maj)) 
                        basis.append(self.GetLabel(listy,n,maj))
                    else: #we are dealing with an end operator
                        basis.append(self.GetLabelPure(listy,en))
                else:
                    basis.append(self.GetLabelPure(listy,en))

            print 'creating:',lab,len(basis) #,basis
            print 'basis',basis


            self.Assemble(basis,lab)  #make pre basis operators

            if(self.POINT and self.SYM): #make post operator basis if required
                listyo,basis,nums,fuggle=self.CondenseSym(listyo,basis,fuggle)
                self.figgle[lab]=fuggle #add the condensed operators 
                self.Assemble(basis,lab+'Post',nums=nums) #save the post operators
        print basStr
        return basStr


    #convert cartesian to spherical polar
    def ConvSph(self,v):
        r=Sqrt(v[0]**2.+v[1]**2.+v[2]**2.)
        theta=ACos(v[2]/r)
        if(v[0]==0 and v[1]==0):
            psi=0
        elif(v[0]==0 and v[1]>0):
            psi=numpy.pi/2
        elif(v[0]==0 and v[1]<0):
            psi=3*numpy.pi/2
        elif(v[1]==0 and v[0]>0):
            psi=0
        elif(v[1]==0 and v[0]<0):
            psi=Pi
        else:
            #psi=numpy.arctan(v[1]/v[0])
            psi=ATan2(v[1],v[0])
        while(psi<0):
            psi+=Pi*2.
        while(psi>Pi*2):
            psi-=Pi*2.
        return r,Fraction(theta/Pi).limit_denominator()*numpy.pi,Fraction(psi/Pi).limit_denominator() *Pi

    #add AX dipolar interactions
    #between set of vectors
    def AddDipolarAX(self,vv,rot='static'):
        dlist=[]
        for i,v in enumerate(vv): #get A->X bond vectors
            r,th,ph=self.ConvSph(v)
            dlist.append((str(i+2),'1',th,ph))
        #co-ordinate v[0] for X is spin '2'
        #co-ordinate v[-1] for X in spin [basisSize]
        for d in dlist:
            self.GetDip('H'+d[0],'C'+d[1],'d',beta=d[2],alpha=d[3],rot=rot)

    #add XX dipolar interactions
    #between set of vectors
    def AddDipolarXX(self,vv,rot='static'):
        ef={}
        ref='efghijkl'
        elist=[]
        

        veco=[]
        for i,v1 in enumerate(vv):
            for j,v2 in enumerate(vv):
                if(j>i): #include all pairs of dipolars
                    #print v1,v2
                    v3=numpy.array(v2)-numpy.array(v1) #dipolar vector

                    veco.append(v3)
                    r,th,ph=self.ConvSph(v3)
                    #print i+2,j+2,r,th,ph
                    key=str(r)[:-1] #performing numerical approximation to condense keys...
                    #print key
                    if(key not in ef.keys()):
                        ef[key]=ref[len(ef.keys())]
                    elist.append((str(i+2),str(j+2),th,ph,ef[key]))

        """
        ong=[]
        for i,vec in enumerate(veco):
            for j,voc in enumerate(veco):
                if(i>j):
                    #print i,j
                    #nvec=(vec[0]**2.+vec[1]**2.+vec[2]**2.)**0.5
                    #nvoc=(voc[0]**2.+voc[1]**2.+voc[2]**2.)**0.5
                    #print numpy.linalg.norm(nvoc),nvoc
                    #print numpy.linalg.norm(nvec),nvec
                    P0=(3.*(numpy.dot(vec,voc)/numpy.linalg.norm(vec)/numpy.linalg.norm(voc))**2.-1)/2.
                    #ong.append(P0)
                    #coss= (3.*((vec[0]*voc[0]+vec[1]*voc[1]+vec[2]*voc[2])/nvec/nvoc)**2.-1)*0.5
                    #ong.append(coss)
                    #ang=ACos(coss)
                    #ong.append(ang)
                    ang=ACos(numpy.dot(vec,voc)/numpy.linalg.norm(vec)/numpy.linalg.norm(voc))/Pi*180.
                    if(ang>90):
                        ang=180-ang

                    P0=(3.*(Cos(ang/180.*Pi))**2.-1)/2.
                    if(ang!=90): #why der fek do we ignore the 90o dudes in tetrahedron?
                        ong.append(P0)
                    #print 'aaa',i,j,ang

        print len(veco)
        print len(ong)

        print 'averageP0:',numpy.average(ong)
        """
        #print ((3*Cos(numpy.average(ong))**2.-1)*0.5)
        #print Fraction((3*Cos(numpy.average(ong))**2.-1)*0.5).limit_denominator()

        #tetrahedron
        #12: 60o   4:
        #3: 90o    1:

        for key,vals in ef.items():
            self.edict[vals]=float(key)
        for e in elist:
            #print e[2]
            self.GetDip('H'+e[0],'H'+e[1],e[4],beta=e[2],alpha=e[3],rot=rot)

        self.CalcCrossCorrVec(vv)




        #sys.exit(100)
    #add CSA for X spin
    def AddXCSA(self,vv,rot='static'):
        dlist=[]
        for i,v in enumerate(vv):
            r,th,ph=self.ConvSph(v)
            dlist.append((str(i+2),'1',th,ph))
        for d in dlist:
            self.GetCSA('H'+d[0],'cX',beta=d[2],alpha=d[3],rot=rot)

    #add CSA for A spin
    def AddACSA(self,vv,rot='static'):
        for i,ta in enumerate(self.baseType):
            if(self.baseType[i]=='C'):
                self.GetCSA('C'+str(i+1),'cA',rot=rot)

                
    #Euler-Rodriguez formula
    #https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula
    #specify polar axis of rotation from an axis defined by a point (kx,ky,kz)
    #and a degrees phi
    def DoRot(self,k,phi):
        k = k / Sqrt(numpy.dot(k, k))
        a=Cos(phi/2.)
        b=-k[0]*Sin(phi/2.)
        c=-k[1]*Sin(phi/2.)
        d=-k[2]*Sin(phi/2.)
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

    #https://www.lfd.uci.edu/~gohlke/code/transformations.py.html
    def DoRef(self,point, normal):
        normal = normal / Sqrt(numpy.dot(normal,normal))
        M = numpy.identity(3)
        M[:3, :3] -= 2.0 * numpy.outer(normal, normal)
        #M[:3, 3] = (2.0 * numpy.dot(point[:3], normal)) * normal
        return M

    def WriteGeometry(self,outfile,vv):
        outy=open(outfile,'w')
        vv=numpy.array(vv)
        for i,v in enumerate(vv):
            #print vv[i,0]
            #print vv[i,1]
            #print vv[i,2]
            outy.write('%i\t%f\t%f\t%f\t\n' % (i,0,0,0))
            outy.write('%i\t%f\t%f\t%f\t\n' % (i,vv[i,0],vv[i,1],vv[i,2]))
            outy.write('\n\n')
        outy.close()

    #try rotation and see if it was symmetric
    def TestSym(self,axis,vv,ang,lab): #apply rotation of phi about axis
        #r=self.DoRot(vv[0],120./180.*numpy.pi)
        r=self.DoRot(axis,ang*1./180.*Pi)
        vnew=numpy.dot(vv,r.transpose())
        #self.WriteGeometry('test1',vnew)
        self.GetCnew(vv,vnew,lab)

    #add identity
    def TestSymIdentity(self,vv,lab): #apply rotation of phi about axis
        self.GetCnew(vv,vv,lab)

    #do inverse
    def TestSymInverse(self,vv,lab):
        sqdist=squareform(pdist(vv))
        #print sqdist
        cnew=(numpy.fabs(sqdist-2.)<1E-6)*1.
        #print cnew
        if( (numpy.sum(cnew,axis=0)==numpy.sum(cnew,axis=1)).all()==0 or numpy.sum(cnew)!=cnew.shape[0]):
            print 'shit. inverse broke!'
            sys.exit(100)
        self.maps[lab]=cnew
        #print self.maps[lab]

    #try reflection about point and normal
    def TestSymRef(self,point,normal,vv,lab): #apply reflection plane
        r=self.DoRef(point,normal)
        vnew=numpy.dot(vv,r.transpose()) #rotate each co-ordinate
        self.WriteGeometry('test1',vnew)
        self.GetCnew(vv,vnew,lab)

    #do improper rotation
    def TestSymImp(self,axis,ang,point,vv,lab):
        r1=self.DoRot(axis,ang*1./180*Pi)
        r2=self.DoRef(point,axis)
        r=numpy.dot(r1,r2)
        vnew=numpy.dot(vv,r.transpose()) #rotate each co-ordinate
        #self.WriteGeometry('test1',vnew)
        self.GetCnew(vv,vnew,lab)

    #test if symmetry operation is valid
    def GetCnew(self,vv,vnew,lab):
        cc=cdist(vv,vnew)
        cnew=(cc<1E-6)*1
        for val in self.maps.values():
            if (cnew==val).all() :
                print 'transform is a repeat:',lab,
                for key,vals in self.maps.items():
                    if( (vals==cnew).all()):
                        print 'matches',key
                #return 
        if( (numpy.sum(cnew,axis=0)==numpy.sum(cnew,axis=1)).all()==0 or numpy.sum(cnew)!=cnew.shape[0]):

            print numpy.sum(cnew)
            print cnew.shape[0]
            print cnew
            print 'shit - numbers dont add up:',lab
            #sys.exit(100)
            #return
        print 'adding',lab #,cnew
        self.maps[lab]=cnew

    #read in the character tables
    def GetCharTable(self,sym):
        self.chars={}
        if(sym=='C2v'):
            #C2v	E	C2 (z)	v(xz)	v(yz)	
            self.charval=1,1,1,1
            self.charOps='E','c2','sigv1','sigv2'
            self.chars['A1']=	+1,	+1,	+1,	+1,	#z	x2, y2, z2	z3, x2z, y2z
            self.chars['A2']=	+1,	+1,	-1,	-1,	#Rz	xy	xyz
            self.chars['B1']=	+1,	-1,	+1,	-1,	#x, Ry	xz	xz2, x3, xy2
            self.chars['B2']=	+1,	-1,	-1,	+1,	#y, Rx	yz	yz2, y3, x2y

        elif(sym=='C3v'):
            #C3v	E	2C3 (z)	3v	
            self.charval=1,2,3
            self.charOps='E','c3','sigv'
            self.chars['A1']=	+1,	+1,	+1, #	z	x2+y2, z2	z3, x(x2-3y2), z(x2+y2)
            self.chars['A2']=	+1,	+1,	-1, #	Rz	-	y(3x2-y2)
            self.chars['E']=	+2,	-1,	0,  #	(x, y) (Rx, Ry)	(x2-y2, xy) (xz, yz)	(xz2, yz2) [xyz, z(x2-y2)] [x(x2+y2), y(x2+y2)]

        elif(sym=='Oh'):
            #Oh	E	8C3	6C2	6C4	3C2 =(C4)2	i	6S4	8S6	3h	6d	
            self.charval=1,8,6,6,3,1,6,8,3,6
            self.charOps='E','c3','c2','c4','c42','i','s4','s6','sigh','sigd'	
            self.chars['A1g']=	+1,	+1,	+1,	+1,	+1,	+1,	+1,	+1,	+1,	+1, #	-	x2+y2+z2	-
            self.chars['A2g']=	+1,	+1,	-1,	-1,	+1,	+1,	-1,	+1,	+1,	-1, #	-	-	-
            self.chars['Eg']=	+2,	-1,	0,	0,	+2,	+2,	0,	-1,	+2,	0, #	-	(2z2-x2-y2, x2-y2)	-
            self.chars['T1g']=	+3,	0,	-1,	+1,	-1,	+3,	+1,	0,	-1,	-1, #	(Rx, Ry, Rz)	-	-
            self.chars['T2g']=	+3,	0,	+1,	-1,	-1,	+3,	-1,	0,	-1,	+1, #	-	(xz, yz, xy)	-
            self.chars['A1u']=	+1,	+1,	+1,	+1,	+1,	-1,	-1,	-1,	-1,	-1, #	-	-	-
            self.chars['A2u']=	+1,	+1,	-1,	-1,	+1,	-1,	+1,	-1,	-1,	+1, #	-	-	xyz
            self.chars['Eu']=	+2,	-1,	0,	0,	+2,	-2,	0,	+1,	-2,	0, #	-	-	-
            self.chars['T1u']=	+3,	0,	-1,	+1,	-1,	-3,	-1,	0,	+1,	+1, #	(x, y, z)	-	(x3, y3, z3) [x(z2+y2), y(z2+x2), z(x2+y2)]
            self.chars['T2u']=	+3,	0,	+1,	-1,	-1,	-3,	+1,	0,	+1,	-1, #	-	-	[x(z2-y2), y(z2-x2), z(x2-y2)]
        elif(sym=='Ih'):
            #Ih	E	12C5	12(C5)2	20C3	15C2	i	12S10	12(S10)3	20S6	15	
            cos=numpy.cos
            pi=numpy.pi
            self.charval=1,12,12,20,15,1,12,12,20,15
            self.charOps='E','c5','c52','c3','c2','i','s101','s102','s6','sigh'	
            self.chars['Ag']=	+1,	+1,	+1,	+1,	+1,	+1,	+1,	+1,	+1,	+1, #	-	x2+y2+z2	-
            self.chars['T1g']=	+3,	-2*cos(4.*pi/5),	-2*cos(2.*pi/5),	0,	-1,	+3,	-2*cos(2.*pi/5),	-2*cos(4.*pi/5),	0,	-1, #	(Rx, Ry, Rz)	-     
            self.chars['T2g']=	+3,	-2*cos(2.*pi/5),	-2*cos(4.*pi/5),	0,	-1,	+3,	-2*cos(4.*pi/5),	-2*cos(2.*pi/5),	0,	-1, #	-	-	-
            self.chars['Gg']=	+4,	-1,	-1,	+1,	0,	+4,	-1,	-1,	+1,	0,  #	-	-	-
            self.chars['Hg']=	+5,	0,	0,	-1,	+1,	+5,	0,	0,	-1,	+1, #	-	[2z2-x2-y2, x2-y2, xy, xz, yz]	-
            self.chars['Au']=	+1,	+1,	+1,	+1,	+1,	-1,	-1,	-1,	-1,	-1, #	-	-	-
            self.chars['T1u']=	+3,	-2*cos(4.*pi/5),	-2*cos(2.*pi/5),	0,	-1,	-3,	+2*cos(2.*pi/5),+2*cos(4.*pi/5),0,+1, # (x, y, z)-[x(z2+y2), y(z2+x2), z(x2+y2)]
            self.chars['T2u']=	+3,	-2*cos(2.*pi/5),	-2*cos(4.*pi/5),	0,	-1,	-3,	+2*cos(4.*pi/5),+2*cos(2.*pi/5),0,+1, # -	-[x3, y3, z3]
            self.chars['Gu']=	+4,	-1,	-1,	+1,	0,	-4,	+1,	+1,	-1,	0, #	-	-	[x(z2-y2), y(z2-x2), z(x2-y2), xyz]
            self.chars['Hu']=	+5,	0,	0,	-1,	+1,	-5,	0,	0,	+1,	-1, #	-	-	-
        elif(sym=='Td'):
            #Td	E	8C3	3C2	6S4	6d	
            self.charval=1,8,3,6,6
            self.charOps='E','c3','c2','s4','sigd'	
            self.chars['A1']=	+1,	+1,	+1,	+1,	+1, #	-	x2+y2+z2	xyz
            self.chars['A2']=	+1,	+1,	+1,	-1,	-1, #	-	-	-
            self.chars['E']=	+2,	-1,	+2,	0,	0, #	-	(2z2-x2-y2, x2-y2)	-
            self.chars['T1']=	+3,	0,	-1,	+1,	-1, #	(Rx, Ry, Rz)	-	[x(z2-y2), y(z2-x2), z(x2-y2)]
            self.chars['T2']=	+3,	0,	-1,	-1,	+1, #	(x, y, z)	(xy, xz, yz)	(x3, y3, z3) [x(z2+y2), y(z2+x2), z(x2+y2)]

        self.charOps=numpy.array(self.charOps)
        self.charval=numpy.array(self.charval)
        for key,vals in self.chars.items():
            self.chars[key]=numpy.array(vals)
        #print self.chars
        #print self.charOps
        #print self.charval

    #return xyz co-ordinates of preprogrammed geometries
    def GetCoords(self,tag):
        vv=[] #place 6 atoms in octahedron
        if(tag=='CH2' or tag=='CHH'): 
            vv.append((1,0,-1))
            vv.append((-1,0,-1))
            vv=numpy.array(vv)
        elif(tag=='CH3' or tag=='CHHH'): 
            #vv.append((1,0,-1/sqrt(2)))
            vv.append((-1,0,-1/sqrt(2)))
            vv.append((0,1,1/sqrt(2)))
            vv.append((0,-1,1/sqrt(2)))
            vv=numpy.array(vv)*(1.5**0.5)/1.5

        elif(tag=='octahedron' or tag=='CH6'):
            vv.append((0,0,1.))
            vv.append((0,0,-1.))
            vv.append((1.,0,0))
            vv.append((-1.,0,0))
            vv.append((0,1.,0))
            vv.append((0,-1.,0))
        elif(tag=='tetrahedron' or tag=='CH4'):
            vv.append((1,0,-1/sqrt(2)))
            vv.append((-1,0,-1/sqrt(2)))
            vv.append((0,1,1/sqrt(2)))
            vv.append((0,-1,1/sqrt(2)))
            vv=numpy.array(vv)*(1.5**0.5)/1.5
        elif(tag=='cube' or tag=='CH8'):
            vv.append((1,1,1))
            vv.append((1,1,-1))
            vv.append((-1,1,1))
            vv.append((-1,1,-1))
            vv.append((1,-1,1))
            vv.append((1,-1,-1))
            vv.append((-1,-1,1))
            vv.append((-1,-1,-1))
            vv=numpy.array(vv)/(3**0.5)
        elif(tag=='icosahedron' or tag=='CH12'):
            #cyclic permuatations of 0 +/-1 +/- phi where
            phi=(1 + sqrt(5))/2.
            vv.append((0,1,phi))
            vv.append((0,1,-phi))
            vv.append((0,-1,phi))
            vv.append((0,-1,-phi))
            vv.append((1,phi,0))
            vv.append((1,-phi,0))
            vv.append((-1,phi,0))
            vv.append((-1,-phi,0))
            vv.append((phi,0,1))
            vv.append((-phi,0,1))
            vv.append((phi,0,-1))
            vv.append((-phi,0,-1))
            vv=numpy.array(vv) /1.90211303259030
        elif(tag=='dodecahedron' or tag=='CH20'):
            phi=(1 + sqrt(5))/2.
            vv.append((phi,phi,phi))
            vv.append((phi,phi,-phi))
            vv.append((phi,-phi,phi))
            vv.append((phi,-phi,-phi))
            vv.append((-phi,phi,phi))
            vv.append((-phi,phi,-phi))
            vv.append((-phi,-phi,phi))
            vv.append((-phi,-phi,-phi))
            vv.append((0,phi**2,1))
            vv.append((0,phi**2,-1))
            vv.append((0,-phi**2,1))
            vv.append((0,-phi**2,-1))
            vv.append((phi**2,1,0))
            vv.append((phi**2,-1,0))
            vv.append((-phi**2,1,0))
            vv.append((-phi**2,-1,0))
            vv.append((1,0,phi**2.))
            vv.append((-1,0,phi**2.))
            vv.append((1,0,-phi**2.))
            vv.append((-1,0,-phi**2.))
            vv=numpy.array(vv)/2.8025170768881473
        else:
            print 'notag here for geometry'
            print tag
            sys.exit(100)

        outy=open('output/eps/v.out','w')
        for i,v in enumerate(vv):
            outy.write('0 0 0\n')
            outy.write('%f\t%f\t%f\t%i\n\n\n\n' % (v[0],v[1],v[2],i+2))
        outy.close()

        gnu=open('output/eps/gnu.gp','w')
        gnu.write('set term eps enh color solid\n')
        gnu.write('set size square\n')
        gnu.write('set ticslevel 0\n')
        gnu.write('set output \'output/eps/%s.eps\'\n' % tag)
        gnu.write('set xlabel \'x\'\n')
        gnu.write('set ylabel \'y\'\n')
        gnu.write('set zlabel \'z\'\n')
        gnu.write('splot ')
        for i,v in enumerate(vv):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'output/eps/v.out\' i %i u 1:2:3 noti w li,\'\' u 1:2:3:4 noti w labels' % i)
        gnu.close()
# DS:
        os.system('gnuplot.exe output/eps/gnu.gp')
        print self.GNUFILE + ' output/esp/gnu.gp'
        return numpy.array(vv)
    
    #work out dictionary (self.maps)
    #full of how symmetry operators map
    #one atom to another.
    def GetSymmetry(self,vv,tag):
        self.maps={} #will contain mappings of each element
        if(tag=='CHH' or tag=='CH2'):
            vCH2=self.GetCoords('CH2')
            self.TestSymIdentity(vCH2,'E')
            #self.TestSymInverse(vCH2,'i')

            self.TestSym((0,0,1),vCH2,180,'c2_1')
            self.TestSymRef((0,0,0),(1,0,0),vCH2,'sigv1_1')
            self.TestSymRef((0,0,0),(0,1,0),vCH2,'sigv2_1')

        elif(tag=='CHHH' or tag=='CH3'):
            vCH3=self.GetCoords('CH3')
            vTet=self.GetCoords('tetrahedron')
            self.TestSymIdentity(vCH3,'E')

            #self.WriteGeometry('test',vCH3)

            self.TestSym(vTet[0],vCH3,120.,'c3_1')
            self.TestSym(vTet[0],vCH3,-120.,'c3_2')

            edges=[] #stripped down from tetrahedron
            edges.append((2,3))
            edges.append((1,2))
            #edges.append((2,0))
            edges.append((1,3))
            #edges.append((1,0))
            #edges.append((3,0))

            for i,edge in enumerate(edges):
                self.TestSymRef((0,0,0),(vTet[edge[0]]-vTet[edge[1]])/2 ,vCH3,'sigv_'+str(i+1))             

        elif(tag=='octahedron' or tag=='cube' or tag=='CH6' or tag=='CH8'):
            #self.WriteGeometry('test',vv)
            #symmetry operations for Oh and Cube
            #transforms and planes should be self-consistent
            #have been defined relatve to input co-ordinates.
            #Oh 
            #8C3  #6C2	#6C4       
            #3C2 =(C4)2	
            #i	
            #6S4	
            #8S6	
            #3h	
            #6d
            vv=self.GetCoords(tag)
            #print vv
            vcube=self.GetCoords('cube')
            voct=self.GetCoords('octahedron')

            self.TestSymIdentity(vv,'E')
            self.TestSymInverse(vv,'i')

            for i,v in enumerate(voct): #one for each vertex of octahedron
                self.TestSym(v,vv,90.,'c4_'+str(i+1))
                self.TestSymImp(v,90.,(0,0,0),vv,'s4_'+str(i+1))
            for j,i in enumerate((1,3,5)): #one for each vertex (no duplicates)
                self.TestSym(voct[i],vv,180.,'c42_'+str(j+1))
                self.TestSymRef((0,0,0),voct[i],vv,'sigh_'+str(j+1))
            for i,v in enumerate(vcube): #one at centre of each face (dual vertex locations)
                self.TestSym(v,vv,120.,'c3_'+str(i+1))
                self.TestSymImp(v,60.,(0,0,0),vv,'s6_'+str(i+1))

            edges=[]
            edges.append((1,1,0))
            edges.append((1,-1,0))
            edges.append((1,0,1))
            edges.append((1,0,-1))
            edges.append((0,1,1))
            edges.append((0,1,-1))
            
            for i,edge in enumerate(edges): #one for each halfway point on edges
                self.TestSym(edge,vv,180.,'c2_'+str(i+1))
                self.TestSymRef((0,0,0),edge,vv,'sigd_'+str(i+1))

        elif(tag=='icosahedron' or tag=='dodecahedron' or tag=='CH20' or tag=='CH12'):

            vicos=self.GetCoords('icosahedron')
            vdodec=self.GetCoords('dodecahedron')

            #self.WriteGeometry('test',vv)
            #12 72o rotations   +108o rotoreflection
            #12 144o rotations  +36o rotoreflection
            #20 120o rotations #one for each face + 60o rotoreflection
            #15 180o rotations #one for each pair of edges + reflections
            self.TestSymIdentity(vv,'E')            
            self.TestSymInverse(vv,'i')
            for i,v in enumerate(vicos): #12 of these.
                self.TestSym(v,vv,72.,'c5_'+str(i+1))
                self.TestSym(v,vv,72.*2,'c52_'+str(i+1))
                self.TestSymImp(v,36.,(0,0,0),vv,'s101_'+str(i+1))
                self.TestSymImp(v,108.,(0,0,0),vv,'s102_'+str(i+1))

            for i,v in enumerate(vdodec): #a C3 on each vertex of dodecahedron
                self.TestSym(v,vv,120.,'c3_'+str(i+1))
                self.TestSymImp(v,60.,(0,0,0),vv,'s6_'+str(i+1))

            edges=[]
            edges.append((0,2)) #1
            edges.append((0,8)) #2
            edges.append((2,8)) #3
            edges.append((8,4)) #4
            edges.append((4,0)) #5
            edges.append((8,5)) #6
            edges.append((5,2)) #7
            edges.append((5,7)) #8
            edges.append((7,2)) #9
            edges.append((9,7)) #10
            edges.append((8,10)) #11
            edges.append((9,2)) #12
            edges.append((0,9)) #13
            edges.append((6,9)) #14
            edges.append((0,6)) #15

            #for i in range(12):
            #    cnt=0
            #    for edge in edges:
            #        if(i in edge):
            #            cnt+=1
            #    print i,cnt
            for i,edge in enumerate(edges): #mid point of edges
                self.TestSym((vicos[edge[0]]+vicos[edge[1]])/2.,vv,180.,'c2_'+str(i+1))
                self.TestSymRef((0,0,0),(vicos[edge[0]]+vicos[edge[1]])/2.,vv,'sigh_'+str(i+1))
        elif(tag=='tetrahedron' or tag=='CH4'):
            #Td	E	8C3	3C2	6S4	6d
            vTet=self.GetCoords('tetrahedron')

            self.WriteGeometry('test',vTet)
            self.TestSymIdentity(vTet,'E')            

            for i,v in enumerate(vTet):                
                self.TestSym(v,vTet,120.,'c3_'+str(i+1))
            for i,v in enumerate(vTet):                
                self.TestSym(v,vTet,-120.,'c3_'+str(i+1+4))

            edges=[]
            edges.append((2,3))
            edges.append((1,2))
            edges.append((2,0))
            edges.append((1,3))
            edges.append((1,0))
            edges.append((3,0))

            for i,edge in enumerate(edges):
                if(i<3):
                    self.TestSym((vTet[edge[0]]+vTet[edge[1]])/2.,vTet,180.,'c2_'+str(i+1))
                self.TestSymImp((vTet[edge[0]]+vTet[edge[1]])/2.,90.,(0,0,0),vTet,'s4_'+str(i+1))
                self.TestSymRef((0,0,0),(vTet[edge[0]]-vTet[edge[1]])/2 ,vTet,'sigd_'+str(i+1))             

    #if test is orthoganol to all reps in new, add it.
    def IsOrth(self,test,new):
        for ne in new:
            if(numpy.sum(ne*test)==0):
                if( (ne*-1==test).al()):
                    return False
                #print 'yes!'
                pass
            else:
                #print 'no'
                return False
        #print 'new guy! adding',test
        new.append(test)
        return True

    #normalise each representation
    def NormRep(self,arr):
        arr=numpy.array(arr)
        norm=numpy.sum((arr*arr),axis=1)**0.5
        for i in range(len(arr)):
            arr[i,:]=arr[i,:]/norm[i]
        return arr

    #take the raw projections together
    #to try and get degenerate representations
    def AddRep(self,degen,raw):
        new=[]
        if(degen==1):
            new.append(raw[0])
            if(len(new)==degen):
                return self.NormRep(new)
        
        for i in range(degen): #take linear combinations
            test=raw[0]+raw[i+1]
            if(self.IsOrth(test,new)):
                #print 'new!'
                pass
            if(len(new)==degen):
                return self.NormRep(new)

        #test=2*raw[3]-raw[1]-raw[2]
        #if(self.IsOrth(test,new)==0):
        #    print 'shit! problem with degeneracy.'
        #    sys.exit(100)

        return self.NormRep(new)

    def NormProj(self,new):
        #orthoganol test...
        for i,ni in enumerate(new):
            for j,nj in enumerate(new):
                if(i!=j):
                    tr=0
                    for key,vals in ni.items():
                        if(key in nj.keys()):
                            tr+=vals+nj[key]
                    #if(tr!=0):
                        #print 'shit. not ortho'
                        #print i,ni,j,nj
        #now normalise.
        for ne in new:
            norm=0
            for key,vals in ne.items():
                norm+=vals**2.
            norm=norm**0.5
            for key,vals in ne.items():
                ne[key]=vals/norm
        return new


    def FormatRepresentation(self,rep):
        st=''
        #print rep
        for i,key in enumerate(rep.keys()):
            frac=Fraction(numpy.fabs(rep[key][0])).limit_denominator()
            if(frac!=0):
                sq=0
                if(len(str(frac.numerator))>3):
                    frac2=Fraction(numpy.fabs((rep[key]**2.)[0])).limit_denominator()
                    if(len(str(frac2.numerator))<len(str(frac.numerator))):
                        frac=frac2
                        sq=1

                if(i!=0):
                    if(rep[key]<0):
                        st+='-'
                    else:
                        st+='+'
                if(i==0):
                    if(rep[key]<0):
                        st+='-'
                if(frac.denominator==1):
                    st+=str(frac.numerator)
                else:
                    st+='('+str(frac.numerator)+'/'+str(frac.denominator)+')'
                if(sq):
                    st+='^(0.5)'
                st+=key
        return st

    #try and find irreducible representation
    #testing all combinations of a+-b
    def DoRepCombs(self,targ,new,raw,sign):
        i=0
        j=0
        while(len(new)!=targ):
            if(i!=j):
                self.TestProj(i,j,raw,sign,new)
            #print len(new),targ
            if(len(new)==targ):
                break
            i+=1
            if(i==len(raw)):
                i=0
                j+=1
            if(j==len(raw)):
                break

        if(len(new)==targ):
            return True
        else:
            return False

    #trying all combinations of three operators
    def DoRepCombsBig(self,targ,new,raw,sign):
        print 'getting big'
        i=0
        j=0
        k=0
        while(len(new)!=targ):
            if(i!=j):
                self.TestProj(k,(i,j),raw,sign,new,big=True)
            #print len(new),targ
            if(len(new)==targ):
                break
            i+=1
            if(i==len(raw)):
                i=0
                j+=1

            if(j==len(raw)):
                j=0
                k+=1
            if(k==len(raw)):
                break

        if(len(new)==targ):
            return True
        else:
            return False

    #trying adding each term on its own
    def DoRepCombsDiag(self,targ,new,raw,sign):
        i=0
        while(len(new)!=targ):
            self.TestProj(i,i,raw,sign,new)
            if(len(new)==targ):
                break
            i+=1
            if(i==len(raw)):
                break
        if(len(new)==targ):
            return True
        else:
            return False

    #trying to see if we can get a new orthonormal irreducible representation
    def TestProj(self,a,b,raw,sign,new,big=False):
        if(big==False):
            success,nj=self.CombineProj(a,b,raw,sign)
        else:
            success,nj=self.CombineProjBig(a,b[0],b[1],raw,sign)
        if(success):
            for ni in new:
                tr=0
                for key,vals in ni.items():
                    if(key in nj.keys()):
                        #tr+=vals+nj[key]
                        tr+=vals*nj[key]
                #print 'fuckeyt',tr,ni,nj
                if(tr!=0):
                    #print 'shit. not ortho'
                    #print ni,nj    
                    success=False
                    break

        if(success):
            new.append(nj)

    #take linear combinations of the raw 
    #projected bases, to find orthoganol combinations
    #aiming to find the target number of representations.
    def AddRepPerm(self,degen,irred,raw):
        irr=irred[1]
        #print irred[0],irred[0][0],round(irred[0][0]),degen
        targ=degen*int(round(irred[0][0])) #number of desired representations.
        new=[]

        if(degen==1):
            #new.append(raw[0])
            #self.TestProj(0,0,raw,'p',new) #good for tetrahedron
            #if(len(new)==targ):
            #    return self.NormProj(new)
            if(self.DoRepCombsDiag(targ,new,raw,'p')):
                return self.NormProj(new)
            #if(self.DoRepCombs(targ,new,raw,'p')):
            #    return self.NormProj(new)

        elif(degen==3): #T
            #self.TestProj(0,1,raw,'p',new)#good for tetrahedron
            #self.TestProj(0,2,raw,'p',new)#good for tetrahedron
            #self.TestProj(0,3,raw,'p',new)#good for tetrahedron

            #if(len(new)==targ):
            #    return self.NormProj(new)

            if(self.DoRepCombsDiag(targ,new,raw,'p')):
                return self.NormProj(new)
            if(self.DoRepCombs(targ,new,raw,'p')):
                return self.NormProj(new)
            if(self.DoRepCombs(targ,new,raw,'m')):
                return self.NormProj(new)

        #if(irr[0]=='E'):
        elif(degen==2 or degen==4):
            if(len(raw)<3):
                print 'shit. edit the E in AddRepPerm'
                sys.exit(100)

            #self.TestProj(0,1,raw,'m',new) #good for tetrahedron
            #self.TestProj(1,2,raw,'m',new) #good for tetrahedron
            #if(len(new)==targ):
            #    return self.NormProj(new)

            if(self.DoRepCombsDiag(targ,new,raw,'p')):
                return self.NormProj(new)
            if(self.DoRepCombs(targ,new,raw,'m')):
                return self.NormProj(new)
            if(self.DoRepCombs(targ,new,raw,'p')):
                return self.NormProj(new)

            new=[] #try a different combination.
            #works sometimes. not sure why!
            if(self.DoRepCombs(targ,new,raw,'p')):
                return self.NormProj(new)
            if(self.DoRepCombs(targ,new,raw,'m')):
                return self.NormProj(new)
            if(self.DoRepCombsDiag(targ,new,raw,'p')):
                return self.NormProj(new)


        elif(degen==5):
            if(self.DoRepCombsDiag(targ,new,raw,'p')):
                return self.NormProj(new)
            if(self.DoRepCombs(targ,new,raw,'m')):
                return self.NormProj(new)
            if(self.DoRepCombs(targ,new,raw,'p')):
                return self.NormProj(new)

        else:

            print 'no entry for', irr
            sys.exit(100)


        if(self.DoRepCombsBig(targ,new,raw,'p')): #getting desperate...
            return self.NormProj(new)


        print 'not enough representations found:',irr,'degen:',degen,'targ:',targ
        #print ne
        for ne in new:
            print ne
        print 'found:',len(new),'target:',targ
        print 'yo',irred
        sys.exit(100)

    #get reducible representation
    def GetReducible(self):
        rep={}
        for key,vals in self.maps.items():
            #print key,vals
            #print key,numpy.trace(vals)
            koi=key.split('_')[0]
            if(koi not in rep.keys()):
                rep[koi]=numpy.trace(vals)
            else:
                rep[koi]+=numpy.trace(vals)
        return rep

    #take a function and reduce it
    def GetReducibleFunc(self,test):
        tast=numpy.array(test)
        rep={}
        for key,vals in self.maps.items():
            cha=0
            for tost in tast:
                tost=numpy.array(list(tost))
                tist=(tost=='a')*1.
                nob=numpy.dot(vals,tist)
                op=numpy.repeat('c',len(tost))
                op[nob==1]='a'  #set 'a's
                op[nob==0]='b'  #set 'b's
                if( (op==tost).all()):
                    cha+=1
            koi=key.split('_')[0]
            if(koi not in rep.keys()):
                rep[koi]=cha
            else:
                rep[koi]+=cha
        #for key,vals in rep.items():
        #    print key,vals
        return rep

    #reduce representation
    #apply reduction formula
    def GetIrreducible(self,rep):
        order=len(self.maps.keys())
        print 'number of symmetry operators:',order
        print rep
        irreds=[]
        for grp,tab in self.chars.items():
            tt=0
            #print grp
            for koi,val in rep.items():
                #print koi,val
                #print grp,tab
                mask=koi==self.charOps
                #print mask,tab[mask],self.charOps[mask]
                tt+=tab[mask]*val
                #print '   ',koi,val,self.charval[mask],tab[mask] ,tt,val,tab[mask]*val
            #print '   ',tt/order
            if(numpy.fabs(tt*1./order)>1E-6):
                irreds.append((tt*1./order,grp))
            #print grp,tab,tt,tt*1./order
            #print numpy.trace(self.maps['E'])
        #print irreds
        return irreds
    def PrintIrreducible(self,irreds,test):
        mult={}
        mult['E']=2
        mult['T']=3
        mult['A']=1
        mult['B']=1
        mult['H']=5
        mult['G']=4

        cnt=0
        st=''
        for i,irred in enumerate(irreds):
            if(i!=0):
                st+='+'
            num=int(round(irred[0][0]))
            if(num==1):
                st+=irred[1]
            else:
                st+=str(num)+irred[1]
            cnt+=num*mult[irred[1][0]]
        #print irreds
        print 'irreducible representation:',st 
        print 'total number of representations:',cnt,'required:',len(test)
        if(cnt!=len(test)):
            print 'shit. wrong number. aborting'
            sys.exit(1)
        return irreds

    def PrintRepresentation(self,reppy,test):
        cnt=0
        for irr,vals in reppy.items():
            print irr,':'
            for val in vals:
                print self.FormatRepresentation(val)
                cnt+=1
        if(len(test)!=cnt):
            print 'Number of representations',cnt,'does not match target:',len(test)
            sys.exit(100)
    def GetRepresentations(self,irreds):
            #now perform projection operator to get final representation
            for irred in irreds:
                ref=0                
                raw=[]
                for ref in range(len(vv)):
                    proj=numpy.zeros((len(vv)))
                    for key,vals in self.maps.items(): #for each symmetry operation...
                        #print vals[0]==1
                        #print numpy.nonzero(vals[0]==1)
                        koi=key.split('_')[0]  #get symmetry operator of map
                        mask=koi==self.charOps #get value from character table - col number
                        proj[numpy.nonzero(vals[ref]==1)[0]]+=self.chars[irred[1]][mask]
                    raw.append(proj)

                    
                degen=int(self.chars[irred[1]][self.charOps=='E'])
                print 'degeneracy:',degen

                new=self.AddRep(degen,raw)
                #are they orthogonal?
                for i in range(len(new)):
                    for j in range(len(new)):
                        if(i!=j):
                            if(numpy.sum(new[i]*new[j])!=0):
                                print 'shit. non degenerate representations'
                                print i,j,numpy.sum(new[i]*new[j])
                                sys.exit(100)
                #print norm
                print irred[1],new


    def GetRepresentationsFunc(self,irreds,test):
        reppy={}
        for irred in irreds:
            ref=0                
            raw=[]
            for tast in test:
                proj={}
                tost=numpy.array(list(tast))
                for key,vals in self.maps.items(): #for each symmetry operation...
                    #print tost
                    
                    tist=(tost=='a')*1.
                    
                    nob=numpy.dot(vals,tist)
                    #print nob
                    op=numpy.repeat('c',len(tost))
                    #op='c','c','c','c' #new placeholder
                    #op=numpy.array(op) #make it a numpy array
                    op[nob==1]='a'  #set 'a's
                    op[nob==0]='b'  #set 'b's
                    oppy=op.tostring()
                        
                    koi=key.split('_')[0]  #get symmetry operator of map
                    mask=koi==self.charOps #get value from character table - col number
                    #proj[numpy.nonzero(vals[ref]==1)[0]]+=self.chars[irred[1]][mask]
                    
                    if(oppy in proj.keys()):
                        proj[oppy]+=self.chars[irred[1]][mask]
                    else:
                        proj[oppy]=self.chars[irred[1]][mask]
                    
                raw.append(proj)

            degen=int(self.chars[irred[1]][self.charOps=='E'])
            #print 'degeneracy:',degen
            #print raw
            new=self.AddRepPerm(degen,irred,raw)
            #print 'adding:',irred[1],new
            reppy[irred[1]]=new
        return reppy

    def CombineProj(self,a,b,raw,sign='p'):
        new={}
        for key,vals in raw[a].items():
            if(key in new):
                new[key]+=vals
            else:
                new[key]=vals.copy()
        #print raw[b].items()
        for key,vals in raw[b].items():
            #print key,vals
            if(sign=='p'):
                v=+1
            else:
                v=-1
            try:
                new[key]+=vals*v
            except:
                new[key]=v*vals.copy()

        #now normalise.
        norm=0
        for key,vals in new.items():
            norm+=vals**2.
        if(norm==0): #normalisation is zero. aborting.
            return False,0
        norm=norm**0.5
        for key,vals in new.items():
            new[key]=vals/norm
        return True,new

    def CombineProjBig(self,a,b,c,raw,sign='p'):
        new={}
        for key,vals in raw[a].items():
            if(key in new):
                new[key]+=2*vals
            else:
                new[key]=2*vals.copy()
        #print raw[b].items()
        for key,vals in raw[b].items():
            try:
                new[key]-=vals
            except:
                new[key]=-vals.copy()
        for key,vals in raw[c].items():
            try:
                new[key]-=vals
            except:
                new[key]=-vals.copy()

        #now normalise.
        norm=0
        for key,vals in new.items():
            norm+=vals**2.
        if(norm==0):
            return False,0
        norm=norm**0.5
        for key,vals in new.items():
            new[key]=vals/norm
        return True,new

    #return coordinates of platonic solids. constant distance centre-point=1
    def GetGeometry(self,tag): 
        vv=self.GetCoords(tag)

        if(self.SYM): #work out symmetry transforms if doing this symbolically
            self.GetSymmetry(vv,tag) 
            if(tag=='cube' or tag=='octahedron'):
                self.GetCharTable('Oh')
            if(tag=='icosahedron' or tag=='dodecahedron'):
                self.GetCharTable('Ih')
            if(tag=='tetrahedron'):
                self.GetCharTable('Td')
            
            #rep=self.GetReducible(test)  #get reduible representation
            #irreds=self.GetIrreducible(rep) #reduce it...
            #reppy=GetRepresentations(irreds) #get representations.


            #test='aaab','aaba','abaa','baaa'
            #test='aabb','abab','abba','baab','baba','bbaa'

            #test='aaaaab','aaaaba','aaabaa','aabaaa','abaaaa','baaaaa'

            """
            for i in range(len(vv)+1): #for each distinct chemical shift
                print
                bas=numpy.concatenate((numpy.repeat('b',len(vv)-i),numpy.repeat('a',i)))
                test=list(more_itertools.distinct_permutations(bas)) #get the permutations
                print 'reducing:',bas,len(test)
                rep=self.GetReducibleFunc(test)
                irreds=self.GetIrreducible(rep) #reduce it...
                self.PrintIrreducible(irreds,test)
                reppy=self.GetRepresentationsFunc(irreds,test)
                self.PrintRepresentation(reppy,test)
            """
            #sys.exit(100)
            #we can now calculate irreducible representation etc.    
            #sys.exit(100)
        return vv

    def FormatStr(self,stry):
        return stry.replace('a','xxaxx').replace('b','xxbxx').replace('xxaxx','{/Symbol a}').replace('xxbxx','{/Symbol b}').replace('z','_z')

    def DoTcPlot(self,tcMin,tcMax,grid,autoRates,outfig,geofig=''):
        lin=numpy.arange(grid)
        tcList=tcMin*10**(numpy.log10(tcMax/tcMin)*lin/(grid-1))
        outy=open('output/eps/outy.out','w')
        for tc in tcList:
            self.tc=tc   #s
            outy.write('%e\t' % (tc))
            for start in autoRates:
                rate=self.CalcRate(start,start,verb='n')
                outy.write('%e\t' % (rate))

            try:
                flemming=self.CalcFlemming()
                for flem in flemming:
                    outy.write('%e\t' % (flem))
            
                tugarinov=self.CalcTugarinov()
                for flem in tugarinov:
                    outy.write('%e\t' % (flem))
            except:
                pass
            outy.write('\n')
        outy.close()


        outy=open('output/eps/gnu.gp','w')
        outy.write('set term post eps color enh solid\n')
        outy.write('set size square\n')
        outy.write('set key left\n')
        outy.write('set logscale \n')
        outy.write('set yrange[*:*]\n')
        outy.write('set xlabel \'{/Symbol t}_c(s)\'\n')
        outy.write('set ylabel \'Rate(s^{-1})\'\n')
        outy.write('set title \'%i %s\'\n' % (self.sfrq,outfig))
        outy.write('set output \'output/eps/%s\'\n' % (outfig))
        outy.write('plot ')
        for i,start in enumerate(autoRates):
            if(i!=0):
                outy.write(',')
            outy.write('\'output/eps/outy.out\' u 1:%i ti \'%s\'w li lw 9 lt %i' % (i+2,self.FormatStr(start),i+1))
        #flemlab='faaa','faab','fabb','fbbb'
        #for i,flem in enumerate(flemlab):
        #    outy.write(',\'outy.out\' u 1:%i ti \'%s\'w li lw 6lt %i' % (len(autoRates)+2+i,flem,i+1))
        #flemlab='taaa','taab','tabb','tbbb'
        #for i,flem in enumerate(flemlab):
        #    outy.write(',\'outy.out\' u 1:%i ti \'%s\'w li lw 1 lt %i' % (len(autoRates)+len(flemlab)+2+i,flem,i+1))
        outy.write('\n')
        outy.close()
        os.system('gnuplot.exe output/eps/gnu.gp')
        if(self.LATEX):
            outy=open(self.LATEXFILE,'a')
            outy.write('\\begin{figure}[h]\n')
            outy.write('\center\n')
            if(geofig!=''):
                outy.write('\\includegraphics[width=0.40\\textwidth]{output/eps/'+geofig+'}\n')
            outy.write('\\includegraphics[width=0.49\\textwidth]{output/eps/'+outfig+'}\n')
            outy.write('\\caption[]{Relaxation rates versus correlation time $\\tau_c$ for spin system %s. macro=%s. %s }\n' % (self.baseTag,self.macro,self.plotTxt))
            outy.write('\end{figure} \n')

    def CalcNOESYCurve(self, rate_matrix, t_len=0.5):
        steps=1000
        t_arr = numpy.linspace(0, t_len, steps)
        
        num = len(rate_matrix[:,0])
        calc_matrix = numpy.zeros((num, num, steps))
        for i in numpy.arange(steps):
            calc_matrix[:,:,i] = expm(-rate_matrix[:,:] * t_arr[i])
        return calc_matrix
   


    def DoRatePlot(self, dist_matrix, rate_matrix, outfig, t_len=0.5):
        steps=1000
        t_arr = numpy.linspace(0, t_len, steps)
        num = rate_matrix.shape[0]
        calc_matrix = numpy.zeros((num, num,steps))
        for i in numpy.arange(steps):
            calc_matrix[:,:,i] = expm(-rate_matrix[:,:] * t_arr[i])

        dist_matrix = dist_matrix.flatten()
        out_matrix = calc_matrix[:,:,400].flatten()
        
        mask = numpy.logical_and(dist_matrix < 10, dist_matrix > 0.1)
        
        dist_matrix = dist_matrix[mask]
        out_matrix = out_matrix[mask]
        
        outy=open('eps/outy.out','w')
        for i in range(len(dist_matrix)):
            outy.write('%f\t' % (dist_matrix[i]))
            outy.write('%f\t' % (out_matrix[i]))
            outy.write('\n')
        outy.write('\n')
        outy.close()

        outy=open('output/eps/gnu.gp','w')
        outy.write('set term post eps color enh solid\n')
        outy.write('set size square\n')
        outy.write('set key left\n')
        outy.write('set yrange[*:*]\n')
        outy.write('set xlabel \'Distance, Angstroms\'\n')
        outy.write('set ylabel \'Magnetization at 200ms\'\n')
        outy.write('set title \'%i %s\'\n' % (self.sfrq,outfig))
        outy.write('set output \'output/eps/%s\'\n' % (outfig))
        outy.write('plot ')
        outy.write('\'output/eps/outy.out\' with points pt 7')
        outy.write('\n')
        outy.close()
        os.system('gnuplot.exe output/eps/gnu.gp')
        
        max_pos = out_matrix[numpy.argmax(out_matrix[:])]
        
        if(self.LATEX):
            outy=open(self.LATEXFILE,'a')
            outy.write('\\begin{figure}[h]\n')
            outy.write('\center\n')
            outy.write('\\includegraphics[width=0.49\\textwidth]{output/eps/'+outfig+'}\n')
            outy.write('\\caption[]{Rate constant against distance for protein yada. $\\tau_c$ for spin system = %s. macro=%s. Max = %.2f s %s }\n' % (self.baseTag,self.macro,max_pos,self.plotTxt))
            outy.write('\end{figure} \n')
            
        

            
    def DoNOESYPlot(self, start, end, t_len, outfig, geofig=''):
        #t_len sets the timeframe we're looking at.
        steps=1000
        t_arr = numpy.linspace(0, t_len, steps)

        r11 = self.CalcRate(start, start, verb='n')
        r22 = self.CalcRate(end, end, verb='n')
        r12 = self.CalcRate(start, end, verb='n')
        r21 = self.CalcRate(end, start, verb='n')

        r_matrix = numpy.array([[r11, r21], [r12, r22]])    
        calc_matrix = numpy.zeros((2,2,steps))
        for i in numpy.arange(steps):
            calc_matrix[:,:,i] = expm(-r_matrix[:,:] * t_arr[i])

        outy=open('output/eps/outy.out','w')
        for i in range(len(t_arr)):
            outy.write('%e\t' % (t_arr[i]))
            outy.write('%e\t' % (calc_matrix[0,0,i]))
            outy.write('%e\t' % (calc_matrix[0,1,i]))
            outy.write('%e\t' % (calc_matrix[1,0,i]))
            outy.write('%e\t' % (calc_matrix[1,1,i]))
            outy.write('\n')
        outy.write('\n')
        outy.close()

        outy=open('output/eps/gnu.gp','w')
        outy.write('set term post eps color enh solid\n')
        outy.write('set size square\n')
        outy.write('set key left\n')
        outy.write('set yrange[*:*]\n')
        outy.write('set xlabel \'t/s\'\n')
        outy.write('set ylabel \'Magnetization\'\n')
        outy.write('set title \'%i %s\'\n' % (self.sfrq,outfig))
        outy.write('set output \'output/eps/%s\'\n' % (outfig))
        outy.write('plot ')
        outy.write('\'output/eps/outy.out\' u 1:%i ti \'%s\'w li lw 9 lt %i' % (2,'11',1))
        outy.write(',')
        outy.write('\'output/eps/outy.out\' u 1:%i ti \'%s\'w li lw 9 lt %i' % (3,'12',2))
        outy.write(',')
        outy.write('\'output/eps/outy.out\' u 1:%i ti \'%s\'w li lw 9 lt %i' % (4,'21',3))
        outy.write(',')
        outy.write('\'output/eps/outy.out\' u 1:%i ti \'%s\'w li lw 9 lt %i' % (5,'22',4))
        outy.write('\n')
        outy.close()
        os.system('gnuplot.exe output/eps/gnu.gp')
        
        max_pos = t_arr[numpy.argmax(calc_matrix[0,0,:])]
        
        if(self.LATEX):
            outy=open(self.LATEXFILE,'a')
            outy.write('\\begin{figure}[h]\n')
            outy.write('\center\n')
            if(geofig!=''):
                outy.write('\\includegraphics[width=0.40\\textwidth]{output/eps/'+geofig+'}\n')
            outy.write('\\includegraphics[width=0.49\\textwidth]{output/eps/'+outfig+'}\n')
            outy.write('\\caption[]{NOESY - Magnetization transfer from methyl to external spin. $\\tau_c$ for spin system = %s. macro=%s. Max = %.2f s %s }\n' % (self.baseTag,self.macro,max_pos,self.plotTxt))
            outy.write('\end{figure} \n')

    def DoBasisPlot(self,testbasis,outfig):

        self.EvalRate(testbasis,testbasis,verb='n')  #get symbolic rates
        vals=numpy.zeros((len(testbasis),len(testbasis)))
    
        outy=open('output/eps/outy.out','w')
        for i,start in enumerate(testbasis):
            for j,finish in enumerate(testbasis):
                vals[i,j]=self.CalcRate(start,finish,verb='n')
                #outy.write('%s\t%s\t%f\t%f\t%f\t%f\t%e\n' % (self.FormatStr(start),self.FormatStr(finish),i,j,i-0.5,j-0.5,vals[i,j]))
                #outy.write('%s\t%s\t%f\t%f\t%f\t%f\t%e\n' % (self.FormatStr(start),self.FormatStr(finish),i,j,i-0.5,j+0.5,vals[i,j]))
                outy.write('%s\t%s\t%f\t%f\t%f\t%f\t%e\n' % (start,finish,i,j,i-0.5,j-0.5,vals[i,j]))
                outy.write('%s\t%s\t%f\t%f\t%f\t%f\t%e\n' % (start,finish,i,j,i-0.5,j+0.5,vals[i,j]))
            outy.write('\n')
            for j,finish in enumerate(testbasis):
                vals[i,j]=self.CalcRate(start,finish,verb='n')
                #outy.write('%s\t%s\t%f\t%f\t%f\t%f\t%e\n' % (self.FormatStr(start),self.FormatStr(finish),i,j,i+0.5,j-0.5,vals[i,j]))
                #outy.write('%s\t%s\t%f\t%f\t%f\t%f\t%e\n' % (self.FormatStr(start),self.FormatStr(finish),i,j,i+0.5,j+0.5,vals[i,j]))
                outy.write('%s\t%s\t%f\t%f\t%f\t%f\t%e\n' % (start,finish,i,j,i+0.5,j-0.5,vals[i,j]))
                outy.write('%s\t%s\t%f\t%f\t%f\t%f\t%e\n' % (start,finish,i,j,i+0.5,j+0.5,vals[i,j]))
            outy.write('\n')
        outy.close()
        outy=open('output/eps/gnu.gp','w')
        outy.write('set term post eps color enh solid\n')
        outy.write('set size square\n')
        outy.write('set key left\n')
        outy.write('unset key\n')
        outy.write('set cblabel \'rate(s-1)\'\n')
        outy.write('set cbrange[-%f:%f]\n' % (numpy.max(vals),numpy.max(vals)))
        outy.write('set palette defined(-1"blue",-0.5"green",0"white",0.5"orange",1"red")\n')
        outy.write('set pm3d map\n')
        outy.write('set title \'tc %.1e tm %.1e %s\'\n' % (self.tc,self.tm,outfig))
        outy.write('set output \'output/eps/%s\'\n' % (outfig))
        #outy.write('plot ')
        #for i,start in enumerate(testbasis):
        #    if(i!=0):
        #        outy.write(',')
        outy.write('splot \'output/eps/outy.out\' u 5:6:7:xticlabels(1):yticlabels(2) lc palette')
        outy.write('\n')
        outy.close()
#         DS:
        os.system('gnuplot.exe output/eps/gnu.gp')
#         os.system('gnuplot eps/gnu.gp')

        if(self.LATEX):
            outy=open(self.LATEXFILE,'a')
            outy.write('\\begin{figure}[h]\n')
            outy.write('\center\n')
            #if(geofig!=''):
            #    outy.write('\\includegraphics[width=0.40\\textwidth]{'+geofig+'}\n')
            extra='. $\\tau_c$= %.2f ns $\\tau_m=$ %.2f ps' % (1E9*self.tc,1E12*self.tm)

            outy.write('\\includegraphics[width=0.49\\textwidth]{output/eps/'+outfig+'}\n')
            outy.write('\\caption[]{Relaxation rates versus correlation time $\\tau_c$ for spin system %s. macro=%s. %s }\n' % (self.baseTag,self.macro,self.plotTxt+extra))
            outy.write('\end{figure} \n')



    ###############################################################################
    #calculate neccessary parameters for calculating numerical relaxation rates
    #numbers saved in self.pars. These are converted to values in self.const
    #this is semi-manual. for a new arrangement, will need a new entry here.
    def SetPars(self,tag):
        self.plotTxt=''
        self.plotTxt+='$\nu_0$ %.2f MHz, ' % (self.sfrq)
        if(tag=='CH3'):
            i = numpy.arange(3)
            v_m1 = self.pars['v_m1']
            sites1 = self.pars['sites1']
            c1_site = numpy.repeat(self.pars['c1_site'][:,numpy.newaxis,:],3,axis=1)
            cA = self.pars['cA']
            cX = self.pars['cX']
            self.CalcConstsDip('d1', 'C', 'H', sites1, c1_site)
            self.CalcConstsDip('e1', 'H', 'H', sites1, sites1[:,(i+1)%3,:])
            self.CalcConstsCSA('cX1', 'H', cX, sites1)
            self.CalcConstsCSA('cA1', 'C', cA, v_m1, typ='static')
            for tag_k in self.Cint:
                for tag_l in self.Cint:
                    self.CalcConsts(tag_k, tag_l)
        elif(tag=='CH3CH3'):
            i = numpy.arange(3)
            v_m1 = self.pars['v_m1']
            v_m2 = self.pars['v_m2']
            v_mm = self.pars['v_mm']
            sites1 = self.pars['sites1']
            sites2 = self.pars['sites2']
            c1_site = numpy.repeat(self.pars['c1_site'][:,numpy.newaxis,:],3,axis=1)
            c2_site = numpy.repeat(self.pars['c2_site'][:,numpy.newaxis,:],3,axis=1)
            cA = self.pars['cA']
            cX = self.pars['cX']
            self.CalcConstsDip('d1', 'C', 'H', sites1, c1_site)
            self.CalcConstsDip('d2', 'C', 'H', sites2, c2_site, transpose=True)
            self.CalcConstsDip('e1', 'H', 'H', sites1, sites1[:,(i+1)%3,:])
            self.CalcConstsDip('e2', 'H', 'H', sites2, sites2[:,(i+1)%3,:], transpose=True)
            self.CalcConstsCSA('cX1', 'H', cX, sites1)
            csa_2 = sites2 - c2_site
            self.CalcConstsCSA('cX2', 'H', cX, csa_2, transpose=True)
            self.CalcConstsCSA('cA1', 'C', cA, v_m1, typ='static')
            self.CalcConstsCSA('cA2', 'C', cA, v_m2, typ='static', transpose=True)
            self.CalcConstsMet('m', 'H', 'H', sites1, sites2)
            
            for tag_k in self.Cint:
                for tag_l in self.Cint:
                    self.CalcConsts(tag_k, tag_l)
#             for key, value in self.consts.iteritems():
#                 print 'TAG = ', key, '\n', self.consts[key]
#                 print '\n'
#             print 'm = ',self.consts['m']
#             print 'self.consts[C(m^2)] = \n',self.consts['C(m^2)']
#             print 'self.consts[S_tm_2(m^2)^2] = \n',self.consts['S_tm_2(m^2)^2']
#             print 'self.consts[S_tm(m^2)^2] = \n',self.consts['S_tm(m^2)^2']
#             print 'self.consts[S_inf(m^2)^2] = \n',self.consts['S_inf(m^2)^2'] 
#             print 'TOTAl = \n',self.consts['S_tm_2(m^2)^2'] + self.consts['S_tm(m^2)^2'] + self.consts['S_inf(m^2)^2'] 
            
#             print 'e1 = ',self.consts['e1']
#             print 'self.consts[C(e1^2)] = \n',self.consts['C(e1^2)']
#             print 'self.consts[1-S(e1^2)^2] = \n',self.consts['(1-S(e1^2)^2)']
#             print 'self.consts[S(e1^2)^2] = \n',self.consts['S(e1^2)^2']
            
#             print 'd1 = ',self.consts['d1']
#             print 'self.consts[C(d1^2)] = \n',self.consts['C(d1^2)']
#             print 'self.consts[1-S(d1^2)^2] = \n',self.consts['(1-S(d1^2)^2)']
#             print 'self.consts[S(d1^2)^2] = \n',self.consts['S(d1^2)^2']
        else:
            print 'Unrecognised tag:',tag
            sys.exit(100)


    def InitLatex(self):
        if(os.path.exists(self.LATEXFILE)):
            os.system('!rm latex.tex')
        print self.LATEXFILE
        outy=open(self.LATEXFILE,'w')
        outy.write('\\documentclass[showkeys,aps,prb,prepreint,amssymb, amsmath,nobibnotes]{revtex4}\n')
        outy.write('\\usepackage{bm,setspace}\n')
        outy.write('\\usepackage{graphicx,graphics,booktabs}\n')
        #outy.write('\\usepackage[showframe=true]{geometry}\n')
        outy.write('\\usepackage{accents}\n')
        outy.write('\\newlength{\dhatheight}\n')
        #outy.write('\\newcommand{\doublehat}[1]{%\settoheight{\dhatheight}{\ensuremath{\hat{#1}}}%\addtolength{\dhatheight}{-0.35ex}%\hat{\vphantom{\rule{1pt}{\dhatheight}}%\smash{\hat{#1}}}}\n')
        outy.write('\\begin{document}\n')
        outy.write('\\section{%s}\n' % self.baseTag);
        outy.close()



    def MacroTable(self,typ,n):


        resMac={} #first take macromolecular limit
        for Pre,Posts in self.resSet.items():
            for Post,ds in Posts.items():
                #print Pre,Post
                if(Pre not in resMac.keys()):
                    resMac[Pre]={}
                if(Post not in resMac[Pre].keys()):
                    resMac[Pre][Post]={}
                for d,js in ds.items():
                    for j,val in js.items():
                        #print 'before:',j,val
                        tick=0
                        if(j=='tm' or j=='J(0,0)' or j=='J(0,wC)'):
                            pass
                        elif(j.split('(')[1].split(',')[0]!='0'):
                            j='tm'
                        else:
                            tick=1

                            
                        if(tick==0):
                            if(d=='f^2C(f^2)(1-S(f^2)^2)'):
                                d='f^2(1-S^2)'
                            if(d=='f^2C(f^2)S(f^2)^2'):
                                d='f^2S^2'

                            #print 'now:',j,val
                            if(Pre not in resMac.keys()):
                                resMac[Pre]={}
                            if(Post not in resMac[Pre].keys()):
                                resMac[Pre][Post]={}
                            if(d not in resMac[Pre][Post].keys()):
                                resMac[Pre][Post][d]={}
                            if(j not in resMac[Pre][Post][d].keys()):
                                resMac[Pre][Post][d][j]=val
                            else:
                                resMac[Pre][Post][d][j]+=val
                                
        #sys.exit(100)

        if(typ=='rot'):
            order=[]
            for PP in 'P0','P1','P2':
                if(PP=='P0'):
                    order.append('d^2%s%s' % (PP,PP))
                    order.append('dcA%s' % PP)
                    order.append('cA^2' )
                    order.append('%s^2' % 'e')
                else:
                    order.append('d^2%s%s' % (PP,PP))
                    #order.append('dcA%s' % PP)

                order.append('dcX%s%s' % (PP,PP))
                if(PP=='P2'):
                    order.append('%scX%s' % ('e',PP))
                order.append('cX^2%s%s' % (PP,PP))

            order.append('f^2S^2')
            order.append('f^2(1-S^2)')
            jord={}
            for PP in 'P0','P1','P2':
                if(PP=='P0'):
                    jord['d^2P0P0']='J(0,0)','J(0,wC)'
                    jord['dcAP0']='J(0,0)','J(0,wC)'
                    jord['cA^2']='J(0,0)','J(0,wC)'            
                    jord['e^2']='J(0,0)','tm'          
                    jord['ecXP0']='J(0,0)',          
                    jord['cX^2P0P0']='J(0,0)',
                    jord['dcXP0P0']='J(0,0)',
                else:
                    for ord in order:
                        if(ord not in jord.keys()):
                            jord[ord]='tm',


            jord['f^2S^2']='J(0,0)',#'tm'       
            jord['f^2(1-S)^2']='tm',

        if(typ=='static'):
            Ps=0
            for Pre,Posts in resMac.items():
                for Post,ds in Posts.items():
                    for d,js in ds.items():
                        if(len(d.split('P0'))>1):
                            Ps=1
                            break

            order=[]

            if(Ps==1): #if any of the rates have P0/P1/P2 still in them...

                for PP in 'P0','P1','P2':
                    if(PP=='P0'):
                        order.append('d^2%s%s' % (PP,PP))
                        order.append('dcA%s' % PP)
                        order.append('cA^2' )
                        order.append('%s^2' % 'e')
                    else:
                        order.append('d^2%s%s' % (PP,PP))
                        #order.append('dcA%s' % PP)

                    order.append('dcX%s%s' % (PP,PP))
                    #if(PP=='P2'):
                    #    order.append('%scX%s' % ('e',PP))
                    order.append('cX^2%s%s' % (PP,PP))

                order.append('f^2S^2')
                #order.append('f^2(1-S^2)')
                jord={}
                for PP in 'P0','P1','P2':
                    jord['d^2'+PP+PP]='J(0,0)','J(0,wC)'
                    jord['dcX'+PP+PP]='J(0,0)',
                    jord['cX^2'+PP+PP]='J(0,0)',
                    
                    #if(PP=='P2'):
                    #    jord['ecX'+PP]='J(0,0)',          
                    
                    if(PP=='P0'):

                        jord['dcAP0']='J(0,0)','J(0,wC)'
                        jord['cA^2']='J(0,0)','J(0,wC)'            
                        jord['e^2']='J(0,0)',         
                        
                        #else:
                        #    for ord in order:
                        #        if(ord not in jord.keys()):
                        #            jord[ord]='tm',
                jord['f^2S^2']='J(0,0)', #'J(0,wC)' #'tm'       
                #jord['f^2(1-S)^2']='tm',

            else:
                order.append('d^2')
                for e in self.edict.keys():
                    order.append('%s^2' % e)
                order.append('dcX')
                for e in self.edict.keys():
                    order.append('%scX' % e)
                order.append('cX^2')
                order.append('f^2S^2')
                #order.append('f^2(1-S^2)')

                jord={}
                jord['d^2']='J(0,0)','J(0,wC)'
                for ord in order:
                    if(ord!='d^2'):
                        jord[ord]='J(0,0)',
                jord['f^2S^2']='J(0,0)', #'J(0,wC)' #'tm'       

        done={}

        outy=open('output/eps/latex.tex','a')
        outy.write('\\begin{table}\n')
        outy.write('\\begin{small}\n')
        #outy.write('\\begin{adjustwidth}{-1cm}{}\n')
        outy.write('\\hspace*{-2cm}\n')
        outy.write('\\begin{tabular}{c|')
        for ord in order:
            for i,jo in enumerate(jord[ord]):
                outy.write('c')
        outy.write('}\n')
        #outy.write('&')
        for ord in order:
            for i,jo in enumerate(jord[ord]):
                outy.write('&$%s$' % (ord.replace('X','_X').replace('A','_A').replace('P0P0','P_0^2').replace('P1P1','P_1^2').replace('P2P2','P_2^2').replace('P0','P_0').replace('P1','P_1').replace('P2','P_2')     ))
        outy.write('\\\\\n')
        #outy.write('&')
        for ord in order:
            for i,jo in enumerate(jord[ord]):
                jo=jo.replace('J(0,0)','\\tau_c').replace('w','\\omega').replace('C','_C').replace('tm','\\tau_m')
                                            
                outy.write('&$%s$' % (jo))
        outy.write('\\\\\n')
        outy.write('\hline\\\\\n')

        for i in range(n+1):
            for Pre,Posts in resMac.items():
                for Post,ds in Posts.items():
                    if(Pre.count('a')==i and Pre.count('H')==0 and Pre.count('h')==0 and Pre==Post):
                        outy.write('%i%i' % (i,n-i))
                        for ord in order:
                            for jo in jord[ord]:
                                outy.write('&')
                                for d,js in ds.items():
                                    if(d==ord):
                                        for j,val, in js.items():
                                            if(j==jo):
                                                outy.write('$%s$' % (self.FormatLatexFrac(val,plus=False)))
                                                del resMac[Pre][Post][d][j]
                                                if(Pre not in done.keys()):
                                                    done[Pre]=[]
                                                done[Pre].append(Post)
                                                #print i,Pre,Post,d,js
                                                #vol=Fraction(val).limit_denominator()
                                                #outy.write('$%s$' % (self.FormatLatexFrac(val)))
                        outy.write('\\\\\n')

        for i in range(n+1):
            for j in range(n+1-i):
                jj=i+j+1
                for Pre,Posts in resMac.items():
                    for Post,ds in Posts.items():
                        if(Pre.count('a')==i and Pre.count('H')==0 and Pre.count('h')==0 and Post.count('H')==0 and Post.count('h')==0 and Post.count('a')==jj):
                            print 'yy',Pre,Post,i,n-i,jj,n-jj
                            if(Pre in done.keys() and Post in done[Pre]):
                                continue
                            else:
                                outy.write('%i%i %i%i' % (i,n-i,jj,n-jj))
                                for ord in order:
                                    for jo in jord[ord]:
                                        outy.write('&')
                                        for d,js in ds.items():
                                            if(d==ord):
                                                for j,val, in js.items():
                                                    if(j==jo):
                                                        #print i,Pre,Post,d,js
                                                        outy.write('$%s$' % (self.FormatLatexFrac(val,plus=False)))
                                                        del resMac[Pre][Post][d][j]
                                                        if(Pre not in done.keys()):
                                                            done[Pre]=[]
                                                        done[Pre].append(Post)
                                outy.write('\\\\\n')


        for test in 'Ha','Hb','CzHz','CzHzHz':
            for Pre,Posts in resMac.items():
                if(Pre==test):
                    Prey=Pre.replace('a','_\\alpha').replace('b','_\\beta').replace('z','_z')
                    outy.write('$%s$' % (Prey))
                    ds=Posts[Pre]
                    #for Post,ds in Posts.items():
                    for ord in order:
                        for jo in jord[ord]:
                            outy.write('&')
                            for d,js in ds.items():
                                if(d==ord):
                                    for j,val, in js.items():
                                        if(j==jo):
                                            outy.write('$%s$' % (self.FormatLatexFrac(val,plus=False)))

                                            del resMac[Pre][Pre][d][j]

                                            if(Pre not in done.keys()):
                                                done[Pre]=[]
                                            done[Pre].append(Pre)
                    outy.write('\\\\\n')




        outy.write('\\end{tabular}\n')
        #outy.write('\\end{adjustwidth}\n')

        outy.write('\\end{small}\n')
        outy.write('\\caption{Macromolecular limit relaxation rates. Left over: ')
        for Pre in done.keys():
            for Post in done[Pre]:
                ds=resMac[Pre][Post]
                for d,js in ds.items():
                    for j,val in js.items():
                        if(val!=0):
                            outy.write(' $'+Pre+Post+d+j+str('%.2f' %val)+'$ ')
        outy.write('}\n\n')
        outy.write('\\end{table}\n')
        outy.close()

            


        """
            outy=open('eps/latex.tex','a')
            outy.write('\\begin{table}\n')
            outy.write('\\begin{small}\n')
            outy.write('\\begin{tabular}{cc|')
            for ord in order:
                for i,jo in enumerate(jord[ord]):
                    outy.write('c')
            outy.write('}\n')
            outy.write('&')
            for ord in order:
                for i,jo in enumerate(jord[ord]):
                    outy.write('&$%s$' % (ord.replace('X','_X').replace('A','_A').replace('P0P0','P_0^2').replace('P1P1','P_1^2').replace('P2P2','P_2^2').replace('P0','P_0').replace('P1','P_1').replace('P2','P_2')     ))
                    #outy.write('&$%s$' % (ord.replace('X','_X').replace('P0P0','P_0^2').replace('P0','P_0')     ))
            outy.write('\\\\\n')
            outy.write('&')
            for ord in order:
                for i,jo in enumerate(jord[ord]):
                    jo=jo.replace('J(0,0)','\\tau_c').replace('w','\\omega').replace('C','_C')
                    outy.write('&$%s$' % (jo))
            outy.write('\\\\\n')
            for i in range(n+1):
                for Pre,Posts in resMac.items():
                    for Post,ds in Posts.items():
                        #print Pre,Pre.count('a'),i,Pre.count('H')
                        if(Pre.count('a')==i and Pre.count('H')==0 and Pre.count('h')==0 and Pre==Post):
                            #print 'hello!',i,n-1,Pre,Post
                            outy.write('%i&%i' % (i,n-i))
                            for ord in order:
                                for jo in jord[ord]:
                                    outy.write('&')
                                    #if(i==0 or i==n):
                                    #    outy.write('0')
                                    
                                    for d,js in ds.items():
                                        if(d==ord):
                                            for j,val, in js.items():
                                                if(j==jo):
                                                    outy.write('$%s$' % (self.FormatLatexFrac(val,plus=False)))
                                                    #vol=Fraction(val).limit_denominator()
                                                    #outy.write('$\\frac{%i}{%i}$' % (vol.numerator,vol.denominator))
                            outy.write('\\\\\n')

            for i in range(n+1):
                for j in range(n+1-i):
                    jj=i+j+1
                    for Pre,Posts in resMac.items():
                        for Post,ds in Posts.items():
                            print Pre,Post,Pre.count('a'),Post.count('a')
                            if(Pre.count('a')==i and Pre.count('H')==0 and Pre.count('h')==0 and Post.count('a')==jj):
                                outy.write('%i%i&%i%i' % (i,n-i,jj,n-jj))
                                for ord in order:
                                    for jo in jord[ord]:
                                        outy.write('&')
                                        for d,js in ds.items():
                                            if(d==ord):
                                                for j,val, in js.items():
                                                    if(j==jo):
                                                        print i,Pre,Post,d,js
                                                        #vol=Fraction(val).limit_denominator()
                                                        outy.write('$%s$' % (self.FormatLatexFrac(val,plus=False)))
                                                        #outy.write('$\\frac{%i}{%i}$' % (vol.numerator,vol.denominator))
                                outy.write('\\\\\n')


            for test in 'Ha','Hb','CzHz','CzHzHz':
                for Pre,Posts in resMac.items():
                    #print Pre,Pre.count('a'),i,Pre.count('H')
                    for Post,ds in Posts.items():
                        if(Pre==test and Post==Pre):
                            Pre=Pre.replace('a','_\\alpha').replace('b','_\\beta').replace('z','_z')
                            outy.write('$%s$&' % (Pre))
                            for ord in order:
                                for jo in jord[ord]:
                                    outy.write('&')
                                    for d,js in ds.items():
                                        if(d==ord):
                                            for j,val, in js.items():
                                                if(j==jo):
                                                    outy.write('$%s$' % (self.FormatLatexFrac(val,plus=False)))
                                                    #vol=Fraction(val).limit_denominator()
                                                    #outy.write('$\\frac{%i}{%i}$' % (vol.numerator,vol.denominator))
                                                    
                                                    #print Pre,Pre,d,js
                            outy.write('\\\\\n')




            outy.write('\\end{tabular}\n')
            outy.write('\\end{small}\n')
            outy.write('\\caption{Macromolecular limit relaxation rates}\n')
            outy.write('\\end{table}\n')
            outy.close()

        """
    def CloseLatex(self):
        if(self.LATEX):
            outy=open(self.LATEXFILE,'a')
            outy.write('\\end{document}\n')
            outy.close()
#             DS:
            os.system('pdflatex ' + self.LATEXFILE)
            os.system('!mv latex.pdf output/pdf/' + self.baseTag+'.pdf')
            os.system('!rm latex.aux latex.log')

    ############################################################

    def ReadCommDictDouble(self):
        self.commDictDouble={}
        inny=open('redfield/refDoubleComm.txt')
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                
                ops=test[0].split(']]')[0]
                c=ops.split('[')[1].split(',')[0]#outer operator
                b=ops.split('[')[2].split(',')[0]#middle operator 
                a=ops.split(',')[2]              #inner operator

                num=test[0].split('=')[1].split('(')[0]  #
                res=test[0].split('(')[1].split(')')[0]  #
                
                if(num[-1]=='i'):
                    if(len(num[:-1].split('/'))>1):
                        val=eval( 'float('+num[:-1].split('/')[0]+')'+ '/' +'float('+ num[:-1].split('/')[1] +')')*complex(0,1)
                    else:
                        val=eval('float('+num[:-1]+')')*complex(0,1)
                else:
                    if(len(num[:-1].split('/'))>1):
                        val=eval( 'float('+num[:].split('/')[0]+')'+ '/' +'float('+ num[:].split('/')[1] +')')
                    else:
                        val=eval('float('+num[:]+')')

                if(a not in self.commDictDouble.keys()):
                    self.commDictDouble[a]={}
                if(b not in self.commDictDouble[a].keys()):
                    self.commDictDouble[a][b]={}
                self.commDictDouble[a][b][c]=val,res,test[0]


    def ReadCommDictSingle(self):
        self.commDictSingle={}
        inny=open('redfield/refSingleComm.txt')
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                ops=test[0].split(']')[0]
                b=ops.split('[')[1].split(',')[0]#middle operator 
                a=ops.split(',')[1]              #inner operator
                num=test[0].split('=')[1].split('(')[0]  #
                res=test[0].split('(')[1].split(')')[0]  #
                if(num[-1]=='i'):
                    if(len(num[:-1].split('/'))>1):
                        val=eval( 'float('+num[:-1].split('/')[0]+')'+ '/' +'float('+ num[:-1].split('/')[1] +')')*complex(0,1)
                    else:
                        val=eval('float('+num[:-1]+')')*complex(0,1)
                else:
                    if(len(num[:-1].split('/'))>1):
                        val=eval( 'float('+num[:].split('/')[0]+')'+ '/' +'float('+ num[:].split('/')[1] +')')
                    else:
                        val=eval('float('+num[:]+')')

                vals=[]
                #vols=[]
                r=''
                vol=val
                for c in res:
                    if(c=='+' or c=='-'):
                        vals.append((vol,numpy.array(list(r))))
                        #vals.append((numpy.array(list(r))))
                        #vols.append(vol)
                        if(c=='+'):
                            vol=val
                        if(c=='-'):
                            vol=-1*val
                        r=''
                    else:
                        r+=c

                vals.append((vol,numpy.array(list(r))))
                #vals.append((numpy.array(list(r))))
                #vols.append(vol)

                if('a' not in b and 'b' not in b): #no need for alpha/beta in outer (hamiltonians)
                    if(a not in self.commDictSingle.keys()):
                        self.commDictSingle[a]={}
                    self.commDictSingle[a][b]=vals

    def TraceSym(self,keys,kois):
        tr=0.
        for ke in keys: #for each operator in key
            k1=ke[1]
            for ko in kois: #for each other operator in koi
                k2=ko[1]    #sum over individual products...
                trA,zz=self.DoTrace(k1,k2)
                if(zz):
                    tr+=trA*ko[0]*ke[0]
        return tr

    def DoTrace(self,k1,k2):
        trA=1
        for ka,kb in zip(k1,k2):
            try:
                trA*=self.TraceDict[ka][kb]
            except:
                return 0,False
            #if(kb in self.TraceDict[ka].keys()):
            #    #print i,k,k2[i],self.TraceDict[k][k2[i]]
            #    trA*=self.TraceDict[ka][kb]
            #else:
            #    return 0,False
        return trA,True

    #operators either C1{xyzab} or E in dictionary
    def GetTrSym(self,key):
        if(key=='E'):
            return 'E'
        else:
            return key[2]
    #read the trace dictionary
    def ReadTraceDict(self):
        self.TraceDict={}
        inny=open('redfield/refSingleTrace.txt')
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                a=self.GetTrSym(test[0].split('|')[0].split('<')[1])
                b=self.GetTrSym(test[0].split('>')[0].split('|')[1])
                c=eval('float('+test[0].split('=')[1]+')')
                if(a not in self.TraceDict.keys()):
                    self.TraceDict[a]={}
                self.TraceDict[a][b]=c
    
    def TestCommDictDouble(self):
        print 'Manually testing all double commutator expressions'
        #first verify commDict
        cnt=0
        for a in self.commDictDouble.keys():
            for b in self.commDictDouble[a].keys():
                for c in self.commDictDouble[a][b].keys():
                    comm1=comm(self.base[self.Expand(b)],self.base[self.Expand(a)])
                    comm2=comm(self.base[self.Expand(c)],comm1)
                    val=self.commDictDouble[a][b][c][0]
                    op=self.commDictDouble[a][b][c][1]
                    if(self.IsZero(comm2-val*self.base[self.Expand(op)])==False):
                        print 'testing'
                        print self.commDict[a][b][c][2]
                        print '['+c+',['+b+',[',a,']]='+str(val)+op
                        print '['+Expand(c)+',['+self.Expand(b)+',[',self.Expand(a),']]='+str(val)+self.Expand(op)
                        print 'shit'
                        print comm2
                        print val*self.base[self.Expand(op)]
                        sys.exit(100)
                    cnt+=1
        print 'Tested ',cnt,'expressions. All correct.'         
   
    def TestCommDictSingle(self):
        #first verify commDictSingle
        print 'Manually testing all single commutator expressions'
        cnt=0
        for a in self.commDictSingle.keys():
            for b in self.commDictSingle[a].keys():
                comm1=comm(self.base[self.Expand(b)],self.base[self.Expand(a)])

                vals=self.commDictSingle[a][b]
                op=''
                for i,val in enumerate(vals):
                    if(i==0):
                        refmat=copy.deepcopy(self.base[self.Expand(val[1])])*val[0]
                        op+=self.Expand(val[1])
                    else:
                        refmat+=self.base[self.Expand(val[1])]*val[0]
                        if(val[0].real<0 or val[0].imag<0):
                            op+='-'+self.Expand(val[1])
                        else:
                            op+='+'+self.Expand(val[1])
                #print
                #print refmat
                #print comm1
                if(self.IsZero(comm1-refmat)==False):
                    print 'testing'
                    print self.commDictSingle[a][b]
                    print '['+b+',[',a,']='+str(val)+op
                    print '['+self.Expand(b)+',[',self.Expand(a),']='+str(val)+'('+op+')'
                    print 'shit'
                    print comm2
                    print val*self.base[self.Expand(op)]
                    sys.exit(100)
                cnt+=1
        print 'Tested ',cnt,'expressions. All correct.'

    #take string operator, Exxyp and return complete one
    # note: basis and number count from right.
    #       the string operator counts from left.
    def Expand(self,a):
        coh=''
        #sys.exit(100)
        for i in range(self.baseSize):
            ii=self.baseSize-i-1
            if(a[ii]!='E'):
                coh+=self.baseType[ii]+str(ii+1)+a[ii]
        return coh

    def ExpandLatex(self,a):
        coh=''
        #sys.exit(100)
        for i in range(self.baseSize):
            ii=self.baseSize-i-1
            if(a[ii]!='E'):
                coh+=self.baseType[ii]+'^'+str(ii+1)+'_'+a[ii]
        return coh.replace('m','-').replace('p','+').replace('a','\\alpha ').replace('b','\\beta ')


    #calculate commutator [dip,key] symbollically
    #NOTE: rate limiting factor here are the join statements
    def CommSym(self,dip,vals,verb='n'):
        #requires outer operator to have only one element and for numerical weighting to be 1.
        #self.baseOp[coh]=(1,op,mask,active,compressedRep),  #save operator representation from makeop
        b=dip[0][3]  #active elements (xx,pm etc). Must be for two spins and single string f dictionary
        res=[] #store results
        #print vals
        #each entry=(val,op,compressRep)
        #verb='y'
        #print vals
        for val,op1,x,y,z in vals: #for each numerical factor+symbolic operator...
            a=op1[dip[0][2]].tostring() #needed for dictionary reference
            try:
                newvals=self.commDictSingle[a][b]
            except:
                continue 
            #if(verb=='y'):
            #    print '[',b,',',a,']=',newvals
            for newval in newvals: #append set of coherences as result in expanded format
                s=op1.copy()
                s[dip[0][2]]=newval[1] #update elements
                res.append((val*newval[0],s,1,1,s.tostring()))
                #value and expanded representation used for commutators
                #compressed representation used for indexing blist.

        #if(verb=='y'):
        #    print 'res:',res
        #    #return to standard format
        #    print '[',self.CommFormat(dip),',',self.CommFormat(vals),']=',self.CommFormat(res)
        #    #print 'num res:',
        #    #print comm(self.base[dip],self.CommVal(vals))

        if(len(res)==0):
            return 0
        else:
            return res

    #format commutator
    def CommFormat(self,vals):
        if(len(vals)==0):
            return 0
        res=''
        for i,val in enumerate(vals):
            num=val[0]
            com=0
            plus=1
            #print num
            if(numpy.fabs(val[0].imag)>0):
                num*=complex(0,1)*-1
                com=1
            if(num.imag>0):
                print 'problem with formatting!'
                sys.exit(100)

            if(num.real<0):
                plus=-1
            frac=Fraction(numpy.fabs(num.real)).limit_denominator()

            if(frac.denominator==1):
                numStr=str(frac.numerator)
            else:
                numStr=str(frac.numerator)+'/'+str(frac.denominator)
            if(com==1):
                numStr+='i'
            if(plus==1):
                res+='+'
            else:
                res+='-'
            res+=numStr+'('+str(self.Expand("".join(val[1])))+')'
        return res

    #work out matrix representation and check
    #symbolic commutator
    def CommVal(self,vals):
        if(len(vals)==0):
            return 0
        if(self.SYM):
            fug=True
        else:
            fug=False
        if(fug):
            self.SYM=False
        for i ,val in enumerate(vals):
            test=self.Expand(val[1])
            #if(test not in self.base.keys()):
            if(i==0):
                refmat=copy.deepcopy(self.MakeOp(test,save=False))*val[0] #create pre operator, if not already there
            else:
                refmat+=self.MakeOp(test,save=False)*val[0]
        if(fug):
            self.SYM=True
        return refmat




#make refSingleComm.txt:
#a file contianing all non-zero single commutators between 
#all operators in shift and cartesian basis
def MakeSingleCommutatorLibrary():
    import basis,copy,scipy,numpy

    inst=CommAux('CH')
    #inst.Sparse=False
    inst.base=GetBasis('CH',Sparse=True)  #get complete basis

    inst.GetDip('H2','C1','d',beta=0,alpha=0) #all dipolar 
    inst.GetCSA('C1','cA',beta=0,alpha=0)     #all dipolar-CSA
    inst.GetCSA('H2','cX',beta=0,alpha=0)     #all dipolar-CSA

    n1=numpy.array(inst.Aint[0])[:,1] #operators 1
    n2=numpy.array(inst.Aint[1])[:,1] #operators 2
    n3=numpy.array(inst.Aint[2])[:,1] #operators 3
    basisDip=numpy.concatenate((n1,n2,n3)) #all possible operators

    print 'making database of sums and differences'
    refbase={} #make a reference database including all sums and differences of operators
    for key,vals in inst.base.items():
        for koi,vols in inst.base.items():
            if(key!=koi):
                tast=vals+vols
                refbase[key+'P'+koi]=tast
                tast=vals-vols
                refbase[key+'M'+koi]=tast

    outy=open('redfield/refSingleComm.txt','w')
    for key,vals in inst.base.items(): #loop over all operators
        for dip in basisDip:  #for all Hamiltonian operators
            comm1=comm(inst.base[dip],vals) #do numerical commutator 1
            if(inst.IsZero(comm1)==0): #if not zero... find its symbolic identity
                frac,dom,resy,baso=compare(comm1,inst.base,maxy=4,even=True)
                if(resy=='none'): #have another go. #try again using the sums and differences database
                    frac,dom,resy,baso=compare(comm1,refbase,maxy=4,even=True)
                if(resy=='none'): #no symbolic operator for the result.
                    print
                    print 'shit'
                    print '[',dip,',',key,']=',frac,dom,resy
                    print resy,baso
                    print
                    sys.exit(100)
                else: #success. print to file
                    op1=inst.GetOpRef(key)
                    op2=inst.GetOpRef(dip)
                    op4=inst.GetOpRef(resy)
                    #print '[',key,',',koi,']=',frac,dom,resy
                    if(frac.denominator==1):
                        num=str(frac.numerator)
                    else:
                        num=str(frac.numerator)+'/'+str(frac.denominator)
                    if(dom=='i'):
                        num+='i'
                    #print
                    line ='['+dip+','+key+']='+num+'('+resy+')'
                    print line
                    line='['+op2+','+op1+']='+num+'('+op4+')'
                    outy.write(line+'\n')
                    print line
                                
                    #test commutator
                    nom=float(frac.numerator)/float(frac.denominator)
                    if(dom=='i'):
                        nom*=complex(0,1)

                    if(resy in inst.base.keys()):
                        if(inst.IsZero(comm1-nom*inst.base[resy])==False):
                            print 'shit'
                            print comm1
                            print nom*inst.base[resy]
                        else:
                            print 'True'
                    if(resy in refbase.keys()):
                        if(inst.IsZero(comm1-nom*refbase[resy])==False):
                            print 'shit'
                            print comm1
                            print nom*inst.base[resy]
                        else:
                            print 'True'
    outy.close()


def TraceSym(key,koi):
    trRef={}
    sz=1
    for i in range(sz):
        pass

#make refSingleTrace.txt:
#contain all single operator trace results
def MakeSingleTraceLibrary():
    import basis,copy,scipy,numpy
    inst=CommAux('C')
    #inst.Sparse=False
    inst.base=GetBasis('C',Sparse=True)  #get complete basis for 1 spin

    inst.GetCSA('C1','cA',beta=0,alpha=0)  #prepare all Hamiltonian operators
    basisDip=numpy.array(inst.Aint[0])[:,1]

    outy=open('redfield/refSingleTrace.txt','w')
    trRef={}
    for key,vals in inst.base.items(): #for all operators in the basis
        for koi,vols in inst.base.items(): #for all other operators in the basis
            tr=inst.Trace(vals,vols)  #calculate trace
            if(tr!=0): #if not zero, save.
                line='<'+key+'|'+koi+'>='+str(tr)
                outy.write(line+'\n')
                s1=inst.GetTrSym(key)
                s2=inst.GetTrSym(koi)
                if(s1 not in trRef.keys()):
                    trRef[s1]={}
                print s1,s2
                trRef[s1][s2]=tr
    outy.close()                    
         
    #test the trace dictionary
    for key,vals in inst.base.items(): #post
        for koi,vols in inst.base.items(): #pre
            s1=inst.GetTrSym(key)
            s2=inst.GetTrSym(koi)
            tr=inst.Trace(vals,vols)  #finish the expression
            trR=0
            if(s1 in trRef.keys()):
                if(s2 in trRef[s1].keys()):
                    trR=trRef[s1][s2] #indexed post/pre
            print '<',key,'|',koi,'>=',tr,trR
            

    #do more testing.
    """
    import basis,copy,scipy,numpy
    inst=CommAux('CH2')
    inst.base=basis.GetBasis('CH2',Sparse=True)  #get complete basis
    inst.GetDip('H2','C1','d',beta=0,alpha=0)
    inst.GetCSA('C1','cA',beta=0,alpha=0)
    inst.GetCSA('H2','cX',beta=0,alpha=0)
    n1=numpy.array(inst.Aint[0])[:,1]
    n2=numpy.array(inst.Aint[1])[:,1]
    n3=numpy.array(inst.Aint[2])[:,1]
    basisDip=numpy.concatenate((n1,n2,n3))
    for key,vals in inst.base.items():
        for koi,vols in inst.base.items():
            tr=basis.TraceSparse(vals,vols)  #finish the expression
            if(tr!=0):
                line='<'+key+'|'+koi+'>='+str(tr)
                print line,inst.TraceSym(key,koi)
    #check the trace dictionary
    for key,vals in inst.base.items():
        for koi,vols in inst.base.items():
            s1=inst.GetTrSym(key)
            s2=inst.GetTrSym(koi)
            tr=basis.TraceSparse(vals,vols)  #finish the expression
            trR=0
            if(s1 in trRef.keys()):
                if(s2 in trRef[s1].keys()):
                    trR=trRef[s1][s2]
            print '<',key,'|',koi,'>=',tr,trR
    """


#test the commutator libraries
def CommutatorTest():
    import basis,numpy
    inst=CommAux('CH')
    inst.base=GetBasis('CH',Sparse=True)  #get complete basis
    inst.GetDip('H2','C1','d',beta=0,alpha=0)
    basisDip=numpy.array(inst.Aint[0])[:,1]
    inst.ReadCommDictDouble()
    inst.TestCommDictDouble()
    inst.ReadCommDictSingle()
    inst.TestCommDictSingle()



def cohOrder(op):
    ps=int(op.count('p'))
    ms=int(op.count('m'))
    return ps-ms

def MakeDoubleCommutatorLibrary():
    import copy,scipy,numpy

    inst=CommAux('CH')
    inst.base=GetBasis('CH',Sparse=True)  #get complete basis

    inst.GetDip('H2','C1','d',beta=0,alpha=0)
    #inst.GetDip('H3','C1','d',beta=0,alpha=0)
    #inst.GetDip('H2','H3','d',beta=0,alpha=0)
    inst.Aint=numpy.array(inst.Aint)
    
    basisDip=inst.Aint[0,:,1]
    #basisDip=numpy.concatenate((basisDip,inst.Aint[1,:,1]))
    #basisDip=numpy.concatenate((basisDip,inst.Aint[2,:,1]))


    print 'making database of sums and differences'
    refbase={}
    #interestingly, not needed for double commutators

    for key,vals in inst.base.items():
        for koi,vols in inst.base.items():
            if(key!=koi):
                tast=vals+vols
                refbase[key+'P'+koi]=tast
                tast=vals-vols
                refbase[key+'M'+koi]=tast
                #print key+'M'+koi
    #sys.exit(100)
    if(1==1):
        outy=open('redfield/refDoubleComm.txt','w')
        for key,vals in inst.base.items():
            for dip in basisDip:
                comm1=comm(inst.base[dip],vals)
                if(inst.IsZero(comm1)==0):
                    for dip2 in basisDip:

                        op2=inst.GetOp(dip).tostring()
                        op3=inst.GetOp(dip2).tostring()

                        c1=cohOrder(op2)
                        c2=cohOrder(op3)
                        
                        if(c1+c2!=0): #condition pi-pj must be zero for tumbling.
                            continue #the double commutator will not occur.

                        comm2=comm(inst.base[dip2],comm1)
                        if(inst.IsZero(comm2)==0):
                            frac,dom,resy,baso=compare(comm2,inst.base,maxy=4,even=True)
                            bas=0
                            if(resy=='none'): #have another go.
                                frac,dom,resy,baso=compare(comm2,refbase,maxy=4,even=True)
                                bas=1
                            if(resy=='none'):
                                print
                                print 'shit'
                                print '[',dip2,',[',dip,',',key,']=',frac,dom,resy
                                print resy,baso
                                sys.exit(100)
                                print
                            else:

                                op1=inst.GetOp(key).tostring()

                                if(bas==1):
                                    test=resy.split('M')
                                    if(len(test)>1):
                                        op4a=inst.GetOp(test[0]).tostring()
                                        op4b=inst.GetOp(test[1]).tostring()
                                        op4=op4a+'-'+op4b
                                        bas=-1.
                                    else:
                                        test=resy.split('P')
                                        op4a=inst.GetOp(test[0]).tostring()
                                        op4b=inst.GetOp(test[1]).tostring()
                                        op4=op4a+'+'+op4b


                                else:
                                    op4=inst.GetOp(resy).tostring()


                                #print '[',key,',',koi,']=',frac,dom,resy
                                if(frac.denominator==1):
                                    num=str(frac.numerator)
                                else:
                                    num=str(frac.numerator)+'/'+str(frac.denominator)
                                if(dom=='i'):
                                    num+='i'
                                #print
                                line ='['+dip2+',['+dip+','+key+']]='+num+'('+resy+')'
                                print line
                                line='['+op3+',['+op2+','+op1+']]='+num+'('+op4+')'
                                outy.write(line+'\n')
                                print line
                                
                                #test commutator

                                nom=float(frac.numerator)/float(frac.denominator)
                                if(dom=='i'):
                                    nom*=complex(0,1)


                                if(resy in inst.base.keys() and inst.IsZero(comm2-nom*inst.base[resy])==False):
                                    print 'shit. error in double comm.'
                                    print comm2
                                    print nom*inst.base[resy]
                                    sys.exit(100)
                                elif(resy in refbase.keys() and inst.IsZero(comm2-nom*refbase[resy])==False):
                                    print 'shit. error in double comm.'
                                    print comm2
                                    print nom*inst.base[resy]
                                    sys.exit(100)

                                cm=inst.base[inst.Expand(op3)]
                                bm=inst.base[inst.Expand(op2)]
                                am=inst.base[inst.Expand(op1)]
                                if(bas!=0):
                                    res=inst.base[test[0]]+bas*inst.base[test[1]]
                                else:
                                    res=inst.base[inst.Expand(op4)]
                                comm2new=comm(cm,comm(bm,am))
                                if(inst.IsZero(comm2new-res*nom)==False):
                                    print comm2
                                    print
                                    print res*nom
                                    print
                                    print bas
                                    print 'double commutator failed'
                                    sys.exit(100)
        outy.close()


def DoCommTest4():
    testing=False  #numerically test symbolic result

    inst=CommAux('CH3')
    inst.base=GetBasis('CH3',Sparse=True)  #get complete basis

    inst.GetDip('H2','C1','d',beta=0,alpha=0)
    inst.GetDip('H3','C1','d',beta=0,alpha=0)
    inst.GetDip('H4','C1','d',beta=0,alpha=0)
    inst.GetDip('H2','H3','d',beta=0,alpha=0)
    inst.GetDip('H2','H4','d',beta=0,alpha=0)
    inst.GetDip('H3','H4','d',beta=0,alpha=0)
    inst.GetCSA('C1','cA',beta=0,alpha=0)

    n1=numpy.array(inst.Aint[0])[:,1]
    n2=numpy.array(inst.Aint[1])[:,1]
    n3=numpy.array(inst.Aint[2])[:,1]
    basisDip=numpy.concatenate((n1,n2,n3))

    comms={}
    #inst.ReadCommDictDouble()
    inst.ReadCommDictSingle()



    for key,vals in inst.base.items():
        inst.MakeOp(key)
        baseops=inst.baseOp[key]

        op1=inst.GetOp(key).tostring()
        if(op1.count('y')>0):
            continue
        if(op1.count('x')>0):
            continue

        for dip in basisDip:
            #comm1=comm(inst.base[dip],vals)
            #res1=inst.CommSym(dip,((1,key),),verb='y') #take operators: return symbolic result
            a=inst.CommSym(inst.baseOp[dip],baseops) #[O2,[O1,pre]]
            if(a==0):
                continue
                        
            for dip2 in basisDip: 

                op2=inst.GetOp(dip).tostring()
                op3=inst.GetOp(dip2).tostring()

                c1=cohOrder(op2)
                c2=cohOrder(op3)
                    
                if(c1+c2!=0): #condition pi-pj must be zero for tumbling.
                    continue #the double commutator will not occur.

                b=inst.CommSym(inst.baseOp[dip2],a) #[O2,[O1,pre]]
                #comm2=basis.comm(inst.base[dip2],comm1)
                #res2=inst.CommSym(dip2,res1,verb='n') #take operators: return symbolic result
                if(b==0):
                    continue
                


                stry=''


                c1=cohOrder(op1)
                cref=cohOrder(b[0][1].tostring())
                for op in b:
                    #print op
                    stry+='('
                    if(op[0].real>0):                    
                        num=Fraction(op[0].real)
                        if(num.denominator==1):
                            stry+='+'+str(num.numerator)
                        else:
                            stry+='+'+str(num.numerator)+'/'+str(num.denominator)                        

                    elif(op[0].real<0):
                        num=Fraction(-1*op[0].real)
                        if(num.denominator==1):
                            stry+='-'+str(num.numerator)
                        else:
                            stry+='-'+str(num.numerator)+'/'+str(num.denominator)                        

                    if(op[0].imag>0):
                        num=Fraction(op[0].imag)
                        if(num.denominator==1):
                            stry+='+'+str(num.numerator)
                        else:
                            stry+='+'+str(num.numerator)+'/'+str(num.denominator)                        
                        stry+='i'
                    elif(op[0].imag<0):
                        num=Fraction(-1*op[0].imag)
                        if(num.denominator==1):
                            stry+='-'+str(num.numerator)
                        else:
                            stry+='-'+str(num.numerator)+'/'+str(num.denominator)                        
                        stry+='i'

                    stry+=')'+op[1].tostring()
                    c2=cohOrder(op[1].tostring())
                    if(c2!=cref):
                        print 'shit! coherence order change'
                        print c2,cref
                        sys.exit(100)
                    

                print '['+op3+',['+op2+','+op1+']=',stry

                comms['['+op3+',['+op2+','+op1+']=']=stry

                if(c1!=cref):
                    print 'Shit! Change in coherence order!'
                    print c1,cref
                    sys.exit(100)
                #if(c1==4):
                #    sys.exit(100)

                if(testing):
                    refmat2=inst.CommVal(b) #calculate corresponding matrix
                    if(inst.IsZero(comm(inst.base[dip2],comm(inst.base[dip],inst.base[key]))-refmat2)==False):
                            print 'shit. commutators do not add up'
                            sys.exit(100)
        print 'Number of unique non-zero double comms:',len(comms.keys())
