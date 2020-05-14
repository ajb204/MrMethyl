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
#         self.DoTestSystem(macro=True,ax=True,xx=True,cx=True,ca=True, sparse=True, sym=True, irreducible=False, num=False)
#         sys.exit(100)
        
        
        self.cross=CommAux('CH3CH3',latex_file='output/eps/cross_latex.tex',sparse=sparse,sym=sym)
        #OPEN LOG FILE FOR WRITING TO:
        print_log=open('output/log/log_Comm','w')
        self.call_func=0
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
        cross_matrix,reflected,d1,d2,c1,h1,c2,h2,self.map_key_nmr,self.map_index_nmr,pairs,mask_matrix=self.peak_file_parser(peak_file_dir,peak_list_dir)
        #Generate the "wobble mask" which screens out guys who are dynamic:
        wobble_mask=self.wobble_masker(cross_matrix,self.map_index_nmr,self.map_key_nmr,c1,h1,c2,h2)
        #Generate data mask:
        self.mask=numpy.logical_and(wobble_mask,mask_matrix)
        #Mask out the diagonals: (I think we want to do this?)
        diag_mask=numpy.ones(self.mask.shape,dtype=bool)
        numpy.fill_diagonal(diag_mask,False) 
        self.peak_mask=numpy.logical_and(self.mask,diag_mask)
        
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
        self.num_samples=4000#SET BY AJB: WAS 2000
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
        for key in self.map_key_nmr:
            xrd_indices.append(self.map_key_xrd[key])
            nmr_indices.append(self.map_key_nmr[key])
        xrd_indices,nmr_indices=numpy.array(xrd_indices),numpy.array(nmr_indices)
        self.ii_xrd,self.jj_xrd=numpy.meshgrid(xrd_indices,xrd_indices,indexing='ij')
        self.ii_nmr,self.jj_nmr=numpy.meshgrid(nmr_indices,nmr_indices,indexing='ij')

        #print self.ii_nmr.shape
        #sys.exit(100)

        #Create and calculate distance matrix:
        self.dist=numpy.zeros(self.ii_nmr.shape)
        dist_matrix=numpy.sqrt(numpy.sum(v_mm*v_mm,axis=2))
        self.dist[self.ii_nmr,self.jj_nmr]=dist_matrix[self.ii_xrd,self.jj_xrd]
        
        #Create and fill the NMR peak matrix:
        self.nmr_peaks=numpy.zeros(self.ii_nmr.shape)
        self.nmr_peaks[self.ii_nmr,self.jj_nmr]=cross_matrix[self.ii_nmr,self.jj_nmr]
        #Create the XRD peak matrix: (This will be populated in the optimization loop, when we do a function call)
        self.xrd_peaks=numpy.zeros(self.ii_nmr.shape)
        
        #INITIALIZATION HERE:
        #Initialize the initial intensities as the diagonal NMR peaks (this should do as a halfway decent proxy)
        d=numpy.arange(len(nmr_indices))
        self.init_intensities=numpy.zeros(self.nmr_peaks.shape)
        self.init_intensities[d,d]=self.nmr_peaks[self.ii_nmr,self.jj_nmr][d,d]*3.
        #Initialize t_m,t_c:
        self.cross.tm=5.E-12
        self.cross.tc=15.E-9
        
        
        
        xinit=self.pack()
        print xinit
        print 'xinit:',xinit
        print 'Initial model run:',self.GetChi2(xinit)
        x0=leastsq(self.chi_func,xinit)
        self.unpack(x0[0])
        xpost=self.pack()
        
       
                
                
        
        print 'xpost:',xpost
        print 'Post model run:',self.GetChi2(xpost)
        i,j=numpy.arange(self.nmr_peaks.shape[0]),numpy.arange(self.nmr_peaks.shape[0])
        ii,jj=numpy.meshgrid(i,j,indexing='ij')
        for nmr,xrd,i,j in zip(self.nmr_peaks[self.peak_mask],self.xrd_peaks[self.peak_mask],ii[self.peak_mask],jj[self.peak_mask]):
            print self.map_index_nmr[i],self.map_index_nmr[j],'%0.5e,%0.5e,%0.5e,%0.5e'%(nmr,xrd,nmr-xrd,(nmr-xrd)/nmr*100.)
        
        x0_jiggled=self.jiggler()
        
        print 'x0_jiggled:',x0_jiggled
        print 'Post jiggled run:',self.GetChi2(x0_jiggled[0])
        print 'pearsonr = ',pearsonr(self.nmr_peaks[self.peak_mask],self.xrd_peaks[self.peak_mask])
        print 'spearmanr = ',spearmanr(self.nmr_peaks[self.peak_mask],self.xrd_peaks[self.peak_mask])
        print "% error sum = \n",numpy.average((self.nmr_peaks[self.peak_mask]-self.xrd_peaks[self.peak_mask])/self.nmr_peaks[self.peak_mask]*100.)
        for nmr,xrd,i,j in zip(self.nmr_peaks[self.peak_mask],self.xrd_peaks[self.peak_mask],ii[self.peak_mask],jj[self.peak_mask]):
            print self.map_index_nmr[i],self.map_index_nmr[j],'%0.5e,%0.5e,%0.5e,%0.5e'%(nmr,xrd,nmr-xrd,(nmr-xrd)/nmr*100.)
        
         #TEST FOR NUMBER OF ITERATIONS NEEDED:
        for num_samples in (250,500,1000,1500,2000,2500,4000,8000):
            self.num_samples=num_samples
            print num_samples
            num=self.xrd_peaks.shape[0]
            data_matrix=numpy.zeros((num,num,10))
            diag_matrix=numpy.zeros((num,10))
            for i in numpy.arange(10):
                xrd_old=copy.deepcopy(self.xrd_peaks)
                start=datetime.now()
                self.CalcModel()
                data_matrix[:,:,i]=self.xrd_peaks
                diag_matrix[:,i]=(self.xrd_peaks.diagonal()-self.nmr_peaks.diagonal())/self.nmr_peaks.diagonal()
                end=datetime.now()
                print end-start
#                 print end-start,numpy.sum((self.xrd_peaks[self.peak_mask]-xrd_old[self.peak_mask])**2),numpy.average((self.xrd_peaks[self.peak_mask]-xrd_old[self.peak_mask])**2)/numpy.average(self.xrd_peaks[self.peak_mask]**2)*100.
            mean_matrix=numpy.average(data_matrix,axis=2)
            std_matrix=numpy.std(data_matrix,axis=2)
            plt.scatter(self.nmr_peaks[self.peak_mask],mean_matrix[self.peak_mask])
            plt.show()
            plt.errorbar(self.nmr_peaks[self.peak_mask],mean_matrix[self.peak_mask],yerr=std_matrix[self.peak_mask],linestyle='None',fmt='o')
            plt.title('Correlation of simulated and experimental cross-peak intensities\nSampling rate = '+str(num_samples))
            plt.text(0.,1.,str(pearsonr(self.nmr_peaks[self.peak_mask],mean_matrix[self.peak_mask])))
            plt.xlabel('Experimental cross-peak intensity')
            plt.ylabel('Simulated cross-peak intensity')
            plt.savefig('figs/correlation_plot'+str(num_samples)+'.svg')
            plt.show()
            plt.errorbar(self.dist[self.peak_mask],mean_matrix[self.peak_mask],yerr=std_matrix[self.peak_mask],linestyle='None',fmt='o')
            plt.title('Simulated cross-peak intensities as a function of C-C separation\nSampling rate = '+str(num_samples))
            plt.xlabel('Distance/Angstroms')
            plt.ylabel('Simulated cross-peak intensity')
            plt.savefig('figs/dist_plot'+str(num_samples)+'.svg')
            plt.show()
            
            
            mean_diag=numpy.average(diag_matrix,axis=1)
            std_diag=numpy.std(diag_matrix,axis=1)
            res_indices=numpy.arange(mean_diag.shape[0])
            plt.bar(res_indices,mean_diag,yerr=std_diag)
            plt.title('Order parameter from diagonal intensities\nSampling rate = '+str(num_samples))
            plt.xlabel('Residue number')
            plt.ylabel('Order parameter from diagonal intensities')
            plt.savefig('figs/order_plot'+str(num_samples)+'.svg')
            plt.show()
            
            
            
            numpy.savetxt('output/mean_matrix'+str(num_samples),mean_matrix)
            numpy.savetxt('output/std_matrix'+str(num_samples),std_matrix)
            numpy.savetxt('output/nmr_matrix'+str(num_samples),self.nmr_peaks)
        
        
        plt.scatter(self.nmr_peaks[self.peak_mask],self.xrd_peaks[self.peak_mask])
        plt.title("Correlation plot - Intensity of simulated against experimental")
        plt.xlabel("Experimental intensities")
        plt.ylabel("Simulated intensities")
        plt.savefig('output/figs/non_sampling_correlation.png',dpi=900)
        plt.show()
        
        plt.scatter(self.dist[self.peak_mask],self.xrd_peaks[self.peak_mask])
        plt.title("Simulated intensity against distance")
        plt.xlabel("Distance/Angstroms")
        plt.ylabel("Simulated intensities")
        plt.savefig('output/figs/non_sampling_sim_dist.png',dpi=900)
        plt.show()
        
        plt.scatter(self.dist[self.peak_mask],self.nmr_peaks[self.peak_mask])
        plt.title("Experimental intensity against distance")
        plt.xlabel("Distance/Angstroms")
        plt.ylabel("Experimental intensities")
        plt.savefig('output/figs/non_sampling_exp_dist.png',dpi=900)
        plt.show()
        
#         self.bootVals=100
#         self.booty(x0_jiggled[0])
        
        
        #Now run a sampled calculation:
        print 'Sampled calculation: '
        print 'x0_jiggled:',x0_jiggled
        print 'Post jiggled run:',self.GetChi2_sampled(x0_jiggled[0])
        print 'pearsonr = ',pearsonr(self.nmr_peaks[self.peak_mask],self.xrd_peaks[self.peak_mask])
        print 'spearmanr = ',spearmanr(self.nmr_peaks[self.peak_mask],self.xrd_peaks[self.peak_mask])
        print "% error sum = \n",numpy.average((self.nmr_peaks[self.peak_mask]-self.xrd_peaks[self.peak_mask])/self.nmr_peaks[self.peak_mask]*100.)
        for nmr,xrd,i,j in zip(self.nmr_peaks[self.peak_mask],self.xrd_peaks[self.peak_mask],ii[self.peak_mask],jj[self.peak_mask]):
            print self.map_index_nmr[i],self.map_index_nmr[j],'%0.5e,%0.5e,%0.5e,%0.5e'%(nmr,xrd,nmr-xrd,(nmr-xrd)/nmr*100.)
        
        plt.scatter(self.nmr_peaks[self.peak_mask],self.xrd_peaks[self.peak_mask])
        plt.title("Correlation plot - Intensity of simulated against experimental")
        plt.xlabel("Experimental intensities")
        plt.ylabel("Simulated intensities")
        plt.savefig('output/figs/sample_correlation.png',dpi=900)
        plt.show()
        
        plt.scatter(self.dist[self.peak_mask],self.xrd_peaks[self.peak_mask])
        plt.title("Simulated intensity against distance")
        plt.xlabel("Distance/Angstroms")
        plt.ylabel("Simulated intensities")
        plt.savefig('output/figs/sampled_sim_dist.png',dpi=900)
        plt.show()
        
        plt.scatter(self.dist[self.peak_mask],self.nmr_peaks[self.peak_mask])
        plt.title("Experimental intensity against distance")
        plt.xlabel("Distance/Angstroms")
        plt.ylabel("Experimental intensities")
        plt.savefig('output/figs/sampled_dist.png',dpi=900)
        plt.show()
        
        for i in numpy.arange(self.nmr_peaks.shape[0]):
            print self.map_index_nmr[i],'%0.5e\t %0.5e'%(self.nmr_peaks[i,i],self.xrd_peaks[i,i])
        
        self.Package=[]
        for rep in numpy.arange(self.num_samples): #calculate sampling sechdules.
            list_methyls=self.methyl_list_generator(self.residues)
            index_list=numpy.array([self.map_key_xrd[x] for x in list_methyls]) #this can be pre-computed
            #i,j are for the expanded matrix. Maps the compacted matrix onto the expanded matrix. 
            ii,jj=numpy.meshgrid(index_list,index_list,indexing='ij')
            #Load up and take matrix exponential:
            self.Package.append((ii,jj)) #numpy method
        
        
        
        sys.exit(100)
        
        
        #NOW SET UP THE LOOP:
        for n in numpy.arange(num_loop):
            self.cross.tm=20.*numpy.random.random()*1.E-12
            self.cross.tc=(5.+10.*numpy.random.random())*1.E-9
            self.init_intensities[d,d]=self.nmr_peaks[self.ii_nmr,self.jj_nmr][d,d]*(1.+4.*numpy.random.random(d.shape))
            x_init=self.pack()
            x0=leastsq(self.chi_func,x_init)
            self.unpack(x0[0])
            chi2=self.GetChi2(x0[0])
            
            i,j=numpy.arange(self.peak_mask.shape[0]),numpy.arange(self.peak_mask.shape[0])
            ii,jj=numpy.meshgrid(i,j,indexing='ij')
            print n,'\t',chi2,'\t','\n',x_init,'\n',x0[0]
            for i,j in zip(ii[self.peak_mask],jj[self.peak_mask]):
                print self.map_index_nmr[i],self.map_index_nmr[j],self.xrd_peaks[i,j],self.nmr_peaks[i,j],100.*(self.nmr_peaks[i,j]-self.xrd_peaks[i,j])/self.nmr_peaks[i,j]
            print '\n'
            if(n==0):
                min_chi=chi2
                min_parameters=x0[0]
            else:
                if(chi2<min_chi):
                    min_chi=chi2
                    min_parameters=x0[0]
        
        print 'min_chi,min_parameters'
        print min_chi,min_parameters
        chi2=self.GetChi2(min_parameters)
        print chi2,'\n',min_parameters
        for i,j in zip(ii[self.peak_mask],jj[self.peak_mask]):
            print self.map_index_nmr[i],self.map_index_nmr[j],self.xrd_peaks[i,j],self.nmr_peaks[i,j],100.*(self.nmr_peaks[i,j]-self.xrd_peaks[i,j])/self.nmr_peaks[i,j]
        sys.exit(100)

        #DEBUGGING:
        print self.nmr_peaks[self.peak_mask]
        #Run the optimization:
        x_init=self.pack()
        print x_init
        print x_init.shape
        print 'xinit:',x_init
        print 'Initial model run:',self.GetChi2(x_init)
        print 'Initial chi2:',self.GetChi2(x_init)
        x0=leastsq(self.chi_func,x_init)
        self.unpack(x0[0])
        print 'Final chi2:  ',self.GetChi2(x0[0])
        print'OPTIMIZED PARAMETERS: \n',x0[0]
        numpy.savetxt('output/log/optimized_parameters',x0[0])


        outy=open('output/outpars.out','w')
        for i in range(len(self.nmr_peaks[self.peak_mask])):
            outy.write('%e\t%e\n' % (self.nmr_peaks[self.peak_mask][i],self.xrd_peaks[self.peak_mask][i]))
        outy.close()
        #print self.nmr_peaks[self.peak_mask]
        #print self.xrd_peaks[self.peak_mask]
    
        print 'R2:',(scipy.stats.pearsonr(self.nmr_peaks[self.peak_mask],self.xrd_peaks[self.peak_mask]))[0]**2.
        print 'aveErr%:',  numpy.average(numpy.fabs(self.nmr_peaks[self.peak_mask]-self.xrd_peaks[self.peak_mask])/numpy.max(self.nmr_peaks[self.peak_mask])*100)

    
    def unpack(self,x):
#         self.cross.tm=x[0,1]
#         self.cross.tc=x[1,1]
#         self.init_intensities=x[:,0]
        
        #Adapt to set tau_m to effectively infinitely
#         self.cross.tm=1.E20
        self.cross.tm=numpy.absolute(x[0])
        self.cross.tc=numpy.absolute(x[1])
        d=numpy.arange(self.nmr_peaks.shape[0])
        #self.init_intensities[d,d][self.mask[d,d]]=copy.deepcopy(x[2:]*1E8)
        self.init_intensities[self.mask[d,d][d],self.mask[d,d][d]]=copy.deepcopy(x[2:])
        
        #print '%.9e\t%.9e\t%.9e' % (x[2],x[3],x[4])
        #print '%.9e\t%.9e\t%.9e' % (self.init_intensities[0,0],self.init_intensities[1,1],self.init_intensities[2,2])

        #print self.init_intensities
    
    def pack(self):
        d=numpy.arange(self.nmr_peaks.shape[0])
        x=numpy.array([self.cross.tm,self.cross.tc])
        x=numpy.append(x,self.init_intensities[d,d][self.mask[d,d]])
        return x
        
    def GetChi2(self,x):
        self.unpack(x)
#         self.CalcModel()
        self.CalcModel_non_sampling()
#         for i,j in zip(self.ii_nmr[self.peak_mask],self.jj_nmr[self.peak_mask]):
#             print self.map_index_nmr[i],self.map_index_nmr[j],self.xrd_peaks[i,j],self.nmr_peaks[i,j],100.*(self.nmr_peaks[i,j]-self.xrd_peaks[i,j])/self.nmr_peaks[i,j]
#         print '\n'
#         print 'xrd_peaks',self.xrd_peaks[self.peak_mask]
#         print 'nmr_peaks',self.nmr_peaks[self.peak_mask]
#         print 'difference',100.*(self.nmr_peaks[self.peak_mask]-self.xrd_peaks[self.peak_mask])/self.nmr_peaks[self.peak_mask]
        return numpy.sum((self.nmr_peaks[self.peak_mask]-self.xrd_peaks[self.peak_mask])**2.)


    def GetChi2_sampled(self,x):
        self.unpack(x)
        self.CalcModel()
#         for i,j in zip(self.ii_nmr[self.peak_mask],self.jj_nmr[self.peak_mask]):
#             print self.map_index_nmr[i],self.map_index_nmr[j],self.xrd_peaks[i,j],self.nmr_peaks[i,j],100.*(self.nmr_peaks[i,j]-self.xrd_peaks[i,j])/self.nmr_peaks[i,j]
#         print '\n'
#         print 'xrd_peaks',self.xrd_peaks[self.peak_mask]
#         print 'nmr_peaks',self.nmr_peaks[self.peak_mask]
#         print 'difference',100.*(self.nmr_peaks[self.peak_mask]-self.xrd_peaks[self.peak_mask])/self.nmr_peaks[self.peak_mask]
        return numpy.sum((self.nmr_peaks[self.peak_mask]-self.xrd_peaks[self.peak_mask])**2.)
    
    def chi_func(self,x):
        self.unpack(x)
#         self.CalcModel()
        self.CalcModel_non_sampling()
        return (self.nmr_peaks[self.peak_mask]-self.xrd_peaks[self.peak_mask]).flatten()
        
        
    def chi_func_boot(self,x):
        self.unpack(x)
#         self.CalcModel()
        self.CalcModel_non_sampling()
        return (self.nmr_peaks[self.peak_mask][self.boot_mask]-self.xrd_peaks[self.peak_mask][self.boot_mask]).flatten()        
    
    
    def CalcModel_non_sampling(self):
        self.call_func+=1
        #Calculate the rates:
        r_cross=self.cross.CalcRate('M1z', 'M2z',verb='n')
        r_auto=self.cross.CalcRate('M1z', 'M1z',verb='n')
        #Create matrices to store results:
        num_methyls=r_cross.shape[0]
        i,j=numpy.arange(num_methyls),numpy.arange(num_methyls)
        ii,jj=numpy.meshgrid(i,j,indexing='ij')
        methyl_mask=numpy.zeros(r_cross.shape)
        methyl_mask,bool_mask=self.methyl_mask_generator(self.residues,methyl_mask)
        rate_matrix=numpy.zeros(r_cross.shape)
        sim_matrix=numpy.zeros(r_cross.shape)
        
        
        mask1=numpy.ones(r_cross.shape)
        mask1[bool_mask]=0.
        #Make adjustments based on populations:
        r_cross=r_cross*methyl_mask
        r_auto=r_auto*methyl_mask
        
        #Construct rate matrix:
        rate_matrix[ii,jj]=r_cross[ii,jj]
        rate_matrix[i,j]=numpy.sum((r_auto[ii,jj]-r_auto[ii,ii]),axis=1)+r_auto[i,j]
        #Take matrix exponential:
        sim_matrix=expm(-rate_matrix[:,:]*self.tau)*methyl_mask
        #Multiply by initial intensities:
        self.xrd_peaks[self.ii_nmr,self.jj_nmr]=sim_matrix[self.ii_xrd,self.jj_xrd]*numpy.absolute(self.init_intensities[self.jj_nmr,self.jj_nmr])
    
    def CalcModel(self):
        self.call_func+=1
        #Calculate the rates:
        r_cross=self.cross.CalcRate('M1z', 'M2z',verb='n')
        r_auto=self.cross.CalcRate('M1z', 'M1z',verb='n')
        #Create matrices to store results:
        num_methyls=r_cross.shape[0]
        
        result_matrix=numpy.zeros((num_methyls,num_methyls,self.num_samples))
        #Calculate intensities by repeating for 2000 combinations:

        for rep in numpy.arange(self.num_samples):
            list_methyls=self.methyl_list_generator(self.residues)
            #print list_methyls
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
        #print result_matrix
            
        #AVERAGE OVER 2000 samples:
        sim_matrix=numpy.average(result_matrix[:,:,:], axis=2)
        
        numpy.set_printoptions(floatmode='unique')
        #print sim_matrix*self.init_intensities[self.jj_nmr,self.jj_nmr]

        #MATCH SIMULATED XRD DATA ONTO NMR INDICES:
        self.xrd_peaks[self.ii_nmr,self.jj_nmr]=sim_matrix[self.ii_xrd,self.jj_xrd]*numpy.absolute(self.init_intensities[self.jj_nmr,self.jj_nmr])
#         print self.xrd_peaks[self.peak_mask]
#         print self.nmr_peaks[self.peak_mask]
    
    
    def color_calc(self,frac):
        r,g,b=0.,0.,0.
#         r=0.5*(-frac+1.)
#         b=0.5*(1.+frac)
        r=frac
        b=1.-frac
        return r,g,b

    def frac_calc(self,nmr,xrd):
        hypo=numpy.sqrt(nmr**2 + xrd**2) #Calculate "hypotenuse" to normalize fraction to -1<x<1
        frac=(nmr-xrd)/hypo
        return frac
    
    def CalcModel_pre_sampled(self):
        #Calculate the rates:
        r_cross=self.cross.CalcRate('M1z', 'M2z',verb='n')
        r_auto=self.cross.CalcRate('M1z', 'M1z',verb='n')
        self.num_methyls=r_cross.shape[0]
        num_samples=len(self.Package)
        r_diff=numpy.transpose(numpy.transpose(r_auto)-r_auto[numpy.diag_indices(self.num_methyls)])#subtract off the diagonal (same as r_auto[ii,jj]-r_auto[ii,ii]?)
        r_cross[numpy.diag_indices(self.num_methyls)]=r_auto[numpy.diag_indices(self.num_methyls)] #do this once
        result_matrix=numpy.zeros((self.num_methyls,self.num_methyls))
        for ii,jj in self.Package:
            rate_matrix=r_cross[ii,jj] #diagonal r_auto already included
            rate_matrix.ravel()[::rate_matrix.shape[1]+1]+=numpy.sum(r_diff[ii,jj],axis=1) #+r_auto[i,j]  #add sum of off diagonals to diagonal
            result_matrix[ii,jj]+=scipy.linalg.expm(-rate_matrix*self.tau)#update matrix
        result_matrix/=num_samples
            
        self.xrd_peaks[self.ii_nmr,self.jj_nmr]=result_matrix[self.ii_xrd,self.jj_xrd]*numpy.absolute(self.init_intensities[self.jj_nmr,self.jj_nmr])
        
    
    def convergence_test(self):
        r_cross=self.cross.CalcRate('M1z', 'M2z',verb='n')
        r_auto=self.cross.CalcRate('M1z', 'M1z',verb='n')
        #Create matrices to store results:
        num_methyls=r_cross.shape[0]
        
        
        num_samples=50000
        result_matrix=numpy.zeros((num_methyls,num_methyls,num_samples))
        for rep in numpy.arange(num_samples):
            list_methyls=self.methyl_list_generator(self.residues)
            #print list_methyls
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
        #print result_matrix
            
        #PREPARE COMPARISON ARRAY:
        compare_matrix=numpy.average(result_matrix[:,:,:], axis=2)
            
        #AVERAGE OVER SAMPLE NUMBER:
        num_arr=numpy.array([1,10,100,250,500,1000,1500,2000,3000,4000,5000,7500,10000,15000,20000,25000,40000,50000])
        compare_list=[]
        for reps in num_arr:
            arr=numpy.average(result_matrix[:,:,0:reps],axis=2)
            chi_2=numpy.sum((compare_matrix-arr)**2)
            compare_list.append(chi_2)
        compare_list=numpy.array(compare_list)
        print compare_list
        print num_arr
        print numpy.sum(compare_matrix**2)
        plt.show()
        plt.scatter(num_arr,compare_list)
        plt.xlabel("Number of samples")
        plt.ylabel("chi**2 - Deviation from converged result")
        plt.title("Convergence dependent on sampling rate")
        plt.savefig('output/figs/convergence_plt',dpi=900)
        plt.show()
        
        numpy.set_printoptions(floatmode='unique')
        #print sim_matrix*self.init_intensities[self.jj_nmr,self.jj_nmr]

        
        
        
    
    def fake(self):
        num=numpy.sum(self.peak_mask)
        self.boot_mask=numpy.random.randint(0,num,num)
    
    def booty(self,xinit):
        print 'Running the bootstrap:'
        self.xinitBack=copy.deepcopy(xinit)
        booty=[]
        for i in range(self.bootVals):
            #print i
            #make fake datset
            self.fake()
            
            x0=leastsq(self.chi_func_boot,self.xinitBack)
            booty.append(x0[0])
        booty=numpy.array(booty)
        
        for i in numpy.arange(booty.shape[1]):
            numpy.savetxt('output/parameter_list/par'+str(i),booty[:,i])
            
        plt.hist(booty[:,0]*1.E12,bins=20)
        plt.title('tau_m histogram plot: ')
        plt.xlabel('tau_m/ps')
        plt.ylabel('Number of hits')
        plt.savefig('output/figs/non_sampling_tau_m.png',dpi=900)
        plt.show()
        
        plt.hist(booty[:,1]*1.E9,bins=20)
        plt.title('tau_c histogram plot: ')
        plt.xlabel('tau_c/ns')
        plt.ylabel('Number of hits')
        plt.savefig('output/figs/non_sampling_tau_c.png',dpi=900)
        plt.show()
        
        print 'tau_m:',xinit[0],numpy.average(booty[:,0]),numpy.std(booty[:,0])
        print 'tau_c:',xinit[1],numpy.average(booty[:,1]),numpy.std(booty[:,1])
    
    def AddHeat(self,pars):
        print 'pars before heat:',pars
        parNew=pars*numpy.random.normal(1.,self.sigma,len(pars))
        print 'pars after heat: ',parNew
        return parNew
    

    def jiggler(self):
        print 'Unleashing the jiggler...'
        self.jiggles=25  #max number of goes before aborting
        self.sigma=0.10   #standard deviation of normal distribution

        xinit=self.pack()
        x0=leastsq(self.chi_func,xinit)

        currParBest=copy.deepcopy(x0[0])
        chi2Best=self.GetChi2(currParBest)

        cnt=0
        while(1==1):
            parCurr=copy.deepcopy(currParBest)
            parCurr=self.AddHeat(parCurr)
            x0=leastsq(self.chi_func,parCurr)
            chi2Curr=self.GetChi2(x0[0])
            print 'best:',chi2Best,'current:',chi2Curr
            if(chi2Curr<chi2Best):
                print '    yay! we have improved! saving.'
                i,j=numpy.arange(self.nmr_peaks.shape[0]),numpy.arange(self.nmr_peaks.shape[0])
                ii,jj=numpy.meshgrid(i,j,indexing='ij')
                for nmr,xrd,i,j in zip(self.nmr_peaks[self.peak_mask],self.xrd_peaks[self.peak_mask],ii[self.peak_mask],jj[self.peak_mask]):
                    print self.map_index_nmr[i],self.map_index_nmr[j],'%0.5e,%0.5e,%0.5e,%0.5e'%(nmr,xrd,nmr-xrd,(nmr-xrd)/nmr*100.)
                self.currParBest=copy.deepcopy(x0[0])
                chi2Best=chi2Curr
                cnt=0
            else:
                print 'no improvment. number of jiggles:',cnt

            cnt+=1
            if(cnt==self.jiggles):
                break
            
    
        x0=leastsq(self.chi_func,currParBest)
        chi2Final=self.GetChi2(x0[0])
        print 'final chi2:',chi2Final
        return x0
        
        
        
        
        
        
        
        
        
    #NOW HAVE ALL THE OTHER RANDOM FUNCTIONS:    
    def wobble_masker(self,cross_matrix,map_index_nmr,map_key_nmr,c1,h1,c2,h2,lim=0.85):
        write_file=open('output/wobble_mask','w')
        compare_file=open('output/compare_file','w')
        
        indices=[]
        for i in numpy.arange(cross_matrix.shape[0]):
            res=map_index_nmr[i]
            res_num=int(res.split()[1])
            if(res[0:3]=='ALA'):
                i=map_key_nmr[res]
                if(c1[i,i]!=0.):
                    indices.append(i)
            elif(res[0:3]=='ILE'):
                p_g=(14.8-c1[i,i])/5.5
                r,g,b=self.color_calc(1.-2.*numpy.sqrt((p_g-0.5)**2))
                write_file.write('set_color color%i = [%.2f,%.2f,%.2f]\n'%(res_num,r,g,b))
                write_file.write('color color%i, resi %i\n'%(res_num,res_num))
                compare_file.write('%i\t%.2f\n'%(res_num,(2.*numpy.sqrt((p_g-0.5)**2))))
                if(c1[i,i]!=0. and (p_g<1.-lim or p_g>lim)):
                    indices.append(i)
            elif(res[0:3]=='LEU'):
                res_a,res_b=res[0:-1]+'A',res[0:-1]+'B'
                i_a,i_b=map_key_nmr[res_a],map_key_nmr[res_b]
                c_a,c_b=c2[i_a,i_a],c2[i_b,i_b]
                p_t=0.5+(c_a-c_b)/10.
                r,g,b=self.color_calc(1.-2.*numpy.sqrt((p_t-0.5)**2))
                write_file.write('set_color color%i = [%.2f,%.2f,%.2f]\n'%(res_num,r,g,b))
                write_file.write('color color%i, resi %i\n'%(res_num,res_num))
                compare_file.write('%i\t%.2f\n'%(res_num,(2.*numpy.sqrt((p_t-0.5)**2))))
                if(c_a!=0 and c_b!=0 and (p_t<1.-lim or p_t>lim)):
                    indices.append(i)
            elif(res[0:3]=='VAL'):
                res_a,res_b=res[0:-1]+'A',res[0:-1]+'B'
                i_a,i_b=map_key_nmr[res_a],map_key_nmr[res_b]
                c_a,c_b=c2[i_a,i_a],c2[i_b,i_b]
                if(c_a!=0. and c_b!=0. and ((c_a>21.0 and c_b>21.0) or (c_a<20. and c_b<20.))):
                    indices.append(i)
            
        write_file.close()
        compare_file.close()
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
                cross_matrix[i2,i1]=reflected[res1][res2]
                reflected_matrix[i1,i2]=reflected[res1][res2]
                d1_matrix[i1,i2]=diag_1[res1][res2]
                d2_matrix[i1,i2]=diag_2[res1][res2]
                c1_matrix[i1,i2],c2_matrix[i1,i2]=c1_shift[res1][res2],c2_shift[res1][res2]
                h1_matrix[i1,i2],h2_matrix[i1,i2]=h1_shift[res1][res2],h2_shift[res1][res2]
                if(cross_matrix[i1,i2]!=0 and reflected_matrix[i1,i2]!=0. and d1_matrix[i1,i2]!=0. and d2_matrix[i1,i2]!=0. and c1_matrix[i1,i2]!=0. and c2_matrix[i1,i2]!=0. and h1_matrix[i1,i2]!=0. and h2_matrix[i1,i2]!=0.):
                    mask_matrix[i1,i2],mask_matrix[i2,i1]=True,True
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
#         print 'METHYL LIST = \n',methyl_list
        return methyl_list

    def methyl_mask_generator(self,residue_list,mask):
        bool_mask=numpy.zeros(mask.shape,dtype=bool)
        indices_VL=[]
        for residue in residue_list:
            if(residue[0:3]=='VAL' or residue[0:3]=='LEU'):
                a,b=residue+' A',residue+' B'
                indices_VL.append(self.map_key_xrd[a])
                indices_VL.append(self.map_key_xrd[b])
                bool_mask[self.map_key_xrd[a],self.map_key_xrd[b]]=True
                bool_mask[self.map_key_xrd[b],self.map_key_xrd[a]]=True
        ii,jj=numpy.ones(mask.shape),numpy.ones(mask.shape)
        indices_VL=numpy.array(indices_VL)
        ii[indices_VL,:],jj[:,indices_VL]=0.5,0.5
        mask=ii*jj
        mask[:,:]=1.
        mask[bool_mask]=0.
        return mask,bool_mask


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



    def DoTestSystem(self,macro=True, ax=True, xx=True, cx=True, ca=True, sparse=True, sym=True, irreducible=False, num=False):
        #First set up the spin system and then cross-correlate:
        print 'Calculating methyl-methyl rates'
        self.cross=CommAux('CH3CH3', latex_file='output/eps/cross_latex.tex', sparse=sparse, sym=sym)
        print_log=open('output/log/log_Comm','w')
        if(ax):
            self.cross.GetDip('H2', 'C1', 'd1', tag ='single')
            self.cross.GetDip('H3', 'C1', 'd1', tag ='single')
            self.cross.GetDip('H4', 'C1', 'd1', tag ='single')
            self.cross.GetDip('H6', 'C5', 'd2', tag ='single')
            self.cross.GetDip('H7', 'C5', 'd2', tag ='single')
            self.cross.GetDip('H8', 'C5', 'd2', tag ='single')
        if(xx):
            self.cross.GetDip('H2', 'H3', 'e1', tag ='single')
            self.cross.GetDip('H2', 'H4', 'e1', tag ='single')
            self.cross.GetDip('H3', 'H4', 'e1', tag ='single')
            self.cross.GetDip('H6', 'H7', 'e2', tag ='single')
            self.cross.GetDip('H6', 'H8', 'e2', tag ='single')
            self.cross.GetDip('H7', 'H8', 'e2', tag ='single')
        if(cx):
            self.cross.GetCSA('H2', 'cX1', tag='single')
            self.cross.GetCSA('H3', 'cX1', tag='single')
            self.cross.GetCSA('H4', 'cX1', tag='single')
            self.cross.GetCSA('H6', 'cX2', tag='single')
            self.cross.GetCSA('H7', 'cX2', tag='single')
            self.cross.GetCSA('H8', 'cX2', tag='single')
        if(ca):
            self.cross.GetCSA('C1', 'cA1', tag='static')
            self.cross.GetCSA('C5', 'cA2', tag='static')

        index = numpy.arange(3)
        for i in (2,3,4):
            i_str = 'H'+str(i)
            for j in (6,7,8):
                j_str = 'H'+str(j)
                self.cross.GetDip(i_str, j_str, 'm', tag='double')

        self.cross.macro=macro
        self.cross.CrossCorrHam()

        self.cross.tm=10E-12
        self.cross.tc=10E-9
        t_point = 400*1E-3
        self.cross.SetFreq(600)  #MHz. Set spectrometer proton frequency
        self.cross.Assemble(('H2z','H3z','H4z'),'M1z')
        self.cross.Assemble(('H6z','H7z','H8z'),'M2z')

        cA=20.
        cX=1.

        #Iterate over all of the pairings:
        #So each methyl in "methyls" will consists of 2 3-vectors - one giving position and one giving orientation
        #These three variables hold - r_single holds auto-relaxation of methyl i on it's own, r_auto holds auto-relaxation of i with j present; r_cross holds cross-relaxation from j to i.


        num=2
        location = numpy.zeros((num, 3))
        orientation = numpy.zeros((num, 3))
        sites = numpy.zeros((num, 3, 3))

    #      methyl = {'loc':c_m, 'orientation':v_m, 'sites':numpy.transpose(numpy.array(sites))}
    #     def CalcMetGeom(c_m, c_a, c_b, r_m=1.09, theta_m=numpy.arccos(1./3.)):
    #     m1,m2=CalcMetGeom([0.,0.,1.],[0.,0.,0.],[1.,0.,0.]),CalcMetGeom([0.,5.,1.],[0.,5.,0.],[1.,5.,0.])
    #     m1,m2=CalcMetGeom([0.,0.,1.],[0.,0.,0.],[1.,0.,0.]),CalcMetGeom([0.,0.,6.],[0.,0.,5.],[1.,0.,5.])
        m1,m2=self.CalcMetGeom([0.,0.,1.],[0.,0.,0.],[1.,0.,0.]),self.CalcMetGeom([5.,0.,1.],[5.,0.,0.],[4.,0.,0.])
        location[0,:],location[1,:]=m1['loc'],m2['loc']
        orientation[0,:],orientation[1,:]=m1['orientation'],m2['orientation']
        sites[0,:,:],sites[1,:,:]=m1['sites'],m2['sites']

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
        self.cross.pars['c2_site']=c_m2
        self.cross.pars['c1_site']=c_m1
        self.cross.pars['sites1']=sites1
        self.cross.pars['v_m1']=v_m1
        self.cross.pars['v_mm']=v_mm
        self.cross.pars['cA']=cA
        self.cross.pars['cX']=cX
        self.cross.pars['v_m2']=v_m2
        self.cross.pars['sites2']=sites2
        self.cross.SetPars('CH3CH3')

        r_cross = self.cross.CalcRate('M1z', 'M2z')
        r_auto = self.cross.CalcRate('M1z', 'M1z')
        rate_matrix=numpy.zeros((num,num))
        i,j=numpy.arange(num),numpy.arange(num)
        ii,jj=numpy.meshgrid(i,j,indexing='ij')
        rate_matrix[ii,jj]=r_cross[ii,jj]
        rate_matrix[i,j]=numpy.sum((r_auto[ii,jj]-r_auto[ii,ii]),axis=1)+r_auto[i,j]
        exp_matrix=expm(-rate_matrix[:,:]*t_point)

        #DEBUGGING:
        numpy.savetxt('output/log/test_cross',r_cross)
        numpy.savetxt('output/log/test_auto',r_auto)
        numpy.savetxt('output/log/test_exp',exp_matrix)

        
mr_noesy('input/correlate','input/ein_noes','input/1EZA.pdb')