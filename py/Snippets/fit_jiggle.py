#!/usr/bin/python

import scipy,numpy,os,sys,copy
from scipy.optimize import leastsq


class fitty():

    def __init__(self):
        print
        print 'I am a fitting program!'
        self.m=70
        self.c=5

        self.bootVals=100
        datsize=10
        self.xvals=numpy.linspace(0,0.1,datsize)
        self.ydat=self.xvals*self.m+self.c + numpy.random.normal(0,1,datsize)
        self.cally=0

        self.m=numpy.random.normal(100,0.1,1)
        self.c=numpy.random.normal(100,0.1,1)

        xinit=self.pack()
        x0=leastsq(self.chifunc,xinit)
        self.unpack(x0[0])
        
        print 'hello!:',self.m
        print 'hello!:',self.c
        print 'calls:',self.cally

        self.jiggler()

        self.booty(x0[0])

        outy=open('outy.out','w')
        for i in range(datsize):
            outy.write('%f\t%f\n' % (self.xvals[i],self.ydat[i]))
        outy.write('\n\n')

        #I am increasing the number of x points to make a nice line.
        self.xvals=numpy.linspace(0,0.1,datsize*10)
        self.calcModel()

        for i in range(datsize*10):
            outy.write('%f\t%f\n' % (self.xvals[i],self.ycalc[i]))
        outy.write('\n\n')
        outy.close()


    def fake(self):
        mask=numpy.random.randint(0,len(self.xvals),len(self.xvals))
        self.xvals=self.xvalsOld[mask]
        self.ydat=self.ydatOld[mask]


    def booty(self,xinit):
        print 'Running the bootstrap:'
        self.xvalsOld=copy.deepcopy(self.xvals)
        self.ydatOld=copy.deepcopy(self.ydat)
        self.xinitBack=copy.deepcopy(xinit)

        booty=[]
        boots=10
        for i in range(self.bootVals):
            #print i
            #make fake datset
            self.fake()
            
            x0=leastsq(self.chifunc,self.xinitBack)
            booty.append(x0[0])
            #fit it
            #make histogram of parameters
        booty=numpy.array(booty)
        
        print ' m:',xinit[0],numpy.average(booty[:,0]),numpy.std(booty[:,0])
        print ' c:',xinit[1],numpy.average(booty[:,1]),numpy.std(booty[:,1])
        
    def pack(self):
        x=[]
        x.append(self.m)
        x.append(self.c)
        return x
    def unpack(self,x):
        self.m=x[0]
        self.c=x[1]

    def GetChi2(self,x):
        self.unpack(x)
        self.calcModel()
        #print 'chi2',numpy.sum((self.ydat-self.ycalc)**2.)        
        return numpy.sum((self.ydat-self.ycalc)**2.)        

    def chifunc(self,x):
        self.unpack(x)
        self.calcModel()
        #print 'chi2',numpy.sum((self.ydat-self.ycalc)**2.)
        return self.ydat-self.ycalc
        
    def calcModel(self):
        self.ycalc=self.xvals*self.m+self.c
        self.cally+=1

    def AddHeat(self,pars):
        print 'pars before heat:',pars
        parNew=pars*numpy.random.normal(1.,self.sigma,len(pars))
        print 'pars after heat: ',parNew
        return parNew
    

    def jiggler(self):
        print 'Unleashing the jiggler...'
        self.jiggles=20  #max number of goes before aborting
        self.sigma=0.1   #standard deviation of normal distribution

        xinit=self.pack()
        x0=leastsq(self.chifunc,xinit)


        currParBest=copy.deepcopy(x0[0])
        chi2Best=self.GetChi2(currParBest)

        cnt=0
        while(1==1):
            parCurr=copy.deepcopy(currParBest)
            parCurr=self.AddHeat(parCurr)
            x0=leastsq(self.chifunc,parCurr)
            chi2Curr=self.GetChi2(x0[0])
            print 'best:',chi2Best,'current:',chi2Curr
            if(chi2Curr<chi2Best):
                print '    yay! we have improved! saving.' 
                self.currParBest=copy.deepcopy(x0[0])
                chi2Best=chi2Curr
                cnt=0
            else:
                print 'no improvment. number of jiggles:',cnt

            cnt+=1
            if(cnt==self.jiggles):
                break
            
    
        x0=leastsq(self.chifunc,currParBest)
        chi2Final=self.GetChi2(x0[0])
        print 'final chi2:',chi2Final
        return x0[0]
inst=fitty()
    


    
