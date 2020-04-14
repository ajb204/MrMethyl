#!/usr/bin/python


import basis,string
import sys
from math import sqrt,fabs
import numpy

from fractions import Fraction

from scipy import mat,zeros

from basis import comm,GetMat,HermConj,compare,IsZero,Trace

#############################################################################

def ElementZero(test):
    tat=0
    for i in range(len(test)):
        for j in range(len(test)):
            if(fabs(test[i,j].real)>0 or fabs(test[i,j]).imag>0):
                tat+=1
    if(tat==0):
        print 'Matrix is zero'
    else:
        print 'Matrix is NOT zero'




def MatConCat(A,B):
    for i in range(len(B)):
        A.append(B[i])

def DoubleComm(b,a,target,basis,label,verb):

    if(verb=='y'):
        opo2=a #inner guy
        opo1=b #outer guy
        print '           testing: [',opo1,'[',opo2,',',target,']]'  #the required double commutator
    test=GetMat(target,basis,label)   #get matrix for test operator

    one=GetMat(a,basis,label) #get the inner commutator matrix
    two=GetMat(b,basis,label) #get the outer commutator matrix

    check1=comm(one,test) #compute the double commutator (inner)
    verb='n' #speeds things up massively! use if need to get identity of commutator

    if(verb=='y'):
        fac,resy,baso=compare(check1,basis,label)
    if(IsZero(check1)==0):
        if(verb=='y'):
            print '           First commutator is zero'
        fac=0;resy='zero';retmat=check1
    else:
        if(verb=='y'):
            if(resy!='none'):
                print    '           First commutator [',opo2,',',target,'] is ',fac,resy
            else:
                print    '           First commutator [',opo2,',',target,'] gives undefined result'
            print '           Performing second commutator ',fac,',[',opo1,',',resy,'] ...'
        check2=comm(two,check1) #perform second commutator (outer)
        if(IsZero(check2)>0): #if the double commutator is not zero...
            if(verb=='y'):
                fac,resy,baso=compare(check2,basis,label)
                print "           Double commutator:",fac,resy
                if(resy=='none'):
                    print 'BASIS NOT FOUND:'
                    print
                    print check2
                    print
                    print test
                    print
            else:
                fac=1;resy=''
        else:
            if(verb=='y'):
                print '           Double commutator is zero'
            fac=0;resy='zero';
        retmat=check2

    return fac,resy,retmat




#work out frequency of operator from its name
def GetFreq(coh):
    Isum=0
    Ssum=0
    if(len(string.split(coh,'I1p'))>1):
        Isum+=1
    if(len(string.split(coh,'I1m'))>1):
        Isum-=1
    if(len(string.split(coh,'I2p'))>1):
        Isum+=1
    if(len(string.split(coh,'I2m'))>1):
        Isum-=1
    if(len(string.split(coh,'I3p'))>1):
        Isum+=1
    if(len(string.split(coh,'I3m'))>1):
        Isum-=1


    if(len(string.split(coh,'S1p'))>1):
        Ssum+=1
    if(len(string.split(coh,'S1m'))>1):
        Ssum-=1
    
    freqstr=''
    if(Isum==0 and Ssum==0):
        freqstr+=str(0)
    if(Isum==1):
        freqstr+='wI'
    if(Isum==-1):
        freqstr+='-wI'        
    if(Ssum==1):
        if(len(freqstr)==0):
            freqstr+='wS'
        else:
            freqstr+='+wS'
    if(Ssum==-1):
        freqstr+='-wS'        

    #correct for symmetry of spectral density function
    if(freqstr=='-wI'):
        freqstr='wI'
    if(freqstr=='-wS'):
        freqstr='wS'
    if(freqstr=='-wI-wS'):
        freqstr='wI+wS'

    if(freqstr=='-wI+wS'):
        freqstr='wI-wS'

    q=Isum+Ssum
    return freqstr,q


#evaluate trigonometric integral
def GetTrig(trig1,trig2):
    #Get the angular integrals

    #dipolar integral:self terms
    if(trig1=='(1 - 3cos^2t)' and trig2=='(1 - 3cos^2t)'):
        Ffac=Fraction(4,5)
    elif(trig1=='sint cost exp(+p)' and trig2=='sint cost exp(+p)'):
        Ffac=Fraction(3,10)
    elif(trig1=='sin^2t exp(+2p)' and trig2=='sin^2t exp(+2p)'):
        Ffac=Fraction(3,10)

    #CSA integral:self term
    elif(trig1=='sin2t exp(+p)' and trig2=='sin2t exp(+p)'):
        Ffac=Fraction(8,15)


    #dipolar/dipolar cross terms
    elif(trig1=='(1 - 3cos^2t)' and trig2=='sint cost exp(+p)'):
        Ffac=0.
    elif(trig1=='sint cost exp(+p)' and trig2=='sin^2t exp(+2p)'):
        Ffac=0.
    elif(trig1=='sin^2t exp(+2p)' and trig2=='sint cost exp(+p)'):
        Ffac=0.


    #CSA/dipolar non-zero cross terms
    elif(trig1=='sin2t exp(+p)' and trig2=='sint cost exp(+p)'):
        Ffac=Fraction(2,5)
    elif(trig1=='sint cost exp(+p)' and trig2=='sin2t exp(+p)'):
        Ffac=Fraction(2,5)

    elif(trig1=='sin^2t exp(+2p)' and trig2=='(1 - 3cos^2t)'):
        Ffac=Fraction(1,5)
    elif(trig1=='(1 - 3cos^2t)' and trig2=='sin^2t exp(+2p)'):
        Ffac=Fraction(1,5)


    #CSA/dipolar zero cross terms
    elif(trig1=='(1 - 3cos^2t)' and trig2=='sin2t exp(+p)'):
        Ffac=0.
    elif(trig1=='sint cost exp(+p)' and trig2=='(1 - 3cos^2t)'):
        Ffac=0.
    elif(trig1=='sin^2t exp(+2p)' and trig2=='sin2t exp(+p)'):
        Ffac=0.
    else:
        print 'Cannot evaluate intrgral:',trig1,trig2
        sys.exit(100)

    Ffac=Ffac*5.  #factor to make J(w)=2/5 tc/(1+tc^2w^2) (inserted into the trig formula)
    return Ffac #return in the 2/5 J form
 

#add to array if test is a new entry
def addifnew(test,array):
    tick=0
    for i in range(len(array)):
        if(array[i]==test):
            tick=1
    if(tick==0):
        array.append(test)
    return array


def FormatSpec(face,nfac,j1):
    if(face>0):
        stry=' +'
    else:
        stry=' -'
    if(face<0):
        face=face*-1
    if(face.denominator==1):
        face=face.numerator
        stry+= '('+str(face.numerator)+') '+nfac+' J('+j1+')'
    else:
        stry+= '('+str(face.numerator)+'/'+str(face.denominator)+') '+nfac+' J('+j1+')'
    return stry


#Analyse Hamiltonian and get Hermitian Conjugates of operators
def AnalHamiltonian(ham):
    hamNew=[]
    for i in range(len(ham)):
        #print i+1,'of',len(ham)
        ref=ham[i][1]
        #conj=compare(HermConj(GetMat(HermConj(ref),basis,label)),basis,label)[1]
        n=''
        for j in range(len(ref)):
            if(ref[j]=='m'):
                n+='p'
            elif(ref[j]=='p'):
                n+='m'
            else:
                n+=ref[j]
        conj=n
        #print ref,conj,n
        hamNew.append((ham[i][0],conj,ham[i][2],ham[i][3]))
    return hamNew

def EvalDoubleComm(Pre,Post,Adipolar,AdipolarA,verb,tag=''):
    res=''   #to store string
    freq=[]  #frequnecy of term i
    coeff=[] #coefficient of term i
    nfac=[]
    for i in range(len(Adipolar)):#for each operator in the hamiltonian

        j1,q=GetFreq(Adipolar[i][1]) #from dipolar operator, work out evolution frequency, and trigonometric integral
        j1b,qb=GetFreq(AdipolarA[i][1]) #from dipolar operator, work out evolution frequency, and trigonometric integral
        Ffac=GetTrig(Adipolar[i][3],AdipolarA[i][3]) #from dipolar operator, work out evolution frequency, and trigonometric integral
        
        if(verb=='y'):
            print i+1,': Evaluating operator',Adipolar[i][1],'at frequency',j1,'coherence order:',q
            print '    with operator',AdipolarA[i][1],'at frequency',j1b,'coherence order:',qb
            print '           Trigonometric integral',Adipolar[i][3],AdipolarA[i][3],':',Ffac

        if(Ffac!=0 and q+qb==0): #proceed if trigonometric integral is non-zero, and apply secular approximation
            fac,resy,baso=DoubleComm(AdipolarA[i][1],Adipolar[i][1],Pre,basis,label,verb) #evaulate [hermconj(A),[A,Pre]]
            if(resy!='zero'):
                prefac=Adipolar[i][0]*AdipolarA[i][0] #take the product of the hamiltonian prefactors...
                tr=Trace(GetMat(Post,basis,label),baso)                      #finish the expression
                trD=Trace(GetMat(Post,basis,label),GetMat(Post,basis,label)) #normalisation
                if(verb=='y'):
                    print '           Prefactor:',Adipolar[i][0]*AdipolarA[i][0]
                    print '           [A*,[A,',Pre,']:',(prefac*fac*Ffac),resy,'J(',j1,')' 
                    print '           Calculating '+str(fac)+'<',Post,'|',resy,'> = ',tr
                    print '           Calculating <',Post,'|',Post,'> = ',trD
                if(fabs(tr)>1E-6):#if trace is non zero...
                    face=Fraction((tr*prefac*Ffac/trD/2*2/5)).limit_denominator()  #the numerical weighting of the term (factor of 2 is from the definition, 2/5 is for form of correlation function)
                    coeff.append((face)) #store coefficient
                    if(verb=='y'):
                        
                        print '           FINAL: (',str(face.numerator),'/',str(face.denominator),') '+Adipolar[i][2]+AdipolarA[i][2]+' J(',j1,')'
                    freq.append(j1)      #store frequency

                    if(Adipolar[i][2]==AdipolarA[i][2]):
                        nfac.append(Adipolar[i][2]+'^2')
                    else:
                        nfac.append((Adipolar[i][2]+AdipolarA[i][2]))
                    res+=FormatSpec(face,nfac[len(nfac)-1],j1)
                else:
                    if(verb=='y'):
                        print '          (trace is zero)'
        else:
            if(verb=='y'):
                print 'Trigonometric integral is zero:',Adipolar[i][3],AdipolarA[i][3]
                print 'or non-secular: sum of coherence:',q+qb
    
    #clean up the expression for printing to screen when finished
    if(len(freq)!=0):
        if(verb=='y'):
            print 'Expression for',Post,'on',Pre,':',res
        sum=[]    #sort out and consolidtate terms
        for i in range(len(freq)):
            sum=addifnew(freq[i],sum)
        fin=''
        for i in range(len(sum)):
            val=0
            for j in range(len(freq)):
                if(freq[j]==sum[i]):
                    val+=coeff[j]
            fin+=FormatSpec(Fraction(val).limit_denominator(),nfac[i],sum[i])
        print 'Simplified for ',Post,'on',Pre,':',fin,'(',tag,')'
    


#produce Hamiltonian lists for guys to cross correlate
#note list includes one on two, and two on one.
def CrossCorrHam(one,two,oneA,twoA,verb='n'):
    if(verb=='y'):
        print 'Testing Hamiltonian inteferrence relaxation'
    Aone=[]
    Atwo=[]

    for j in range(len(one)):
        if(verb=='y'):
            print
            print j+1,'of ',len(one),':Testing:',one[j][1]
            print
        for ii in range(len(two)):
            Aone.append(one[j])
            Atwo.append(two[ii])
            #Aone.append(oneA[j])
            #Atwo.append(twoA[ii])
            #Aone=addifnew(one[j],Aone)
            #Atwo=addifnew(two[ii],Atwo)
            #Aone=addifnew(oneA[j],Aone)
            #Atwo=addifnew(twoA[ii],Atwo)
    #cnt=0
    #for i in range(len(Aone)):
    #    for j in range(len(Aone)):
    #        if(i!=j):
    #            if(Aone[i]==Aone[j]):
    #                if(Atwo[i]==Atwo[j]):
    #                    print 'Shit!'
    #                    print i,Aone[i],Atwo[i]
    #                    print j,Aone[j],Atwo[j]
    #                    print
    #                    cnt+=1
                        

    #print len(Aone),cnt
    #sys.exit(100)
    return Aone,Atwo



############################################################
# START FROM HERE!

#d is gammaIgammaS / 4pi e0 r^6
#csa is (parallel-permendicular) * gamma * B0)
# J(w) is tc/(1+wc^2tc^2) - Kay type - numerical constants are as stated

baseType='CH'  #specify basis, C,CH,CH2,CH3

basis,label=basis.GetBasis(baseType)  #get basis
        
#############################################################################

#Dipolar Hamiltonian
#1. Hamiltonian (1/4 for zero quantum terms
#2. factor of 1/2 to take operator from 2IaSb from IaSb
#3. factor of 3 takes Allard's definition for d

Adipolar=[]
Adipolar.append((Fraction(1,2),'I1zS1z','d','(1 - 3cos^2t)')) 
Adipolar.append((Fraction(1,4*2),'I1pS1m','d','(1 - 3cos^2t)'))
Adipolar.append((Fraction(1,4*2),'I1mS1p','d','(1 - 3cos^2t)'))
Adipolar.append((Fraction(1,2),'I1zS1p','d','sint cost exp(+p)'))
Adipolar.append((Fraction(1,2),'I1zS1m','d','sint cost exp(+p)'))
Adipolar.append((Fraction(1,2),'I1pS1z','d','sint cost exp(+p)'))
Adipolar.append((Fraction(1,2),'I1mS1z','d','sint cost exp(+p)'))
Adipolar.append((Fraction(1,2),'I1pS1p','d','sin^2t exp(+2p)'))
Adipolar.append((Fraction(1,2),'I1mS1m','d','sin^2t exp(+2p)'))

AcsaI=[] 
AcsaI.append(( (sqrt(5)/3),'I1z','cI','(1 - 3cos^2t)')) 
AcsaI.append(( (sqrt(5)/4),'I1p','cI','sin2t exp(+p)'))
AcsaI.append(( (sqrt(5)/4),'I1m','cI','sin2t exp(+p)'))

AcsaS=[] 
AcsaS.append(( (sqrt(5)/3),'S1z','cS','(1 - 3cos^2t)')) 
AcsaS.append(( (sqrt(5)/4),'S1p','cS','sin2t exp(+p)'))
AcsaS.append(( (sqrt(5)/4),'S1m','cS','sin2t exp(+p)'))



AcsaI=[] 
AcsaI.append(( (1/3.),'I1z','cI','(1 - 3cos^2t)')) 
AcsaI.append(( (1/4.),'I1p','cI','sin2t exp(+p)'))
AcsaI.append(( (1/4.),'I1m','cI','sin2t exp(+p)'))

AcsaS=[] 
AcsaS.append(( (1/3.),'S1z','cS','(1 - 3cos^2t)')) 
AcsaS.append(( (1/4.),'S1p','cS','sin2t exp(+p)'))
AcsaS.append(( (1/4.),'S1m','cS','sin2t exp(+p)'))



AdipolarA = AnalHamiltonian(Adipolar) #get hermitian conjugates of hamiltonian operators
AcsaIA    = AnalHamiltonian(AcsaI)    #get hermitian conjugates of hamiltonian operators
AcsaSA    = AnalHamiltonian(AcsaS)    #get hermitian conjugates of hamiltonian operators
AdipcsaI,AdipcsaIA=CrossCorrHam(Adipolar,AcsaIA,AdipolarA,AcsaI)
AdipcsaS,AdipcsaSA=CrossCorrHam(Adipolar,AcsaSA,AdipolarA,AcsaS)
AcsaIcsaS,AcsaIcsaSA=CrossCorrHam(AcsaI,AcsaSA,AcsaIA,AcsaS)


verb='n'

print 'Comparisons to Keeler'
EvalDoubleComm('I1z','I1z',Adipolar,AdipolarA,verb,'Iz/Iz dipolar (eq19)')  #correct
EvalDoubleComm('S1z','S1z',Adipolar,AdipolarA,verb,'Sz/Sz dipolar (eq19)')  #correct
EvalDoubleComm('I1z','S1z',Adipolar,AdipolarA,verb,'Sz/Iz cross dipolar (eq19)')  #correct
EvalDoubleComm('I1zS1z','I1zS1z',Adipolar,AdipolarA,verb,'IzSz dipolar (page 6.23)')  #correct
EvalDoubleComm('I1x','I1x',Adipolar,AdipolarA,verb,'Ix dipolar (page 6.24)')  #correct
EvalDoubleComm('I1pS1z','I1pS1z',Adipolar,AdipolarA,verb,'I+Sz dipolar (page 6.27)')  #correct
EvalDoubleComm('I1pS1p','I1pS1p',Adipolar,AdipolarA,verb,'I+S+ dipolar (page 6.28)')  #correct 
EvalDoubleComm('I1pS1m','I1pS1m',Adipolar,AdipolarA,verb,'I+S- dipolar (page 6.28)')  #correct 
EvalDoubleComm('I1z','I1z',AcsaI,AcsaIA,verb,'Iz csa (page 6.30)')  #correct 
EvalDoubleComm('I1x','I1x',AcsaI,AcsaIA,verb,'Ix csa (page 6.30)')  #correct 
EvalDoubleComm('I1z','I1zS1z',AdipcsaI,AdipcsaIA,'n','dipolar/CSA_I')  #correct
EvalDoubleComm('I1pS1z','I1p',AdipcsaI,AdipcsaIA,'n','dipolar/CSA_I (page 6.33)')  #correct - note small error on 6.33, not on next page

print
EvalDoubleComm('I1xS1z','I1x',AdipcsaI,AdipcsaIA,'n','dipolar/CSA_I (page 6.33)')  #correct 
EvalDoubleComm('S1x','I1zS1x',AdipcsaS,AdipcsaSA,'n','dipolar/CSA_I (page 6.33)')  #correct 
EvalDoubleComm('I1xS1x','I1yS1y',Adipolar,AdipolarA,'n','dipolar MQ transfer')  #correct 
EvalDoubleComm('I1yS1x','I1xS1y',Adipolar,AdipolarA,'n','dipolar MQ transfer')  #correct 

print
EvalDoubleComm('I1aS1x','I1aS1x',AdipcsaS,AdipcsaSA,'n','IaSx diploar/csaS (AntiTrosy)')  #correct 
EvalDoubleComm('I1aS1x','I1bS1x',AdipcsaS,AdipcsaSA,'n','IaSxIbSx diploar/csaS (AntiTrosy)')  #will be zero as orthoganol
EvalDoubleComm('I1bS1x','I1bS1x',AdipcsaS,AdipcsaSA,'n','IbSx dipolar csaS (Trosy)')  #correct 

print
EvalDoubleComm('I1aS1x','I1aS1x',AdipcsaS,AdipcsaSA,'n','IaSx diploar/csaS (AntiTrosy)')  #correct 
EvalDoubleComm('I1bS1x','I1bS1x',AdipcsaS,AdipcsaSA,'n','IbSx dipolar csaS (Trosy)')  #correct
EvalDoubleComm('I1aS1x','I1aS1x',Adipolar,AdipolarA,'n','IaSx dipolar (AntiTrosy)')  #correct 
EvalDoubleComm('I1bS1x','I1bS1x',Adipolar,AdipolarA,'n','IbSx dipolar (Trosy)')  #correct 

print 'end'
EvalDoubleComm('I1xS1a','I1x',AcsaI,AcsaIA,'n','CSA_I')  #correct 
EvalDoubleComm('I1xS1b','I1x',AcsaI,AcsaIA,'n','CSA_I')  #correct 
EvalDoubleComm('I1xS1a','I1x',Adipolar,AdipolarA,'n','CSA_I')  #correct 
EvalDoubleComm('I1xS1b','I1x',Adipolar,AdipolarA,'n','CSA_I')  #correct 




Anew=[]
MatConCat(Anew,Adipolar)
#MatConCat(Anew,AcsaS)
#MatConCat(Anew,AcsaI)

AnewA=[]
MatConCat(AnewA,AdipolarA)
#MatConCat(AnewA,AcsaSA)
#MatConCat(AnewA,AcsaIA)

#a test of the new cross correlations...
Anewnew,AnewnewA=CrossCorrHam(Anew,AnewA,AnewA,Anew)



EvalDoubleComm('I1bS1x','I1bS1x',Anewnew,AnewnewA,'n','IbSx dipolar (Trosy)')  #correct 
EvalDoubleComm('I1bS1x','I1bS1x',Adipolar,AdipolarA,'n','IbSx dipolar (Trosy)')  #correct 


#for k in ('I1p','I1x','I1y','S1x','S1z','I1zS1z','I1xS1z','I1xS1x','I1pS1m','I1pS1p'):



sys.exit(100)


#################################
print 'doing calc for CH3'

import basis
baseType='CH3'  #specify basis, C,CH,CH2,CH3

basis,label=basis.GetBasis(baseType)  #get basis


Adipolar1=[]
Adipolar1.append((Fraction(1,2),'I1zS1z','d','(1 - 3cos^2t)')) 
Adipolar1.append((Fraction(1,4*2),'I1pS1m','d','(1 - 3cos^2t)'))
Adipolar1.append((Fraction(1,4*2),'I1mS1p','d','(1 - 3cos^2t)'))
Adipolar1.append((Fraction(1,2),'I1zS1p','d','sint cost exp(+p)'))
Adipolar1.append((Fraction(1,2),'I1zS1m','d','sint cost exp(+p)'))
Adipolar1.append((Fraction(1,2),'I1pS1z','d','sint cost exp(+p)'))
Adipolar1.append((Fraction(1,2),'I1mS1z','d','sint cost exp(+p)'))
Adipolar1.append((Fraction(1,2),'I1pS1p','d','sin^2t exp(+2p)'))
Adipolar1.append((Fraction(1,2),'I1mS1m','d','sin^2t exp(+2p)'))

Adipolar2=[]
Adipolar2.append((Fraction(1,2),'I2zS1z','d','(1 - 3cos^2t)')) 
Adipolar2.append((Fraction(1,4*2),'I2pS1m','d','(1 - 3cos^2t)'))
Adipolar2.append((Fraction(1,4*2),'I2mS1p','d','(1 - 3cos^2t)'))
Adipolar2.append((Fraction(1,2),'I2zS1p','d','sint cost exp(+p)'))
Adipolar2.append((Fraction(1,2),'I2zS1m','d','sint cost exp(+p)'))
Adipolar2.append((Fraction(1,2),'I2pS1z','d','sint cost exp(+p)'))
Adipolar2.append((Fraction(1,2),'I2mS1z','d','sint cost exp(+p)'))
Adipolar2.append((Fraction(1,2),'I2pS1p','d','sin^2t exp(+2p)'))
Adipolar2.append((Fraction(1,2),'I2mS1m','d','sin^2t exp(+2p)'))

Adipolar3=[]
Adipolar3.append((Fraction(1,2),'I3zS1z','d','(1 - 3cos^2t)')) 
Adipolar3.append((Fraction(1,4*2),'I3pS1m','d','(1 - 3cos^2t)'))
Adipolar3.append((Fraction(1,4*2),'I3mS1p','d','(1 - 3cos^2t)'))
Adipolar3.append((Fraction(1,2),'I3zS1p','d','sint cost exp(+p)'))
Adipolar3.append((Fraction(1,2),'I3zS1m','d','sint cost exp(+p)'))
Adipolar3.append((Fraction(1,2),'I3pS1z','d','sint cost exp(+p)'))
Adipolar3.append((Fraction(1,2),'I3mS1z','d','sint cost exp(+p)'))
Adipolar3.append((Fraction(1,2),'I3pS1p','d','sin^2t exp(+2p)'))
Adipolar3.append((Fraction(1,2),'I3mS1m','d','sin^2t exp(+2p)'))

Adipolar1A = AnalHamiltonian(Adipolar1) #get hermitian conjugates of hamiltonian operators
Adipolar2A = AnalHamiltonian(Adipolar2) #get hermitian conjugates of hamiltonian operators
Adipolar3A = AnalHamiltonian(Adipolar3) #get hermitian conjugates of hamiltonian operators]



Adipolar23=[]
MatConCat(Adipolar23,Adipolar2)
MatConCat(Adipolar23,Adipolar3)

Adipolar12=[]
MatConCat(Adipolar12,Adipolar1)
MatConCat(Adipolar12,Adipolar2)

Adipolar13=[]
MatConCat(Adipolar13,Adipolar1)
MatConCat(Adipolar13,Adipolar3)


Adipolar12A = AnalHamiltonian(Adipolar12) #get hermitian conjugates of hamiltonian operators
Adipolar23A = AnalHamiltonian(Adipolar23) #get hermitian conjugates of hamiltonian operators
Adipolar13A = AnalHamiltonian(Adipolar13) #get hermitian conjugates of hamiltonian operators

Adipdip1_23,Adipdip1_23A=CrossCorrHam(Adipolar1,Adipolar23A,Adipolar1A,Adipolar23)
Adipdip2_13,Adipdip2_13A=CrossCorrHam(Adipolar2,Adipolar13A,Adipolar2A,Adipolar13)
Adipdip3_12,Adipdip3_12A=CrossCorrHam(Adipolar3,Adipolar12A,Adipolar3A,Adipolar12)


Adipolar=[]
AdipolarA=[]
MatConCat(Adipolar,Adipolar1)
MatConCat(Adipolar,Adipolar2)
MatConCat(Adipolar,Adipolar3)
MatConCat(Adipolar,Adipdip1_23)
MatConCat(Adipolar,Adipdip2_13)
MatConCat(Adipolar,Adipdip3_12)

MatConCat(AdipolarA,Adipolar1A)
MatConCat(AdipolarA,Adipolar2A)
MatConCat(AdipolarA,Adipolar3A)
MatConCat(AdipolarA,Adipdip1_23A)
MatConCat(AdipolarA,Adipdip2_13A)
MatConCat(AdipolarA,Adipdip3_12A)


Adipdip,Adipdip=CrossCorrHam(Adipolar,AdipolarA,AdipolarA,Adipolar) #proceed to cross correlate



AcsaI=[] 
AcsaI.append(( (sqrt(5)/3),'I1z','cI','(1 - 3cos^2t)')) 
AcsaI.append(( (sqrt(5)/4),'I1p','cI','sin2t exp(+p)'))
AcsaI.append(( (sqrt(5)/4),'I1m','cI','sin2t exp(+p)'))

AcsaS=[] 
AcsaS.append(( (sqrt(5)/3),'S1z','cS','(1 - 3cos^2t)')) 
AcsaS.append(( (sqrt(5)/4),'S1p','cS','sin2t exp(+p)'))
AcsaS.append(( (sqrt(5)/4),'S1m','cS','sin2t exp(+p)'))



AcsaI=[] 
AcsaI.append(( (1/3.),'I1z','cI','(1 - 3cos^2t)')) 
AcsaI.append(( (1/4.),'I1p','cI','sin2t exp(+p)'))
AcsaI.append(( (1/4.),'I1m','cI','sin2t exp(+p)'))

AcsaS=[] 
AcsaS.append(( (1/3.),'S1z','cS','(1 - 3cos^2t)')) 
AcsaS.append(( (1/4.),'S1p','cS','sin2t exp(+p)'))
AcsaS.append(( (1/4.),'S1m','cS','sin2t exp(+p)'))

#print label
#print len(label),basis[0],basis[0].shape

print len(Adipolar)

print 'Getting Hermitian Conjugates'
#AdipolarA = AnalHamiltonian(Adipolar) #get hermitian conjugates of hamiltonian operators
AcsaIA    = AnalHamiltonian(AcsaI)    #get hermitian conjugates of hamiltonian operators
AcsaSA    = AnalHamiltonian(AcsaS)    #get hermitian conjugates of hamiltonian operators
AdipcsaI,AdipcsaIA=CrossCorrHam(Adipolar,AcsaIA,AdipolarA,AcsaI)

AdipNew=[]
AdipNewA=[]

MatConCat(AdipNew,Adipolar1)
MatConCat(AdipNew,Adipolar2)
MatConCat(AdipNew,Adipolar3)

MatConCat(AdipNewA,Adipolar1A)
MatConCat(AdipNewA,Adipolar2A)
MatConCat(AdipNewA,Adipolar3A)


AdipcsaS,AdipcsaSA=CrossCorrHam(AdipNew,AcsaSA,AdipNewA,AcsaS)




AcsaIcsaS,AcsaIcsaSA=CrossCorrHam(AcsaI,AcsaSA,AcsaIA,AcsaS)

#Adipdip,AdipdipA=CrossCorrHam(Adipolar,AdipolarA,AdipolarA,Adipolar)



EvalDoubleComm('I3aI2aI1aS1p','I3aI2aI1aS1p',Adipolar,AdipolarA,'n','C+HaHaHa dipolar L')  #correct
EvalDoubleComm('I3aI2aI1aS1p','I3aI2aI1aS1p',AcsaS,AcsaSA,'n','C+HaHaHa csa')  #correct
EvalDoubleComm('I3aI2aI1aS1p','I3aI2aI1aS1p',AdipcsaS,AdipcsaSA,'n','C+HaHaHa csa/dipolar')  #correct

#EvalDoubleComm('I3aI2aI1aS1p','I3aI2aI1aS1p',Adipdip,AdipdipA,'n','C+HaHaHa dipolar L')  #correct

EvalDoubleComm('I3bI2aI1aS1p','I3bI2aI1aS1p',Adipolar,AdipolarA,'n','C+HaHaHa dipolar L')  #correct
EvalDoubleComm('I3bI2aI1aS1p','I3bI2aI1aS1p',AcsaS,AcsaSA,'n','C+HaHaHa dipolar L')  #correct
EvalDoubleComm('I3bI2aI1aS1p','I3bI2aI1aS1p',AdipcsaS,AdipcsaSA,'n','C+HaHaHa dipolar L')  #correct

EvalDoubleComm('I3bI2bI1aS1p','I3bI2bI1aS1p',Adipolar,AdipolarA,'n','C+HaHaHa dipolar L')  #correct
EvalDoubleComm('I3bI2bI1aS1p','I3bI2bI1aS1p',AcsaS,AcsaSA,'n','C+HaHaHa dipolar L')  #correct
EvalDoubleComm('I3bI2bI1aS1p','I3bI2bI1aS1p',AdipcsaS,AdipcsaSA,'n','C+HaHaHa dipolar L')  #correct


EvalDoubleComm('I3bI2bI1bS1p','I3bI2bI1bS1p',Adipolar,AdipolarA,'n','C+HaHaHa dipolar L')  #correct
EvalDoubleComm('I3bI2bI1bS1p','I3bI2bI1bS1p',AcsaS,AcsaSA,'n','C+HaHaHa dipolar L')  #correct
EvalDoubleComm('I3bI2bI1bS1p','I3bI2bI1bS1p',AdipcsaS,AdipcsaSA,'n','C+HaHaHa dipolar L')  #correct


#EvalDoubleComm('I3bI2aI1aS1p','I3bI2aI1aS1p',Adipdip,AdipdipA,'n','C+HaHaHa dipolar L')  #correct

sys.exit(100)


EvalDoubleComm('I3aI2aI1aS1p','I3aI2aI1aS1p',Adipolar,AdipolarA,'n','C+HaHaHa dipolar L')  #correct






print

EvalDoubleComm('I3bI1bS1p','I3bI1bS1p',Adipolar,AdipolarA,'n','C+HbHbHb dipolar L')  #correct
EvalDoubleComm('I3bI1bS1p','I3bI1bS1p',AcsaS,AcsaSA,'n','C+HbHbHb csa')  #correct
EvalDoubleComm('I3bI1bS1p','I3bI1bS1p',AdipcsaS,AdipcsaSA,'n','C+HbHbHb csa/dipolar')  #correct

print

EvalDoubleComm('I1bS1p','I1bS1p',Adipolar,AdipolarA,'y','C+HbHbHb dipolar L')  #correct

EvalDoubleComm('I1bS1p','I1bS1p',AcsaS,AcsaSA,'n','C+HbHbHb csa')  #correct
EvalDoubleComm('I1bS1p','I1bS1p',AdipcsaS,AdipcsaSA,'n','C+HbHbHb csa/dipolar')  #correct

print

EvalDoubleComm('I3bI1bS1p','I3bI1bS1p',Adipolar,AdipolarA,'n','C+HbHbHb dipolar L')  #correct
EvalDoubleComm('I3bI1bS1p','I3bI1bS1p',AcsaS,AcsaSA,'n','C+HbHbHb csa')  #correct
EvalDoubleComm('I3bI2bI1bS1p','I3bI2bI1bS1p',AdipcsaS,AdipcsaSA,'n','C+HbHbHb csa/dipolar')  #correct

print

EvalDoubleComm('I3aI2aI1bS1p','I3aI2aI1bS1p',Adipolar,AdipolarA,'n','C+HaHaHb dipolar')  #correct
EvalDoubleComm('I3aI2aI1bS1p','I3aI2aI1bS1p',AcsaS,AcsaSA,'n','C+HaHaHb csa')  #correct

EvalDoubleComm('I3aI2bI1bS1p','I3aI2bI1bS1p',Adipolar,AdipolarA,'n','C+HaHbHb dipolar')  #correct

EvalDoubleComm('I2bI1bS1p','I2bI1bS1p',Adipolar,AdipolarA,'n','C+HbHb dipolar')  #correct
EvalDoubleComm('I1zS1p','I1zS1p',Adipolar,AdipolarA,'n','C+ dipolar')  #correct
EvalDoubleComm('I2zI1zS1p','I2zI1zS1p',Adipolar,AdipolarA,'n','C+ dipolar')  #correct
EvalDoubleComm('I3zI2zI1zS1p','I3zI2zI1zS1p',Adipolar,AdipolarA,'n','C+ dipolar')  #correct

EvalDoubleComm('I2bS1p','I2bS1p',Adipolar,AdipolarA,'n','C+ dipolar')  #correct
EvalDoubleComm('I3bS1p','I3bS1p',Adipolar,AdipolarA,'n','C+ dipolar')  #correct
EvalDoubleComm('S1p','S1p',Adipolar,AdipolarA,'n','C+ dipolar')  #correct


######################



sys.exit(100)

#NOTE  - spectral density functions
#keeler - Real[int -inf+inf exp(-t/tc) exp(-i w t) dt ] = 2tc/(1+tc^2w^2)
#palmer - Real[int -inf+inf 1/5 exp(-t/tc) exp(-i w t) dt ] = 2/5 tc/(1+tc^2w^2)
#Allard takes same as Palmer

#NOTE - CSA Hamiltonian
#need to check. Allard and Keeler differ by a root(5). Not sure why.
#AX IS COMPLETE! We can recapitulate the results of Allard.



#compute CH3 relaxation


for k in ('I1p','I1x','I1y','S1x','S1z','I1zS1z','I1xS1z','I1xS1x','I1pS1m','I1pS1p'):
    print
    print 'Deriving relaxation rate for',k
    for i in range(len(label)):#try to relax with all operators in basis

        if(len(string.split(label[i],'a'))<2 and len(string.split(label[i],'b'))<2 and len(string.split(label[i],'p'))<2 and len(string.split(label[i],'m'))<2 ):
            EvalDoubleComm(k,label[i],Adipolar,AdipolarA,verb,'dipolar')  #get expression. Same as in Palmer's book
#            print 'Considering chemical shift anisotropy on I'
            EvalDoubleComm(k,label[i],AcsaI,AcsaIA,verb,'CSA_I')  #get expression. Same as in Palmer's book
#            print 'Considering chemical shift anisotropy on S'
            EvalDoubleComm(k,label[i],AcsaS,AcsaSA,verb,'CSA_S')  #get expression. Same as in Palmer's book

 #           print 'Considering cross-correlation between CSA and dipolar:'
            for j in range(len(AcsaIA)):
                Atest=[]
                for ii in range(len(Adipolar)):
                    Atest.append(AcsaIA[j])
                EvalDoubleComm(k,label[i],Adipolar,Atest,verb,'dipolar/CSA_I')  #get expression. Same as in Palmer's book
    print




