#!/usr/bin/python


import sys,numpy,copy,scipy
from math import sqrt
from scipy import mat,zeros,kron
from scipy.sparse import kron as kronSparse
from fractions import Fraction

#Creates bases and contains functions for performing mathematical operations

#compute the outer product of a and b
def OuterProd(a,b):
    #return kron(b,a)
    return numpy.kron(b,a)
def OuterProdSparse(a,b):
    return kronSparse(b,a)
    #return numpy.kron(b,a)

########################################################################################
# First define the starting basis

def GetSpin1():
    #define cartesian basis for spin 1/2
    E1=mat(zeros((3,3)),dtype=complex)
    E1[0,0]=complex(1.,0)
    E1[1,1]=complex(1.,0)
    E1[2,2]=complex(1.,0)

    Ix1=mat(zeros((3,3)),dtype=complex)
    Ix1[0,1]=complex(1.,0.)/numpy.sqrt(2.)
    Ix1[1,0]=complex(1.,0.)/numpy.sqrt(2.)
    Ix1[1,2]=complex(1.,0.)/numpy.sqrt(2.)
    Ix1[2,1]=complex(1.,0.)/numpy.sqrt(2.)


    Iy1=mat(zeros((3,3)),dtype=complex)
    Iy1[0,1]=complex(0.,-1.)/numpy.sqrt(2.)
    Iy1[1,0]=complex(0.,1.)/numpy.sqrt(2.)
    Iy1[1,2]=complex(0.,-1)/numpy.sqrt(2.)
    Iy1[2,1]=complex(0.,1.)/numpy.sqrt(2.)

    Iz1=mat(zeros((3,3)),dtype=complex)
    Iz1[0,0]=complex(1.,0.)#*numpy.sqrt(2.)
    Iz1[2,2]=complex(-1.,0.)#*numpy.sqrt(2.)

    """
    Ip1=mat(zeros((3,3)),dtype=complex)
    Ip1[0,1]=complex(1.,0.)*numpy.sqrt(2.)
    Ip1[1,2]=complex(1.,0.)*numpy.sqrt(2.)

    Im1=mat(zeros((3,3)),dtype=complex)
    Im1[1,0]=complex(1.,0.)*numpy.sqrt(2.)
    Im1[2,1]=complex(1.,0.)*numpy.sqrt(2.)
    """

    return E1,Ix1,Iy1,Iz1


def GetSpinHalf():
    #define cartesian basis for spin 1/2
    E1=mat(zeros((2,2)),dtype=complex)
    E1[0,0]=complex(1,0)
    E1[1,1]=complex(1,0)

    Ix1=mat(zeros((2,2)),dtype=complex)
    Ix1[0,1]=complex(0.5,0.)
    Ix1[1,0]=complex(0.5,0.)

    Iy1=mat(zeros((2,2)),dtype=complex)
    Iy1[0,1]=complex(0.,-0.5)
    Iy1[1,0]=complex(0.,0.5)

    Iz1=mat(zeros((2,2)),dtype=complex)
    Iz1[0,0]=complex(0.5,0.)
    Iz1[1,1]=complex(-0.5,0.)

    return E1,Ix1,Iy1,Iz1

def Get1SpinBasis(cart_flg,shift_flg,spin=1/2):
    if(spin==1/2):
        E1,Ix1,Iy1,Iz1=GetSpinHalf()
    elif(spin==1):
        E1,Ix1,Iy1,Iz1=GetSpin1()
    else:
        print 'Spin not recognised:',spin
        sys.exit(100)
    
    base={}

    if(cart_flg=='y'):
        base['E']=E1
        base['C1x']=Ix1
        base['C1y']=Iy1
        base['C1z']=Iz1

    if(shift_flg=='y'):
        #define the shift basis
        Ip1=(Ix1+Iy1*complex(0,1))
        Im1=(Ix1-Iy1*complex(0,1))
        Ia1=(E1/2.+Iz1)
        Ib1=(E1/2.-Iz1)

        #basis2 - shift
        base['C1p']=Ip1
        base['C1m']=Im1
        base['C1a']=Ia1
        base['C1b']=Ib1

        if(spin==1):
            #print Iy1
            Iy2=Iy1.copy()
            Iy2[1,2]*=-1
            Iy2[1,0]*=-1
            Ix2=Ix1.copy()
            Ix2[1,2]*=-1
            Ix2[1,0]*=-1

            Ip2=(Ix2+Iy2*complex(0,1))
            Im2=(Ix2-Iy2*complex(0,1))
            
            #It1=zeros_like(Ix1)
            It1=mat(zeros((3,3)),dtype=complex)
            It1[2,0]=1.

            base['C1t']=It1            
            #base['C1u']=Iy2            
            base['C1v']=Ip2            
            base['C1w']=Im2            
            #print Ip2
            #print Im2
            #sys.exit(100)
    return base


#
def MakeLab(spin,ilab,basekey,newbasekey):
    if(basekey=='E'):
        if(newbasekey=='E'):
            return 'E'
    lab=''
    if(basekey!='E'):
        lab+=spin+str(ilab)+basekey[2]
    if(newbasekey!='E'):
        lab+=newbasekey #otherwise add other label
    return lab

#add a new element to basis
def AddNewX(base,Newbase,spin,ilab):
    baseNew={}
    basekeys=base.keys()
    newbasekeys=Newbase.keys()
    for i in range(len(basekeys)):
        for j in range(len(newbasekeys)):
            lab=MakeLab(spin,ilab,basekeys[i],newbasekeys[j])
            baseNew[lab]=OuterProd(base[basekeys[i]],Newbase[newbasekeys[j]])
    return baseNew

#initialise spin 1 and 1/2
def InitBasis(cart_flg='y',shift_flg='y',sparse=True):
    baseHalf=Get1SpinBasis(cart_flg,shift_flg,spin=1/2)
    base1=Get1SpinBasis(cart_flg,shift_flg,spin=1)
    if(sparse):
        from scipy.sparse import csr_matrix
        for key,vals in base1.items():
            base1[key]=csr_matrix(vals)
        for key,vals in baseHalf.items():
            baseHalf[key]=csr_matrix(vals)
    return baseHalf,base1


def GetBasis(tag,cart_flg='y',shift_flg='y',Sparse=False):

    baseHalf=Get1SpinBasis(cart_flg,shift_flg,spin=1/2)
    base1=Get1SpinBasis(cart_flg,shift_flg,spin=1)

    #old tests:
    #print comm(base1['C1x'],base1['C1y'])    
    #print comm(baseHalf['C1x'],baseHalf['C1y'])    
    #print numpy.dot(base1['C1x'],base1['C1x'])
    #print numpy.dot(baseHalf['C1x'],baseHalf['C1x'])

    test='n'
    if(test=='y'):
        CartesianTest1Spin(baseHalf)    
        CartesianTest1Spin(base1)    

    #print base1['C1m']
    #print base1['C1y']
    #print base1['C1z']
    
    if(len(tag)==1):
        if(tag[0]=='D'): #if just want single spin basis
            return base1
        else:
            return baseHalf
    
    Newbase=copy.deepcopy(baseHalf) #start with a spin 1/2.

    #Newbase1=copy.deepcopy(base1)
    if(tag=='C2H2'):
        Newbase=AddNewX(base,Newbase,2,spin='S')  
        Newbase=AddNewX(base,Newbase,1,spin='I')
        Newbase=AddNewX(base,Newbase,2,spin='I')     
        return NewBasis,NewLabel

    try:
        Cint=int(tag.split('C')[1][0])
    except:
        if(len(tag.split('C'))>1):
            Cint=1
        else:
            Cint=0
    try:
        Hint=int(tag.split('H')[1][0])
    except:
        if(len(tag.split('H'))>1):
            Hint=1
        else:
            Hint=0
    try:
        Dint=int(tag.split('D')[1][0])
    except:
        if(len(tag.split('D'))>1):
            Dint=1
        else:
            Dint=0
    sys.stdout.write('Constructing C%iH%iD%i\n' % (Cint,Hint,Dint))

    j=2
    if(Cint>1):
        for i in range(Cint-1):
            Newbase=AddNewX(baseHalf,Newbase,'C',j)        
            j+=1
    for i in range(Hint):
        print 'Adding H...'
        Newbase=AddNewX(baseHalf,Newbase,'H',j)        
        print len(Newbase.keys())
        j+=1
    for i in range(Dint):
        Newbase=AddNewX(base1,Newbase,'D',j)        
        j+=1


    if(test=='y' and j==1):
        CartesianTest2Spin(Newbase)

    print 'done'
    if(Sparse):
        from scipy.sparse import csr_matrix
        for key,vals in Newbase.items():
            Newbase[key]=csr_matrix(vals)




    return Newbase







#################################################################



def CartesianTest1Spin(base):
    print '1Spin:'
    print 'Cartesian commutators:'
    TestComm('C1x','C1y',complex(0,1),'C1z',base)
    TestComm('C1z','C1x',complex(0,1),'C1y',base)
    TestComm('C1y','C1z',complex(0,1),'C1x',base)
    TestComm('C1y','C1x',complex(0,-1),'C1z',base)
    TestComm('C1x','C1z',complex(0,-1),'C1y',base)
    TestComm('C1z','C1y',complex(0,-1),'C1x',base)
    print 'Shift commutators'
    TestComm('C1p','C1z',-1,'C1p',base)
    TestComm('C1z','C1p',1,'C1p',base)
    TestComm('C1m','C1z',1,'C1m',base)
    TestComm('C1z','C1m',-1,'C1m',base)
    TestComm('C1p','C1m',2,'C1z',base)
    TestMult('C1x','C1y',complex(0,0.5),'C1z',base)
    TestMult('C1p','C1m',1,'C1a',base)
    TestMult('C1m','C1p',1,'C1b',base)
    TestMult('C1b','C1m',1,'C1m',base)
    TestMult('C1m','C1b',0,'C1m',base)
    TestMult('C1a','C1p',1,'C1p',base)
    TestMult('C1p','C1a',0,'C1p',base)


def CartesianTest2Spin(base):
    print '2Spin:'
    TestMult('I1p','I1a',0,'I1p',base)
    TestComm('I1zC1z','I1yC1z',complex(0,-1),'I1x',base)
    TestComm('I1zC1z','I1x',complex(0,1),'I1yC1z',base)
    TestComm('I1xC1y','I1zC1z',0,'E',base)
    TestComm('I1zC1z','I1p',1,'I1pC1z',base)
    TestCommSum('I1pC1m','I1mC1p',4,'I1z',-1,'C1z',base)
    TestCommSum('I1pC1p','I1mC1m',4,'I1z',+1,'C1z',base)
    TestComm('I1pC1m','I1z',-1,'I1pC1m',base)
    TestComm('I1pC1z','I1mC1z',2,'I1z',base)
    TestComm('I1mC1m','I1z',1,'I1mC1m',base)


#evaluate [1,2-fac3]==? if zero, commutator is correct. Otherwise, report does not equal
def TestComm(target1,target2,target3fac,target3,base):
    test=comm(base[target1],base[target2] )-target3fac*base[target3]
    #print comm(base[target1],base[target2])+base[target3]*complex(0,1)
    #print test
    if(IsZero(test)):
        print '[',target1,',',target2,']=',target3fac,target3
    else:
        print 'warning: [',target1,',',target2,']!=',target3fac,target3

def TestCommSum(target1,target2,target3fac,target3,target4fac,target4,base):
    test=comm(base[target1],base[target2] )-target3fac*(base[target3]+target4fac*base[target4])
    if((test==base[target3]*0).sum()==base[target3].size):
        print '[',target1,',',target2,']=',target3fac,'(',target3,target4fac,target4,')'
    else:
        print '[',target1,',',target2,']!=',target3fac,'(',target3,target4fac,target4,')'
        print test
def TestMult(target1,target2,target3fac,target3,base):
    test=(base[target1]*base[target2] )-target3fac*base[target3]
    if((test==base[target3]*0).sum()==base[target3].size):
        print '',target1,'*',target2,'=',target3fac,target3
    else:
        print '',target1,'*',target2,'!=',target3fac,target3
        print test



##################################################################
# Math functions

#return the commutator of a and b
#[a,b]=a*b-b*a
#sparse functions
def comm(a,b):
    #return numpy.dot(a,b)-numpy.dot(b,a) 
    return a*b-b*a

#def commSparse(a,b):
#    return a*b-b*a numpy.dot(a,b)-numpy.dot(b,a) 


#search through the labels to find the appropriate basis matrix
#def GetMat(a,base):
#    for i in range(len(label)):
#        if(a==label[i]):
#            return basis[i]
#    sys.stdout.write('could not find matrix %s\n' %a)
#    sys.exit(100)
    
#return hermitian conjugate
def HermConj(a):
    #b=numpy.transpose(a).conjugate()
    return a.transpose().conjugate()
    #return compare(b,base)[1]

def classMat(a):
    amaxr = a.real.max()  # maximum real part
    amaxi = a.imag.max()  # maximum imaginary part
    aminr = a.real.min()  # maximum real part
    amini = a.imag.min()  # maximum imaginary part
    if(amaxi==0 and amini==0):        #max val is real
        vtp='r'
        if(numpy.fabs(amaxr)>=numpy.fabs(aminr)):
            val=amaxr
        else:
            val=aminr
    elif(amaxr==0 and aminr==0):        #max val is imaginary
        vtp='i'
        if(numpy.fabs(amaxi)>=numpy.fabs(amini)):
            val=amaxi
        else:
            val=amini
    else:
        vtp='c'
        normMin=aminr**2.+amini**2.
        normMax=amaxr**2.+amaxi**2.
        if(normMax>normMin):
            val=normMax**0.5
        else:
            val=-1.*normMin**0.5
        pass
    return val,vtp

def compMethod1(a,base):
    verb='n'
    aval,atp=classMat(a)
    for key,vals in base.items():
        #if(key=='H2p'):
        #    verb='y'
        #else:
        #    verb='n'

        vval,vtp=classMat(vals)
        
        if(verb=='y'):
            print 'comparing;'
            print key,aval,atp,vval,vtp
            print a.real.max()  # maximum real part
            print a.imag.max()  # maximum imaginary part
            print a.real.min()  # maximum real part
            print a.imag.min()  # maximum imaginary part        
            print vals.real.max()  # maximum real part
            print vals.imag.max()  # maximum imaginary part
            print vals.real.min()  # maximum real part
            print vals.imag.min()  # maximum imaginary part        
            print aval/vval
        num=aval/vval
        numT='1'
        if(atp==vtp):
            numV=num
        elif( (atp=='r' and vtp=='i') or (atp=='i' and vtp=='r')):
            numV=num*complex(0,1)
            numT='i'
        else:
            #we're complex
            numV=num
            #print num,atp,vtp
            pass

        if(verb=='y'):
            print a
            print vals*num
        #print vals
        #print scipy.equal(a,vals*numV).all()
        if((a!=vals*numV).nnz==0):

        #if(scipy.equal(a,vals*numV).all()==True):
            num=Fraction(num).limit_denominator() 
            return num,numT,key,vals                   
        if((a!=vals*numV*-1).nnz==0):
        #if(scipy.equal(a,-1*vals*numV).all()==True):
            num=Fraction(-1*num).limit_denominator() 
            return num,numT,key,vals                   
    return 0,'shit','none',a



def compMethod2(a,base,maxy,even):
    for j in range(maxy):
        if(even==True and j%2==0 and j!=0): #skip odd js >=3
            continue
        for key,vals in base.items():
            for k in 1,-1:
                if((a!=1.*k/(j+1)*vals).nnz==0):
                #if(scipy.equal(a,1.*k/(j+1)*vals).all()==True):
                    ans=key
                    baso=vals
                    return Fraction(k,(j+1)),1,ans,baso                    
                if((a!=1.*k*(j+1)*vals).nnz==0):
                #if(scipy.equal(a,-1.*k/(j+1)*vals).all()==True):
                    ans=key
                    baso=vals
                    return Fraction(k*(j+1),1),1,ans,baso                    
                if((a!=1.*k*(j+1)*vals*complex(0,1)).nnz==0):
                #if(scipy.equal(a,-k*(j+1)*vals*complex(0,1)).all()==True):
                    ans=key
                    baso=vals
                    return Fraction(k*(j+1),1),'i',ans,baso                    
                if((a!=1.*k/(j+1)*vals*complex(0,1)).nnz==0):
                #if(scipy.equal(a,-1.*k/(j+1)*vals*complex(0,1)).all()==True):
                    ans=key
                    baso=vals
                    return Fraction(k,(j+1)),'i',ans,baso                    

    #sys.stdout.write('           Found no exact answer.\n')
    return 0,'shit','none',a

#compare a given matrix to the basis and look for a symbolic match
def compare(a,base,maxy=20,even=False):

    a1,b1,c1,d1=compMethod1(a,base)
    #a2,b2,c2,d2=compMethod2(a,base,maxy,even)
    
    #print 'new',a1,b1,c1,d1
    #print 'old',a2,b2,c2,d2
    
    """
    if(a1!=a2):
        print
        print 'different results'
        print a1,b1,c1
        print a2,b2,c2
        print d1*a1
        print d2*a2
        #sys.exit(100)
    """
    return a1,b1,c1,d1
        

#compare two matrices to find if they are the same
def evalComp(a,b):
    c=a-b
    return numpy.sum(numpy.abs(c))


#return the magnitude of the matrix. Complex magnitude of each element taken.
def IsZero(a):
    return not numpy.any(a)
def IsZeroSparse(a):
    return scipy.sparse.csr_matrix.count_nonzero(a)==0


#return the trace of b*a
#  <a|b> = tr (a^* b)
def Trace(a,b):
    c=a.transpose().conj() * b
    tr= scipy.trace(c)
    #tr= c.diagonal().sum() #
    return (float(tr.real))

def TraceSparse(a,b):
    c=a.transpose().conj() * b
    #tr= scipy.trace(c)
    tr= c.diagonal().sum() #
    return (float(tr.real))
