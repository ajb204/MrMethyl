#!/usr/bin/python

import basis,string
import sys
from math import sqrt,fabs
import numpy

from fractions import Fraction

from scipy import mat,zeros

from basis import comm,GetMat,HermConj,compare,IsZero,Trace,GetBasis,TestComm
#from commute import compare

basis,label=GetBasis('C2H2')
#print label
#a=GetMat('I1x',basis,label)
#b=GetMat('I2zI1y',basis,label)
#c=GetMat('I2zI1z',basis,label)

#print (a*b-b*a).shape

#test=comm(GetMat('I2zI1z',basis,label),GetMat('I1x',basis,label)-(complex(0,1))*GetMat('I2zI1y',basis,label))
#print test

def ApplyProdOp(term,ham,coeff):

    NewTerm=[]
    for i in range(len(term)):
        rho=term[i][0]
        pre=term[i][1]

        one=GetMat(rho,basis,label)
        test=GetMat(ham,basis,label)
        check1=comm(one,test)
        if(numpy.sum(numpy.abs(check1))==0):
            fac=0
            #print rho,'---',coeff,ham,'--->',rho
            NewTerm.append((rho,pre))

        else:
            fac,resy,baso=compare(check1,basis,label,maxy=1)
            #print rho,'---',coeff,ham,'--->',rho,'cos(',coeff,')','+',str(fac)+'*'+resy,'sin(',coeff,')'
            if(coeff=='pi'):
                NewTerm.append((rho,'-1*('+pre+')'))
            elif(coeff=='pi/2'):
                NewTerm.append((resy,'+'+str(fac)+'*('+pre+')'))
            else:
                NewTerm.append((rho,'+cos('+coeff+')*('+pre+')'))
                NewTerm.append((resy,'+'+str(fac)+'*sin('+coeff+')*('+pre+')'))

    return NewTerm

def ShowTerm(NewTerm):
    for term in NewTerm:
        print term

    

#sort out cos/sin pi/pi0.5 terms, and any 1/4Js specified
def Simplify(terms,Js):
    NewTerm=[]
    print 
    print 'Simplifying:'
    print 'starting with ',len(terms),'terms'
    for term in terms:
        pre=term[1]
        tig=0
        if(len(pre.split('sin(pi)'))>1):
           tig=1
        if(len(pre.split('cos(pi/2)'))>1):
           tig=1

        for i in range(len(Js)):
            if(len(pre.split('cos(pi'+Js[i]+'t)'))>1):
                tig=1
            pre=pre.replace('sin(pi'+Js[i]+'t)*','1*')
        pre=pre.replace('cos(pi)*','1*')
        pre=pre.replace('sin(pi/2)*','1*')
        pre=pre.replace('+-','-')
        if(tig==0):
            try:
                NewTerm.append((term[0],str(eval(pre))))
            except:
                NewTerm.append((term[0],pre))

    print 'ending with ',len(NewTerm),'terms'
    return NewTerm

def FreePrecess(time,NewTerm,StrongC,StrongH,crossCH):
    NewTerm=ApplyProdOp(NewTerm,'I1zS1z','piJc1h1'+time)
    NewTerm=ApplyProdOp(NewTerm,'I2zS2z','piJc2h2'+time)

    if(crossCH=='y'):
        NewTerm=ApplyProdOp(NewTerm,'I1zS2z','piJc1h2'+time)
        NewTerm=ApplyProdOp(NewTerm,'I2zS1z','piJc2h1'+time)

    if(StrongC[0]==1):
        NewTerm=ApplyProdOp(NewTerm,'S2zS1z','piJcc'+time)
        if(StrongC[1]=='y'):
            NewTerm=ApplyProdOp(NewTerm,'S2xS1x','piJcc'+time)
            NewTerm=ApplyProdOp(NewTerm,'S2yS1y','piJcc'+time)
    if(StrongH[0]==1):
        NewTerm=ApplyProdOp(NewTerm,'I2zI1z','piJhh'+time)
        if(StrongH=='y'):
            NewTerm=ApplyProdOp(NewTerm,'I2xI1x','piJhh'+time)
            NewTerm=ApplyProdOp(NewTerm,'I2yI1y','piJhh'+time)
    return NewTerm



def RunHMQC(start,Js,strongC,strongH,crossHC):
    print 'Starting HMQC:'
    ShowTerm(start)
    print
    NewTerm=start
    
    NewTerm=ApplyProdOp(NewTerm,'I1x','pi/2')
    NewTerm=ApplyProdOp(NewTerm,'I2x','pi/2')

    print 'FreePrecession1:'
    NewTerm=FreePrecess('t',NewTerm,strongC,strongH,crossHC)
    print 'Terms:',len(NewTerm)

    print 'Applying pulses:'
    NewTerm=ApplyProdOp(NewTerm,'S1x','pi/2')
    NewTerm=ApplyProdOp(NewTerm,'S2x','pi/2')
    NewTerm=ApplyProdOp(NewTerm,'I1x','pi')
    NewTerm=ApplyProdOp(NewTerm,'I2x','pi')
    NewTerm=ApplyProdOp(NewTerm,'S1x','pi/2')
    NewTerm=ApplyProdOp(NewTerm,'S2x','pi/2')
    print 'Terms:',len(NewTerm)
    
    print 'FreePrecession2:'
    NewTerm=FreePrecess('t',NewTerm,strongC,strongH,crossHC)
    print 'Terms:',len(NewTerm)    
    
    NewTerm=Simplify(NewTerm,Js)
    
    ShowTerm(NewTerm)
    print 'Terms:',len(NewTerm)
    return NewTerm


def Inept(time,NewTerm,strongC,strongH,crossHC):
    print 'Inept FreePrecession1:'
    NewTerm=FreePrecess(time,NewTerm,strongC,strongH,crossHC)


    NewTerm=ApplyProdOp(NewTerm,'S1x','pi')
    NewTerm=ApplyProdOp(NewTerm,'S2x','pi')
    NewTerm=ApplyProdOp(NewTerm,'I1x','pi')
    NewTerm=ApplyProdOp(NewTerm,'I2x','pi')
    print 'Inept FreePrecession2:'
    NewTerm=FreePrecess(time,NewTerm,strongC,strongH,crossHC)
    print 'Terms:',len(NewTerm)
    return NewTerm

def Purge(terms):
    NewTerm=[]
    for term in terms:
        tig=0
        if(len(term[0].split('x'))>1):
            tig=1
        if(len(term[0].split('y'))>1):
            tig=1
        if(tig==0):
            NewTerm.append(term)


    return NewTerm

def RunHSQC(start,Js,strongC,strongH,crossHC):
    print 'Starting HSQC:'
    ShowTerm(start)
    print
    NewTerm=start
    
    NewTerm=ApplyProdOp(NewTerm,'I1x','pi/2')
    NewTerm=ApplyProdOp(NewTerm,'I2x','pi/2')

    NewTerm=Inept('t/2',NewTerm,strongC,strongH,crossHC) #first inept


    print 'Applying pulses:'
    NewTerm=ApplyProdOp(NewTerm,'I1y','pi/2')
    NewTerm=ApplyProdOp(NewTerm,'I2y','pi/2')

    NewTerm=Purge(NewTerm)

    NewTerm=ApplyProdOp(NewTerm,'S1x','pi/2')
    NewTerm=ApplyProdOp(NewTerm,'S2x','pi/2')


    NewTerm=ApplyProdOp(NewTerm,'I1x','pi')
    NewTerm=ApplyProdOp(NewTerm,'I2x','pi')



    NewTerm=ApplyProdOp(NewTerm,'S1x','pi/2')
    NewTerm=ApplyProdOp(NewTerm,'S2x','pi/2')

    NewTerm=Purge(NewTerm)

    NewTerm=ApplyProdOp(NewTerm,'I1x','pi/2')
    NewTerm=ApplyProdOp(NewTerm,'I2x','pi/2')

    NewTerm=Inept('t/2',NewTerm,strongC,strongH,crossHC) #final inept
    
    NewTerm=Simplify(NewTerm,Js)
    
    ShowTerm(NewTerm)
    print 'Terms:',len(NewTerm)
    return NewTerm

#print NewTerm
#TestComm('S1x','S1y',complex(0,1),'S1z',basis,label)
#TestComm('S1z','S1x',complex(0,1),'S1y',basis,label)
##TestComm('S1y','S1z',complex(0,1),'S1x',basis,label)
##TestComm('S1x','S2zS1y',complex(0,1),'S2zS1z',basis,label)
#TestComm('S1y','S2zS1x',complex(0,-1),'S2zS1z',basis,label)
#TestComm('S1x','S2yS1z',complex(0,-1),'S2yS1y',basis,label)
#TestComm('I1x','I1yS1z',complex(0,1),'I1zS1z',basis,label)
#TestComm('I2x','I2yS1z',complex(0,1),'I2zS1z',basis,label)
#TestComm('I1y','I1xS1z',complex(0,-1),'I1zS1z',basis,label)
#TestComm('I2y','I2xS1z',complex(0,-1),'I2zS1z',basis,label)
#TestComm('I1y','I1xS1z',complex(0,1),'I1zS1z',basis,label)
#TestComm('I2y','I2xS1z',complex(0,1),'I2zS1z',basis,label)

#a,b,c= compare(complex(-1,0)*GetMat('I2y',basis,label),basis,label,maxy=1)


        #fac,resy,baso=compare(check1,basis,label)

def SageOperator(op):
    tag=''
    if(len(op)>=3):
        tag+=op[0:3]
    if(len(op)>=6):
        tag+='*'+op[3:6]
    if(len(op)>=9):
        tag+='*'+op[6:9]
    if(len(op)>=12):
        tag+='*'+op[9:12]
    #print 'rawTag',tag,op
    if(tag[0:2]!='I2'):
        tag='I2E*'+tag
    if(tag[4:6]!='I1'):
        if(tag[-1:]!="*"):
            tag+='*'
        tag=tag[:4]+'I1E*'+tag[4:]
    if(tag[8:10]!='S2'):
        if(tag[-1:]!="*"):
            tag+='*'
        tag=tag[:8]+'S2E*'+tag[8:]
    if(tag[12:14]!='S1'):
        if(tag[-1:]!="*"):
            tag+='*'
        tag+='S1E'

    if(len(tag)==16):
        tag=tag[:15]
    return tag


def MakeSage(NewTerm):


    

    NewTerm=numpy.array(NewTerm)
    uni,args,cnts=numpy.unique( NewTerm[:,0] ,return_inverse=True,return_counts=True)
    print 'Unique terms:',len(uni)
    print uni


    #print NewTerm[args==0,0]
    print
    print
    print 'b=0'
    for i in range(len(uni)):
        sage=''
        for j in range(cnts[i]):
            op=NewTerm[args==i,0][j]
            va=NewTerm[args==i,1][j]
            op=SageOperator(op)
            sage+='+'+op+'*('+va.replace('t','*t').replace('pi','pi*')+')'
        print 'b'+str(i)+'= (matrix(SR,[str(\''+sage+'\'),])[0][0]).expand()'#.full_simplify()'
    #print 'veccy=b.subs(Jc2h2=Jc1h1).subs(t=1/(2*Jc1h1)).expand().factor().simplify_trig()'
    #print 'print veccy'
#.factor()#.full_simplify()
#print veccy


Js=('Jc1h1','Jc2h2',)
Js=[]



start=[]
start.append(('I1z','1'))
start.append(('I2z','1'))

NewTerm=RunHMQC(start,Js,(1,'y'),(1,'y'),'y')

#NewTerm=RunHSQC(start,Js,(1,'n'),(1,'n'),'y')
MakeSage(NewTerm)


