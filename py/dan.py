#!/usr/bin/python


import random,numpy,sys,datetime
from numba import njit,prange


def methyl_list_generator(residue_list):
    #EXAMPLES OF RESIDUE NAMES:
    #'ILE   57', 'LEU  123', 'VAL  176'
    methyl_list = []
    for residue in residue_list:
        if(residue[0:3]=='ILE'):
            methyl_list.append(residue+' A')
        if(residue[0:3]=='VAL' or residue[0:3]=='LEU'):
            ab = random.choice([' A', ' B'])
            methyl_list.append(residue+ab)
    return methyl_list



#Create fake ILV list.

resNum=100
tests=10000
intys=numpy.random.randint(0,3,resNum)  #random numbers 0,1,2, resNum times

#get residue number from random number generator list.
resys=numpy.random.randint(0,500,1000) #random numbers between 0-499, 1000 times
resys=numpy.unique(resys)    #make these unique...
resys=sorted(resys)[0:resNum]  #take the first resnum
resys=numpy.array(resys).astype('str')

#setup array of methyl types (strings)
methyl_typ=intys.copy().astype('str')
iles=intys==0  #find the iles
leus=intys==1  #find the leus
vals=intys==2  #find the vals
maskLV=numpy.logical_or(leus,vals) #get LVs
sz=int(numpy.sum(maskLV*1.))
print 'LVs:',sz

methyl_typ[iles]='ILE'
methyl_typ[leus]='LEU'
methyl_typ[vals]='VAL'
#synthesise residue type and number into one list
methyl_list= [i + j for i, j in zip(methyl_typ,resys)]  

Anew=methyl_typ.copy()
 
print 'synthesised methyl list:'
print methyl_list



AorB=intys.copy()  #will be an array with only As(0) or Bs(1)
AorB[iles]=0  #set the iles to A

now=datetime.datetime.now()
for i in range(tests):
    res=methyl_list_generator(methyl_list)
now2=datetime.datetime.now()

print res
print 'Old method: ',now2-now





def methyl_list_generator2():
    LVab=numpy.random.randint(0,2,sz)
    AorB[maskLV]=LVab  #input random As or Bs for LVs
    Anew[AorB==0]='A'
    Anew[AorB==1]='B'
    #return map(str.__add__,methyl_list,Anew)
    return [i + j for i, j in zip(methyl_list,Anew)] 


def methyl_list_generator3():
    LVab=numpy.random.randint(0,2,sz)
    #AorBinstance=AorB.copy()
    #AorBinstance[maskLV]=LVab  #input random As or Bs for LVs
    #return AorBinstance
    #AorBinstance=AorB.copy()
    AorB[maskLV]=LVab  #input random As or Bs for LVs
    return AorB
    #Anew=methyl_typ.copy()
    #Anew[AorBinstance==0]='A'
    #Anew[AorBinstance==1]='B'
    #return [i + j for i, j in zip(methyl_list,Anew)] 

@njit
def methyl_list_generator3a(sz,AorB,maskLV,tests):
    for i in range(tests):
        LVab=numpy.random.randint(0,2,sz)
        AorB[maskLV]=LVab  #input random As or Bs for LVs
        spec=AorB.copy()

#method using strings
now=datetime.datetime.now()
for i in range(tests):
    res=methyl_list_generator2()
now2=datetime.datetime.now()

print res
print methyl_list[numpy.arange(res]
print 'New method: ',now2-now



#return list with 0s or 1s if you want A or B
now=datetime.datetime.now()
for i in range(tests):
    res=methyl_list_generator3()
now2=datetime.datetime.now()

print res
print methyl_list[res]
print 'New method2:',now2-now



#method using strings
methyl_list_generator3a(sz,AorB,maskLV,0) #compile the jit...
now=datetime.datetime.now()
methyl_list_generator3a(sz,AorB,maskLV,tests)
now2=datetime.datetime.now()

print 'New method3: ',now2-now

