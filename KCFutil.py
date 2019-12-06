# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 11:04:26 2018

@author: kcf
"""
#########################################################################################################
#                    Peptide Class + some functions and support classes
#                          written by Kilian Conde-Frieboes
#                                     11.06.2019
# For question come by and ask or write to kcf@novonordisk.com
#########################################################################################################







import re
import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
from math import sqrt


#########################################################################################################
def factorial(n):
    """
    n should be a integer n>=0
    It calculates n!=n*(n-1)*(n-2)... *(n-(n-1))
    Limit seems to be the number of recursion possible, on Win32 n=2970
    factorial calls itself with n-1 down to n=1    
    n!=n*factorial(n-1)
    0! = 1 by definition
    """
    if type(n)!=int:
        print ("only int allowed\n")
        return
    if n==0:
        return 1
    if n<0:
        print (" Only n>=0 allowed \nPrinting abs(n)!")
        n=-n
    m=n*factorial(n-1)
    return m
###########################################################################################################




############################################################################################################
class molecule:
    '''
    Molecule Class
    extract_count() and parse_sf() convert an entered sumformula string (CH4) into numbers elements('C')=1
    calc_peaks() and isopattern calculate the isotope peak abundance. Code was org. written ~15 years ago in C/C++ 
    and converted now into python
    '''
  
    def __init__(self,sumformula='CH3Cl'):
        
        self.elements = {
            'C':0,
            'H':0,
            'N':0,
            'O':0,
            'Cl':0,
            'S':0
        }
        self.sumformula=sumformula
        self.parse_sf(self.sumformula)
        self.mw=0.0
        self.monomass=0.0
        self.isotopepeak=[1.0]*20 # array to safe the isotope pattern      
        for ele in self.elements:
            self.mw+=element_mass_avg[ele]*self.elements[ele]
            self.monomass+=element_mass_mono[ele][0]*self.elements[ele]
#        self.isopattern() # subroutine to calc the isotope pattern
        
    def info(self):
        self.isopattern() # subroutine to calc the isotope pattern
        x=np.arange(20)
        for i in range (20):
            x[i]+=self.monomass
        plt.bar(x, height=self.isotopepeak)
        plt.xlabel('m/z')
        plt.ylabel('rel. Abundance')
        plt.title('Isotope pattern')
        info='MW {:8.4f}\nMono isotopic Mass {:8.4f}'.format(self.mw, self.monomass)
#        plt.savefig('BayesianRidge trainset PCA', dpi=300,bbox_inches='tight', pad_inches=0.2)
        plt.gcf().text(1.0,0.75,info, fontsize=12)
        plt.show()
        return
        
    def extract_count(self, sumformula='CH4',element='C'):
        'extracting the number of a certain element in a sum formula'
        position=sumformula.find(element)
        if len(sumformula)>position+len(element):
            if sumformula[position+len(element)].isalpha():
                self.elements[element]=1
                re.sub(element,'',sumformula)
            else:
                m=re.findall(element+'(\d+)',sumformula)
                self.elements[element]=int(m[0])
                re.sub(element+m[0],'',sumformula)
        else:
            self.elements[element]=1
            re.sub(element,'',sumformula)
        return sumformula;
    
    def parse_sf(self,sumformula): #double lettered elements on top. 
        if 'Cl' in sumformula:
            sumformula=self.extract_count(sumformula,'Cl')
        if 'C' in sumformula:
            sumformula=self.extract_count(sumformula,'C')
        if 'H' in sumformula:
            sumformula=self.extract_count(sumformula,'H')
        if 'N' in sumformula:
            sumformula=self.extract_count(sumformula,'N')
        if 'O' in sumformula:
            sumformula=self.extract_count(sumformula,'O')   
        if 'S' in sumformula:
            sumformula=self.extract_count(sumformula,'S')

            
        
    
        
    def calcpeaks(self,array1=[],array2=[]):               #subroutine for isotope patter calc
        matrix=[[0.0]*20 for i in range(20)]
        for i in range (20):
            for j in range (20):
                if i+j<20:
                    matrix[i][i+j]=array1[i]*array2[j]
                    
        for i in range (20):
            array1[i]=0.0
            for j in range (20):
                array1[i]+=matrix[j][i]
        
    def isopattern(self):                                   #subroutine for isotope patter calc
        elementpeaks = {
            'C':0,
            'H':0,
            'N':0,
            'O':0,
            'Cl':0,
            'S':0
        }
        for ele in self.elements:
            elementpeaks[ele]= [0.0]*20
        dummy = [0.0]*20
        for ele in self.elements:
            if self.elements[ele]>0:
                for i in range (20):
                    if (i<4):
                        elementpeaks[ele][i]=element_ab[ele][i]/element_ab[ele][0]
                        dummy[i]=element_ab[ele][i]/element_ab[ele][0]
                    else:
                        elementpeaks[ele][i]=0.0
                        dummy[i]=0.0
            for i in range(2,self.elements[ele]+1):
                self.calcpeaks(elementpeaks[ele],dummy) #subroutine here
        for i in range(20):
            self.isotopepeak[i]=0.0
        self.isotopepeak[0]=1.0
        for ele in self.elements:
            if self.elements[ele]>0:
                self.calcpeaks(self.isotopepeak,elementpeaks[ele]) #subroutine here
        maxpeak=0.0                                                         # Finding the max isotope and norm it to 1.0
        for i in range(20):                                                 #
            maxpeak=max(maxpeak,self.isotopepeak[i])                        #
        for i in range(20):                                                 #
            self.isotopepeak[i]=self.isotopepeak[i]/maxpeak                 #
                

                
   
            
######################################################################################################################################################
class peptide(molecule):
    '''
    peptide class skeleton
    '''
    def __init__(self,sequence='KILIAN', amide=False, acetyl=False, cyclic=False, label='', qda = 0.0, exist=False):
        self.sequence = sequence #list of amino acids objects
        self.length=len(sequence)
        self.elements = {
            'C':0,
            'H':0,
            'N':0,
            'O':0,
            'Cl':0,
            'S':0
        }
        self.sumformula=''
        self.mw=0.0
        self.monomass=0.0
        self.acetyl=acetyl
        self.amide=amide
        self.cyclic=cyclic
        self.label=label
        self.biol = []
        self.biolpre = 0.0
        self.qda = qda        #mass measured
        self.exist=exist
        self.zscale3=[]
        self.zscale5=[]
        self.zscale5mat=[]
        if self.cyclic==False:
            if self.amide==False:
                self.elements['O']+=1
                self.elements['H']+=1
            else:
                self.elements['N']+=1
                self.elements['H']+=2
            if self.acetyl==False:
                self.elements['H']+=1
            else:
                self.elements['0']+=1
                self.elements['H']+=3
                self.elements['C']+=2
        self.parse_sequence(sequence)
        self.parse_sf(self.sumformula)
        self.isotopepeak=[1.0]*20 # array to safe the isotope pattern      
        for ele in self.elements:
            self.mw+=element_mass_avg[ele]*self.elements[ele]
            self.monomass+=element_mass_mono[ele][0]*self.elements[ele]
        #self.isopattern() # subroutine to calc the isotope pattern
        for aa in self.sequence:
            self.zscale5mat.append(aminoacids[aa][2])
            for i in range(5):
                if i<3:
                    self.zscale3.append(aminoacids[aa][2][i])
                self.zscale5.append(aminoacids[aa][2][i])
    
    def parse_sequence(self,sequence):
        for i in range(len(sequence)):
            #self.sequence.append(aminoacids[sequence[i]])
            for ele in aminoacids[sequence[i]][1]:
                self.elements[ele]+=aminoacids[sequence[i]][1][ele]
        sum=''        
        for ele in self.elements:
            if self.elements[ele]==0:
                continue
            else:
                sum+=ele
            if self.elements[ele]>1:
                sum+=str(self.elements[ele])
        self.sumformula=sum
        
    def get_zscale(self, number = 0, zscale = 0): #aa number from N-terminal by convention
        aa = self.sequence[number]
        
        return (aminoacids[aa][2][zscale])
    
            
        
######################################################################################################################################################

class aminoacid:
    '''
    amino acid class skeleton
    '''
    def __init__(self,AA='A'):
        self.type=AA
        
######################################################################################################################################################
element_mass_avg = {
    'C':12.011,
    'H':1.00794,
    'Cl':35.4527,
    'O':15.9994,
    'N':14.00674,
    'S':32.059
}

element_mass_mono = {
    'C':[12.0],
    'H':[1.007825],
    'Cl':[34.968853],
    'O':[15.994915],
    'N':[14.003074],
    'S':[31.972071]
}
element_ab = {
    'C':[0.989, 0.011, 0.0, 0.0],
    'H':[0.99985, 0.00015, 0.0, 0.0],
    'Cl':[0.7577, 0.0, 0.2423, 0.0],
    'O':[0.99762, 0.00038, 0.002, 0.0],
    'N':[0.99634, 0.00366, 0.0, 0.0],
    'S':[0.9502, 0.0075, 0.0421, 0.0002]
}
aminoacids = {
    'A':['C3H5NO',{'C':3,'H':5,'N':1,'O':1,'S':0},[0.24,-2.32,0.6,-0.14,1.3]],
    'C':['C3H5NOS',{'C':3,'H':5,'N':1,'O':1,'S':1},[0.84,-1.67,3.71,0.18,-2.65]],
    'D':['C4H5NO3',{'C':4,'H':5,'N':1,'O':3,'S':0},[3.98,0.93,1.93,-2.46,0.75]],
    'E':['C5H7NO3',{'C':5,'H':7,'N':1,'O':3,'S':0},[3.11,0.26,-0.11,-3.04,-0.25]],
    'F':['C9H9NO',{'C':9,'H':9,'N':1,'O':1,'S':0},[-4.22,1.94,1.06,0.54,-0.62]],
    'G':['C2H3NO',{'C':2,'H':3,'N':1,'O':1,'S':0},[2.05,-4.06,0.36,-0.82,-0.38]],
    'H':['C6H7N3O',{'C':6,'H':7,'N':3,'O':1,'S':0},[2.47,1.95,0.26,3.9,0.09]],
    'I':['C6H11NO',{'C':6,'H':11,'N':1,'O':1,'S':0},[-3.89,-1.73,-1.71,-0.84,0.26]],
    'K':['C6H12N2O',{'C':6,'H':12,'N':2,'O':1,'S':0},[2.29,0.89,-2.49,1.49,0.31]],
    'L':['C6H11NO',{'C':6,'H':11,'N':1,'O':1,'S':0},[-4.28,-1.3,-1.49,-0.72, 0.84]],
    'M':['C5H9NOS',{'C':5,'H':9,'N':1,'O':1,'S':1},[-2.85,-0.22,0.47,1.94,-0.98]],
    'N':['C4H6N2O2',{'C':4,'H':6,'N':2,'O':2,'S':0},[3.05,1.62,1.04,-1.15,1.61]],
    'P':['C5H7NO',{'C':5,'H':7,'N':1,'O':1,'S':0},[-1.66,0.27,1.84,0.7,2.0]],
    'Q':['C5H8N2O2',{'C':5,'H':8,'N':2,'O':2,'S':0},[1.75,0.5,-1.44,-1.34,0.66]],
    'R':['C6H12N4O',{'C':6,'H':12,'N':4,'O':1,'S':0},[3.52,2.5,-3.5,1.99,-0.17]],
    'S':['C3H5NO2',{'C':3,'H':5,'N':1,'O':2,'S':0},[2.39,-1.07,1.15,-1.39,0.67]],
    'T':['C4H7NO2',{'C':4,'H':7,'N':1,'O':2,'S':0},[0.75,-2.18,-1.12,-1.46,-0.4]],
    'V':['C5H9NO',{'C':5,'H':9,'N':1,'O':1,'S':0},[-2.59,-2.64,-1.54,-0.85,-0.02]],
    'W':['C11H10N2O',{'C':11,'H':10,'N':2,'O':1,'S':0},[-4.36,3.94,0.59,3.44,-1.59]],
    'Y':['C9H9NO2',{'C':9,'H':9,'N':1,'O':2,'S':0},[-2.54,2.44,0.43,0.04,-1.47]],
#    'X':['C4H7NO',{'C':4,'H':7,'N':1,'O':1,'S':0},[-1.33, -2.8, -0.61, -0.55, 0.4]] #Aib
}
####################################################################################################################################################
def ala_scan(peptidelist, pep):
    cyclic=pep.cyclic
    amide=pep.amide
    acetyl=pep.acetyl
    for i in range (pep.length):
        seq=pep.sequence[:i]+'A'+pep.sequence[i+1:]
        newpeptide = peptide(sequence=seq,cyclic=cyclic, amide=amide, acetyl=acetyl)
        peptidelist.append(newpeptide)
    
def con_mut(peptidelist, pep):
    #basic amino acids  K R H
    #aromatic           F Y W
    #acidic             E D
    #hydrophilic        S T N Q
    #hydrophobic        F Y W I L V M
    #subroutine to change the input sequence pep using conservative mutations
    cyclic=pep.cyclic
    amide=pep.amide
    acetyl=pep.acetyl
    for i in range (pep.length):
        if pep.sequence[i] in 'KRH':
            for aa in 'KRH':
                if aa!=pep.sequence[i]:
                    seq=pep.sequence[:i]+aa+pep.sequence[i+1:]
                    newpeptide = peptide(sequence=seq,cyclic=cyclic, amide=amide, acetyl=acetyl)
                    peptidelist.append(newpeptide)
            continue
        if pep.sequence[i] in 'FYW':
            for aa in 'FYW':
                if aa!=pep.sequence[i]:
                    seq=pep.sequence[:i]+aa+pep.sequence[i+1:]
                    newpeptide = peptide(sequence=seq,cyclic=cyclic, amide=amide, acetyl=acetyl)
                    peptidelist.append(newpeptide)
            continue
        if pep.sequence[i] in 'ED':
            for aa in 'ED':
                if aa!=pep.sequence[i]:
                    seq=pep.sequence[:i]+aa+pep.sequence[i+1:]
                    newpeptide = peptide(sequence=seq,cyclic=cyclic, amide=amide, acetyl=acetyl)
                    peptidelist.append(newpeptide)
            continue
        if pep.sequence[i] in 'STNQ':
            for aa in 'STNQ':
                if aa!=pep.sequence[i]:
                    seq=pep.sequence[:i]+aa+pep.sequence[i+1:]
                    newpeptide = peptide(sequence=seq,cyclic=cyclic, amide=amide, acetyl=acetyl)
                    peptidelist.append(newpeptide)
            continue
        if pep.sequence[i] in 'FYWILVM':
            for aa in 'FYWILVM':
                if aa!=pep.sequence[i]:
                    seq=pep.sequence[:i]+aa+pep.sequence[i+1:]
                    newpeptide = peptide(sequence=seq,cyclic=cyclic, amide=amide, acetyl=acetyl)
                    peptidelist.append(newpeptide)
            continue

def charge_scan(peptidelist, pep, basic='K',acid='E'):
    cyclic=pep.cyclic
    amide=pep.amide
    acetyl=pep.acetyl
    for i in range (pep.length):
        seq=pep.sequence[:i]+basic+pep.sequence[i+1:]
        newpeptide = peptide(sequence=seq,cyclic=cyclic, amide=amide, acetyl=acetyl)
        peptidelist.append(newpeptide)
        seq=pep.sequence[:i]+acid+pep.sequence[i+1:]
        newpeptide = peptide(sequence=seq,cyclic=cyclic, amide=amide, acetyl=acetyl)
        peptidelist.append(newpeptide)

def saltbridge_scan(peptidelist, pep, basic='K',acid='E', ra=0, dis=4):
    if ra==0:
        ra=pep.lenght
        
    cyclic=pep.cyclic
    amide=pep.amide
    acetyl=pep.acetyl
    for i in range (ra-dis):
        seq=pep.sequence[:i]+basic+pep.sequence[i+1:]
        seq=seq[:i+dis]+acid+seq[i+dis+1:]
        newpeptide = peptide(sequence=seq,cyclic=cyclic, amide=amide, acetyl=acetyl)
        peptidelist.append(newpeptide)
        seq=pep.sequence[:i]+acid+pep.sequence[i+1:]
        seq=seq[:i+dis]+basic+seq[i+dis+1:]
        newpeptide = peptide(sequence=seq,cyclic=cyclic, amide=amide, acetyl=acetyl)
        peptidelist.append(newpeptide)





def random_walk(peptidelist, pep, n=2, m=10, exclude_C=True): #n= number of mutations, m=number of peptides
    cyclic=pep.cyclic
    amide=pep.amide
    acetyl=pep.acetyl
    for i in range(m):
        for j in range(n):
            mut=random.randint(0,len(pep.sequence)-1) #pythons randint goes from x to z both included, numpy excludes z
            possible_aa = list(aminoacids.keys())
            if exclude_C == True:
                possible_aa.remove('C')
            random_aa = random.choice(possible_aa)
            #print (possible_aa)
            seq=pep.sequence[:mut]+random_aa+pep.sequence[mut+1:]
        newpeptide = peptide(sequence=seq,cyclic=cyclic, amide=amide, acetyl=acetyl)
        peptidelist.append(newpeptide)
        
def random_walk2(peptidelist, pep, n=2, m=10, exclude_C=True, old_pep_dict={}): 
    #n= number of mutations, m=number of peptides
    # this routine checks for duplicates and creates a list with m unique sequences
    # neither this nor random_walk is working properly
    # seems to work now, 12.Nov. 2019
    # will also check whether the new peptides already exist in a dict
    # now with up n mutation, adding a randint(1,n) statement
    
    cyclic=pep.cyclic                                           
    amide=pep.amide
    acetyl=pep.acetyl                                            
    listofstrings = {}                                           
    
    listofstrings[pep.sequence]=1
    while len(listofstrings)<m+1:
        seq=pep.sequence
        number_of_mut=random.randint(1,n)
        for j in range(number_of_mut):
            mut=random.randint(0,len(pep.sequence)-1) #pythons randint goes from x to z both included, numpy excludes z
            possible_aa = list(aminoacids.keys())
            if exclude_C == True:
                possible_aa.remove('C')
            random_aa = random.choice(possible_aa)
            #print (pep.sequence[:mut]+'\x1b[6;30;42m'+random_aa+'\x1b[0m'+pep.sequence[mut+1:])
            seq=seq[:mut]+random_aa+seq[mut+1:]
        if seq in old_pep_dict.keys():
            continue
        listofstrings[seq]=1
    listofstrings.pop(pep.sequence)
    for i in listofstrings:        
        newpeptide = peptide(sequence=i,cyclic=cyclic, amide=amide, acetyl=acetyl)
        peptidelist.append(newpeptide)


#this routine checks for duplicates and creates a list with m unique sequences out of the top x peptides
def random_walk3(peptidedict, top=10, n=2, m=10, exclude_C=True, feat=0):
    possible_aa = list(aminoacids.keys())
    if exclude_C == True:
        possible_aa.remove('C')
    peptidelist=[]
    for peptides in peptidedict:
        peptidelist.append(peptidedict[peptides])
    peptidelist.sort(key=lambda x: x.biol[feat], reverse=False)
    for i in range (top):
        k=0
        cyclic=peptidelist[i].cyclic
        amide=peptidelist[i].amide
        acetyl=peptidelist[i].acetyl
        while k <m:
            seq=peptidelist[i].sequence
            for j in range(n):
                mut=random.randint(0,len(seq)-1) #pythons randint goes from x to z both included, numpy excludes z
                random_aa = random.choice(possible_aa)
                seq=seq[:mut]+random_aa+seq[mut+1:]
            if seq in peptidedict:
                continue
            else:
                k=k+1
                newpeptide=peptide(sequence=seq,cyclic=cyclic, amide=amide, acetyl=acetyl)
                peptidedict[seq]=newpeptide
            
    return

def write_lib(peplist, filename='test.txt', row1='.space', row2='.space', n=80):
    seq_list= []
    for k in range(n//80):
        seq_list.append('; Plate '+str(k+1)) 
        for i in range (8):
            seq_list.append(row1)
            seq_list.append(row2)
            for j in range (10):
                blanked_seq=''.join(peplist[i*10+j+k*80].sequence[m]+' 'for m in range(len(peplist[i*10+j+k*80].sequence)-1))
                blanked_seq=blanked_seq+peplist[i*10+j+k*80].sequence[len(peplist[i*10+j+k*80].sequence)-1]
                seq_list.append(blanked_seq)   
    newfile=pd.DataFrame(data=seq_list)  
    newfile.to_csv(path_or_buf=filename, header=False, index=False, sep=',')
    
    
#*************************************************************************************************************************************
def distance(peptide1,peptide2):
    if type(peptide1)!=peptide:
        print ('only peptides allowed!\n')
        return
    if peptide1.length!=peptide2.length:
        print ('Peptides have to have the same length!\n')
        return
    sqdis=0
    for i in range (peptide1.length):
        for j in range (5):
            sqdis=sqdis+(peptide1.get_zscale(i,j)-peptide2.get_zscale(i,j))**2
    return sqrt(sqdis)
#*************************************************************************************************************************************