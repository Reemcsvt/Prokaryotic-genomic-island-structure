##############################################################################################
# Copyright and terms of use (DO NOT REMOVE):
# The code is made freely available for non-commercial uses only, provided that the copyright 
# header in each file not be removed, and suitable citation(s) (see below) be made for papers 
# published based on the code.
#
# The copyright of the code is retained by the authors.  By downloading/using this code you
# agree to all the terms stated above.
#
#   Reem Aldaihani, Lenwood S. Heath 
#   "Connecting Genomic Islands across Prokaryotic and Phage Genomes via Protein Families" 
#   Department of Computer Science, Virginia Tech.
#   Blacksburg, VA, 2022. 
#
# Copyright (c) 2022, Reem Aldaihani All rights reserved.
##############################################################################################

from collections import Counter
from itertools import combinations 
import itertools
Gen=[]
Spe=[]
AccN=[]
with open('Genus_Species_AccN_New_3.txt', 'r') as f: #Species names and accession number
    for line in f:
       Gen.append(line.split(None, 1)[0])
       Spe.append(line.split(None, 2)[1])
       AccN.append(line.split(None, 3)[2])
print('Done read Genus species Acc')     
################################################################################
dataset=[]
dataset_=[]
count=0
with open('Apriori_New_5.txt', 'r') as f: #Apriori_New_5.txt  Ap.txt
  for line in f:
      dataset.append([])
      dataset_.append([])
      PF=int(line.split(None, 2)[1])
      for i in range (2,PF+2):
            dataset[count].append(line.split(None, i+1)[i])   
            dataset_[count].append(line.split(None, i+1)[i]+'$')              
      count=count+1   
          
print('Done read data set')   
################################################################################
GIs=[]
PFs=[]
PFs_All=[]
Indx=-1
with open('New_GIsANDtheirProteins_New_5.txt', 'r') as f: #New_GIsANDtheirProteins_New_5.txt  GIs_PFs.txt
  for line in f:
    A=line.split(None, 1)[0]
    if '$' not in A:
        GIs.append(A)
        PFs.append([])
        PFs_All.append([]) 
        Indx=Indx+1
    else:
        PFs[Indx].append(A.split('$')[1])
        PFs_All[Indx].append(A)

################################################################################
Wanted_PFs=open('Apri_Patterns_Occ_New_5.txt', 'w')

Indx=-1
Position=0
for Pattern in dataset:    
    #Save pattern
    Geno_Islands=0
    Indx=Indx+1
    print(Indx)
    Pattern_Size=len(Pattern)
                           
    for i in PFs:
       if len(set(Pattern).intersection(set(i)))==len(Pattern):                    
             Geno_Islands=Geno_Islands+1 
        
    #if Geno_Islands>10:  
    Wanted_PFs.write(str(Pattern_Size))
    for i in Pattern:
       Wanted_PFs.write(' '+i)           
    Wanted_PFs.write(' '+str(Geno_Islands)+'\n')     

print('Done G_S_PFs_Apri_InOrder_Top_New')    














        