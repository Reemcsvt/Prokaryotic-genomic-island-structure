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
################################################################################
dataset=[]
dataset_=[]
count=0
Pattern_N=[]
Apriori_AllComp=open('Apriori_AllCombinations_New_5.txt', 'w')

with open('Apriori_New_5.txt', 'r') as f: #Apriori_New_5.txt  Ap.txt
  for line in f:
      Temp=[]
      PF=int(line.split(None, 2)[1])
      for i in range (2,PF+2):
            Temp.append(line.split(None, i+1)[i])        
      c = list(itertools.permutations(Temp, len(Temp)))     
      
      Pattern_N.append(count)
      for i in c:
          dataset.append([])
          dataset_.append([])
          for j in i:
            dataset[count].append(j) 
            dataset_[count].append(j+'$')
            Apriori_AllComp.write(j+' ')   
                      
          Apriori_AllComp.write('\n')            
          count=count+1
          
print('Done read data set')   
################################################################################
GIs=[]
PFs=[]
PFs_All=[]
Indx=-1
with open('New_GIsANDtheirProteins_New_5.txt', 'r') as f: #New_GIsANDtheirProteins_New_5.txt  GIs_PFs_Apri.txt
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
Wanted_PFs=open('Apriori_AllCombinations_Occ_New_5.txt', 'w')
Wanted_PFs_WithIndx=open('Apriori_AllCombinations_Occ_Index_New_5.txt', 'w')

Indx=-1
Position=0
for Pattern in dataset:    
    #Save pattern
    Geno_Islands=0
    Indx=Indx+1
    print(Indx)
    Pattern_Size=len(Pattern)
            
    PFs_Indx=-1               
    for i in PFs:
       PFs_Indx=PFs_Indx+1
       if len(set(Pattern).intersection(set(i)))==len(Pattern):                    
          idx=0
          for z in i:
           if idx == Pattern_Size:
             break
           else:  
            if Pattern[idx] ==z:
               idx=idx+1 
          
          if idx== Pattern_Size:             
             Geno_Islands=Geno_Islands+1 
        
    if Geno_Islands>10:
      if Position<len(Pattern_N):
        #print('Position'+str(Position))
        if Indx == Pattern_N[Position]:
           Wanted_PFs_WithIndx.write('***************************\n')
           Position=Position+1
  
      Wanted_PFs.write(str(Pattern_Size))
      Wanted_PFs_WithIndx.write(str(Pattern_Size))
      for i in Pattern:
         Wanted_PFs.write(' '+i) 
         Wanted_PFs_WithIndx.write(' '+i)             
      Wanted_PFs.write(' '+str(Geno_Islands)+'\n') 
      Wanted_PFs_WithIndx.write(' '+str(Geno_Islands)+'\n')    

print('Done Apriori AllCombinations Occurance')    




 









        
