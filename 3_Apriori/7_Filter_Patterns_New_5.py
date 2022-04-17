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

Genus=['Bacillus', 'Streptococcus', 'Xanthomonas', 'Streptomyces', 'Bordetella', 'Klebsiella', 'Burkholderia', 'Pseudomonas', 'Salmonella', 'Escherichia']
Order=['Bacillales', 'Lactobacillales', 'Xanthomonadales', 'Actinomycetales', 'Burkholderiales', 'Enterobacterales', 'Burkholderiales', 'Pseudomonadales' , 'Enterobacterales' , 'Enterobacterales']

#Genus=[]
#Order=[]
#with open('Order_Genus_New_5.txt', 'r') as f: #
#    for line in f:
#       Order.append(line.split(None, 2)[1])
#       Genus.append(line.split(None, 3)[2])
#################################################################################
SpeG=[]
AccN=[]
with open('Genus_Species_AccN_New_3.txt', 'r') as f: #Species names and accession number
    for line in f:
       SpeG.append(line.split(None, 1)[0])
       AccN.append(line.split(None, 3)[2])
print('Done read species')
#################################################################################
Patterns=[]
Pattern_Occ=[]
count=-1
with open('Apriori_AllCombinations_Occ_New_5.txt', 'r') as f: #Apriori_AllCombinations_Occ_New_5.txt  Ap_AllComp.txt
  for line in f:
      PFs=int(line.split(None, 1)[0])
      if PFs >2 and PFs<7:
        Occ=int(line.split(None, PFs+2)[PFs+1])
        if Occ>00:
          Patterns.append([])
          count=count+1
          for i in range (1,PFs+1):
                Patterns[count].append(line.split(None, i+1)[i])       
          Pattern_Occ.append(Occ) 
print('Done read Apriori_AllCombinations')                
#################################################################################
GIs=[]
PFs=[]
#PFs_All=[]
Indx=-1
with open('New_GIsANDtheirProteins_New_5.txt', 'r') as f: #New_GIsANDtheirProteins_New_5.txt  GIs_PFs.txt
  for line in f:
    A=line.split(None, 1)[0]
    if '$' not in A:
        GIs.append(A)
        PFs.append([])
        #PFs_All.append([]) 
        Indx=Indx+1
    else:
        PFs[Indx].append(A.split('$')[1])
        #PFs_All[Indx].append(A)
print('Done read GIsANDtheirProteins')        
#################################################################################
FILE_PFs='Wanted_GIs_Patterns_Genus_New_5.txt'
Wanted_PFs=open(FILE_PFs, 'w')
Pattern_GIs=[]
Pattern_GIs_Matrix=[]
Pattern_Genus_Matrix=[]
count=-1
for Pattern in Patterns:
    #Save pattern
    count=count+1
    Pattern_GIs_Matrix.append([])
    Pattern_Genus_Matrix.append([])
    Pattern_GIs=[]
    Geno_Islands=0
    Indx=Indx+1        
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
             Pattern_GIs.append(GIs[PFs_Indx])
             Pattern_GIs_Matrix[count].append(GIs[PFs_Indx])
             G=GIs[PFs_Indx].split('_B_')[0] 
             Pattern_Genus_Matrix[count].append(SpeG[AccN.index(G)])             
             idx=0
    Wanted_PFs.write(str(Pattern_Size)+' '+str(Pattern)+' '+str(Pattern_Occ[count])+'\n')
    Wanted_PFs.write(str(len(Pattern_GIs))+' '+str(Pattern_GIs)+'\n') 
    Wanted_PFs.write(str(len(Pattern_Genus_Matrix[count]))+' '+str(Pattern_Genus_Matrix[count])+'\n')          
Wanted_PFs.close()
print('Done read Wanted')    
##############################################################################
from collections import Counter

class geeks:  
    def __init__(self, name, roll):  
        self.name = name  
        self.roll = roll 
PatternsGenusMatrix=open('Patterns_Genus_Matrix_New_5.txt', 'w')
PatternsGenusMatrix_New_5=open('Patterns_OrderGenus_Matrix_New_5.txt', 'w')

for i in range(0,len(Patterns)):
    Patterns_GM=[]
    Patterns_Count=[]
    for j in Genus:
       W=Pattern_Genus_Matrix[i].count(j) 
       Patterns_Count.append(W)
       if W>0:
          Patterns_GM.append(Order[Genus.index(j)])
    
    ProList=[]
    List_Pro_Freq=Counter(Patterns_GM)
    for key,value in sorted(List_Pro_Freq.items()):
        ProList.append( geeks(key, value) ) 
    
    if len(ProList)>=7: 
       PatternsGenusMatrix.write(str(len(Patterns[i]))+' '+str(Patterns[i])+' ')
       PatternsGenusMatrix_New_5.write(str(len(Patterns[i]))+' '+str(Patterns[i])+'\n')
       for k in Patterns_Count:
           PatternsGenusMatrix.write(str(k)+' ')
       PatternsGenusMatrix.write('\n')     
        
        
       
       
       
       
        
