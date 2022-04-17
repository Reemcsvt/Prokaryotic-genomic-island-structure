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
PFs1=[]
with open('New_ProteinFamilies_5.txt', 'r') as f: # New_ProteinFamilies_5.txt
    for line in f:
      PFs1.append(line.split(None, 1)[0]) 
#PFs1=['XFP_N', 'DUF4277', 'NinG', 'DUF2528', 'Phage_CII', 'Phage_Nu1', 'Phage_antitermQ']
print('Done read protein families')
################################################################################
GIs=[]
with open('New_GenomicIslands_New_5.txt', 'r') as f: # New_GenomicIslands_New_5.txt
    for line in f:
      GIs.append(line.split(None, 1)[0]) 
#GIs=['NZ_CP022699.1_B_8', 'NZ_CP022699.1_B_6', 'NC_016894.1_B_12', 'NC_014248.1_B_16', 'NC_009925.1_B_13', 'NZ_CP011389.1_B_8']
print('Done read protein families')
################################################################################
#Read Matrix
Matrix=[]
count=0
with open('New_InputMatrix_New_5.txt', 'r') as f: #New_InputMatrix_New_5.txt  a.txt
  for line in f:    
     if count>0:
       print('Count= '+str(count))
       PF=0
       Matrix.append([])
       for i in line:
         if i=='1':
            Matrix[count-1].append(PFs1[PF])
         
         if i=='1' or i=='0':   
            PF=PF+1          
     count=count+1
print('Done read proteins from matrix')         
################################################################################
dataset=[]
count=0
with open('Apriori_New_5.txt', 'r') as f: #Apriori_New_5.txt  Ap.txt
  for line in f:
      dataset.append([])
      PF=int(line.split(None, 2)[1])
      for i in range (2,PF+2):
            dataset[count].append(line.split(None, i+1)[i])              
      count=count+1
print('Done read data set')    
################################################################################
Genus=[]
Species=[]
count=0
#g=[]
for i in dataset:
   indx=-1
   visit=0
   Species.append([])
   Genus.append([])
   #g.append([])
   for j in Matrix:
      indx=indx+1      
      if len(set(i).intersection(set(j)))==len(i):
        visit=1
        G=GIs[indx].split('_B_')[0] 
        Species[count].append(Spe[AccN.index(G)])
        Genus[count].append(Gen[AccN.index(G)])
        #g[count].append(GIs[indx])
        
   if visit==1:
      count=count+1
print('Done finding species')              
###################################################################################        
L=[]
Sp=[]
count=0
for i in Species:
    Sp.append([])    
    S=Counter(i)
    L.append(len(S))
    for key,value in sorted(S.items()):
        #print(key+' '+str(value)) 
        Sp[count].append(key)
    count=count+1       
###################################################################################
Ln=[]
Gens=[]
count=0
for i in Genus:
    Gens.append([])    
    S=Counter(i)
    Ln.append(len(S))
    for key,value in sorted(S.items()):
        #print(key+' '+str(value)) 
        Gens[count].append(key)
    count=count+1       
###################################################################################
AS=open('Apri_Genus_Species_New_5.txt','w')  
count=0 
with open('Apriori_New_5.txt', 'r') as f: #Apriori_New_5.txt  Ap.txt
  for line in f:
     AS.write(line.rstrip()+' '+str(L[count])+' ')
     for i in Sp[count]:
        AS.write(i+'&')
        
     AS.write(' '+str(Ln[count])+' ') 
     for i in Gens[count]:
        AS.write(i+'&')     
       
     AS.write('\n')   
     count=count+1
   











                   