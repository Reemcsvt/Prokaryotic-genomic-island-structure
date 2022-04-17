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


ProS=[]
ProF=[]
lines_seen = []
E_values = []
EVs=[]
ProteinFamilies=[]
ProteinSeqS=[]

with open('HmmscanResults.txt', 'r') as f:   #HmmscanResults.txt  HmmS.txt
    for line in f:
        ProteinSeqS.append(line.split(None, 1)[0])
        EVs.append(float(line.split(None, 5)[4]))
        ProteinFamilies.append(line.split(None, 3)[2]) 
print('Done read HmmscanResults')        
    
INDEX=-1
for j in ProteinSeqS:   
        INDEX=INDEX+1
        print('ProteinSeq = '+str(INDEX))
         
        ProteinSeq = j
        EV = EVs[INDEX]

        if ProteinSeq in lines_seen: 
              In=lines_seen.index(ProteinSeq)
              E=E_values[In]
              if E>EV:
                 E_values[In]=EV
                 ProteinFamily = ProteinFamilies[INDEX]         
                 ProS[In]=ProteinSeq
                 ProF[In]=ProteinFamily
        else:
              lines_seen.append(ProteinSeq)  
              E_values.append(EV)            
              ProteinFamily = ProteinFamilies[INDEX]         
              ProS.append(ProteinSeq)
              ProF.append(ProteinFamily)
f.close()   
print('Done save HmmscanResults')
#########################################################################
HmmscanResults_New=open('HmmscanResults_New.txt', 'w')
s=0
for i in ProS:
   HmmscanResults_New.write(i+' '+ProF[s]+'\n')
   s=s+1

#########################################################################
Totals=[]
Ribosomals=[]
GIs=[]

import os
from Bio import SeqIO
os.listdir(path='.')
def getListOfFiles(dirName):
    listOfFile = os.listdir(dirName)
    allFiles = list()
    for entry in listOfFile: 
        fullPath = os.path.join(dirName, entry)
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(entry)
        else:
            if entry[-2:]=='fa':
               allFiles.append(entry)                
    return allFiles
dirName = '/All_Proteins';
listOfFiles1 = getListOfFiles(dirName)


#listOfFiles1=['NZ_CP022699.1_B_8.fa', 'NC_016894.1_B_12.fa', 'NZ_CP022699.1_B_12.fa' ,'NZ_CP022699.1_B_6.fa', 'NZ_CP022699.1_B_9.fa']

for i in listOfFiles1:
   i=i[:-3]
   GIs.append(i)    
print('Done read GIs')     

########################################################################

print(len(GIs))

for i in range(0,len(GIs)+1):
    Totals.append(0)
    Ribosomals.append(0)  

count=-1
for i in ProS:
        count=count+1
        s=i
        l=len(s.split(')_'))-1
        GI=s.split(')_')[l]
        Index_GI=GIs.index(GI)
        
        if 'Ribosomal' in ProF[count]:
           Ribosomals[Index_GI]=Ribosomals[Index_GI]+1
        Totals[Index_GI]=Totals[Index_GI]+1   
        
########################################################################        
GIs_PFs_Ribosomals=open('GIs_PFs_Ribosomals.txt','w') 
count=-1
for i in GIs:
   count=count+1
   GIs_PFs_Ribosomals.write(i+' '+str(Totals[count])+' '+str(Ribosomals[count])+'\n')   
   
########################################################################  
GIs_WRibosomals=open('GIs_WRibosomals_G50.txt','w')  
count=-1
for i in GIs:
   count=count+1
   if Totals[count] !=0:
     if(Ribosomals[count]/Totals[count])*100  > 50:
        GIs_WRibosomals.write(i+'\n')  
   
 
 
 
 
 
             
