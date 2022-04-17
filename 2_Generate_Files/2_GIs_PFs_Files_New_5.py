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

GIsWR=[]
with open('GIs_WRibosomals_G50.txt', 'r') as f: #R.txt  GIs_WRibosomals_G50.txt
    for line in f:
      GIsWR.append(line.split(None, 1)[0]) 
#GIsWR=['NZ_CP022699.1_B_8', 'NC_016894.1_B_12']
#############################################################################################
#                             Get all proteins families from the file                       #
#############################################################################################
ProFamily=[]
from collections import Counter
wordlist=[]

with open('HmmscanResults_New.txt', 'r') as f:   #HmmscanResults_New.txt  teto.txt
    for line in f:
        
        ProteinSeq=line.split(None, 1)[0]
        ProteinFamily =line.split(None, 2)[1]
        
        l=len(ProteinSeq.split(')_'))-1
        GI=ProteinSeq.split(')_')[l] 
        
        if GI not in GIsWR:
          wordlist.append(ProteinFamily)
          ProFamily.append(GI+'_'+ProteinSeq+'$'+ProteinFamily)

print('Done read HmmscanResults')        
f.close()   

#############################################################################################
#                              Save protein families in a list                              #
#############################################################################################
PFL=open('New_ProteinFamilies_5.txt','w')   
class geeks:  
    def __init__(self, name, roll):  
        self.name = name  
        self.roll = roll 
ProList=[]
List_Pro_Freq=Counter(wordlist)
for key,value in sorted(List_Pro_Freq.items()):
    ProList.append( geeks(key, value) ) 
for i in ProList:  
   PFL.write(i.name+'\n')   
##############################################################################################
##                                 Get the names of the GIs                                  #
##############################################################################################
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

#listOfFiles1=['NZ_CP022699.1_B_8.fa', 'NC_016894.1_B_12.fa', 'NZ_CP022699.1_B_12.fa','NZ_CP022699.1_B_6.fa', 'NZ_CP022699.1_B_9.fa']

listOfFiles=[]
for i in listOfFiles1:
   visit=0
   i=i[:-3]
   for j in GIsWR:
     if i == j:
       visit=1
   if visit==0:   
      listOfFiles.append(i)    
print('Done read GIs')      
############################################################################################        
GIsANDProteins=open("New_GIsANDtheirProteins_5.txt", "w") 
GIsANDProFreq=open("New_GIsANDProFreq_5.txt", "w")
Er=open("Er.txt", "w")
class Families:  
    def __init__(self):
        self.GI=0
        self.PFs = []
    def addPF(self, PF):
        self.PFs.append(PF)                
        
ProFamList = [] 
ProFamListNew = [] 
Proteins_seen = set() 
GIII=0
countt=1
GIsFreq=[]
ProFamily.sort()
print('length='+str(len(ProFamily)))
print('length='+str(len(listOfFiles)))
for i in listOfFiles:  
          #i=i[:-3]   
          print(str(GIII)+' '+i) 
          GIII=GIII+1           
          ProFamList1 = Families() 
          ProFamList2 = Families() 
          ProFamList1.GI=i
          ProFamList2.GI=i
          GIsANDProteins.write(ProFamList1.GI+'\n')
          inc=0
          count=0
          visit=0
          for k in ProFamily:
              if i+'$' in k:
                 visit=1
                 ProFamList1.addPF(k) #Add Seq and PF
                 ProFamList2.addPF(k.split('$')[1])    #(ProFamilyP[count])  #Add PF
                 if len(k.split('$'))>2:
                    Er.write(k+'\n')
                 GIsANDProteins.write(ProFamList1.PFs[inc]+'\n')
                 inc=inc+1
              else:
                 if visit==1:
                    break
              count=count+1
          
          GIsANDProFreq.write(ProFamList1.GI+' '+str(len(ProFamList1.PFs))+'\n')  
          GIsFreq.append(len(ProFamList1.PFs))          
          ProFamList.append(ProFamList1) 
          ProFamListNew.append(ProFamList2)                                 
print('Done build the FamList')  #Done get the Protein Sequences for each GI

#print(GIsFreq)
c=0
for i in GIsFreq:
  if i ==0:
    c=c+1
    

#############################################################################################
#                                 Build the Bimax matrix                                    #
#############################################################################################
Bimax=open("New_InputMatrix_5.txt", "w")
GIs=open("New_GenomicIslands_5.txt", "w")  
Number_Of_GIs=len(listOfFiles)
Number_Of_ProteinsFamilies=len(ProList)
print(str(Number_Of_GIs)+'\n')
print(str(Number_Of_ProteinsFamilies)+'\n')

Bimax.write(str(Number_Of_GIs-c)+' '+str(Number_Of_ProteinsFamilies)+' '+'0 0'+'\n')

Row=1
countt=1
indx=-1
for i in ProFamListNew: # each GI
  indx=indx+1
  if GIsFreq[indx]>0:
          GIs.write(i.GI+'\n')
          print('Row='+str(Row))
          Row=Row+1
          for j in ProList:
                   Has_PF=0 
                   for k in i.PFs:
                        if j.name == k:
                           Bimax.write(' 1')
                           Has_PF=1
                           break
                   if Has_PF==0:
                        Bimax.write(' 0')  
          Bimax.write('\n')           



          