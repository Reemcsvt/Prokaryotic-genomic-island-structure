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
from io import StringIO
import os
import Bio
from Bio import SeqIO
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, Scorer
from Bio.Phylo.NewickIO import Parser
import matplotlib.pyplot as plt
import itertools
from itertools import combinations
from difflib import SequenceMatcher

Path='/ExpressoDB/'
Filtering=Path+'Patterns_Filtering_Information/'  
Patterns_Proteins_Path=Path+'Patterns_Blast_Proteins_Information/'
Output=Path+'Viruses/Output/'
####################################################################################################################
Patterns_Proteins=open(Patterns_Proteins_Path+"Patterns_Proteins_New_5.txt", 'w')
Patterns_Proteins_Unique=open(Patterns_Proteins_Path+"Patterns_Proteins_Unique_New_5.txt", 'w')
Patterns_Proteins_GIs=open(Patterns_Proteins_Path+"Patterns_Proteins_GIs_New_5.txt", 'w')
Patterns_Proteins_Blast=open(Patterns_Proteins_Path+"Patterns_Proteins_Blast_New_5.txt", 'w')
RatioFile=open(Patterns_Proteins_Path+"Ratio_New_5.txt", 'w')
####################################################################################################################
File_Index=0
with open(Filtering+'Patterns_Filtering_Information_New_5.txt', 'r') as f: #Test_.txt  #Patterns_Filtering_Information_GIs_Sorted_New_5.txt
 for line_Pattern in f:
    #Terminase_1&Phage_capsid&Phage_portal& Bacillus_phage_vB_BtS_BMBtp16 KT372714.1 2 NZ_CP015176.1,NZ_CP016163.1 NZ_CP015176.1_B_10,NZ_CP016163.1_B_18 0.0,9.01e-145 4 Firmicutes,Terrabacteria_group,Bacilli,Bacteria,
    File_Index=File_Index+1
    print(File_Index)
    Pattern=line_Pattern.split(None, 1)[0]; PFs=Pattern.split('&'); virus_name=line_Pattern.split(None, 2)[1]; virus=line_Pattern.split(None, 3)[2]
    GIs=line_Pattern.split(None, 6)[5]; GIs_List_=GIs.split(','); GIs_List=[]
    for i in GIs_List_:  GIs_List.append(i+'$')    
    for Pattern_PF in PFs:     
      if Pattern_PF != '': 
         Proteins_PFs=open('GIs_'+Pattern_PF+'_Proteins.txt', 'w')
         Pro=[]; P_GIs=[]; seen1=[]; seen=[]; Flage=0
         Protein=[Pattern_PF+'$']; Wanted_Proteins=[]; Wanted_Proteins_Unique=[]; 
         with open('New_GIsANDtheirProteins_New_5.txt', 'r') as f: 
            for line in f:
               A=line.split(None, 1)[0]
               if '$' not in A:
                 if A+'$' in GIs_List: 
                     GI=A; Flage=1;  seen1=[]
                 else: Flage=0            
               else:
                 if Flage==1:
                   PF=A.split('$')[1]
                   if PF+'$' in Protein:
                     if PF+'$' not in seen1:
                        seen1.append(PF+'$') 
                        P=A.split('|')[0];
                        M=P[P.rfind(GI)+len(GI)+1:len(P)]                          
                        Wanted_Proteins.append(M) 
                        if M+'$' not in seen:
                           Pro.append(M); P_GIs.append(GI+',')
                           seen.append(M+'$') 
                           Wanted_Proteins_Unique.append(M) 
                           Proteins_PFs.write(GI+' '+A+'\n')
                        else:
                           P_GIs[Pro.index(M)]+=GI+','   
                        Flage=0                  
         Proteins_PFs.close()
         os.remove('GIs_'+Pattern_PF+'_Proteins.txt')
         #################################
         #Now I have List of all GIs that belong to the protein family Pattern_PF             
         for k in Pro:
           Patterns_Proteins_GIs.write(str(File_Index)+' '+Pattern+' '+line_Pattern.split(None,2)[1]+' '+virus+' '+Pattern_PF+' '+k+' '+(P_GIs[Pro.index(k)])[0:-1]+'\n')
         ##################################
         for j in range(0,6):
            Patterns_Proteins.write(line_Pattern.split(None, j+1)[j]+' ')
            Patterns_Proteins_Unique.write(line_Pattern.split(None, j+1)[j]+' ')
         #All
         for j in range(0,len(Wanted_Proteins)-1):  Patterns_Proteins.write(Wanted_Proteins[j]+',')
         #Unique
         for j in range(0,len(Wanted_Proteins_Unique)-1): Patterns_Proteins_Unique.write(Wanted_Proteins_Unique[j]+','); 
         ###   
         Patterns_Proteins.write(Wanted_Proteins[len(Wanted_Proteins)-1]+' ')   
         Patterns_Proteins.write(line_Pattern.split(None, 7)[6]+' '); Patterns_Proteins.write(line_Pattern.split(None, 8)[7]+' ')
         Patterns_Proteins.write(line_Pattern.split(None, 9)[8]+'\n')
         ###
         Patterns_Proteins_Unique.write(Wanted_Proteins_Unique[len(Wanted_Proteins_Unique)-1]+' ')   
         Patterns_Proteins_Unique.write(line_Pattern.split(None, 7)[6]+' '); Patterns_Proteins_Unique.write(line_Pattern.split(None, 8)[7]+' ')
         Patterns_Proteins_Unique.write(line_Pattern.split(None, 9)[8]+'\n')
         #################################
         #Get Viral protein
         for seq_record in SeqIO.parse(Output+virus+'_Output_New.txt', "fasta"): 
              if Pattern_PF+'$' in seq_record.description:
                 SD=seq_record.description;
                 virus_protein=SD[SD.rfind('$'+virus+'$')+len('$'+virus+'$'):SD.find(' |')]
                 break
         #################################         
         Wanted_Proteins_Unique.append(virus_protein)
         result=list(combinations(Wanted_Proteins_Unique, 2))
         #print(list(result))
         
         #print(str(File_Index)+' '+Pattern+' '+virus+' '+Pattern_PF)
         
         Remo1=[]; Remo2=[]
         File1=open('File1.fa', 'w');  File2=open('File2.fa', 'w');
         IDx=-1
         for FL in Wanted_Proteins_Unique:                
           IDx+=1
           File1.write(FL+'\n'); Remo1.append(FL)
           if IDx<len(Wanted_Proteins_Unique)-1:
              File2.write(FL+'\n'); Remo2.append(FL)
         File1.close(); File2.close()
         
         S=('blastp -query File1.fa -subject File2.fa -out Blast_Result.txt -outfmt "6 qseqid sseqid  pident evalue  staxid ssciname"')
         os.system(S)
                                                                              
         with open('Blast_Result.txt', 'r') as f: 
            for line in f:
               Protein1=line.split(None, 1)[0]; Protein2=line.split(None, 2)[1];
               FL_Q=Protein1.split('|');   FL_Sub=Protein2.split('|');
               Iden=line.split(None, 3)[2]; Evalue=line.split(None, 4)[3];  taxID=line.split(None, 5)[4]; 
               ssciname=line[line.rfind(taxID)+len(taxID):line.find('\n')]
               if FL_Q[len(FL_Q)-2] != FL_Sub[len(FL_Sub)-2]:
                  Patterns_Proteins_Blast.write(str(File_Index)+' '+Pattern+' '+virus_name+' '+virus+' '+Pattern_PF+' '+FL_Q[len(FL_Q)-2]+' '+FL_Sub[len(FL_Sub)-2]+' '+Iden+' '+Evalue+' '+taxID+' '+(ssciname.lstrip()).replace(" ", "_")+'\n')
                  
                  if virus_protein == FL_Q[len(FL_Q)-2]: # Query and virus
                     Ratio=SequenceMatcher(None, virus_name, ssciname.lstrip()).ratio();                      
                     if Ratio == 0:
                        RatioFile.write(str(File_Index)+' '+Pattern+' '+virus_name+' '+virus+' '+Pattern_PF+' '+FL_Q[len(FL_Q)-2]+' '+FL_Sub[len(FL_Sub)-2]+' '+Iden+' '+Evalue+' '+taxID+' '+(ssciname.lstrip()).replace(" ", "_")+' '+str(Ratio)+'\n')
         os.remove('File1.fa'); os.remove('File2.fa'); os.remove('Blast_Result.txt');
         
         #################################
 
