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

from io import StringIO
import Bio
from Bio import SeqIO
import sys
from Bio import Entrez
Entrez.email = 'reem.aldaihani@gmail.com'
from difflib import SequenceMatcher
###################################################
Path='/ExpressoDB/'
Proteins_Blast=Path+'Patterns_Blast_Proteins_Information/'
HGT_Path=Path+'Patterns_HGT_Species_Information/'
Output=Path+'Viruses/Output/'
####################################################################################################
GIs_ID=[]; GIs_Pattern=[]; GIs_Virus=[]; GIs_PF=[]; GIs_Protein=[]; GIs_GIs=[]
with open(Proteins_Blast+'Patterns_Proteins_GIs_New_5.txt', 'r') as f: #GIs_test.txt   Patterns_Proteins_GIs_New_5.txt
  for line in f:
#1052 Terminase_1&Phage_capsid&Phage_portal& Pseudomonas_phage_PAJU2 AP009624.1 Terminase_1 WP_048940025.1 NZ_CP012077.1_B_13,NZ_CP024172.1_B_18
      GIs_ID.append(line.split(None, 1)[0]);    GIs_Pattern.append(line.split(None, 2)[1]);    GIs_Virus.append(line.split(None, 4)[3]);
      GIs_PF.append(line.split(None, 5)[4]);    GIs_Protein.append(line.split(None, 6)[5]);
      GIs_=line.split(None, 7)[6];  GIs_GIs.append(GIs_.split(','))
####################################################################################################
Lineage_Acc=[]; Lineage_TaxId=[]; Lineage_TaxId_Parent=[]
with open('Bacteria_TaxIDs_Lineages_New_5.txt', 'r') as f: 
  for line in f:
#NC_000853.1 243274 123 Bacteria,Thermotogae,Thermotogae,Thermotogales,Thermotogaceae,Thermotoga,Thermotoga_maritima,Thermotoga_maritima_MSB8
     Lineage_Acc.append(line.split(None, 1)[0]); Lineage_TaxId.append(line.split(None, 2)[1]); Lineage_TaxId_Parent.append(line.split(None, 3)[2]); 
####################################################################################################
HGT=open(HGT_Path+'Patterns_HGT_Virus_New_5.txt', 'w')
index=0; 
with open(Proteins_Blast+'Patterns_Proteins_Blast_Update_New_5.txt', 'r') as f:  #Test_  Patterns_Proteins_Blast_Update_New_5.txt
  for line in f:
#1055 Terminase_1&Phage_capsid&Phage_portal& uncultured_Caudovirales_phage MF417925.1 Phage_capsid ASN71451.1 WP_063676008.1 92.396 0.0 1428 Bacillus _thuringiensis
      index+=1; print(index)
      a=line.split(None, 1)[0]+line.split(None, 5)[4] #row index and protein family
      HGT_Index=line.split(None, 1)[0];  HGT_Pattern=line.split(None, 2)[1]; HGT_PF=line.split(None, 5)[4]; 
      virus=line.split(None, 4)[3]; virus_Text=line.split(None, 3)[2];  taxID=line.split(None, 10)[9];
      if index ==1:     
         HGT_line=[]; HGT_Protein1=[]; HGT_Protein2=[]; HGT_Iden=[]; HGT_Evalue=[]; HGT_Tax=[]; HGT_Ssname=[];    
         Prev=a; Prev_virus=virus; Prev_virus_Text=virus_Text; Prev_index=HGT_Index; Prev_pattern=HGT_Pattern; Prev_PF=HGT_PF
         HGT_Protein1.append(line.split(None, 6)[5]); HGT_Protein2.append(line.split(None, 7)[6])
         HGT_Iden.append(line.split(None, 8)[7]); HGT_Evalue.append(line.split(None, 9)[8]); 
         HGT_Tax.append(taxID); HGT_Ssname.append(line.split(None, 11)[10])   
         HGT_line.append(line);      
      else:
         if a == Prev:
            HGT_Protein1.append(line.split(None, 6)[5]); HGT_Protein2.append(line.split(None, 7)[6])
            HGT_Iden.append(line.split(None, 8)[7]); HGT_Evalue.append(line.split(None, 9)[8]); 
            HGT_Tax.append(taxID); HGT_Ssname.append(line.split(None, 11)[10])   
            HGT_line.append(line);      
         else:
            #Save the group  #Get Viral protein            
            for seq_record in SeqIO.parse(Output+Prev_virus+'_Output_New.txt', "fasta"): 
                if Prev_PF+'$' in seq_record.description:
                   SD=seq_record.description;
                   virus_protein=SD[SD.rfind('$'+Prev_virus+'$')+len('$'+Prev_virus+'$'):SD.find(' |')]
                   break          
            ###################        
            IDX=-1
            for i in HGT_Protein1:
              IDX+=1
              if i == virus_protein:
                if float(HGT_Iden[IDX]) >= 85:  
                  Ratio=SequenceMatcher(None, Prev_virus_Text, HGT_Ssname[IDX][1:]).ratio();
                  if Ratio<0.4: 
                      Tax_ID=HGT_Tax[IDX] #taxid of second protein
                      #Get GIs of this second protein                  
                      Protein_GIs=''; GIs_Index=-1; Flage=0
                      for j in GIs_ID: #Protein and their GIS
                         GIs_Index+=1
                         if j == Prev_index and GIs_Pattern[GIs_Index]==Prev_pattern and GIs_PF[GIs_Index]==Prev_PF and GIs_Protein[GIs_Index]==HGT_Protein2[IDX] and GIs_Virus[GIs_Index]==Prev_virus:
                              k_index=-1; GI_seen=[]
                              for k in GIs_GIs[GIs_Index]:
                                 k_index+=1
                                 GI=k.split('_B_')[0];
                                 GI_TaxID=Lineage_TaxId[Lineage_Acc.index(GI)]
                                 GI_TaxID_Parent=Lineage_TaxId_Parent[Lineage_Acc.index(GI)]
                                 if (Tax_ID==GI_TaxID or Tax_ID==GI_TaxID_Parent) and (GI not in GI_seen):    
                                    Protein_GIs=Protein_GIs+k+','                    
                                    Flage=1; GI_seen.append(GI)
                         if Flage==1: 
                            HGT.write(Prev_index+' '+Prev_pattern+' '+Prev_virus_Text+' '+Prev_virus+' '+Prev_PF+' '+HGT_Protein1[IDX]+' '+HGT_Protein2[IDX]+' '+HGT_Iden[IDX]+' '+HGT_Evalue[IDX]+' '+HGT_Tax[IDX]+' '+HGT_Ssname[IDX][1:]+' '+Protein_GIs+' '+str(Ratio)+'\n')
                            break         
            
            HGT_line=[]; HGT_Protein1=[]; HGT_Protein2=[]; HGT_Iden=[]; HGT_Evalue=[]; HGT_Tax=[]; HGT_Ssname=[];
            Prev=a; Prev_virus=virus; Prev_virus_Text=virus_Text; Prev_index=HGT_Index; Prev_pattern=HGT_Pattern; Prev_PF=HGT_PF
            HGT_line.append(line); HGT_Protein1.append(line.split(None, 6)[5]); HGT_Protein2.append(line.split(None, 7)[6])
            HGT_Iden.append(line.split(None, 8)[7]); HGT_Evalue.append(line.split(None, 9)[8]); 
            HGT_Tax.append(taxID); HGT_Ssname.append(line.split(None, 11)[10])        
  
######################################################################
#Save the group  #Get Viral protein            
for seq_record in SeqIO.parse(Output+Prev_virus+'_Output_New.txt', "fasta"): 
    if Prev_PF+'$' in seq_record.description:
       SD=seq_record.description;
       virus_protein=SD[SD.rfind('$'+Prev_virus+'$')+len('$'+Prev_virus+'$'):SD.find(' |')]
       break          
###################        
IDX=-1
for i in HGT_Protein1:
  IDX+=1
  if i == virus_protein:
    if float(HGT_Iden[IDX]) >= 85:  
      Ratio=SequenceMatcher(None, Prev_virus_Text, HGT_Ssname[IDX][1:]).ratio();
      if Ratio<0.4: 
          Tax_ID=HGT_Tax[IDX] #taxid of second protein
          #Get GIs of this second protein                  
          Protein_GIs=''; GIs_Index=-1; Flage=0
          for j in GIs_ID: #Protein and their GIS
             GIs_Index+=1
             if j == Prev_index and GIs_Pattern[GIs_Index]==Prev_pattern and GIs_PF[GIs_Index]==Prev_PF and GIs_Protein[GIs_Index]==HGT_Protein2[IDX] and GIs_Virus[GIs_Index]==Prev_virus:
                  k_index=-1; GI_seen=[]
                  for k in GIs_GIs[GIs_Index]:
                     k_index+=1
                     GI=k.split('_B_')[0];
                     GI_TaxID=Lineage_TaxId[Lineage_Acc.index(GI)]
                     GI_TaxID_Parent=Lineage_TaxId_Parent[Lineage_Acc.index(GI)]
                     if (Tax_ID==GI_TaxID or Tax_ID==GI_TaxID_Parent) and (GI not in GI_seen):    
                        Protein_GIs=Protein_GIs+k+','                    
                        Flage=1; GI_seen.append(GI)
             if Flage==1: 
                HGT.write(Prev_index+' '+Prev_pattern+' '+Prev_virus_Text+' '+Prev_virus+' '+Prev_PF+' '+HGT_Protein1[IDX]+' '+HGT_Protein2[IDX]+' '+HGT_Iden[IDX]+' '+HGT_Evalue[IDX]+' '+HGT_Tax[IDX]+' '+HGT_Ssname[IDX][1:]+' '+Protein_GIs+' '+str(Ratio)+'\n')
                break         

HGT.close()
###########################################################################################################################
###########################################################################################################################
####################################Get virus rows that have all Pattern PFs###############################################
###########################################################################################################################
###########################################################################################################################

from collections import Counter
HGT_PFs=open(HGT_Path+'Patterns_HGT_Virus_PFs_New_5.txt', 'w'); index=0;
with open(HGT_Path+'Patterns_HGT_Virus_New_5.txt', 'r') as f:  #Test_HGT_Virus.txt  Patterns_HGT_Virus_New_5.txt
  for line in f:
#1 Phage_capsid&Phage_portal&Terminase_1& Erysipelothrix_phage_phi1605 MF172979.1 Phage_portal ASD51098.1 WP_041721144.1 98.082 0.0 208226 Alkaliphilus_metalliredigens NC_009633.1_B_19, 0.14285714285714285
      index+=1; print(index)
      a=line.split(None, 1)[0]    
      if index==1:
         prev_index=a; HGTPFs_line=[]; HGTPFs_PF=[];  Prev_Pattern_len=len((line.split(None, 2)[1]).split('&'))-1   
         HGTPFs_line.append(line); HGTPFs_PF.append(line.split(None, 5)[4]); 
      else:
         if a == prev_index:
            HGTPFs_line.append(line); HGTPFs_PF.append(line.split(None, 5)[4]); 
         else:
            C=len(Counter(HGTPFs_PF));  
            if C==Prev_Pattern_len:
               for i in HGTPFs_line:
                  HGT_PFs.write(i) 
            ###
            prev_index=a; HGTPFs_line=[]; HGTPFs_PF=[]; Prev_Pattern_len=len((line.split(None, 2)[1]).split('&'))-1   
            HGTPFs_line.append(line); HGTPFs_PF.append(line.split(None, 5)[4]); 
########################               
C=len(Counter(HGTPFs_PF));  
if C==Prev_Pattern_len:
   for i in HGTPFs_line:
      HGT_PFs.write(i) 



    