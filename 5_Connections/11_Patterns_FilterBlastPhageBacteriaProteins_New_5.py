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

from Bio import Entrez
Entrez.email = 'reem.aldaihani@gmail.com'
from difflib import SequenceMatcher
###################################################
Path='/ExpressoDB/'
Proteins_Blast=Path+'Patterns_Blast_Proteins_Information/'
####################################################
#Get only Species Rows
Proteins_Blast_Update=open(Proteins_Blast+'Patterns_Proteins_Blast_Update_New_5.txt', 'w')
#Ranks=open(Proteins_Blast+'Ranks.txt', 'w');  Error1=open('Erorero1.txt', 'w'); Error2=open('Erorero2.txt', 'w')
index=0; seen1=[]; seen2=[]
TaxList1=['species group$','genus$', 'family$', 'clade$', 'order$', 'class$', 'phylum$', 'kingdom$', 'superkingdom$','$']
TaxList2=['species$','sub species$', 'strain$', 'isolate$']
with open(Proteins_Blast+'Patterns_Proteins_Blast_New_5.txt', 'r') as f: #Blast_Test.txt  #Patterns_Proteins_Blast_New_5.txt
  for line in f:
#1055 Terminase_1&Phage_capsid&Phage_portal& uncultured_Caudovirales_phage MF417925.1 Phage_capsid ASN71451.1 WP_063676008.1 92.396 0.0 1428 Bacillus thuringiensis
     index+=1; print(index)
     Tax=line.split(None, 10)[9]  
     if Tax!='N/A':    
         if Tax+'$' not in seen1:           
           search = Entrez.esummary(id = Tax, db = "taxonomy", retmode = "xml")
           a=Entrez.read(search)
           Rank=a[0]['Rank']
           seen1.append(Tax+'$'); seen2.append(Rank)
           #Ranks.write(Tax+' '+Rank+'\n')
         else:
           Rank=seen2[seen1.index(Tax+'$')]
    
         if (Rank+'$' not in TaxList1) and  (Rank+'$' in TaxList2):
             Proteins_Blast_Update.write(line)             
             #Error1.write(Rank+'\n') 
         #else:
             #Error2.write(Rank+'\n')    
Proteins_Blast_Update.close() 
print('Done Ranking')