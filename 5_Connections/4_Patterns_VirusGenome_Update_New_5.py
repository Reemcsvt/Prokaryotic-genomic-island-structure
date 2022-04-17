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

import Bio
from Bio import SeqIO
from Bio import Entrez
Entrez.email = '1##############'
from urllib.error import HTTPError

Source='/Patterns_Table_Information/'
with open('Patterns_Update_Information_New_5.txt', 'r') as f: #
  count=0
  for line_Pattern in f:
      count=count+1;  print('Pattern'+' '+str(count)); dataset=[];  dataset.append([]);  
      PF=int(line_Pattern.split(None, 1)[0])
      for i in range (1,PF+1):
            dataset[0].append(line_Pattern.split(None, i+1)[i])                      
      Pattern_Print=''; Pattern=''; 
      for p in dataset[0]: Pattern_Print=Pattern_Print+p+'_';   Pattern=Pattern+p+'&';   
      Pattern_File=str(PF)+"_"+line_Pattern.split(None, PF+2)[PF+1]+"_"+Pattern_Print+'Local'
      Patterns_Information_New=Source+'PatternInfo_'+Pattern_File+'.txt'
      #######################################################################################
      Source1='/Patterns_Update_Information/'
      Patterns_Update_New_5=open(Source1+"Patterns_Update_Information_"+Pattern_File+".txt", "w")              
      V_G=[]; Bacteria=[]; Evalue=[]
      with open(Patterns_Information_New, 'r') as f: 
          for line in f:
            #a&b& NZ_CP021326.1 NZ_CP021326.1_B_17 Acinetobacter_baumannii Moraxella_phage_Mcat7 1 AKI27325.1 KR093631.1 3.96e-99 
            
            A1=line.split(None, 1)[0];  A2=line.split(None, 2)[1]; A3=line.split(None, 3)[2]; A4=line.split(None, 4)[3];
            A5=line.split(None, 5)[4];  A6=line.split(None, 6)[5]; A8=line.split(None, 7)[6]; A9=line.split(None, 9)[8];
            LINE1=A1+' '+A2+' '+A3+' '+A4+' '+A5+' '+A6+' '+A8+' '
            Genome=line.split(None, 8)[7]
            if 'Error' not in Genome:
               Patterns_Update_New_5.write(line)
#            else:
#              ProteinAcc=A8
#              try:
#                  handle = Entrez.efetch(db="protein",id=ProteinAcc , retmode="xml",retmax=1000000)#, rettype="gb", retmode="text") 
#                  record = Entrez.read(handle);  recordi_index=-1
#                  for recordi in record:
#                      recordi_index=recordi_index+1
#                      ORG=recordi['GBSeq_source-db']
#                      if 'accession' in ORG:
#                          GGG=ORG[ORG.rfind('accession')+len('accession')+1:len(ORG)]   
#                          if len(GGG)>10:    
#                             A=ORG.split('accession ')[1]
#                             Patterns_Update_New_5.write(LINE1+' '+A.split(';')[0]+' '+A9+'\n')####
#                             print('Good1')
#                          else:   
#                             Patterns_Update_New_5.write(LINE1+' '+GGG+' '+A9+'\n')##### 
#                             print('Good2')
#                      elif 'locus' in  ORG:
#                          GGG=ORG[ORG.rfind('locus')+len('locus')+1:len(ORG)]    
#                          if len(GGG)>10:    
#                             A=ORG.split('locus ')[1]
#                             Patterns_Update_New_5.write(LINE1+' '+A.split(';')[0]+' '+A9+'\n') ##### 
#                             print('Good3')
#                          else:   
#                             Patterns_Update_New_5.write(LINE1+' '+GGG+' '+A9+'\n') ####   
#                             print('Good4')              
#              except HTTPError as err: print('Error1')
#              except ValueError: print('Error2')
#              except KeyError as e: print('Error3')













                  