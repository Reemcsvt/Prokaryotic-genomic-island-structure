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
from Bio import SearchIO
import os
######################################################################################
########################################################################################
Source="/Blast_Viral_Table/"
with open('1_Patterns_New_5.txt', 'r') as f: #Apri_Genus_Species_Order_Patterns_Top100_Print.txt  Top_Pattterns.txt
  count=0
  for line_Pattern in f:
      count=count+1;  print('Pattern'+' '+str(count))
      dataset=[]; dataset_=[];  dataset.append([]);   dataset_.append([])
      PF=int(line_Pattern.split(None, 1)[0])
      for i in range (1,PF+1):
            dataset[0].append(line_Pattern.split(None, i+1)[i])   
            dataset_[0].append(line_Pattern.split(None, i+1)[i]+'$')  

      Pattern_Print=''
      for p in dataset[0]: Pattern_Print=Pattern_Print+p+'_'     
      Pattern_File=str(PF)+"_"+line_Pattern.split(None, PF+2)[PF+1]+"_"+Pattern_Print+'Local'                              
      #Pattern_File=str(PF)+"_"+line_Pattern.split(None, PF+2)[PF+1]+"_Local"
      print('Pattern_File '+Pattern_File)
      ################### Read all GIs and their PFs ##################################
      GIs=[]; PFs=[]; PFs_All=[]; Indx=-1
      with open('New_GIsANDtheirProteins_New_5.txt', 'r') as f: #New_GIsANDtheirProteins_New_5.txt  GIs_PFs.txt
        for line in f:
          A=line.split(None, 1)[0]
          if '$' not in A:
              GIs.append(A); PFs.append([]); PFs_All.append([]); Indx=Indx+1
          else:
              PFs[Indx].append(A.split('$')[1])
              PFs_All[Indx].append(A)             
      #################Read all GIs and their PFs with the current pattern#############
      FILE_PFs=Source+'Wanted_PFs/'+'Wanted_PFs_'+Pattern_File+'.txt'      
      Wanted_PFs=open(FILE_PFs, 'w');    Indx=-1
      for Pattern in dataset:
          #print(Pattern) 
          Geno_Islands=0;  Indx=Indx+1;  Pattern_Size=len(Pattern);  PFs_Indx=-1               
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
                   Wanted_PFs.write(GIs[PFs_Indx]+' '+str(Pattern_Size)+'\n')
                   Geno_Islands=Geno_Islands+1;   idx=0
                   for z in zip(i,PFs_All[PFs_Indx]):
                     if idx == Pattern_Size:
                       break
                     else:  
                      if Pattern[idx] ==z[0]:
                         Wanted_PFs.write(z[1]+'\n');   idx=idx+1                                                          
          #print(str(Pattern)+' '+str(Geno_Islands))            
      Wanted_PFs.close()
      print('Done Reading GIs that have the this pattern')   
      ##############################################################################
      directory = "GIs_Proteins_New_"+Pattern_File
      parent_dir = "/GIs_Proteins/"
      path = os.path.join(parent_dir, directory)        
      if not os.path.exists(path):
         os.mkdir(path)     
      indx=0; B=0
      with open(FILE_PFs, 'r') as f: #Wanted_PFs.txt
        for line in f:
          if B==0:
              A=line.split(None, 1)[0];    B=int(line.split(None, 2)[1])
              FILE=path+"/"+A+".fa"
              GI_Proteins=open(FILE, 'w')
              indx=indx+1; GenomicIsland=A  
          else:
              P=line.split(None, 1)[0]
              WP=P[P.rfind('WP')+len('WP'):P.find('|')]
              M='WP'+WP        
              I1=P[P.rfind('_(')+len('_('):P.rfind('..')]
              I2=P[P.rfind('..')+len('..'):P.rfind(')_')]              
              for seq_record in SeqIO.parse("/All_Proteins/"+GenomicIsland+".fa", "fasta"):                           
                  if  I1 in seq_record.description  and I2 in seq_record.description:  
                     GI_Proteins.write('>'+seq_record.description+'\n')
                     GI_Proteins.write(str(seq_record.seq)+'\n')
                     B=B-1; break      
      ######################################################################################
      # run BLAST         
      GIs_Blast=[]; 
      with open(FILE_PFs, 'r') as f: #Wanted_PFs_New_5.txt W_PFs_New_5.txt
        for line in f:
          A=line.split(None, 1)[0]
          if '$' not in A:
              GIs_Blast.append(A+'.fa')
      counter_GIs=1;  L_GIs=len(GIs_Blast)        
      for i in GIs_Blast:   
          print('Pattern No.= '+str(count)+'  GI='+i+'    '+str(counter_GIs)+' out of '+str(L_GIs)); counter_GIs=counter_GIs+1
          S = 'blastp -db Blast_DB_3_25  -query '+path+"/"+i+' -out '+ Source+'Blast_Output_New/'+Pattern_File+'_'+i[:-3]+'_blastP.txt  -evalue 10e-10 -outfmt "6 qseqid sseqid  pident evalue"'      #
          os.system(S)
      ###########################################################################################  
