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
import time
###########################################################
Spe=[]; AccN=[]
with open('Species_AccN_New_3.txt', 'r') as f: #Species names and accession number
    for line in f:
       Spe.append(line.split(None, 1)[0]); AccN.append(line.split(None, 2)[1])
print('Done read species')
###########################################################
Output_Source='/Patterns_Table_Information/'
with open('Patterns_Table_Information_New_5.txt', 'r') as f: #
  count=0
  for line_Pattern in f:
      #if count % 100 == 0: time.sleep(5)
      count=count+1;  print('Pattern'+' '+str(count)); dataset=[];  dataset.append([]);  
      PF=int(line_Pattern.split(None, 1)[0])
      for i in range (1,PF+1):
            dataset[0].append(line_Pattern.split(None, i+1)[i])                      
      Pattern_Print=''; Pattern=''; 
      for p in dataset[0]: Pattern_Print=Pattern_Print+p+'_';   Pattern=Pattern+p+'&';   
      Pattern_File=str(PF)+"_"+line_Pattern.split(None, PF+2)[PF+1]+"_"+Pattern_Print+'Local'
      Pattern_Info=open(Output_Source+'PatternInfo_'+Pattern_File+'.txt','w')
      
      Source="/Blast_Viral_Table/"
      File_Table=Source+'Evalues_Accs/'+'Evalues_Accs_'+Pattern_File+'.txt'
      #################################
      Total=0
      with open(File_Table, 'r') as f: 
          for line in f:
            Total=Total+1
      f.close()
      #################################
      index=0; GI=[]; V=[]; PNo=[]; PAcc=[]; PEvalue=[]; G=[]; Species=[]; Genome=[];  index_Entrez=0
      with open(File_Table, 'r') as f: 
          for line in f:
           Protein_Accession=line.split(None, 5)[4]
           if  Protein_Accession!= 'WP_015980308.1':
            index=index+1; 
            if index==1:   ProteinAcc=Protein_Accession;
            else:          ProteinAcc=ProteinAcc+','+Protein_Accession;   
            GI.append(line.split(None, 2)[1]); GIsland=line.split(None, 2)[1]; G.append(GIsland.split('_B_')[0]); Geno=GIsland.split('_B_')[0]
            Species.append(Spe[AccN.index(Geno)]); V.append(line.split(None, 3)[2]);  PNo.append(line.split(None, 4)[3]); 
            PAcc.append(Protein_Accession); PEvalue.append(line.split(None, 6)[5]); 
            if index==10000:
              index_Entrez=index_Entrez+1
              if index_Entrez ==10: time.sleep(10)
              try:
                  handle = Entrez.efetch(db="protein",id=ProteinAcc , retmode="xml",retmax=1000000)#, rettype="gb", retmode="text") 
                  print("Done Entrez efetch ="+str(index_Entrez))
                  record = Entrez.read(handle);  recordi_index=-1
                  for recordi in record:
                      recordi_index=recordi_index+1
                      ORG=recordi['GBSeq_source-db']
                      if 'accession' in ORG:
                          GGG=ORG[ORG.rfind('accession')+len('accession')+1:len(ORG)]   
                          if len(GGG)>10:    
                             A=ORG.split('accession ')[1]
                             Genome.append(A.split(';')[0])  
                          else:   Genome.append(GGG) 
                      elif 'locus' in  ORG:
                          GGG=ORG[ORG.rfind('locus')+len('locus')+1:len(ORG)]    
                          if len(GGG)>10:    
                             A=ORG.split('locus ')[1]
                             Genome.append(A.split(';')[0]) 
                          else:   Genome.append(GGG)   
                      else:
                          Genome.append('Error')               
              except HTTPError as err: Genome.append('Error')
              except ValueError: Genome.append('Error')
              except KeyError as e: Genome.append('Error')
              index=0
      ################
      if index>0:
        try:
            handle = Entrez.efetch(db="protein",id=ProteinAcc , retmode="xml",retmax=1000000)#, rettype="gb", retmode="text") 
            record = Entrez.read(handle);  recordi_index=-1
            for recordi in record:
                recordi_index=recordi_index+1
                ORG=recordi['GBSeq_source-db']
                if 'accession' in ORG:
                    GGG=ORG[ORG.rfind('accession')+len('accession')+1:len(ORG)]   
                    if len(GGG)>10:    
                       A=ORG.split('accession ')[1]
                       Genome.append(A.split(';')[0])  
                    else:   Genome.append(GGG) 
                elif 'locus' in  ORG:
                    GGG=ORG[ORG.rfind('locus')+len('locus')+1:len(ORG)]    
                    if len(GGG)>10:    
                       A=ORG.split('locus ')[1]
                       Genome.append(A.split(';')[0]) 
                    else:   Genome.append(GGG)   
                else:
                    Genome.append('Error')               
        except HTTPError as err: Genome.append('Error')
        except ValueError: Genome.append('Error')
        except KeyError as e: Genome.append('Error')
      
                         
      #################################                               
      INDX=-1
      for i in G:  
        INDX=INDX+1; print(INDX); 
        print(Pattern+' '+G[INDX]+' '+GI[INDX]+' '+Species[INDX]+' '+V[INDX]+' '+PNo[INDX]+' '+PAcc[INDX]+' '+Genome[INDX]+' '+PEvalue[INDX])    
        Pattern_Info.write(Pattern+' '+G[INDX]+' '+GI[INDX]+' '+Species[INDX]+' '+V[INDX]+' '+PNo[INDX]+' '+PAcc[INDX]+' '+Genome[INDX]+' '+PEvalue[INDX]+'\n')
            
            
            
            
            
            
            
            
            
            
            
            
