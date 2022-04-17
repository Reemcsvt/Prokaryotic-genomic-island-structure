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

#Input: Files in folder Patterns_Information_New
#Output: Patterns_Information_New_5
#Function: For each virus and its genome compute occurrence i.e. how many GIs have proteins belong to this genome viral
########################################################################
from collections import Counter
Source='/Patterns_ViralGenomes_Information/'
with open('Patterns_Occurance_Information_New_5.txt', 'r') as f: #
  count=0
  for line_Pattern in f:
      count=count+1;  print('Pattern'+' '+str(count)); dataset=[];  dataset.append([]);  
      PF=int(line_Pattern.split(None, 1)[0])
      for i in range (1,PF+1):
            dataset[0].append(line_Pattern.split(None, i+1)[i])                       
      Pattern_Print=''; Pattern=''; 
      for p in dataset[0]: Pattern_Print=Pattern_Print+p+'_';   Pattern=Pattern+p+'&';   
      Pattern_File=str(PF)+"_"+line_Pattern.split(None, PF+2)[PF+1]+"_"+Pattern_Print+'Local'
      Patterns_ViralGenomes_Information=Source+'Patterns_ViralGenomes_Information_'+Pattern_File+'.txt'
      #######################################################################################
      Source1='/Patterns_Occurance_Information/'
      Patterns_Occurance_Information=open(Source1+"Patterns_Occurance_Information_GI_"+Pattern_File+".txt", "w")              
      V_G=[]; Bacteria=[]; Evalue=[]; Bacteria_GI=[]
      with open(Patterns_ViralGenomes_Information, 'r') as f: 
          for line in f:
            #a&b& NZ_CP021326.1 NZ_CP021326.1_B_17 Acinetobacter_baumannii Moraxella_phage_Mcat7 KR093631.1 3.96e-99 
            Bacteria.append(line.split(None, 2)[1]);  Bacteria_GI.append(line.split(None, 3)[2]);  Evalue.append(line.split(None, 7)[6])          
            V_G.append(line.split(None, 5)[4]+'&'+line.split(None, 6)[5])  
          Virus_Genome=Counter(V_G); Viruses_Genomes=[]; Viruses_Genomes_Occurance=[]
          for key,value in sorted(Virus_Genome.items()):
              Viruses_Genomes.append(key) 
              Viruses_Genomes_Occurance.append(value) 
          
          #print(Viruses_Genomes)
          #print(Viruses_Genomes_Occurance)
          
          Bacteria_Print=[]; Bacteria_GI_Print=[]; Evalue_Print=[];
          for i in Viruses_Genomes:
             Bacteria_Print.append(''); Bacteria_GI_Print.append('');  Evalue_Print.append('')   
          
          IDX=-1
          for i in V_G:
             IDX=IDX+1
             index=Viruses_Genomes.index(i)
             if Bacteria_Print[index]=='':
                 Bacteria_Print[index]=Bacteria[IDX]
                 Bacteria_GI_Print[index]=Bacteria_GI[IDX]
                 Evalue_Print[index]= Evalue[IDX]             
             else:
                 Bacteria_Print[index]=Bacteria_Print[index]+','+Bacteria[IDX]
                 Bacteria_GI_Print[index]=Bacteria_GI_Print[index]+','+Bacteria_GI[IDX]
                 Evalue_Print[index]= Evalue_Print[index]+','+Evalue[IDX]
          
          IDX=-1
          for i in Viruses_Genomes:
               IDX=IDX+1;        Percent=(Viruses_Genomes_Occurance[IDX]/len(V_G))*100           
               Patterns_Occurance_Information.write(Pattern+' '+i.split('&')[0]+' '+i.split('&')[1]+' ')
               Patterns_Occurance_Information.write(str(Viruses_Genomes_Occurance[IDX])+' '+str(Percent)+' ') 
               Patterns_Occurance_Information.write(Bacteria_Print[IDX]+' '+Bacteria_GI_Print[IDX]+' '+Evalue_Print[IDX]+'\n')
              
  





