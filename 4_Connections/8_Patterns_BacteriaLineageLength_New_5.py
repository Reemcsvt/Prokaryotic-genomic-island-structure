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

####################################################################################
Source='/Patterns_Lineage_Information/'                
InterstingA=[]
Intersting_ViralGenomes=[] 
Source1='/Patterns_Intersting_Information/' 
Patterns_InterstingA_Information=open(Source1+"Patterns_Intersting_Information.txt", "w") 
Patterns_Intersting_ViralGenomes_Information=open(Source1+"Patterns_Intersting_ViralGenomes_Information.txt", "w")       
count=0;
with open('Patterns_Intersting_Information_New_5.txt', 'r') as f: #
  for line_Pattern in f:
      count=count+1;  print('Pattern'+' '+str(count)); dataset=[];  dataset.append([]);  
      PF=int(line_Pattern.split(None, 1)[0])
      for i in range (1,PF+1):
            dataset[0].append(line_Pattern.split(None, i+1)[i])                      
      Pattern_Print=''; Pattern=''; 
      for p in dataset[0]: Pattern_Print=Pattern_Print+p+'_';   Pattern=Pattern+p+'&';   
      Pattern_File=str(PF)+"_"+line_Pattern.split(None, PF+2)[PF+1]+"_"+Pattern_Print+'Local'
      Patterns_Lineage_Information=Source+'Patterns_Lineage_Information_'+Pattern_File+'.txt'
      #######################################################################################
      Patterns_Intersting_Information=open(Source1+"Patterns_Intersting_Information_"+Pattern_File+".txt", "w")
      Intersting=[]
      with open(Patterns_Lineage_Information, 'r') as f: 
          for line in f:
#a&b& Acinetobacter_phage_fEg-Aba01 MT344103.1 2 NZ_CP020586.1,NC_000907.1 NZ_CP020586.1_22,NC_000907.1_12 2.99e-68,2.99e-68 3 Gammaproteobacteria,Bacteria,Proteobacteria,
            A=int(line.split(None, 8)[7]); 
            Viral=line.split(None, 2)[1] ; VG=line.split(None, 3)[2]           
            if A<5:
               Intersting.append((line,A));  InterstingA.append((line,A)); 
               Intersting_ViralGenomes.append((Viral, VG,A))
      Intersting.sort(key = lambda x: x[1])
      for i in Intersting:
         Patterns_Intersting_Information.write(i[0])
      
InterstingA.sort(key = lambda x: x[1])
for i in InterstingA:
   Patterns_InterstingA_Information.write(i[0])

Intersting_ViralGenomes.sort(key = lambda x: x[2])
for i in Intersting_ViralGenomes:
   Patterns_Intersting_ViralGenomes_Information.write(i[0]+' '+i[1]+' '+str(i[2])+'\n')
             