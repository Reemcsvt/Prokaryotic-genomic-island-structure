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
from matplotlib import pyplot as plt

Source='/Patterns_Intersting_Information/'                
Source1='/Patterns_Filtering_Information/' 

Evalues=[]; #Histogram_List=[]
with open(Source+'Patterns_Intersting_Information.txt', 'r') as f: #
  for line_Pattern in f:
      E=line_Pattern.split(None, 7)[6]; Pattern_Evalues=E.split(',')
      for i in Pattern_Evalues: 
         Evalues.append(float(i));   

#Interesting_Evalues=open('Interesting_Evalues.txt','w')
#Evalues_Counter=Counter(Evalues);  
#for key,value in sorted(Evalues_Counter.items()):
#    Interesting_Evalues.write(str(key)+' '+str(value)+'\n') 
    
###########################################################################
Evalues=[]; Phage=[]
Patterns_Filtering_Information=open(Source1+'Patterns_Filtering_Information.txt','w')
Patterns_Filtering_Information_New=open(Source1+'Patterns_Filtering_Information_New.txt','w')
Patterns_Filtering_Information_New_5=open(Source1+'Patterns_Filtering_Information_New_5.txt','w')
with open(Source+'Patterns_Intersting_Information.txt', 'r') as f: #        Patterns_Intersting_Information.txt    test.txt
  for line_Pattern in f:
# Terminase_4&Phage_capsid&HNH& Streptococcus_phage_phi29854 MT311967.1 8 NZ_CP030880.1,NZ_CP030880.1,NZ_CP013457.1,NC_019042.1,NC_017582.1,NZ_CP013046.1,NZ_LR134488.1,NZ_LR134488.1 NZ_CP030880.1_B_9,NZ_CP030880.1_B_10,NZ_CP013457.1_B_49,NC_019042.1_B_8,NC_017582.1_B_2,NZ_CP013046.1_B_25,NZ_LR134488.1_B_10,NZ_LR134488.1_B_11 4.3e-13,4.3e-13,2.18e-11,1.56e-81,1.55e-72,6.06e-11,6.82e-12,6.82e-12 1 Bacteria,
      Bacteria_Acc=line_Pattern.split(None, 5)[4];   Pattern_Bacteria=Bacteria_Acc.split(',');
      Bacteria_GIs=line_Pattern.split(None, 6)[5];   Pattern_GIs=Bacteria_GIs.split(',');
      E=line_Pattern.split(None, 7)[6]; Pattern_Evalues=E.split(','); Flag=1 
      
      ##Remove if At least one is not ok
      for i in Pattern_Evalues: 
         if float(i)>=10e-100: #I need only less than 10e-100            
            Flag=0        
      if Flag==1:
         Patterns_Filtering_Information.write(line_Pattern)            
         #Phage.append(line_Pattern.split(None, 3)[2])   
      
      Flag=1
      ##Other Cases
      for i in Pattern_Evalues: 
         if float(i)<10e-100:             
            Flag=0        
      if Flag==0:
         Phage.append(line_Pattern.split(None, 3)[2])
         
         #Take whole line
         Patterns_Filtering_Information_New.write(line_Pattern)            
           
         #Take only Bacteria with <10e-100
         Patterns_Filtering_Information_New_5.write(line_Pattern.split(None, 1)[0]+' '+line_Pattern.split(None, 2)[1]+' '+line_Pattern.split(None, 3)[2]+' ')      
         line_index=-1; B_Acc=''; GIs_No=''; Es=''; Bacteria_count=0
         for i in Pattern_Evalues: 
             line_index+=1
             if float(i)<10e-100: 
               Bacteria_count+=1
               if B_Acc =='':
                B_Acc=Pattern_Bacteria[line_index];    GIs_No=Pattern_GIs[line_index];   Es=i               
               else: 
                B_Acc=B_Acc+','+Pattern_Bacteria[line_index]
                GIs_No=GIs_No+','+Pattern_GIs[line_index]
                Es=Es+','+i
         Patterns_Filtering_Information_New_5.write(str(Bacteria_count)+' '+B_Acc+' '+GIs_No+' '+Es+' '+line_Pattern.split(None, 8)[7]+' '+line_Pattern.split(None, 9)[8]+'\n')      
          
             
############################################################################         
            
Bacteriophages=open(Source1+'Bacteriophages.txt','w')
Phage_Counter=Counter(Phage);  
for key,value in sorted(Phage_Counter.items()):
    Bacteriophages.write(str(key)+'\n') 
           



         