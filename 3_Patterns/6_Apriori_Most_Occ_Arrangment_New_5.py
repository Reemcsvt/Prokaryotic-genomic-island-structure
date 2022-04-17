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

Indx=0
dataset=[]
Temp=[]
Occ=[]
Pattern=[]
Pattern_Occ=[]
count=1
Most_Occ_Arrangment=open('Apriori_Most_Occ_Arrangment_New_5.txt', 'w')

with open('Apriori_AllCombinations_Occ_New_5.txt', 'r') as f: #Apriori_AllCombinations_Occ_New_5.txt  Apri_AllComp_O.txt
  for line in f:
      print(count)
      count=count+1
      
      Pattern_Size=int(line.split(None, 1)[0])      
      for i in range (1,Pattern_Size+1):
            dataset.append(line.split(None, i+1)[i])   
      Occ=int(line.split(None, Pattern_Size+2)[Pattern_Size+1])
      #print(str(dataset))
      #print(str(Occ))
      if Indx==0:
         Pattern.append(dataset)
         Pattern_Occ.append(Occ)
         Temp=dataset
         dataset=[]
      else:
        if len(set(dataset).intersection(set(Temp)))==len(dataset):                  
           Pattern.append(dataset)
           Pattern_Occ.append(Occ)
           Temp=dataset
           dataset=[]
        else:
         #Find Max    
         #print(str(Pattern_Occ)) 
         Max=max(Pattern_Occ)  
         index=Pattern_Occ.index(Max)
         
         Most_Occ_Arrangment.write(str(len(Pattern[index])))
         for i in Pattern[index]:
            Most_Occ_Arrangment.write(' '+i)
         Most_Occ_Arrangment.write(' '+str(Pattern_Occ[index])+'\n')
         
         #clear parameters 
         Pattern_Occ=[]
         Pattern=[]          
         Pattern.append(dataset)
         Pattern_Occ.append(Occ)                  
         Temp=dataset
         dataset=[]             
      Indx=Indx+1   

Max=max(Pattern_Occ)  
index=Pattern_Occ.index(Max)
Most_Occ_Arrangment.write(str(len(Pattern[index])))
for i in Pattern[index]:
  Most_Occ_Arrangment.write(' '+i)
Most_Occ_Arrangment.write(' '+str(Pattern_Occ[index])+'\n')

      