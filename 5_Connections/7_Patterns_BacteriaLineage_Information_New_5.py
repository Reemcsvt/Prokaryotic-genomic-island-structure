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

Bacteria_Accession=[]; Bacteria_Lineage=[]; index=-1
with open('Bacterias_Lineages_New_5.txt', 'r') as f: 
    for line in f:
       index=index+1; Acc=line.split(None, 1)[0]
       Bacteria_Accession.append(Acc);  Bacteria_Lineage.append([])
       A=line.split(None, 3)[2]; B=A.split(',')
       for i in B:  Bacteria_Lineage[index].append(i)        
print('Done read lineages')
####################################################################################
#Use GIs' GI accession instead of Bacteria name
Source='/Patterns_Occurance_Information/'
with open('Patterns_Lineage_Information_New_5.txt', 'r') as f: #
  count=0
  for line_Pattern in f:
      count=count+1;  print('Pattern'+' '+str(count)); dataset=[];  dataset.append([]);  
      PF=int(line_Pattern.split(None, 1)[0])
      for i in range (1,PF+1):
            dataset[0].append(line_Pattern.split(None, i+1)[i])                      
      Pattern_Print=''; Pattern=''; 
      for p in dataset[0]: Pattern_Print=Pattern_Print+p+'_';   Pattern=Pattern+p+'&';   
      Pattern_File=str(PF)+"_"+line_Pattern.split(None, PF+2)[PF+1]+"_"+Pattern_Print+'Local'
      Patterns_Occurance_Information=Source+'Patterns_Occurance_Information_'+Pattern_File+'.txt'
      #######################################################################################
      Source1='/Patterns_Lineage_Information/'
      Patterns_Lineage_Information=open(Source1+"Patterns_Lineage_Information_"+Pattern_File+".txt", "w")             
      V_G=[]; Bacteria=[]; Evalue=[]
      with open(Patterns_Occurance_Information, 'r') as f: 
          for line in f:
      #a&b& Acinetobacter_phage_fEg-Aba01 MT344103.1 2 NZ_CP020586.1,NZ_CP021326.1 NZ_CP020586.1_22,NZ_CP021326.1_4 2.99e-68,2.99e-68
            A=line.split(None, 5)[4]; B=A.split(',')
            Lineages=[]; Lineages_index=-1
            for i in B:
               BAI=Bacteria_Accession.index(i); Lineages.append([]);  Lineages_index=Lineages_index+1
               for j in Bacteria_Lineage[BAI]:
                  Lineages[Lineages_index].append(j)
            I=set.intersection(*[set(x) for x in Lineages])            
            Level='';
            for i in I: Level=Level+i+',';
            #write
            Patterns_Lineage_Information.write(line.rstrip('\n')+' '+str(len(I))+' '+Level+'\n')
            #if len(I)<3:
            #   Patterns_Lineages_Information_New_5_Intersting.write(line.rstrip('\n')+' '+str(len(I))+' '+Level+'\n')
            