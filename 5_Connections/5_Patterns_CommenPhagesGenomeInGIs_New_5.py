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

#Input: Files in folder Update
#Output: Viral genomes
#Function: Find common genome for each virus in each GI
########################################################################
Source='/Patterns_Update_Information/'
with open('Patterns_ViralGenomes_Information_New_5.txt', 'r') as f: #
  count=0
  for line_Pattern in f:
      count=count+1;  print('Pattern'+' '+str(count)); dataset=[];  dataset.append([]);  
      PF=int(line_Pattern.split(None, 1)[0])
      for i in range (1,PF+1):
            dataset[0].append(line_Pattern.split(None, i+1)[i])                       

      Pattern_Print=''; Pattern=''; 
      for p in dataset[0]: Pattern_Print=Pattern_Print+p+'_';   Pattern=Pattern+p+'&';   
      Pattern_File=str(PF)+"_"+line_Pattern.split(None, PF+2)[PF+1]+"_"+Pattern_Print+'Local'
      Patterns_Update_Information=Source+'Patterns_Update_Information_'+Pattern_File+'.txt'
      #######################################################################################
      Source1='/Patterns_ViralGenomes_Information/'
      Patterns_ViralGenomes_Information=open(Source1+"Patterns_ViralGenomes_Information_"+Pattern_File+".txt", "w")              
      Viral_Name=[]; Protein_No=[];  Protein_Geno=[]; Protein_Evalue=[];  index=-1 
      with open(Patterns_Update_Information, 'r') as f: 
          for line in f:
            #Terminase_4&Phage_capsid&HNH& NZ_CP031126.1 NZ_CP031126.1_B_24 Bacillus_licheniformis Bacillus_phage_Carmel_SA 1 ARW58537.1 KY963371.1 9.91e-12
            index=index+1; Flage=0
            G=line.split(None, 2)[1]; GI=line.split(None, 3)[2]; B=line.split(None, 4)[3];  V=line.split(None, 5)[4] 
            PNo=int(line.split(None, 6)[5]); PGeno=line.split(None, 8)[7]; PEvalue=float(line.split(None, 9)[8]);             
            
            if index==0:
              Pattern_G=G; Pattern_GI=GI; Bacteria=B; ProNum=PNo; 
              Viral_Name.append(V); Protein_No.append(PNo);  Protein_Geno.append(PGeno); Protein_Evalue.append(PEvalue); 
            else:
              if PNo<ProNum: 
                 Print_Bacteria_Line= Pattern+' '+Pattern_G+' '+Pattern_GI+' '+Bacteria+' '                 
                 #WRITE######################################################################                 
                 output = []                 
                 for x in Protein_No: 
                    if x not in output:  output.append(x)  
                 if len(output)==PF:
                    Genomes_Intersection=[]; RANGE=len(output)
                    for i in range (0,RANGE):
                       Genomes_Intersection.append([])                                         
                    index_genomes=-1; 
                    for i in Protein_No:                     
                       index_genomes=index_genomes+1;   Genomes_Intersection[i-1].append(Protein_Geno[index_genomes])  
                    I=set.intersection(*[set(x) for x in Genomes_Intersection]);                           
                    if len(I)>0:
                       E_G=[]
                       for i in I: 
                          index_Evalue=-1; Evalues_Genomes=[]; 
                          for j in Protein_Geno:
                             index_Evalue=index_Evalue+1
                             if i == j:  Evalues_Genomes.append(Protein_Evalue[index_Evalue])
                          E_G.append(max(Evalues_Genomes)) 
                       index_write=-1
                       for i in I: 
                          Virus=Viral_Name[Protein_Geno.index(i)];  index_write=index_write+1  
                          Patterns_ViralGenomes_Information.write(Print_Bacteria_Line+Virus+' '+str(i)+' '+str(E_G[index_write])+'\n')  
                 #Initialization
                 Viral_Name=[]; Protein_No=[]; Protein_Geno=[]; Protein_Evalue=[]; Pattern_G=G; Pattern_GI=GI; Bacteria=B; ProNum=PNo; 
                 Viral_Name.append(V); Protein_No.append(PNo);  Protein_Geno.append(PGeno); Protein_Evalue.append(PEvalue);                
                 ##############################################################################                          
              else: 
                    ProNum=PNo; Viral_Name.append(V); Protein_No.append(PNo);  Protein_Geno.append(PGeno); Protein_Evalue.append(PEvalue);         
##################################################################################################################
      Print_Bacteria_Line= Pattern+' '+Pattern_G+' '+Pattern_GI+' '+Bacteria+' '
      output = []
      for x in Protein_No: 
         if x not in output:  output.append(x)  
      if len(output)==PF:
         Genomes_Intersection=[]; RANGE=len(output)
         for i in range (0,RANGE):
            Genomes_Intersection.append([])                                        
         index_genomes=-1; 
         for i in Protein_No:                     
            index_genomes=index_genomes+1;   Genomes_Intersection[i-1].append(Protein_Geno[index_genomes])  
         I=set.intersection(*[set(x) for x in Genomes_Intersection]);         
         if len(I)>0:
            E_G=[]
            for i in I: 
               index_Evalue=-1; Evalues_Genomes=[]; 
               for j in Protein_Geno:
                  index_Evalue=index_Evalue+1
                  if i == j:  Evalues_Genomes.append(Protein_Evalue[index_Evalue])
               E_G.append(max(Evalues_Genomes)) 
            index_write=-1
            for i in I: 
               Virus=Viral_Name[Protein_Geno.index(i)];  index_write=index_write+1  
               Patterns_ViralGenomes_Information.write(Print_Bacteria_Line+Virus+' '+str(i)+' '+str(E_G[index_write])+'\n')       

