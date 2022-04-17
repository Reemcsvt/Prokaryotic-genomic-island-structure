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

import math
from collections import Counter
######################################################################################
Spe=[]; AccN=[]
with open('Species_AccN_New_3.txt', 'r') as f: #Species names and accession number
    for line in f:
       Spe.append(line.split(None, 1)[0]); AccN.append(line.split(None, 2)[1])
print('Done read species')
#######################################################################################
Pro_Acc=[]; Pro_OrgName=[]
with open('Accession_Organisms_New_5', 'r') as f: #Organisms names and accession number
    for line in f:
       Pro_Acc.append(line.split(None, 1)[0]); Pro_OrgName.append(line.split(None, 2)[1])
print('Done read organisms')
########################################################################################
Source="/Blast_Viral_Table/"
with open('1_Patterns_New_5.txt', 'r') as f: #
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
      FILE_PFs=Source+'Wanted_PFs/'+'Wanted_PFs_'+Pattern_File+'.txt'  
      Evalues_Accs=open(Source+'Evalues_Accs/'+'Evalues_Accs_'+Pattern_File+'.txt','w')
      
      Pattern_Print=''
      for p in dataset[0]: Pattern_Print=Pattern_Print+p+'&'     
      ##############################################################################
      #Viruses1=open(Source+'Viruses/'+'VIRUSES_'+Pattern_File+'.txt', 'w')  
      #Matches1=open(Source+'Viruses/Viruses1_'+Pattern_File+'.txt', 'w') 
      ###############################################################################   
      GIs_Blast=[]; Species=[];  Genomes_dataset=[]; Pattern_GIs=[]
      with open(FILE_PFs, 'r') as f: 
        for line in f:
          A=line.split(None, 1)[0]
          if '$' not in A:
              GIs_Blast.append(A+'.fa'); Pattern_GIs.append(A);  G=A.split('_B_')[0];   
              Genomes_dataset.append(G);   Species.append(Spe[AccN.index(G)])           
      ###################################################################################                  
      M=[]; IDX=-1; GIs_Blast_Counter=0
      for i in GIs_Blast:       
          GIs_Blast_Counter=GIs_Blast_Counter+1; print(str(GIs_Blast_Counter)+'/'+str(len(GIs_Blast))+': '+i)   
          # Write
          A=[]; Protein_Acc=[]; Protein_Evalue=[]; index=-1;  #Matches1.write(i+'\n')      
          with open(Source+'Blast_Output_New/'+Pattern_File+'_'+i[:-3]+'_blastP.txt', 'r') as f: #
             for line in f:                                
               Query=line.split(None, 1)[0];  Accessoin=line.split(None, 2)[1];   Evalue=line.split(None, 4)[3]
               Organism=Pro_OrgName[Pro_Acc.index(Accessoin)] 
               if index==-1: 
                   index=index+1;   A.append([]); Protein_Acc.append([]); Protein_Evalue.append([])                                              
               else:
                 if Query !=  Q:             
                    index=index+1;     A.append([]);   Protein_Acc.append([]); Protein_Evalue.append([])     
               A[index].append(Organism); Protein_Acc[index].append(Accessoin); Protein_Evalue[index].append(Evalue)     
               Q=Query     
          #################################  
          Flag=0 
          if len(A)==PF:
            for A1 in A:
              if len(A1) ==0:
                 Flag=1 
          else: Flag=1
            
          if Flag==1:   A=[]; Protein_Acc=[];  Protein_Evalue=[];                       
          #############################  Intersection between proteins within the GI ###############   
          I=[]; IDX=IDX+1;  M.append([])  #GI viruses
          if len(A)>0:                          
            I=set.intersection(*[set(x) for x in A])         
            if len(I) > 0:
              for e in I:
                 M[IDX].append(e)  #List of viruses in this GI    
              ############################# Save Evalues and Acc for proteins that belong to the intersected viruses ##################
            for k in I: #common viruses
             A_index=-1
             Protein_count=0
             for h in A:  #lists of viruses where each list represent viruses in a protein // h is a list of viruses in a protein
               A_index=A_index+1               
               Protein_count=Protein_count+1
               if k in h:
                 H_index=-1 
                 for j in h:
                   H_index=H_index+1
                   if j==k:
                    Evalues_Accs.write(Pattern_Print+' '+i[:-3]+' '+k+' '+str(Protein_count)+' '+Protein_Acc[A_index][H_index]+' '+Protein_Evalue[A_index][H_index]+'\n')
                             
          #Viruses1.write(str(I)+'\n')
      #Viruses1.close();  Matches1.close();  
      ###########################################################################################                             
      ###########################################################################################  
      ########################## Intersection between GIs #######################################             
      #print('GIs Viruses= '+str(len(M)))                     
      Dataset_Viruses=[]; Dataset_Species=[]; M_New=[]; indx=-1; Index=-1
      Pattern_GIs_Viruses=open(Source+'Evalues_Accs/'+'Pattern_GIs_Viruses_'+Pattern_File+'.txt', 'w')  
      for i in M: #List of viruses in GIs
         Index=Index+1
         if len(i)>0:
            indx=indx+1
            M_New.append([])
            Dataset_Species.append(Species[Index])
            for j in i:
               Dataset_Viruses.append(j)
               M_New[indx].append(j+'$')
               Pattern_GIs_Viruses.write(Pattern_Print+' '+str(Pattern_GIs[Index])+' '+Species[Index]+' '+j+'\n')  
      #print(len(M_New))
      Pattern_GIs_Viruses.close();
      #############################################################
      Pattern_Viruses=open(Source+'Viruses/'+'Pattern_Viruses_'+Pattern_File+'.txt', 'w')  
      S=Counter(Dataset_Viruses)      
      for key,value in sorted(S.items(), key=lambda item: item[1]):
           IDX=-1
           Pattern_Viruses.write(key+' '+str(len(GIs_Blast))+' '+str(len(M_New))+' '+str(value)+' '+str((value/len(M_New))*100)+' ')
           Seen=[]
           for i in M_New:
            IDX=IDX+1
            if key+'$' in i:
              if Dataset_Species[IDX]+'$' not in Seen:
                 Pattern_Viruses.write(Dataset_Species[IDX]+'&')
                 Seen.append(Dataset_Species[IDX]+'$')    
           Pattern_Viruses.write('\n')         
      #############################################################################################
