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
Path='/ExpressoDB/'
HGT_Path=Path+'Patterns_HGT_Species_Information/'
HGT_Table=Path+'Patterns_HGT_Species_Table_Information/'

HGT_Table=open(HGT_Table+'Patterns_HGT_Virus_Table_New_5.txt', 'w'); index=0;
with open(HGT_Path+'Patterns_HGT_Virus_PFs_New_5.txt', 'r') as f:  #HGT_Table_Test.txt Patterns_HGT_Virus_PFs_New_5.txt
  for line in f:
#1 Phage_capsid&Phage_portal&Terminase_1& Erysipelothrix_phage_phi1605 MF172979.1 Phage_portal ASD51098.1 WP_041721144.1 98.082 0.0 208226 Alkaliphilus_metalliredigens NC_009633.1_B_19, 0.14285714285714285
      index+=1; print(index)
      a=line.split(None, 1)[0]; PF=line.split(None, 5)[4];
      if index==1:
         prev_index=a; HGTTable_Genomes=[]; HGTTable_Genomes_GIs=[]; Genomes=[]; Genomes_GIs=[]; seen=[]; HGTTable_PF=line.split(None, 5)[4];
         PFs=[]; PFs.append(PF); Prev_Pattern_len=len((line.split(None, 2)[1]).split('&'))-1;
         GIs=line.split(None, 12)[11];  HGTTable_GIs=GIs.split(',');
         for k in HGTTable_GIs:
           if k!='':   
              G=k.split('_B_')[0]; 
              if G+'$' not in seen: 
                 Genomes.append(G); seen.append(G+'$'); Genomes_GIs.append(k)                
         Prev_FirstPart=prev_index+' '+line.split(None, 2)[1]+' '+line.split(None, 3)[2]+' '+line.split(None, 4)[3]
      ########################
      else:
         if a == prev_index:
             PFs.append(PF)
             if PF==HGTTable_PF:
                GIs=line.split(None, 12)[11];  HGTTable_GIs=GIs.split(',');
                for k in HGTTable_GIs:
                   if k!='':   
                      G=k.split('_B_')[0]; 
                      if G+'$' not in seen: 
                         Genomes.append(G); seen.append(G+'$');  Genomes_GIs.append(k)       
             else:
                HGTTable_Genomes.append(Genomes); HGTTable_Genomes_GIs.append(Genomes_GIs)
                HGTTable_PF=line.split(None, 5)[4]; Genomes=[]; Genomes_GIs=[]; seen=[];
                GIs=line.split(None, 12)[11];  HGTTable_GIs=GIs.split(',');
                for k in HGTTable_GIs:
                   if k!='':   
                      G=k.split('_B_')[0]; 
                      if G+'$' not in seen: 
                         Genomes.append(G); seen.append(G+'$'); Genomes_GIs.append(k)                
         ##########
         else:             
             C=len(Counter(PFs));
             HGTTable_Genomes.append(Genomes); HGTTable_Genomes_GIs.append(Genomes_GIs)
             HGTTable_Genomes_Wanted=set(HGTTable_Genomes[0]).intersection(*HGTTable_Genomes)
             if len( HGTTable_Genomes_Wanted)>0 and C==Prev_Pattern_len:
                HGT_Table.write(Prev_FirstPart+' ')   
                w_index=0; GIs_Genomes=''
                for w in HGTTable_Genomes_Wanted: 
                  w_index+=1
                  if w_index<len(HGTTable_Genomes_Wanted):
                        HGT_Table.write(w+',')
                  else: HGT_Table.write(w)
                  ###
                  for g in HGTTable_Genomes_GIs:
                    for j in g:
                       GenomeGI=j.split('_B_')[0];
                       if GenomeGI==w:
                          GIs_Genomes=GIs_Genomes+','+j
                          break
                  GIs_Genomes=GIs_Genomes+'#'        
                HGT_Table.write(' '+GIs_Genomes+'\n')          
             #########
             prev_index=a; HGTTable_Genomes=[]; HGTTable_Genomes_GIs=[]; Genomes=[]; Genomes_GIs=[]; seen=[]; HGTTable_PF=line.split(None, 5)[4];
             PFs=[]; PFs.append(PF); Prev_Pattern_len=len((line.split(None, 2)[1]).split('&'))-1;
             GIs=line.split(None, 12)[11];  HGTTable_GIs=GIs.split(',');
             for k in HGTTable_GIs:
               if k!='':   
                  G=k.split('_B_')[0]; 
                  if G+'$' not in seen: 
                     Genomes.append(G); seen.append(G+'$'); Genomes_GIs.append(k)                
             Prev_FirstPart=prev_index+' '+line.split(None, 2)[1]+' '+line.split(None, 3)[2]+' '+line.split(None, 4)[3]
                 
###########################################################

C=len(Counter(PFs));
HGTTable_Genomes.append(Genomes); HGTTable_Genomes_GIs.append(Genomes_GIs)
HGTTable_Genomes_Wanted=set(HGTTable_Genomes[0]).intersection(*HGTTable_Genomes)
if len( HGTTable_Genomes_Wanted)>0 and C==Prev_Pattern_len:
  HGT_Table.write(Prev_FirstPart+' ')   
  w_index=0; GIs_Genomes=''
  for w in HGTTable_Genomes_Wanted: 
    w_index+=1
    if w_index<len(HGTTable_Genomes_Wanted):
          HGT_Table.write(w+',')
    else: HGT_Table.write(w)
    ###
    for g in HGTTable_Genomes_GIs:
      for j in g:
         GenomeGI=j.split('_B_')[0];
         if GenomeGI==w:
            GIs_Genomes=GIs_Genomes+','+j
            break
    GIs_Genomes=GIs_Genomes+'#'        
  HGT_Table.write(' '+GIs_Genomes+'\n')          
