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
import itertools

Path='/ExpressoDB/'
HGT_Path=Path+'Patterns_HGT_Species_Information/'
HGT_Table_Path=Path+'Patterns_HGT_Species_Table_Information/'

HGT_Table=open(HGT_Table_Path+'Patterns_HGT_Bacteria_Table_New_5.txt', 'w'); index=0;
with open(HGT_Path+'Patterns_HGT_Bacteria_PFs_New_5.txt', 'r') as f:  #HGT_Bacteria_Table_Test.txt  Patterns_HGT_Bacteria_PFs_New_5.txt
  for line in f:
#1 Phage_capsid&Phage_portal&Terminase_1& Erysipelothrix_phage_phi1605 MF172979.1 Phage_capsid WP_038604150.1 WP_013888252.1 96.305 0.0 191610 258224 Corynebacterium_atypicum Corynebacterium_resistens NZ_CP008944.1_B_2, NC_015673.1_B_7,
      index+=1; print(index)
      a=line.split(None, 1)[0]; PF=line.split(None, 5)[4];
      if index==1:
         prev_index=a; HGTTable_Genomes=[]; HGTTable_Genomes_GIs=[]; Genomes1=[]; Genomes2=[]; Genomes1GIs=[]; Genomes2GIs=[]; 
         seen1=[]; seen2=[]; HGTTable_PF=line.split(None, 5)[4];
         PFs=[]; PFs.append(PF); Prev_Pattern_len=len((line.split(None, 2)[1]).split('&'))-1;
         GIs1=line.split(None, 14)[13];  HGTTable_GIs1=GIs1.split(',');
         GIs2=line.split(None, 15)[14];  HGTTable_GIs2=GIs2.split(',');
         for k in HGTTable_GIs1:
           if k!='':   
              G=k.split('_B_')[0]; 
              if G+'$' not in seen1: 
                 Genomes1.append(G); seen1.append(G+'$'); Genomes1GIs.append(k)                
         for k in HGTTable_GIs2:
           if k!='':   
              G=k.split('_B_')[0]; 
              if G+'$' not in seen2: 
                 Genomes2.append(G); seen2.append(G+'$'); Genomes2GIs.append(k)                               
         #########
         Genomes_product=[];     Genomes_product.append(Genomes1);  Genomes_product.append(Genomes2); 
         Genomes_product_GIs=[]; Genomes_product_GIs.append(Genomes1GIs); Genomes_product_GIs.append(Genomes2GIs);
         result = list(itertools.product(*Genomes_product));       Pattern_GIs_Combination=[]
         resultGIs = list(itertools.product(*Genomes_product_GIs));    Pattern_GIs_Combination_GIs=[]
         IDX=-1
         for k in result:
           IDX+=1
           if k[0]!=k[1]:
              Pattern_GIs_Combination.append(k[0]+'$'+k[1]); Pattern_GIs_Combination.append(k[1]+'$'+k[0])                 
              Pattern_GIs_Combination_GIs.append(resultGIs[IDX][0]+'$'+resultGIs[IDX][1]) 
              Pattern_GIs_Combination_GIs.append(resultGIs[IDX][1]+'$'+resultGIs[IDX][0])
         Prev_FirstPart=prev_index+' '+line.split(None, 2)[1]+' '+line.split(None, 3)[2]+' '+line.split(None, 4)[3]
                          
      ########################
      else:
           if a == prev_index:
               PFs.append(PF)
               if PF==HGTTable_PF:
                  GIs1=line.split(None, 14)[13];  HGTTable_GIs1=GIs1.split(',');
                  GIs2=line.split(None, 15)[14];  HGTTable_GIs2=GIs2.split(',');
                  for k in HGTTable_GIs1:
                     if k!='':   
                        G=k.split('_B_')[0]; 
                        if G+'$' not in seen1: 
                           Genomes1.append(G); seen1.append(G+'$'); Genomes1GIs.append(k)                
                  for k in HGTTable_GIs2:
                     if k!='':   
                        G=k.split('_B_')[0]; 
                        if G+'$' not in seen2: 
                           Genomes2.append(G); seen2.append(G+'$'); Genomes2GIs.append(k)                               
                  #########
                  Genomes_product=[];     Genomes_product.append(Genomes1);  Genomes_product.append(Genomes2); 
                  Genomes_product_GIs=[]; Genomes_product_GIs.append(Genomes1GIs); Genomes_product_GIs.append(Genomes2GIs);
                  result = list(itertools.product(*Genomes_product));       Pattern_GIs_Combination=[]
                  resultGIs = list(itertools.product(*Genomes_product_GIs));    Pattern_GIs_Combination_GIs=[]
                  IDX=-1
                  for k in result:
                     IDX+=1
                     if k[0]!=k[1]:
                        Pattern_GIs_Combination.append(k[0]+'$'+k[1]); Pattern_GIs_Combination.append(k[1]+'$'+k[0])                 
                        Pattern_GIs_Combination_GIs.append(resultGIs[IDX][0]+'$'+resultGIs[IDX][1]) 
                        Pattern_GIs_Combination_GIs.append(resultGIs[IDX][1]+'$'+resultGIs[IDX][0])
               ##################################                              
               else:
                  HGTTable_Genomes.append(Pattern_GIs_Combination);               
                  HGTTable_Genomes_GIs.append(Pattern_GIs_Combination_GIs)
                  HGTTable_PF=line.split(None, 5)[4]; Genomes1=[]; Genomes2=[]; Genomes1GIs=[]; Genomes2GIs=[]; seen1=[]; seen2=[];
                  GIs1=line.split(None, 14)[13];  HGTTable_GIs1=GIs1.split(',');
                  GIs2=line.split(None, 15)[14];  HGTTable_GIs2=GIs2.split(',');

                  for k in HGTTable_GIs1:
                     if k!='':   
                        G=k.split('_B_')[0]; 
                        if G+'$' not in seen1: 
                           Genomes1.append(G); seen1.append(G+'$');  Genomes1GIs.append(k)                             
                  for k in HGTTable_GIs2:
                     if k!='':   
                        G=k.split('_B_')[0]; 
                        if G+'$' not in seen2: 
                           Genomes2.append(G); seen2.append(G+'$');  Genomes2GIs.append(k)                              
                   
                  #########
                  Genomes_product=[];     Genomes_product.append(Genomes1);  Genomes_product.append(Genomes2); 
                  Genomes_product_GIs=[]; Genomes_product_GIs.append(Genomes1GIs); Genomes_product_GIs.append(Genomes2GIs);
                  result = list(itertools.product(*Genomes_product));       Pattern_GIs_Combination=[]
                  resultGIs = list(itertools.product(*Genomes_product_GIs));    Pattern_GIs_Combination_GIs=[]
                  IDX=-1
                  for k in result:
                     IDX+=1
                     if k[0]!=k[1]:
                        Pattern_GIs_Combination.append(k[0]+'$'+k[1]); Pattern_GIs_Combination.append(k[1]+'$'+k[0])                 
                        Pattern_GIs_Combination_GIs.append(resultGIs[IDX][0]+'$'+resultGIs[IDX][1]) 
                        Pattern_GIs_Combination_GIs.append(resultGIs[IDX][1]+'$'+resultGIs[IDX][0])
           ##########
           else: 
               C=len(Counter(PFs));                                                                                        
               HGTTable_Genomes.append(Pattern_GIs_Combination); 
               HGTTable_Genomes_GIs.append(Pattern_GIs_Combination_GIs)                                         
               HGTTable_Genomes_Wanted=set(HGTTable_Genomes[0]).intersection(*HGTTable_Genomes)
                              
               if len(HGTTable_Genomes_Wanted)>0 and C==Prev_Pattern_len:
                 HGT_Table.write(Prev_FirstPart+' ')   
                 
                 w_index=0; print_seen=[]; Flage=0
                 for w in HGTTable_Genomes_Wanted: 
                   w_index+=1; GIs_Genomes=''
                   y=w.split('$'); y1=y[0]; y2=y[1];                      
                   if (y1+'$'+y2 not in print_seen) and (y2+'$'+y1 not in print_seen):
                      print_seen.append(y1+'$'+y2); print_seen.append(y2+'$'+y1) 
                      if w_index<len(HGTTable_Genomes_Wanted):
                             HGT_Table.write(w+','); Flage=1
                      else: HGT_Table.write(w);  
                        
                      for g in HGTTable_Genomes_GIs:
                         for j in g:
                            Genome_GI=j.split('$');
                            GenomeGI0=Genome_GI[0].split('_B_')[0]; GenomeGI1=Genome_GI[1].split('_B_')[0];
                            if GenomeGI0==y1 and GenomeGI1==y2:
                               GIs_Genomes=GIs_Genomes+','+Genome_GI[0]+'$'+Genome_GI[1]
                               break
                         GIs_Genomes=GIs_Genomes+'#'      
                      HGT_Table.write('%'+GIs_Genomes+'@')          
                  
                 if Flage==1:
                    HGT_Table.write('\n')
                                                                      
               #############################################
               prev_index=a; HGTTable_Genomes=[]; HGTTable_Genomes_GIs=[]; Genomes1=[]; Genomes2=[]; Genomes1GIs=[]; Genomes2GIs=[]; 
               seen1=[]; seen2=[]; HGTTable_PF=line.split(None, 5)[4];
               PFs=[]; PFs.append(PF); Prev_Pattern_len=len((line.split(None, 2)[1]).split('&'))-1;
               GIs1=line.split(None, 14)[13];  HGTTable_GIs1=GIs1.split(',');
               GIs2=line.split(None, 15)[14];  HGTTable_GIs2=GIs2.split(',');
               for k in HGTTable_GIs1:
                 if k!='':   
                    G=k.split('_B_')[0]; 
                    if G+'$' not in seen1: 
                       Genomes1.append(G); seen1.append(G+'$'); Genomes1GIs.append(k)                
               for k in HGTTable_GIs2:
                 if k!='':   
                    G=k.split('_B_')[0]; 
                    if G+'$' not in seen2: 
                       Genomes2.append(G); seen2.append(G+'$'); Genomes2GIs.append(k)                               
               #########
               Genomes_product=[];     Genomes_product.append(Genomes1);  Genomes_product.append(Genomes2); 
               Genomes_product_GIs=[]; Genomes_product_GIs.append(Genomes1GIs); Genomes_product_GIs.append(Genomes2GIs);
               result = list(itertools.product(*Genomes_product));       Pattern_GIs_Combination=[]
               resultGIs = list(itertools.product(*Genomes_product_GIs));    Pattern_GIs_Combination_GIs=[]
               IDX=-1
               for k in result:
                 IDX+=1
                 if k[0]!=k[1]:
                    Pattern_GIs_Combination.append(k[0]+'$'+k[1]); Pattern_GIs_Combination.append(k[1]+'$'+k[0])                 
                    Pattern_GIs_Combination_GIs.append(resultGIs[IDX][0]+'$'+resultGIs[IDX][1]) 
                    Pattern_GIs_Combination_GIs.append(resultGIs[IDX][1]+'$'+resultGIs[IDX][0])
               Prev_FirstPart=prev_index+' '+line.split(None, 2)[1]+' '+line.split(None, 3)[2]+' '+line.split(None, 4)[3]                 
###########################################################

C=len(Counter(PFs));                                                                                        
HGTTable_Genomes.append(Pattern_GIs_Combination); 
HGTTable_Genomes_GIs.append(Pattern_GIs_Combination_GIs)                                         
HGTTable_Genomes_Wanted=set(HGTTable_Genomes[0]).intersection(*HGTTable_Genomes)
              
if len(HGTTable_Genomes_Wanted)>0 and C==Prev_Pattern_len:
  HGT_Table.write(Prev_FirstPart+' ')   
  
  w_index=0; print_seen=[]; Flage=0
  for w in HGTTable_Genomes_Wanted: 
    w_index+=1; GIs_Genomes=''
    y=w.split('$'); y1=y[0]; y2=y[1];                      
    if (y1+'$'+y2 not in print_seen) and (y2+'$'+y1 not in print_seen):
       print_seen.append(y1+'$'+y2); print_seen.append(y2+'$'+y1) 
       if w_index<len(HGTTable_Genomes_Wanted):
              HGT_Table.write(w+','); Flage=1
       else: HGT_Table.write(w);  
        
       for g in HGTTable_Genomes_GIs:
          for j in g:
             Genome_GI=j.split('$');
             GenomeGI0=Genome_GI[0].split('_B_')[0]; GenomeGI1=Genome_GI[1].split('_B_')[0];
             if GenomeGI0==y1 and GenomeGI1==y2:
                GIs_Genomes=GIs_Genomes+','+Genome_GI[0]+'$'+Genome_GI[1]
                break
          GIs_Genomes=GIs_Genomes+'#'      
       HGT_Table.write('%'+GIs_Genomes+'@')          
  
  if Flage==1:
     HGT_Table.write('\n')