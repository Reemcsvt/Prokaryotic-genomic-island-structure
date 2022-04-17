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

Path='/ExpressoDB/'
HGT_Table=Path+'Patterns_HGT_Species_Table_Information/'
####################################################
Merge_Table=open(HGT_Table+'Patterns_Merge_Table_New_5.txt', 'w')
Merge_Table_Related=open(HGT_Table+'Patterns_Merge_Table_Related_New_5.txt', 'w')
Bacteria_Table_line=[]; Bacteria_Table_Index=[];  Bacteria_Table_Pattern=[];  Bacteria_Table_VirusT=[]; 
Bacteria_Table_VirusID=[]; Bacteria_Table_GIs=[]; GIs_seen=[]
with open(HGT_Table+'Patterns_HGT_Bacteria_Table_New_5.txt', 'r') as f: #Bacteria_Test.txt   Patterns_HGT_Bacteria_Table_New_5.txt
  for line in f:
#644 Terminase_1&Phage_portal&Phage_capsid& Klebsiella_phage_ST16-OXA48phi5.3 MK416014.1 NZ_CP014762.1$NZ_LS992183.1,%,NZ_CP014762.1_B_27$NZ_LS992183.1_B_28#,NZ_CP014762.1_B_27$NZ_LS992183.1_B_28#,NZ_CP014762.1_B_27$NZ_LS992183.1_B_28#@NZ_CP020820.1$NZ_CP014762.1,%,NZ_CP020820.1_B_8$NZ_CP014762.1_B_27#,NZ_CP020820.1_B_8$NZ_CP014762.1_B_27#,NZ_CP020820.1_B_8$NZ_CP014762.1_B_27#@NZ_LS992183.1$NZ_CP020820.1,%,NZ_LS992183.1_B_28$NZ_CP020820.1_B_8#,NZ_LS992183.1_B_28$NZ_CP020820.1_B_8#,NZ_LS992183.1_B_28$NZ_CP020820.1_B_8#@

      Bacteria_Table_line.append(line)
      Bacteria_Table_Index.append(line.split(None, 1)[0]);  Bacteria_Table_Pattern.append(line.split(None, 2)[1]);  
      Bacteria_Table_VirusT.append(line.split(None, 3)[2]); Bacteria_Table_VirusID.append(line.split(None, 4)[3]); 
      Bacteria_Table_GIs.append(line.split(None, 5)[4]); seen=[]  
      GIs1=(line.split(None, 5)[4]).split('@'); 
      for i in GIs1: # ['a$b,%,a_B_27$b_B_28#,a_B_27$b_B_28#,a_B_27$b_B_28#', 'c$a,%,c_B_8$a_B_27#,c_B_8$a_B_27#,c_B_8$a_B_27#', 'b$c,%,b_B_28$c_B_8#,b_B_28$c_B_8#,b_B_28$c_B_8#', '']
          GIs2=i.split(',%,'); #'a$b,%,a_B_27$b_B_28#,a_B_27$b_B_28#,a_B_27$b_B_28#'
          GIs3=GIs2[0].split('$') #a$b
          if len(GIs3) ==2:
             seen.append(GIs3[0]+'$'); seen.append(GIs3[1]+'$');
      GIs_seen.append(seen)     

index=0
with open(HGT_Table+'Patterns_HGT_Virus_Table_New_5.txt', 'r') as f: #Virus_Test.txt   Patterns_HGT_Virus_Table_New_5.txt
  for line in f:
#384 Phage_capsid&Phage_portal&Terminase_1& Vibrio_phage_Va_90-11-287_p41_Ba35 MK672800.1 NZ_CP021980.1,NZ_CP016095.1 
#  ,NZ_CP021980.1_B_19,NZ_CP021980.1_B_19,NZ_CP021980.1_B_19#,NZ_CP016095.1_B_15,NZ_CP016095.1_B_15,NZ_CP016095.1_B_15#
      index+=1; print(index)
      Virus_Table_line=line
      Virus_Table_index=line.split(None, 1)[0]; Virus_Table_Pattern=line.split(None, 2)[1]; Virus_Table_Virus_Text=line.split(None, 3)[2]; 
      Virus_Table_Virus_ID=line.split(None, 4)[3];  Virus_Table_Genomes=line.split(None, 5)[4]; Virus_Table_GIs=line.split(None, 6)[5];
      
      Bacteria_Index=-1
      for i in Bacteria_Table_Index:
         Bacteria_Index+=1
         if Virus_Table_index==i:
            Merge_Table.write(Virus_Table_line+Bacteria_Table_line[Bacteria_Index])
            #####
            GIs=Virus_Table_Genomes.split(',');
            for j in GIs:
              if j !='':
                if j+'$' in GIs_seen[Bacteria_Index]:
                   Merge_Table_Related.write(Virus_Table_line+Bacteria_Table_line[Bacteria_Index])
                   break
            break
            






