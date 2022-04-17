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

################################################################################
#Read PFs from Get_PFs_New.txt and save them in two lists
PFs1=[]
with open('New_ProteinFamilies_5.txt', 'r') as f: # New_ProteinFamilies_5.txt
    for line in f:
      PFs1.append(line.split(None, 1)[0]) 

#PFs1=['XFP_N', 'DUF4277', 'NinG', 'DUF2528', 'Phage_CII', 'Phage_Nu1', 'Phage_antitermQ']
print('Done read protein families')
################################################################################
#Read Matrix
count=0
PFsInMatrix=open('PFs_InMatrix_G1.txt', 'w')
with open('New_InputMatrix_New_5.txt', 'r') as f: #New_InputMatrix_5.txt  a.txt
  for line in f:
     print('Count= '+str(count))
     if count>0:
       PF=0
       counter = line.count('1') 
       if counter>1:
        for i in line:
         if i=='1':
            PFsInMatrix.write(PFs1[PF]+' ')         
         if i=='1' or i=='0':   
            PF=PF+1   
        PFsInMatrix.write('\n')  
     count=count+1               
     
print('Done read proteins from matrix')





                   