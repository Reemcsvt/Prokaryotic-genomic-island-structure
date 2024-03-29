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

import os
import Bio
from Bio import SeqIO
from collections import Counter

PATH='/Genomes_Proteins_Files_New_5/'
Output='/Output/'
BlastOutput='/BlastOutput/'
Input=PATH

Bacteria_List=[]; index=0
with open(PATH+'Relations_Genomes.txt', 'r') as f: #Relations_Genomes.txt
  for line in f:  
      index+=1 
      #if index <=939:  #500
      Bacteria_List.append(line.split(None, 1)[0])

for i in Bacteria_List:
    S=('hmmsearch  --cut_nc  --tblout   '+BlastOutput+i+'_BlastOutput.fasta'+'     Pfam-A.hmm    '+Input+i+'.txt')
    os.system(S)    ??? K     ??? ??� ??? ???     ???    ?? ??? ??? 0   ???                 ??  ??? ??�     ??  ??� ??� ??? ??? - ?   ??? ???      ?? ??? ??? ??? ???     ??? ??? ???     ??? ??? ??? ???    ??� ??? ??? ??? ??? ??? ???   ???C  ??? ??? ??? ??? ??? ?     ??? ???    ???  ?? ??? ??? ??? ???     ?? ??? ??? @   ???  ??     ??? ??? ???     ??� ??� ??� ?           ??� ??? ???         ??? ??? ??? ???         ??? ??� ?           ??�     ??� ??? ?   ??? ??� ?   ???   ? ? �  ?? ??? (   0   ???   ????      0   (   ??? F     ??? ???  ?? ???     ??? ???  ?? ??� ??� ???        ??? ?   ??  ??� ??� ???   ? ?'  ##########################################################################################
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

import os
import Bio
from Bio import SeqIO
from collections import Counter

PATH='/Genomes_Proteins_Files_New_5/'
Output='/Output/'
BlastOutput='/BlastOutput/'
Input=PATH

Bacteria_List=[]; index=0
with open(PATH+'Relations_Genomes.txt', 'r') as f: #Relations_Genomes.txt
  for line in f:  
      index+=1 
      #if index <=939:  #500
      Bacteria_List.append(line.split(None, 1)[0])

for i in Bacteria_List:
    S=('hmmsearch  --cut_nc  --tblout   '+BlastOutput+i+'_BlastOutput.fasta'+'     Pfam-A.hmm    '+Input+i+'.txt')
    os.system(S)    ??? K     ??? ??� ??? ???     ???    ?? ??? ??? 0   ???                 ??  ?