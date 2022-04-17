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

#############################################################################################
# Integrate the hmmscan results in one file after removing the start and the end of each file #
#############################################################################################

Result = open("HmmscanResults_New_5.txt", "w")
File = ['PFsHmmscan1', 'PFsHmmscan2', 'PFsHmmscan3', 'PFsHmmscan4','PFsHmmscan5','PFsHmmscan6','PFsHmmscan7','PFsHmmscan8','PFsHmmscan9','PFsHmmscan10']
for x in File:
    count=0
    with open('/OutputFiles/'+x+'.faa') as a:
        for lineee in a:
            count += 1      
    count=count-3
    with open('/OutputFiles/'+x+'.faa') as f:
        for _ in range(3):
            next(f)
        for line in f:
            if count>10:
                   Result.write(line)
            count=count-1  
Result.close()
 
      
