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

import pandas as pd
from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import apriori

dataset=[]
count=0
with open('PFs_InMatrix_G1.txt', 'r') as f: #PFs_InMatrix_G1.txt  AAA.txt
  for line in f:      
      PFs=len(line.split())
      if PFs>1:
       dataset.append([])
       for i in range (0,PFs):
            dataset[count].append(line.split(None, i+1)[i])              
       count=count+1
print('Done read data set')      
#############################################################################
te = TransactionEncoder()
print(te)
print('Done TransactionEncoder') 

te_ary = te.fit(dataset).transform(dataset)
print(te_ary)
print('Done fit transform')

df = pd.DataFrame(te_ary, columns=te.columns_)
print('DataFrame')
#print(df)

#############################################################################

Apriori_New=open('Apriori_New_5.txt', 'w')  #Apriori_New_5.txt Example_Presentation.txt

frequent_itemsets = apriori(df, min_support=0.006, use_colnames=True)
frequent_itemsets['length'] = frequent_itemsets['itemsets'].apply(lambda x: len(x))
print(frequent_itemsets)

count=0
for i,j,k in zip(frequent_itemsets.support, list(frequent_itemsets.itemsets), frequent_itemsets.length):
  if k>1:
   print(count)
   count=count+1
   #print(str(i)+' '+str(set(j))+' '+str(k))
   #Apriori_New.write(str(i)+' '+str(k)+' '+str(set(j))+'\n')
   Apriori_New.write(str(i)+' '+str(k))
   for m in j:
      Apriori_New.write(' '+str(set(j)))
   Apriori_New.write('\n')    
   
   
   