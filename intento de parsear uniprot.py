# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 01:13:04 2020

@author: plasp
"""

from Bio import SeqIO
import urllib

handle = urllib.request.urlopen("http://www.uniprot.org/uniprot/P12345.xml")
record = SeqIO.read(handle, "uniprot-xml")
#print(record)
print(record.name)


#%%


from Bio import SeqIO
uniprot = SeqIO.index("uniprot_sprot.dat", "swiss")
for acc in ["P33487", "P19801", "P13689", "Q8JZQ5", "Q9TRC7"]:
    print(uniprot.get_raw(acc))
    
#%%

import pandas as pd

data = pd.read_csv("Pellegrini_sup2.csv")