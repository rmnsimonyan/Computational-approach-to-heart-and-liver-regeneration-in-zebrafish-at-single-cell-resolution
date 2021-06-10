# Transform the table from the previous sh script to normal expression table

import pandas as pd
import math

table = pd.read_csv("final_matrix.csv", header=None)

# the list of unique UMI's
cells = table[0].unique()

# table[1].unique() is the list of unique gene names, we create a table with those names as rownames/indecies and then fill the table
# with the data from each cell one by one in the for loop

fintable = pd.DataFrame(index=table[1].unique())

for i in range(0,len(cells)):
    print(i)
    # Initialize the snptable with the first sample's snp table
    newcell = table[table[0] == cells[i]]
    newcell.index = newcell[1]
    newcell = newcell.drop([0,1,2,3], axis=1)
    newcell.columns = [cells[i]]
    newcell.sort_index(inplace=True)
    exptable = pd.merge(pd.DataFrame(index=table[1].unique()), newcell, left_index=True, right_index=True, how="outer")
    
    fintable = pd.concat([fintable,exptable], axis=1)
    

fintable.to_csv("exptable.csv", index = True)
