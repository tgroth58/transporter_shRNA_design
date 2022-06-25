import os
import pandas as pd
from functools import reduce

libDir='shRNA_targets'

#Read dataframes of shRNA designs for each gene, return as list
shRNA_dfs=[]
for f in [x for x in os.listdir(libDir) if 'shRNA_designs.txt' in x]:
    #Get file prefix:
    pref=f.split('_')[0]
    #Read dataframe:
    df=pd.read_csv('/'.join([libDir,f]),sep='\t',index_col=0)
    #df=pd.read_csv('/'.join([libDir,f]),sep='\t',index_col=0).dropna(subset=['efficiency'])
    df['geneName']=pref
    shRNA_dfs.append(df)

#Concatenate all DFs vertically
total_lib=pd.concat(shRNA_dfs,axis=0)

#Custom filtering:
dgMean=-(33+28)/2
new_df_list=[]
# If there are less than 25 "efficient" shRNAs found by the design tool:
#   Sort by:
#     1. diffdG: the difference between an shRNA's dG and the "mean" dG, defined by the variable "dgMean"
#     2. ddG: the delta delta G value.  The threshold >2, but we sorted in descending order to pick at least 25 shRNAs per gene.

for _,df in total_lib.groupby(['geneName']):
    num_efficient=df.dropna().shape[0]
    if num_efficient>=25:
        df=df.dropna()
    #Sort by diffdG and ddG
    df['diffdG']=[abs(x-dgMean) for x in df['dG']]
    df=df.sort_values(['diffdG','ddG'],ascending=[True,False])
    new_df_list.append(df.head(25))

total_lib_top=pd.concat(new_df_list,axis=0)
#END OF FILTERING

#Read in shRNAs targeting random DNA sequences:
random_shRNAs=pd.read_csv('total_random_seqs_designs.txt',sep='\t',index_col=0).dropna(subset=['efficiency'])
random_shRNAs['geneName']=['_'.join(['random',str(x)]) for x in range(random_shRNAs.shape[0])]
total_lib_top=pd.concat([total_lib_top,random_shRNAs])

#Process sequences into constructs:
hairpin='GTTAATATTCATAGC'
trdict={'A':'T','C':'G','G':'C','T':'A'}
u_trdict={'U':'T'}
revcomp_hairpin=hairpin.translate(trdict)[::-1]
fwd_pref='CTAGC';fwd_suff='TTTTTG'
rev_pref='AATTCAAAAA';

total_lib_top['Forward Primer']=total_lib_top.apply(lambda x:fwd_pref+x['sense'].upper()+hairpin+x['antisense']+fwd_suff.translate(u_trdict),axis=1)
total_lib_top['Reverse Primer']=total_lib_top.apply(lambda x:rev_pref+x['sense'].upper()+revcomp_hairpin+x['antisense'].translate(u_trdict),axis=1)
#Write out CSV file:
total_lib_top.to_csv('total_shRNA_constructs.csv',sep='\t')
