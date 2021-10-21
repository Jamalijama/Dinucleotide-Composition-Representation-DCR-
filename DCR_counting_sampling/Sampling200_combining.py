# -*- coding: utf-8 -*-
"""
Created on Sat Aug  7 18:23:18 2021

@author: Jing Li, Small steps make changes. dnt_seq@163.com
"""

import pandas as pd
import numpy as np
import os

pd.set_option('display.max_columns',1000)
pd.set_option('display.max_colwidth',1000)
pd.set_option('display.max_rows',50000)
pd.set_option('display.width', None)  # 设置字符显示宽度
np.set_printoptions(threshold=10000)


amino_table = ['I', 'D', 'M', 'H', 'E', 'W', 'R', 'L', 'Y', 'Q', 'G', 'A', 'S', 'P', 'C', 'T', 'V', 'F', 'N', 'K']
nt_table = ['t','c', 'a', 'g']
dnt_table = [nt1+nt2 for nt1 in nt_table for nt2 in nt_table]
dnts_table = [dnt1+dnt2 for dnt1 in dnt_table for dnt2 in dnt_table]
codon_table = [nt1+nt2+nt3 for nt1 in nt_table for nt2 in nt_table for nt3 in nt_table]
codon_table1 = codon_table.copy()
codon_table1.remove('taa')
codon_table1.remove('tag')
codon_table1.remove('tga')
codonpair_table = [codon0 + codon1 for codon0 in codon_table1 for codon1 in codon_table1]
dnt_category = ['n12', 'n23','n31']
dntpair_category = ['n12m12', 'n23m23','n31m31','n12n31','n23m12','n31m23']
dnts_cols_list = ['Freq_'+ dnt + '_' + dnt_cat for dnt_cat in dnt_category for dnt in dnt_table]
dntpair_cols_list = ['Freq_'+ dnts + '_' + dnts_cat for dnts_cat in dntpair_category for dnts in dnts_table]
codon_cols_list = ['Freq_'+ codon for codon in codon_table]
codonpair_cols_list = ['Freq_'+ codonpair for codonpair in codonpair_table]
amino_cols_list = ['Freq_'+ amino for amino in amino_table]
# full_cols_list0 = dnts_cols_list + dntpair_cols_list + codon_cols_list + codonpair_cols_list + amino_cols_list
dcr_set_list = dnts_cols_list + dntpair_cols_list + codon_cols_list + codonpair_cols_list + amino_cols_list
print (len(dcr_set_list))
# print ((dcr_set_list))
dcr_set_list2 = dcr_set_list + ['Gene']

for file in os.listdir('./'):
    if file.startswith('df_full_DCR')&file.endswith('RdRp.csv'):
        df_DCR = pd.read_csv (file,index_col = 'accession')
        df_index = pd.read_csv (file[:-8] + 'Glyp_sampling200.csv',index_col = 'accession')
        print (df_DCR.shape)
        print (df_index.shape)
        index_list = df_index.index.tolist()
        # print (index_list)
        df_DCR_sampling2 = df_DCR.loc[index_list,:]
        df_DCR_sampling2['Gene'] = 200*['RdRp']
        # df_DCR_sampling2.to_csv(file[:-4] + '_sampling200.csv')
        print (df_DCR.shape)

df_all = pd.DataFrame()
for file in os.listdir('./'):
    if file.startswith('df_full_DCR')&file.endswith('_sampling200.csv'):
        print (file)
        gene_name = file[-20:-15]
        print (gene_name)
        virus_name = file[21:28]
        print(virus_name)
        
        df_DCR = pd.read_csv (file,index_col = 'accession')
        print (df_DCR.shape)
        accession_list = df_DCR.index.tolist()
        label_list = [virus_name + gene_name + '_' + accession for accession in accession_list]
        # print (label_list)
        df_DCR1 = df_DCR[dcr_set_list2]
        df_DCR1['Label'] = label_list
        df_all = pd.concat((df_all,df_DCR1),axis = 0)
        print (df_all.shape)
df_all.to_csv('df_full_DCR_counting_sampling200_RdRp_Gp_all_restsssss.csv')
