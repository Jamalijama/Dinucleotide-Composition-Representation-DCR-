# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 11:06:36 2021

@author: Jing Li, Small steps make changes. dnt_seq@163.com
"""
import pandas as pd
import numpy as np
from Counting_full_nucleotides import count_DCR
import os

#################### 

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

full_cols_list0 = dnts_cols_list + dntpair_cols_list + codon_cols_list + codonpair_cols_list + amino_cols_list
full_cols_list1 = [col[5:] for col in full_cols_list0]

#################### 

path = r"../Data/"
file_list = os.listdir(path)

target_cols = ['orf_GP', 'orf_RdRp']
for file_obj in file_list:

    if file_obj.startswith ('df0_CoVs_SARS2_excluded_ORFs_final_cleaned') & file_obj.endswith('.csv'):
        print (file_obj)
        csv_data = pd.read_csv(file_obj)
        print (csv_data.shape)
        seq_num = csv_data.shape[0]
        id_list = csv_data.loc[:, 'accession'].tolist()

        for target_col in target_cols:
            print (target_col)
            dnt_count_list = []
            df_NCR = pd.DataFrame ()
            array_count = np.zeros(shape = (seq_num,1))
            orf_seq_list = csv_data.loc[:, target_col].tolist()
            print (len(orf_seq_list))
            for seq_i in range(seq_num):
                my_seq = orf_seq_list[seq_i].lower()        #################### code should be editted if sliding window needed
                seq_len = (len(my_seq))
                freq_dnt = count_DCR(my_seq)
                dnt_count_list.append(freq_dnt)
            # print (len(dnt_count_list))
    
            # print(len(id_list))
            array_freq_dnt = np.array(dnt_count_list)
            print (array_freq_dnt.shape)
            
            df_NCR['accession'] = id_list
            df_NCR[full_cols_list0] = array_freq_dnt
            df_NCR = df_NCR.set_index(['accession'])

            df_NCR.to_csv ('df_full_DCR_counting' + file_obj[3:-4] + '_' + target_col + '.csv')

###########################################################################
###########################################################################
###########################################################################
