# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 09:09:02 2021

@author: Jing Li, Small steps make changes. dnt_seq@163.com
"""

import numpy as np
import pandas as pd
import time
import random
from sklearn.utils import shuffle
from Bio.Seq import Seq
from Bio import SeqIO

#####################################################################
#####################################################################
orf_list = ['ORF1ab', 'S','ORF3a','E','M','ORF6','ORF7a','ORF7b','ORF8','N','ORF10'] 
nt_check_list = []
for item in orf_list:
    nt_check_list.append('nt_check_' + item)
print (nt_check_list)
nt_check_list2 = ['nt_check_ORF1ab', 'nt_check_S', 'nt_check_E', 'nt_check_M', 'nt_check_N']
 
template_file = 'SARS2_template.gb'
orf_loc_dict = {}
for record in SeqIO.parse(template_file, 'genbank'):
    seq_id = record.id
    for feature in record.features:
        feature_type = feature.type
        if feature_type == 'CDS':
            location_feature = feature.location
            gene_name0 = feature.qualifiers.get('gene')[0]
            # print (gene_name0)
            orf_loc_dict[gene_name0] = location_feature
# print (dictionary)
# for orf in orf_list:
    # print (orf_loc_dict[orf])
#####################################################################
#####################################################################


fasta_file = '../Data/2021-07-10_masked.fasta'
orfs_all_list = []
check_all_list = []

number_check2 = []
for record in SeqIO.parse(fasta_file, 'fasta'):
    number_check2.append((record.id))
    if len(number_check2)<=1000:
        orfs_record_list = []
        check_record_list = []
        for orf in orf_list:
            orf_loc =  orf_loc_dict[orf]
            target_orf = str(orf_loc.extract(record.seq).upper())
            seqlen = len(target_orf)
            nt_count = (target_orf.count('T') + target_orf.count('C') + target_orf.count('A') + target_orf.count('G'))
            nt_freq_check = nt_count / seqlen
            orfs_record_list.append(target_orf)
            check_record_list.append(nt_freq_check)
        orfs_all_list.append(orfs_record_list)
        check_all_list.append(check_record_list)

df_fasta = pd.DataFrame()
df_fasta[orf_list] = orfs_all_list
df_fasta[nt_check_list] = check_all_list
print (df_fasta.shape)
for check_i in nt_check_list2:
    df_fasta = df_fasta[df_fasta[check_i]>=0.998]
    print (df_fasta.shape)
df_fasta.to_csv ('parsed_orfs__2021-07-10_masked_ntchecked_test.csv')
