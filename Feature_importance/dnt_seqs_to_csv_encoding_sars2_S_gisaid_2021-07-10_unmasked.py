# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 19:05:10 2021

@author: Jing Li, Small steps make changes. dnt_seq@163.com
"""
import numpy as np
import pandas as pd
import time
import random
from sklearn.utils import shuffle
from Bio.Seq import Seq
from Bio import SeqIO

pd.set_option('display.max_columns',1000)
pd.set_option('display.max_colwidth',1000)
pd.set_option('display.max_rows',50000)
pd.set_option('display.width', None)  # 设置字符显示宽度
np.set_printoptions(threshold=np.inf) #################一定加上这一条，否则打印或打印输出不全
NT_list = ['T','C','A','G','-']
DNT_list = []
for NT1 in NT_list:
    for NT2 in NT_list:
        DNT = NT1+NT2
        DNT_list.append(DNT)
# print (DNT_list)
dnt_list = list('defh-iklm-npqr-svwy------')
# print (dnt_list)
dict_dnt = dict(zip(DNT_list, dnt_list))
print (dict_dnt)

# fw1 = open ('dnt_seqs_sars2_S_gisaid_2021-07-10_unmasked_ntchecked_1pred.fasta','w')
# fw2 = open ('dnt_seqs_sars2_S_gisaid_2021-07-10_unmasked_ntchecked_2pred.fasta','w')
fw1 = open ('dnt_seqs_sars2_S_gisaid_2021-07-10_unmasked_ntchecked_1pred.csv','w')
fw2 = open ('dnt_seqs_sars2_S_gisaid_2021-07-10_unmasked_ntchecked_2pred.csv','w')

pred1_file = 'df_Pred_concat_info_3dCNN_full_DCR_sars2_Gp_gisaid_2021-07-10_unmasked_ntchecked_v10.4_2_1predict.csv'
pred2_file = 'df_Pred_concat_info_3dCNN_full_DCR_sars2_Gp_gisaid_2021-07-10_unmasked_ntchecked_v10.4_2_2predict.csv'
pred1_csv = pd.read_csv (pred1_file)
pred2_csv = pd.read_csv (pred1_file)
pred1_id_list = pred1_csv['Accession'].tolist()
pred2_id_list = pred2_csv['Accession'].tolist()

gisaid_unmask = 'parsed_sars2_orfs_gisaid_2021-07-10_unmasked_ntchecked.csv'
df_gisaid_unmask = pd.read_csv (gisaid_unmask)
seq_num = df_gisaid_unmask.shape[0]

dnt_S_list_all1 = []
dnt_S_list_all2 = []

for seq_i in range(seq_num):
    seq_id = df_gisaid_unmask.loc[seq_i,'seq_id']
    seq_S = df_gisaid_unmask.loc[seq_i,'S'].upper()
    seqlen = len (seq_S)
    # print (seqlen)
    dnt_S = ''
    for i in range (0,seqlen-1,1):
        DNT = seq_S[i:i+2]
        if DNT in DNT_list:
            DNT = DNT
        else:
            DNT = '--'
        dnt_S = dnt_S+dict_dnt[DNT]
    dnt_S_list = list(dnt_S)

    if seq_id in pred1_id_list:
        dnt_S_list_all1.append(dnt_S_list)

    elif seq_id in pred2_id_list:
        dnt_S_list_all2.append(dnt_S_list)

    if seq_i%10000==0:
        print (seq_i)

df_dnt_S1 = pd.DataFrame(dnt_S_list_all1,columns = range(seqlen-1))
df_dnt_S2 = pd.DataFrame(dnt_S_list_all2,columns = range(seqlen-1))
df_dnt_S1.to_csv(fw1)
df_dnt_S2.to_csv(fw2)

fw1.close()
fw2.close()
