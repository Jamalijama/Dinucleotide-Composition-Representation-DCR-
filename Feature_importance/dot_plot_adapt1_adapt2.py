# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 20:47:08 2021

@author: Jing Li, Small steps make changes. dnt_seq@163.com
"""

import numpy as np
import pandas as pd

pd.set_option('display.max_columns',1000)
pd.set_option('display.max_colwidth',1000)
pd.set_option('display.max_rows',50000)
pd.set_option('display.width', None)  # 设置字符显示宽度
np.set_printoptions(threshold=np.inf) #################一定加上这一条，否则打印或打印输出不全

adaptive1 = 'dnt_count_matrix_sars2_S_gisaid_2021-07-10_unmasked_ntchecked_1pred_all.csv'
adaptive2 = 'dnt_count_matrix_sars2_S_gisaid_2021-07-10_unmasked_ntchecked_2pred_noVOC.csv'

csv_adapt1 = pd.read_csv (adaptive1, index_col = 'index')
csv_adapt1 = (csv_adapt1 - csv_adapt1.min()) / (csv_adapt1.max() - csv_adapt1.min())
csv_adapt2 = pd.read_csv (adaptive2, index_col = 'index')
csv_adapt2 = (csv_adapt2 - csv_adapt2.min()) / (csv_adapt2.max() - csv_adapt2.min())
print (csv_adapt1.shape)
dot_list = []

# dot = np.dot(np.array(csv_adapt1), np.array(csv_adapt1))
# print (dot)
for i in range(csv_adapt1.shape[1]):
    vector_adapt1 = np.array(csv_adapt1.iloc[:,i])
    vector_adapt2 = np.array(csv_adapt2.iloc[:,i])
    dot_product = vector_adapt1.dot(vector_adapt2)
    dot_list.append(dot_product)

print (len(dot_list))

df_features = pd.DataFrame (dot_list)
df_features ['Sites'] = list(range(1,3822,1))
print (df_features.shape)
df_features.to_csv ('dnt_count_matrix_sars2_S_gisaid_2021-07-10_adapt1_adapt2_no_voc.csv')