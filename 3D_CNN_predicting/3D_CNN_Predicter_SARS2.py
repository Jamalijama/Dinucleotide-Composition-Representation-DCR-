# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 17:06:06 2021

@author: Jing Li, Small steps make changes. dnt_seq@163.com
"""
import pandas as pd
import numpy as np
import torch
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F
import os

nt_table = ['t','c', 'a', 'g']
dnt_table = [nt1+nt2 for nt1 in nt_table for nt2 in nt_table]
dnts_table = [dnt1+dnt2 for dnt1 in dnt_table for dnt2 in dnt_table]
dntpair_category = ['n12m12', 'n23m23','n31m31','n12n31','n23m12','n31m23']
dntpair_cols_list = ['Freq_'+ dnts + '_' + dnts_cat for dnts_cat in dntpair_category for dnts in dnts_table]

model_name = './model_name.txt'

class CNN (nn.Module):
    def __init__ (self):
        super (CNN, self).__init__()
        self.conv1 = nn.Sequential ( 
            nn.Conv3d ( in_channels = 1,
                        out_channels = 8,
                        kernel_size = (1, 3, 3), # kernel only for 2d data
                        stride =(1,1,1),
                        padding = (0,1,1),
                        bias = True
                        ),            # 
            nn.ReLU (),
            nn.AvgPool3d (kernel_size = (1,2,2)) 
        )
        self.conv2 = nn.Sequential ( 
            nn.Conv3d ( in_channels = 8,
                        out_channels = 16,
                        kernel_size = (1, 3, 3),# kernel only for 2d data
                        stride =(1,1,1),
                        padding = (0,1,1),
                        bias = True
                        ),            # 
            nn.ReLU (),
            nn.AvgPool3d (kernel_size = (1, 2, 2)) # Max or Avg
        )
        self.conv3 = nn.Sequential (
            nn.Conv3d ( in_channels = 16,
                        out_channels = 32,
                        kernel_size = (1, 3, 3),# kernel only for 2d data
                        stride =(1,1,1),
                        bias = True,
                        padding = (0,1,1)
                        ),
            nn.ReLU (),
            nn.AvgPool3d (kernel_size = (1, 2, 2)) #MaxPool3d
        )
        
        self.fc1 = nn.Linear (768, 192)  # adding this step, too slow
        # self.fc2 = nn.Linear (192, 48)  # adding this step, too slow
        self.fc3 = nn.Linear (192, 3)

    def forward (self, x):  # x is the a matrix, however, only lowcase is suggeted to use
        x = self.conv1(x)
        x = self.conv2(x)
        x = self.conv3(x)

        x = x.view (x.size(0), -1) # flat x, similar to reshape of numpy
        # return x.shape

        x = F.sigmoid (self.fc1(x))  # to activate x
        # x = F.sigmoid (self.fc2(x))  # to activate x
        output = self.fc3(x)
        return output

cnn = CNN()
cnn = torch.load(model_name)
if_use_gpu = 1
if if_use_gpu:
    cnn = cnn.cuda()

for file in os.listdir('./'):
    if file.startswith ('TEST_sars2_orfs_unmasked_ntchecked_')&file.endswith('split_S.csv'):
        print (file)
        df_DCR_counting = pd.read_csv (file)
        # print (df_DCR_counting.columns)
        print (df_DCR_counting.shape)
        accession_list = df_DCR_counting.loc[:,'accession'].tolist()
        df_DCR = df_DCR_counting.loc[:,dntpair_cols_list]
        SARS2_DCR_array = np.array(df_DCR)
        SARS2_num = SARS2_DCR_array.shape[0]
        SARS2_DCR_array2 = SARS2_DCR_array.reshape(SARS2_num,1,6,16,16)
        print (SARS2_DCR_array2.shape)
        SARS2_DCR_tensor = torch.tensor(SARS2_DCR_array2)
        SARS2_DCR_tensor = SARS2_DCR_tensor.to(torch.float32)

        SARS2_pred = cnn (Variable(SARS2_DCR_tensor))
        SARS2_pred_array = SARS2_pred.data.numpy()
        SARS2_labels = torch.max (SARS2_pred, 1)[1].data

        res = pd.DataFrame ({'Accession': accession_list,'pred_label': SARS2_labels})
        print (res.shape)
        res [['Score_0','Score_1','Score_2']] = SARS2_pred_array
        print (res.shape)
        file_name = 'df_Pred_3dCNN_' + file[3:-4]  + '.csv'
        res.to_csv (file_name, index = False)