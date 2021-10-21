import logomaker
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


'''
def __init__(self,
             df,
             color_scheme=None,
             font_name='sans',
             stack_order='big_on_top',
             center_values=False,
             baseline_width=0.5,
             flip_below=True,
             shade_below=0.0,
             fade_below=0.0,
             fade_probabilities=False,
             vpad=0.0,
             vsep=0.0,
             alpha=1.0,
             show_spines=None,
             ax=None,
             zorder=0,
             figsize=(10, 2.5),  ## usually x / y == 4
             **kwargs):
'''

# df_dominant_sites = pd.read_csv ('dominant_differential_dnts_for_adapt1_all_adapt2_noVOC.csv')
dominant_sites = (pd.read_csv ('dominant_differential_dnts_for_adapt1_all_adapt2_noVOC.csv')).loc[:,'index'].tolist()

all_list = list(range(0,3821,1))
sorted_list = dominant_sites + list(set(all_list).difference(set(dominant_sites)))
# print (sorted_list[:5])
# df_adapt = pd.read_csv('bayes_by_month/bayes_adapt_%s_except_lambda.csv' % m, index_col=0)
# df_inadapt = pd.read_csv('bayes_by_month/bayes_nonadapt_%s_except_lambda.csv' % m, index_col=0)
df_adapt = pd.read_csv('dnt_count_matrix_sars2_S_gisaid_2021-07-10_unmasked_ntchecked_1pred_all_T.csv')
df_inadapt = pd.read_csv('dnt_count_matrix_sars2_S_gisaid_2021-07-10_unmasked_ntchecked_2pred_noVOC_T.csv')

df_adapt['index'] = df_adapt['index'].astype('category')
df_adapt['index'].cat.reorder_categories(sorted_list, inplace=True)
df_adapt.sort_values('index', inplace=True)

df_inadapt['index'] = df_inadapt['index'].astype('category')
df_inadapt['index'].cat.reorder_categories(sorted_list, inplace=True)
df_inadapt.sort_values('index', inplace=True)


print (df_adapt.shape, df_inadapt.shape)

df_adapt = df_adapt.iloc[:,1:].replace(-1, 0)
df_inadapt = df_inadapt.iloc[:,1:].replace(-1, 0)
df_sort = df_adapt - df_inadapt
print (df_sort.shape)
sort = []
for i in range(25):
    line = df_sort.iloc[i]
    sort.append(line.max())
# df_sort_new = pd.DataFrame({'sort': sort}, index=df_adapt.index)
df_sort_new = pd.DataFrame({'sort': sort}, index=sorted_list[:25])

df_sort_new.sort_values(by='sort', ascending=False).to_csv('sort_diff_dot_product_1pred_all_2pred_noVOC.csv')
df_sort_new = df_sort_new.sort_values(by='sort', ascending=False)
print(df_sort_new)
# df_adapt = pd.concat([df_adapt, df_sort_new], axis=1)
# df_adapt = df_adapt.sort_values(by='sort', ascending=False)
# df_adapt = df_adapt.drop('sort', 1)
# df_adapt.to_csv('sort_bayes_adapt_%s_lambda_rm_x_0813.csv' % m)

# df_inadapt = pd.concat([df_inadapt, df_sort_new], axis=1)
# df_inadapt = df_inadapt.sort_values(by='sort', ascending=False)
# df_inadapt = df_inadapt.drop('sort', 1)
# df_inadapt.to_csv('sort_bayes_nonadapt_%s_lambda_rm_x_0813.csv' % m)
# df_inadapt = -df_inadapt

df_adapt = df_adapt.iloc[:25]
df_adapt = df_adapt.replace(-1, 0)
print(df_adapt.index)
df_adapt.index = range(25)
df_inadapt = df_inadapt.iloc[:25]
df_inadapt = df_inadapt.replace(1, 0)
df_inadapt.index = range(25)
print(len(df_adapt.columns))


logo = logomaker.Logo(df_adapt, color_scheme='NajafabadiEtAl2017', font_name='Stencil Std', vpad=.1, width=.8) # ,figsize=(30, 2.5)
logo.style_spines(visible=False)
logo.style_spines(spines=['left', 'bottom'], visible=True)
logo.style_xticks(rotation=90)
logo.ax.set_ylabel('Counting_number')
logo.ax.set_xticks(df_adapt.index)
# plt.show()
plt.savefig('logo_plot_SARS2_1pred_all.png', format='png', dpi=600)

logo_inadapt = logomaker.Logo(df_inadapt, color_scheme='NajafabadiEtAl2017', font_name='Stencil Std', vpad=.1,width=.8)
logo_inadapt.style_spines(visible=False)
logo_inadapt.style_spines(spines=['left', 'top', 'bottom'], visible=True)
logo_inadapt.style_xticks(rotation=90)
logo_inadapt.ax.set_xticks(df_inadapt.index)
# plt.show()
plt.savefig('logo_plot_SARS2_2pred_noVOC.png', format='png', dpi=600)
