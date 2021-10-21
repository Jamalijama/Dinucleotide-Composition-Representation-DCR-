# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 19:14:42 2021

@author: Jing Li, Small steps make changes. dnt_seq@163.com
"""

import os
import pandas as pd
from Bio import SeqIO

# orf_list = ['ORF1ab', 'S','E','M','N'] 
orf_list = ['ORF1ab', 'S','ORF3a','E','M','ORF6','ORF7a','ORF7b','ORF8','N','ORF10'] 
col_list0 = ['accession', 'organism','mol_type', 'isolate', 'isolation_source',\
             'host', 'country', 'collection_date', 'taxonomy0', 'taxonomy1',\
             'taxonomy2', 'taxonomy3', 'taxonomy4','taxonomy5', 'sequence','seqlen']
ORF1ab_list = ['seq_ORF1ab','seqlen_ORF1ab','product_ORF1ab',' protein_id_ORF1ab',' db_xref_ORF1ab',' translation_ORF1ab']
S_list = ['seq_S','seqlen_S','product_S',' protein_id_S',' db_xref_S',' translation_S']
ORF3a_list = ['seq_ORF3a','seqlen_ORF3a','product_ORF3a',' protein_id_ORF3a',' db_xref_ORF3a',' translation_ORF3a']
E_list = ['seq_E','seqlen_E','product_E',' protein_id_E',' db_xref_E',' translation_E']
M_list = ['seq_M','seqlen_M','product_M',' protein_id_M',' db_xref_M',' translation_M']
ORF6_list = ['seq_ORF6','seqlen_ORF6','product_ORF6',' protein_id_ORF6',' db_xref_ORF6',' translation_ORF6']
ORF7a_list = ['seq_ORF7a','seqlen_ORF7a','product_ORF7a',' protein_id_ORF7a',' db_xref_ORF7a',' translation_ORF7a']
ORF7b_list = ['seq_ORF7b','seqlen_ORF7b','product_ORF7b',' protein_id_ORF7b',' db_xref_ORF7b',' translation_ORF7b']
ORF8_list = ['seq_ORF8','seqlen_ORF8','product_ORF8',' protein_id_ORF8',' db_xref_ORF8',' translation_ORF8']
N_list = ['seq_N','seqlen_N','product_N',' protein_id_N',' db_xref_N',' translation_N']
ORF10_list = ['seq_ORF10','seqlen_ORF10','product_ORF10',' protein_id_ORF10',' db_xref_ORF10',' translation_ORF10']
col_list = col_list0 + ORF1ab_list + S_list + ORF3a_list + E_list + M_list\
            + ORF6_list + ORF7a_list + ORF7b_list + ORF8_list + N_list + ORF10_list
# col_list = col_list0 + ORF1ab_list + S_list + E_list + M_list + N_list
print (len(col_list))


for file in os.listdir('./'):

    if file.startswith('data1_SARS2_')&file.endswith('.gb'):
        print (file)

        target_id_list = []

        if file.startswith('data1_SARS2_')&file.endswith('.gb'):
            for record in SeqIO.parse(file, 'genbank'):
                full_orfs_list = []
                seq_id = record.id
                location_list = []
                for feature in record.features:
                    feature_type = feature.type
                    if feature_type == 'CDS':
                        location_feature = feature.location
                        gene_name0 = feature.qualifiers.get('gene')
                        if gene_name0:
                            gene_name = gene_name0[0]
                        else:
                            gene_name = gene_name0
                        full_orfs_list.append(gene_name)
                        if set(orf_list) == set(full_orfs_list):
                            target_id_list.append(seq_id)
        print (len(target_id_list))

        virus_item_list = []
        index_list = []
        for record in SeqIO.parse(file, 'genbank'):
            item_list = []
            seq_id = record.id

            if seq_id in target_id_list:
                index_list.append(seq_id)
                accession = record.annotations['accessions'][0] #record.annotations['accessions'] is a list 
                taxonomy_list = record.annotations['taxonomy'][-6:] #annotations['taxonomy'] is a list of 6 eles, for some, less than 6 eles
                n = 0
                repeat_control_list = []                

                for feature in record.features:
                    feature_type = feature.type
    
                    if feature_type == 'source':
                        strain_organism = feature.qualifiers.get('organism')

                        if strain_organism:
                            strain_organism = strain_organism[0]
                        else:
                            strain_organism = strain_organism
        
                        strain_mol_type = feature.qualifiers.get('mol_type')

                        if strain_mol_type:
                            strain_mol_type = strain_mol_type[0]
                        else:
                            strain_mol_type = strain_mol_type
        
                        strain_isolate = feature.qualifiers.get('isolate')

                        if strain_isolate:
                            strain_isolate = strain_isolate [0]
                        else:
                            strain_isolate = strain_isolate
        
                        strain_isolation_source0 = feature.qualifiers.get('isolation_source')

                        if strain_isolation_source0:
                            strain_isolation_source = strain_isolation_source0[0]
                        else:
                            strain_isolation_source = strain_isolation_source0
    
                        strain_host0 = feature.qualifiers.get('host')

                        if strain_host0:
                            strain_host = strain_host0[0]
                        else:
                            strain_host = strain_host0
    
                        strain_country0 = feature.qualifiers.get('country')

                        if strain_country0:
                            strain_country = strain_country0[0]
                        else:
                            strain_country = strain_country0
    
                        strain_collection_date = feature.qualifiers.get('collection_date')

                        if strain_collection_date:
                            strain_collection_date = strain_collection_date[0]
                        else:
                            strain_collection_date = strain_collection_date
    
                        location_feature = feature.location
                        sequence = str(location_feature.extract(record.seq).upper())
                        target_seqlen = len(sequence)

                        item_list.append(accession)
                        item_list.append(strain_organism)
                        item_list.append(strain_mol_type)
                        item_list.append(strain_isolate)
                        item_list.append(strain_isolation_source)
                        item_list.append(strain_host)
                        item_list.append(strain_country)
                        item_list.append(strain_collection_date)
                        item_list = item_list + taxonomy_list
                        item_list.append(sequence)
                        item_list.append(target_seqlen)

                    if feature_type == 'CDS':
                        location_feature = feature.location
                        sequence0 = str(location_feature.extract(record.seq).upper())
                        sequencelen = len(sequence0)
                        gene_name0 = feature.qualifiers.get('gene')
                        if gene_name0:
                            gene_name = gene_name0[0]
                        else:
                            gene_name = gene_name0
                        if (gene_name not in repeat_control_list)&(gene_name in orf_list):
                            repeat_control_list.append(gene_name)
                            my_seq = sequence0
                            my_seqlen = sequencelen
                            product0 = feature.qualifiers.get('product')
                            if product0:
                                my_product = product0[0]
                            else:
                                my_product = product0
                            protein_id0 = feature.qualifiers.get('protein_id')
                            if protein_id0:
                                my_protein_id = protein_id0[0]
                            else:
                                my_protein_id = protein_id0
                            db_xref0 = feature.qualifiers.get('db_xref')
                            if db_xref0:
                                my_db_xref = db_xref0[0]
                            else:
                                my_db_xref = db_xref0
                            translation0 = feature.qualifiers.get('translation')
                            if translation0:
                                my_translation = translation0[0]
                            else:
                                my_translation = translation0

                            item_list.append(my_seq)
                            item_list.append(my_seqlen)
                            item_list.append(my_product)
                            item_list.append(my_protein_id)
                            item_list.append(my_db_xref)
                            item_list.append(my_translation)

                virus_item_list.append(item_list)
                n += 1

        print(len(virus_item_list))

        data=pd.DataFrame(virus_item_list,columns=col_list, index=index_list)
        print (data.shape)
        data_check = data

        fw1=open(file[:-3] + '_parsed_ORFs.csv', "w",newline='')
        data_check.to_csv(fw1,sep=',', header=True, index=True)

        fw2=open(file[:-3] + '_parsed_details.txt', "w")
        print ('col_list:','\n', col_list,'\n','data shape original: ','\n', data.shape, '\n', 'data post nt check: ','\n', data_check.shape, file = fw2)

        fw1.close()
        fw2.close()