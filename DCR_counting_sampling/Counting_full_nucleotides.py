# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 13:29:27 2021

@author: Jing Li, Small steps make changes. dnt_seq@163.com
"""

import numpy as np

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

def cds_translator(seq):
    amino_list = ['F','F',\
                  'L','L','L','L','L','L',\
                  'S','S','S','S','S','S',\
                  'Y','Y',\
                  '*','*',\
                  'C','C',\
                  '*',\
                  'W',\
                  'P','P','P','P',\
                  'H','H',\
                  'Q','Q',\
                  'R','R','R','R','R','R',\
                  'I','I','I',\
                  'M',\
                  'T','T','T','T',\
                  'N','N',\
                  'K','K',\
                  'V','V','V','V',\
                  'A','A','A','A',\
                  'D','D',\
                  'E','E',\
                  'G','G','G','G']
    codon_table1 = ['ttt', 'ttc',\
                    'tta', 'ttg','ctt', 'ctc', 'cta', 'ctg',\
                    'tct', 'tcc', 'tca', 'tcg', 'agt', 'agc',\
                    'tat', 'tac',\
                    'taa', 'tag',\
                    'tgt', 'tgc',\
                    'tga',\
                    'tgg',\
                    'cct', 'ccc', 'cca', 'ccg',\
                    'cat', 'cac',\
                    'caa', 'cag',\
                    'cgt', 'cgc', 'cga', 'cgg','aga', 'agg',\
                    'att', 'atc', 'ata',\
                    'atg',\
                    'act', 'acc', 'aca', 'acg',\
                    'aat', 'aac',\
                    'aaa', 'aag',\
                    'gtt', 'gtc', 'gta', 'gtg',\
                    'gct', 'gcc', 'gca', 'gcg',\
                    'gat', 'gac',\
                    'gaa', 'gag',\
                    'ggt', 'ggc', 'gga', 'ggg']
    
    # print (len(amino_list))
    # print (len(codon_table1))
    codon_dict = dict(zip(codon_table1, amino_list))
    # print (codon_dict)
    seqlen = len(seq)
    # if seqlen % 3 != 0:
    #     return 'It is not a cds sequence!'
    # else:
        # seq = seq.lower()
    Protein = ''
    for codon_i in range(0,seqlen,3):
        codon = seq[codon_i:codon_i+3]
        if codon in codon_table1:
            Protein = Protein+codon_dict[codon]
        # else:
        #     Protein = Protein+'?'
    return Protein[:-1]
def count_DCR (cds, dnt_table = dnt_table, dnts_table = dnts_table,\
               codon_table = codon_table, codonpair_table = codonpair_table,\
               amino_table = amino_table, Num0 = 256, Num1 = 64,Num2 = 16, \
               Num3 = 3721, Num4 = 20, freq = True):
    seq_len = len(cds)
    cds = cds.lower()

#################################################   DNT counting
    count_N12 = np.zeros(Num2)
    for i in range(0, seq_len-3, 3):
        cut = cds[i:i+2]
        if cut in dnt_table:
            dnt = dnt_table.index(cut)
            count_N12[dnt] += 1

    count_N23 = np.zeros(Num2)
    for i in range(0, seq_len-3, 3):
        cut = cds[i+1:i+3]
        if cut in dnt_table:
            dnt = dnt_table.index(cut)
            count_N23[dnt] += 1

    count_N31 = np.zeros(Num2)
    for i in range(0, seq_len-3, 3):
        cut = cds[i+2:i+4]
        if cut in dnt_table:
            dnt = dnt_table.index(cut)
            count_N31[dnt] += 1

#################################################   DNTpair counting
    count_N12M12 = np.zeros(Num0)
    for i in range(0, seq_len-6, 3):
        cut1 = cds[i:i+2]
        cut2 = cds[i+3:i+5]
        cut = cut1 + cut2
        if cut in dnts_table:
            dnts = dnts_table.index(cut)
            count_N12M12[dnts] += 1

    count_N23M23 = np.zeros(Num0)
    for i in range(0, seq_len-6, 3):
        cut1 = cds[i+1:i+3]
        cut2 = cds[i+4:i+6]
        cut = cut1 + cut2
        if cut in dnts_table:
            dnts = dnts_table.index(cut)
            count_N23M23[dnts] += 1

    count_N31M31 = np.zeros(Num0)
    for i in range(0, seq_len-6, 3):
        cut1 = cds[i+2:i+4]
        cut2 = cds[i+5:i+7]
        cut = cut1 + cut2
        if cut in dnts_table:
            dnts = dnts_table.index(cut)
            count_N31M31[dnts] += 1

    count_N12N31 = np.zeros(Num0)
    for i in range(0, seq_len-6, 3):
        # cut1 = cds[i:i + 2]
        # cut2 = cds[i+2:i + 4]
        cut = cds[i:i+4]
        if cut in dnts_table:
            dnts = dnts_table.index(cut)
            count_N12N31[dnts] += 1

    count_N23M12 = np.zeros(Num0)
    for i in range(0, seq_len-6, 3):
        cut = cds[i+1:i+5]
        if cut in dnts_table:
            dnts = dnts_table.index(cut)
            count_N23M12[dnts] += 1

    count_N31M23 = np.zeros(Num0)
    for i in range(0, seq_len-6, 3):
        cut = cds[i+2:i+6]
        if cut in dnts_table:
            dnts = dnts_table.index(cut)
            count_N31M23[dnts] += 1

#################################################   codon counting
    count_codon = np.zeros(Num1)
    for i in range(0, seq_len-3, 3):
        cut = cds[i:i+3]
        if cut in codon_table:
            codons = codon_table.index(cut)
            count_codon[codons] += 1

#################################################   codonpair counting
    count_codonpair = np.zeros(Num3)
    for i in range(0, seq_len-9, 3):
        cut = cds[i:i+6]
        if cut in codonpair_table:
            codonpair = codonpair_table.index(cut)
            count_codonpair[codonpair] += 1

#################################################   amino acid counting
    count_amino = np.zeros(Num4)
    protein = cds_translator(cds.lower())
    seqlen = len(protein)
    for i in range(0, seqlen, 1):
        cut = protein[i:i+1]
        if cut in amino_table:
            amino = amino_table.index(cut)
            count_amino[amino] += 1

    count_dnt = np.hstack((count_N12, count_N23, count_N31))
    count_dnts = np.hstack((count_N12M12, count_N23M23, count_N31M31, count_N12N31, count_N23M12, count_N31M23))

    '''
    ## counting twice for each nt, and twice for each dnt.
    ## Therefore, the freq for a dnt = count_dnts*1536 / (seq_len-1)/2
    '''

    if freq:
        return np.hstack ((count_dnt*16*3/seq_len, (count_dnts*256*3 / (seq_len-1)), \
                           count_codon*64*3/seq_len, count_codonpair*3721*3/(seq_len-3), \
                           count_amino*20*3/seq_len))
    else:
        return np.hstack ((count_dnt, count_dnts, count_codon, count_codonpair, count_amino))

################################################## test 
################################################## test 

# Seq_S = 'ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACTCAGGACTTGTTCTTACCTTTCTTTTCCAATGTTACTTGGTTCCATGCTATACATGTCTCTGGGACCAATGGTACTAAGAGGTTTGATAACCCTGTCCTACCATTTAATGATGGTGTTTATTTTGCTTCCACTGAGAAGTCTAACATAATAAGAGGCTGGATTTTTGGTACTACTTTAGATTCGAAGACCCAGTCCCTACTTATTGTTAATAACGCTACTAATGTTGTTATTAAAGTCTGTGAATTTCAATTTTGTAATGATCCATTTTTGGGTGTTTATTACCACAAAAACAACAAAAGTTGGATGGAAAGTGAGTTCAGAGTTTATTCTAGTGCGAATAATTGCACTTTTGAATATGTCTCTCAGCCTTTTCTTATGGACCTTGAAGGAAAACAGGGTAATTTCAAAAATCTTAGGGAATTTGTGTTTAAGAATATTGATGGTTATTTTAAAATATATTCTAAGCACACGCCTATTAATTTAGTGCGTGATCTCCCTCAGGGTTTTTCGGCTTTAGAACCATTGGTAGATTTGCCAATAGGTATTAACATCACTAGGTTTCAAACTTTACTTGCTTTACATAGAAGTTATTTGACTCCTGGTGATTCTTCTTCAGGTTGGACAGCTGGTGCTGCAGCTTATTATGTGGGTTATCTTCAACCTAGGACTTTTCTATTAAAATATAATGAAAATGGAACCATTACAGATGCTGTAGACTGTGCACTTGACCCTCTCTCAGAAACAAAGTGTACGTTGAAATCCTTCACTGTAGAAAAAGGAATCTATCAAACTTCTAACTTTAGAGTCCAACCAACAGAATCTATTGTTAGATTTCCTAATATTACAAACTTGTGCCCTTTTGGTGAAGTTTTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGCTGATTATTCTGTCCTATATAATTCCGCATCATTTTCCACTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAGATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAATTCTAACAATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAAAAAGTCTACTAATTTGGTTAAAAACAAATGTGTCAATTTCAACTTCAATGGTTTAACAGGCACAGGTGTTCTTACTGAGTCTAACAAAAAGTTTCTGCCTTTCCAACAATTTGGCAGAGACATTGCTGACACTACTGATGCTGTCCGTGATCCACAGACACTTGAGATTCTTGACATTACACCATGTTCTTTTGGTGGTGTCAGTGTTATAACACCAGGAACAAATACTTCTAACCAGGTTGCTGTTCTTTATCAGGATGTTAACTGCACAGAAGTCCCTGTTGCTATTCATGCAGATCAACTTACTCCTACTTGGCGTGTTTATTCTACAGGTTCTAATGTTTTTCAAACACGTGCAGGCTGTTTAATAGGGGCTGAACATGTCAACAACTCATATGAGTGTGACATACCCATTGGTGCAGGTATATGCGCTAGTTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGCTAGTCAATCCATCATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTACTCTAATAACTCTATTGCCATACCCACAAATTTTACTATTAGTGTTACCACAGAAATTCTACCAGTGTCTATGACCAAGACATCAGTAGATTGTACAATGTACATTTGTGGTGATTCAACTGAATGCAGCAATCTTTTGTTGCAATATGGCAGTTTTTGTACACAATTAAACCGTGCTTTAACTGGAATAGCTGTTGAACAAGACAAAAACACCCAAGAAGTTTTTGCACAAGTCAAACAAATTTACAAAACACCACCAATTAAAGATTTTGGTGGTTTTAATTTTTCACAAATATTACCAGATCCATCAAAACCAAGCAAGAGGTCATTTATTGAAGATCTACTTTTCAACAAAGTGACACTTGCAGATGCTGGCTTCATCAAACAATATGGTGATTGCCTTGGTGATATTGCTGCTAGAGACCTCATTTGTGCACAAAAGTTTAACGGCCTTACTGTTTTGCCACCTTTGCTCACAGATGAAATGATTGCTCAATACACTTCTGCACTGTTAGCGGGTACAATCACTTCTGGTTGGACCTTTGGTGCAGGTGCTGCATTACAAATACCATTTGCTATGCAAATGGCTTATAGGTTTAATGGTATTGGAGTTACACAGAATGTTCTCTATGAGAACCAAAAATTGATTGCCAACCAATTTAATAGTGCTATTGGCAAAATTCAAGACTCACTTTCTTCCACAGCAAGTGCACTTGGAAAACTTCAAGATGTGGTCAACCAAAATGCACAAGCTTTAAACACGCTTGTTAAACAACTTAGCTCCAATTTTGGTGCAATTTCAAGTGTTTTAAATGATATCCTTTCACGTCTTGACAAAGTTGAGGCTGAAGTGCAAATTGATAGGTTGATCACAGGCAGACTTCAAAGTTTGCAGACATATGTGACTCAACAATTAATTAGAGCTGCAGAAATCAGAGCTTCTGCTAATCTTGCTGCTACTAAAATGTCAGAGTGTGTACTTGGACAATCAAAAAGAGTTGATTTTTGTGGAAAGGGCTATCATCTTATGTCCTTCCCTCAGTCAGCACCTCATGGTGTAGTCTTCTTGCATGTGACTTATGTCCCTGCACAAGAAAAGAACTTCACAACTGCTCCTGCCATTTGTCATGATGGAAAAGCACACTTTCCTCGTGAAGGTGTCTTTGTTTCAAATGGCACACACTGGTTTGTAACACAAAGGAATTTTTATGAACCACAAATCATTACTACAGACAACACATTTGTGTCTGGTAACTGTGATGTTGTAATAGGAATTGTCAACAACACAGTTTATGATCCTTTGCAACCTGAATTAGACTCATTCAAGGAGGAGTTAGATAAATATTTTAAGAATCATACATCACCAGATGTTGATTTAGGTGACATCTCTGGCATTAATGCTTCAGTTGTAAACATTCAAAAAGAAATTGACCGCCTCAATGAGGTTGCCAAGAATTTAAATGAATCTCTCATCGATCTCCAAGAACTTGGAAAGTATGAGCAGTATATAAAATGGCCATGGTACATTTGGCTAGGTTTTATAGCTGGCTTGATTGCCATAGTAATGGTGACAATTATGCTTTGCTGTATGACCAGTTGCTGTAGTTGTCTCAAGGGCTGTTGTTCTTGTGGATCCTGCTGCAAATTTGATGAAGACGACTCTGAGCCAGTGCTCAAAGGAGTCAAATTACATTACACATAA'
# test_count_DCR_S = count_DCR(Seq_S,freq = True)
# print (test_count_DCR_S.shape)
################################################## test 
################################################## test
