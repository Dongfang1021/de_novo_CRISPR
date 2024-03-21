#!/usr/bin/env python
# coding: utf-8



import pandas as pd
import collections
import argparse
import re
from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.Restriction import BsaI
from Bio.Restriction import RestrictionBatch

parser = argparse.ArgumentParser()
parser.add_argument('--guide_file', help='gRNA file path')
argv = vars(parser.parse_args())

if argv['guide_file'] == None:
    	raise Exception ('You should provide guide file!')
else:
	guide_file=argv['guide_file'].strip()

output_path = "/".join(guide_file.strip().split('/')[:-1])
sample = "_".join(guide_file.strip().split('/')[-1].split('-')[:2])



exon = pd.read_csv('human_exon_table.xls', sep='\t')
guide = pd.read_csv(guide_file, sep='\t')
# guide['percentage'] = guide.apply(lambda x: x['position']/x['cds_len'], axis=1)
guide['gene_id'] = guide['Target Gene ID']
guide['gene_id'] = guide['gene_id'].astype(str)
guide['sgRNA Cut Position (1-based)'] = guide['sgRNA Cut Position (1-based)'].fillna(0)
guide['sgRNA Cut Position (1-based)'] = guide['sgRNA Cut Position (1-based)'].astype(int)
guide_num = guide.groupby('gene_id')['Off-Target Rank'].apply(list).reset_index(name='guide_number')
guide_num['guide_number'] = guide_num['guide_number'].apply(len)
guide = guide.merge(guide_num, on='gene_id', how='left')
guide = guide[guide['Target Cut %']>= 5]
guide = guide[guide['Target Cut %']<= 65]



designs = (((BsaI,), 'TGTGGAAAGGACGAAACACCG{}GTTTAAGAGCTATGCTGGAAGGAGA'),)

def design_test(n):
    def _test(seq):
        as_design = Seq(design.format(seq), IUPACAmbiguousDNA)
        is_cut = any(enzyme.search(as_design) for enzyme in enzymes)
        return not (pattern.search(seq) or is_cut)
    pattern = re.compile('(^[Gg]{5})|([Tt]{4})')
    enzymes, design = designs[n]
    return _test


guide['cloneable'] = guide['sgRNA Sequence'].apply(lambda x: design_test(0)(x))
guide = guide[guide['cloneable']==True]
def getGrafType(seq):
    """ check if a guide fulfills the criteria described in Graf et al, Cell Reports 2019
    returns "tt" or "ggc" depending on the motif found or None if none of them were found.
    """
    # guide ends with TTC or TTT
    #seq = seq.upper()
    #if seq.endswith("TTC") or seq.endswith("TTT"):
        #return "tt"

    seq = seq.upper()

    if seq.endswith("TTC") or seq.endswith("TTT"):
        return "tt"

    # the last 4 nucleotides contain only T and C and more than 2 Ts
    suffix = seq[-4:]
    if set(suffix)==set(["T", "C"]) and suffix.count("T")>=2:
        return "tt"

    # the last four positions contain "TT" and at last one more T or C 
    if "TT" in suffix and (suffix.count("T")>=3 or suffix.count("C")>=1):
        return "tt"

    # the guide ends with [AGT]GCC
    if seq.endswith("GCC") and suffix[-4] in ["A", "G", "T"]:
        return "ggc"

    # the guide ends with GCCT
    if seq.endswith("GCCT"):
        return "ggc"

    return "OK"

def gcContent(seq):
    """return GC content as a float"""
    c = 0
    for x in seq:
        if x in ['G', "C"]:
            c += 1
    return (float(c)/len(seq))


guide['getGrafType'] = guide['sgRNA Sequence'].apply(lambda x: getGrafType(x))
guide['grafType'] = guide['getGrafType'].apply(lambda x: 0 if x == 'tt' or x == 'ggc' else 1)
guide = guide[guide['grafType'] == 1]
guide['gcContent'] = guide['sgRNA Sequence'].apply(lambda x: gcContent(x))
guide = guide[guide['gcContent']<= 0.75]
guide = guide[guide['gcContent']>= 0.25]





exon_id = list(exon['CDS_id'])
gene_id = list(exon['gene_id'])
exon_gene_dict = dict(zip(exon_id, gene_id))



exon_dict = {}
exon_chr = list(exon['chr'])
exon_start = list(exon['start'])
exon_end = list(exon['end'])
exon_id = list(exon['CDS_id'])
for x, y, z, n in zip(exon_chr, exon_start, exon_end, exon_id):
    exon_location = str(x)+"_"+str(int(y))+'_'+str(int(z))
    exon_name = n
    exon_dict.setdefault(exon_location, []).append(exon_name)



def target_gene_list(gRNA_chr_start_end, dict_gene):
    blast_chr = gRNA_chr_start_end[0]
    start = gRNA_chr_start_end[1]
    list_gene = []
    for i in dict_gene:
        key_name = i
        i = i.strip().split('_')
        if str(blast_chr) == str("_".join(i[:2])) and (start <= max(int(i[-2]), int(i[-1]))) and ( start >= min(int(i[-2]), int(i[-1]))):
            list_gene.extend(dict_gene[key_name])
    return list(set(list_gene))


guide['perfect_position_list'] = guide[['Reference Sequence','sgRNA Cut Position (1-based)']].astype(str).agg('_'.join, axis=1)

#guide['perfect_position_list'] = guide['perfect_position_list'].apply(list)
# formated_guide = guide[guide['perfect_match_count'] >= 1]
formated_guide = guide


new_list = list(formated_guide['perfect_position_list'])
# print(new_list)
# new_list = [x[2:-2] for x in new_list]
# print(new_list)


gene_list = []
blast_gene_count = []

for i in new_list:
    i_list = i.split("', '")
    gene_list_name = []
    for each in i_list:
        chr_pos = "_".join(each.split('_')[:2])
        guide_pos = each.split('_')[-1]
        position = [chr_pos, int(guide_pos)]
        gene_list_name.extend(target_gene_list(position, exon_dict))
    target_gene_number = len(gene_list_name)
    gene_list.append(gene_list_name)
    blast_gene_count.append(target_gene_number)


formated_guide['blast_CDS_target'] = gene_list
formated_guide['target_CDS_number'] = blast_gene_count



def exon_to_gene(exon_list):
    gene_list = [exon_gene_dict[x] for x in exon_list]
    gene_set = list(set(gene_list))
    return gene_set
        

formated_guide['blast_gene_target'] = formated_guide['blast_CDS_target'].apply(lambda x: exon_to_gene(x))
formated_guide['target_gene_number'] = formated_guide['blast_gene_target'].apply(len)
data = formated_guide.fillna('exon_exon_junction')
#data = data[data['blast_CDS_target'] != 'exon_exon_junction']
data['transcript_id'] = data['Target Transcript']


exon['gene_id'] = exon['gene_id'].astype(str)
exon_gene = exon.groupby('gene_id')['CDS_id'].apply(list).reset_index(name='CDS_id_list')
#print(exon_gene[exon_gene['gene_id']=='127550'])
exon_gene['CDS_id_list'] = exon_gene['CDS_id_list'].apply(lambda x:list(set(x)))
transcript_gene = exon[['gene_id', 'CDS_id', 'transcript_id']].drop_duplicates()
transcript_gene = transcript_gene.merge(exon_gene, on='gene_id', how="left")




data = data.merge(transcript_gene, on='transcript_id', how='left')
#data = data[data['target_CDS_number']!=0]
#data = data.merge(exon_gene, on='gene_id', how='left')

def exon_to_transcript(exon_list, dict_exon_transcript):
    formated_transcript_list = []
    if len(exon_list)>0:
        for exon in exon_list:
            formated_transcript_list.extend(dict_exon_transcript[exon])
        return list(set(formated_transcript_list))


exon_id = list(exon['CDS_id'])
CDS_id = list(exon['CDS_id'])
exon_transcript_dict = {}
for x in range(len(exon_id)):
    if exon_id[x] not in exon_transcript_dict:
        exon_transcript_dict.setdefault(exon_id[x], []).append(CDS_id[x])
    else:
        exon_transcript_dict[exon_id[x]].append(CDS_id[x])
        

data['formated_CDS_list'] = data['blast_CDS_target'].apply(lambda x: exon_to_transcript(x, exon_transcript_dict) )
def fraction(CDS_id_list, formated_transcript_list):
    if not CDS_id_list and not formated_transcript_list:
        a = len(set(formated_transcript_list) & set(CDS_id_list))
        b = len(set(CDS_id_list))
        return float(a)/float(b)
    else:
        return 1




data['fraction'] = data.apply(lambda x: fraction(x['CDS_id_list'], x['formated_CDS_list']), axis=1)


def cellecta_score(off_rank, on_score, fraction, guide_number):
    return ((guide_number-off_rank+1)*on_score*fraction)/guide_number

data['cellecta_score'] = data.apply(lambda x: cellecta_score(x['Off-Target Rank'], x['On-Target Efficacy Score'], x['fraction'], x['guide_number']), axis=1)
data = data.sort_values(by=['Target Gene ID','cellecta_score'], ascending=False)
data.to_csv(output_path+'/'+sample+"_selected_guide.xls", sep='\t', index=False)
# data.to_csv(sample+"_selected_guide.xls", sep='\t', index=False)

