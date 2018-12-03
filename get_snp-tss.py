#Kaushik Srinivasan, used to isolate SNP and TSS regions

#imports
import pandas as pd
import numpy as np

#definitions to get genetic strings from table given a row
def get_gene_bp(row):
    return bp_tss[bp_tss[0] == row[1]].values[0][1]

def get_gene_nbp(row):
    return nbp_tss[nbp_tss[0] == row[1]].values[0][1]

#read data files
bp_snp = pd.read_csv('gene-datasets/bp_SNP_4.csv', header=None)
bp_tss = pd.read_csv('gene-datasets/bp_TSS.csv', header=None)

nbp_snp = pd.read_csv('gene-datasets/nbp_SNP_4.csv', header=None)
nbp_tss = pd.read_csv('gene-datasets/nbp_TSS.csv', header=None)

#Applying series to get all associatetd genetic code
bp_snp[2] = bp_snp.apply(get_gene_bp, axis=1)
nbp_snp[2] = nbp_snp.apply(get_gene_nbp, axis=1)

#isolating rows
bp_snp_tss = bp_snp[[0, 2]]
nbp_snp_tss = nbp_snp[[0, 2]]

#writing to file
bp_snp_tss.to_csv('gene-datasets/bp_snp-tss.csv', header=False, index=False)
nbp_snp_tss.to_csv('gene-datasets/nbp_snp-tss.csv', header=False, index=False)
