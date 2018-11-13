import numpy as np
import pandas as pd
import os
import sys

datadir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL/data_new/pairs_ciseQTL/input'
outdir = '/work-zfs/abattle4/heyuan/forKaushik/DL_project'
gene_annotation = pd.read_csv('/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt',sep='\t', index_col=0)



Y = pd.read_csv('%s/v7_eQTL_ciseQTL_slope_ownPairs.txt' % datadir, sep='\t', index_col = [0, 1], header=None)
Pairs = np.array(Y.index)
genes = [p[0] for p in Pairs]
SNPs  = [p[1] for p in Pairs]
Y.index = pd.MultiIndex.from_arrays([genes, SNPs])

tissues = pd.read_csv('tissues.txt', header=None)
Y.columns = np.array(tissues[0])

Pairs = np.array(Y.index)

Brain_Pairs = Y[[t for t in Y.columns if ('Brain' in t) and ('Brain_Cere' not in t)]]
Non_Brain_Pairs = Y[[t for t in Y.columns if ('Brain' not in t)]]



brain_idx = np.where(np.sum(Brain_Pairs != 0, axis=1)  >= (Brain_Pairs.shape[1]/3))[0]
brain_idx2 = np.where(np.sum(Non_Brain_Pairs != 0, axis=1)  == 0)[0]
brain_pairs = Pairs[np.intersect1d(brain_idx, brain_idx2)]

print '%d brain pairs ' % len(brain_pairs)

fn = open('%s/Brain_pairs.txt' % outdir, 'w')
fn.write('SNP.ID\tSNP.chr\tSNP.location\tgene.ID\tgene.chr\tgene.start_pos\tgene.end_pos\tgene.name\n')
for p in brain_pairs:
    gene = np.array(gene_annotation.loc[p[1]])
    fn.write('%s\tchr%s\t%s\t%s\t%s\t%d\t%d\t%s\n' % (p[0], p[0].split('_')[0], p[0].split('_')[1], p[1], gene[0], gene[1], gene[2], gene[4]))
fn.close()




non_brain_idx = np.where(np.sum(Non_Brain_Pairs != 0, axis=1)  >= (Non_Brain_Pairs.shape[1]/2))[0]
non_brain_idx2 = np.where(np.sum(Brain_Pairs != 0, axis=1)  == 0)[0]
non_brain_pairs = Pairs[np.intersect1d(non_brain_idx, non_brain_idx2)]

print '%d non-brain pairs ' % len(non_brain_pairs)

fn = open('%s/Non_Brain_pairs.txt' % outdir, 'w')
fn.write('SNP.ID\tSNP.chr\tSNP.location\tgene.ID\tgene.chr\tgene.start_pos\tgene.end_pos\tgene.name\n')
for p in non_brain_pairs:
    gene = np.array(gene_annotation.loc[p[1]])
    fn.write('%s\tchr%s\t%s\t%s\t%s\t%d\t%d\t%s\n' % (p[0], p[0].split('_')[0], p[0].split('_')[1], p[1], gene[0], gene[1], gene[2], gene[4]))
fn.close()
