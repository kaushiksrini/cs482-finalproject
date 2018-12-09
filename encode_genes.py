import pandas as pd
import numpy as np
import sys

#get a genetic string and then
def get_kmers(string, k):
    return [string[start:start+k] for start in range(len(string) - k + 1)]

#gets the kmer along with the kmer dict and returns the encoding of it
def kmer_encode(kmer_dict, kmer):
    if kmer in kmer_dict:
        return kmer_dict[kmer]
    else:
        kmer_dict[kmer] = len(kmer_dict)
        return kmer_dict[kmer]

def read_dataframe(dataframe, column, k, dict):
    all_col = []
    for col in dataframe[column].values:
        col_list = get_kmers(col, k)
        encode = [kmer_encode(dict, kmer) for kmer in col_list]
        all_col.append(encode)
    return all_col

def write_file(encode_dict, filename):
    with open(filename, "w") as myfile:
        for key in encode_dict.keys():
            myfile.write("%s, %s\n" % (key, encode_dict[key]))

#begin of main program
bp = pd.read_csv("gene-datasets/bp_snp-tss.csv", names=("SNP", "TSS"))
nbp = pd.read_csv("gene-datasets/nbp_snp-tss.csv", names=("SNP", "TSS"))
kmer_size = int(sys.argv[1])
encode_dict_k = {}

bp_all_snp = read_dataframe(bp, 'SNP', kmer_size, encode_dict_k)
bp_all_tss = read_dataframe(bp, 'TSS', kmer_size, encode_dict_k)

nbp_all_snp = read_dataframe(nbp, 'SNP', kmer_size, encode_dict_k)
nbp_all_tss = read_dataframe(nbp, 'TSS', kmer_size, encode_dict_k)

all_list_bp = [bp_all_snp[i] + bp_all_tss[i] for i in range(len(bp_all_tss))]
all_list_nbp = [nbp_all_snp[i] + nbp_all_tss[i] for i in range(len(nbp_all_tss))]

np_bp_list = np.asarray(all_list_bp)
np_nbp_list = np.asarray(all_list_nbp)

np.savetxt("gene-datasets/bp_snp-tss-encoded-{}.csv".format(str(kmer_size)), np_bp_list, delimiter=",", fmt='%d')
np.savetxt("gene-datasets/nbp_snp-tss-encoded-{}.csv".format(str(kmer_size)), np_nbp_list, delimiter=",", fmt='%d')

write_file(encode_dict_k, "gene-datasets/kmers-{}-encoded.csv".format(str(kmer_size)))
