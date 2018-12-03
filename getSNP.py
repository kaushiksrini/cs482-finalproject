import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from urllib.request import urlopen
import re
import csv
import sys

def getSNPSequence(chrom_number, startLoc, endLoc):
    url = "http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment={}:{},{}.html".format(chrom_number, startLoc, endLoc)
    page = urlopen(url)
    soup = BeautifulSoup(page, 'html.parser')
    geneCode = re.sub("\n|\r|\t", "", soup.dna.text[1:-1])
    return geneCode

def getSNP(dataName, listName, cutSize):
     for i in range(dataName.shape[0]):
            genome = getSNPSequence(dataName.iloc[i]['SNP.chr'], str(int(dataName.iloc[i]['SNP.location']) - cutSize), str(int(dataName.iloc[i]['SNP.location']) + cutSize))
            listName.append([genome, dataName.iloc[i]['gene.name']])
            if i % 50 == 0:
                print("Gene", i)

def writeFile(SNPList, filename):
    with open(filename, 'w') as f:
        for item in SNPList:
            f.write("%s,%s\n"%(item[0],item[1]))

_cutSize = int(sys.argv[1])

data_nbp = pd.read_csv("data/Non_Brain_pairs.txt", sep = '\t')
#data_bp = pd.read_csv("data/Brain_pairs.txt", sep = '\t')

#print("getting brain pairs")
#bp_SNP = []
#getSNP(data_bp, bp_SNP, _cutSize)

print("getting non-brain pairs")
nbp_SNP = []
getSNP(data_nbp, nbp_SNP, _cutSize)

print("writing to files")
#outfile_bp = "gene-datasets/bp_SNP_{}.csv".format(str(_cutSize))
outfile_nbp = "gene-datasets/nbp_SNP_{}.csv".format(str(_cutSize))

#writeFile(bp_SNP, outfile_bp)
writeFile(nbp_SNP, outfile_nbp)
