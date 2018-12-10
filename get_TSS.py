import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from urllib.request import urlopen
import re
import csv

def getSequence(chrom_number, startLoc, endLoc, geneName, geneDict):
    if geneName in geneDict:
        return
    else:
        url = "http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment={}:{},{}.html".format(chrom_number, startLoc, endLoc)
        page = urlopen(url)
        soup = BeautifulSoup(page, 'html.parser')
        geneCode = re.sub("\n|\r|\t", "", soup.dna.text[1:-1])
        geneDict[geneName] = geneCode
        return

def getTSS(varName, dictType, cutSize):
    for i in range(varName.shape[0]):
        getSequence(varName.iloc[i]['gene.chr'], str(int(varName.iloc[i]['gene.start_pos']) - cutSize), str(int(varName.iloc[i]['gene.start_pos']) + cutSize), varName.iloc[i]['gene.name'], dictType)
        if i % 50 == 0:
            print("Gene", i)

def writeFile(geneDict, filename):
    with open(filename, 'w') as f:
        for key in geneDict.keys():
            f.write("%s,%s\n"%(key,geneDict[key]))




data_nbp = pd.read_csv("data/Non_Brain_pairs.txt", sep = '\t')
data_bp = pd.read_csv("data/Brain_pairs.txt", sep = '\t')

bp_geneTSS = {}
nbp_geneTSS = {}

getTSS(data_bp, bp_geneTSS, 49)
writeFile(bp_geneTSS, 'gene-datasets/bp_TSS.csv')

getTSS(data_nbp, nbp_geneTSS, 49)
writeFile(nbp_geneTSS, 'gene-datasets/nbp_TSS.csv')
