import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from urllib.request import urlopen
import re
import csv

# read gene coordinate data from files
data_nbp = pd.read_csv("data/Non_Brain_pairs.txt", sep = '\t')
data_bp = pd.read_csv("data/Brain_pairs.txt", sep = '\t')

# adds a genetic sequence given chromosome number, start location and end location.
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

#for each value in table, it gets the genes
def getGenes(varName, dictType):
    for i in range(varName.shape[0]):
        getSequence(varName.iloc[i]['gene.chr'], varName.iloc[i]['gene.start_pos'], varName.iloc[i]['gene.end_pos'], varName.iloc[i]['gene.name'], dictType)
        if i % 50 == 0:
            print("Gene", i)

#write to the file for a given filename
def writeFile(geneDict, filename):
    with open(filename, 'w') as f:
        for key in geneDict.keys():
            f.write("%s,%s\n"%(key,geneDict[key]))

bp_geneDict = {}
nbp_geneDict = {}

getGenes(data_nbp, nbp_geneDict)
getGenes(data_bp, bp_geneDict)

writeFile(bp_geneDict, 'bp_genes.csv')
writeFile(nbp_geneDict, 'nbp_genes.csv')
