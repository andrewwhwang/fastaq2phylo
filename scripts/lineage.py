from argparse import ArgumentParser
from Bio import Entrez
from ete3 import NCBITaxa
import sys
ncbi = NCBITaxa()
Entrez.email = "A.N.Other@example.com"
parser = ArgumentParser()
parser.add_argument('-file', help="txt file to be parsed")
parser.add_argument('-dbType', help="type of database")
parser.add_argument('-filenum', help="output name")
args = parser.parse_args()

filename = args.file
dbType = args.dbType
filenum  = args.filenum

# filename = 'C:/Users/Andrew.Hwang/Desktop/fastaq2phylo/output/blastout.txt'
# dbType = 'nt'
# filenum  = "0"

memory = {}
writeLines = []
with open(filename, 'r') as f:
    for line in f:
        line_arr = line.split("\t")
        ID=line_arr[1]
        pos=int(round(100*((int(line_arr[2])+int(line_arr[3]))/2.0) / int(line_arr[4])))
        if dbType == 'nt':
            lineage = ncbi.get_lineage(ID)
            names = ncbi.get_taxid_translator(lineage)
            strList = [str(names[taxid]) for taxid in lineage]
            writeLines.append( '; '.join(strList) + "::" + str(pos))
        elif dbType == 'viruses' or dbType == 'blood':
            if not memory.has_key(ID):
                handle = Entrez.efetch(db='nucleotide', id=ID, retmode="xml")
                records = Entrez.read(handle)
                lineage = records[0]['GBSeq_taxonomy']
                organism = records[0]['GBSeq_organism']
                lineage = "; ".join(lineage.split("; ")[:-1] + [organism])
                memory[ID] = lineage
            writeLines.append(memory[ID] + "::" + str(pos))
        else:
            writeLines.append(line_arr[5][:-1] + "::" + str(pos))

with open('output/lineage.'+filenum+'.txt', 'w') as f:
    f.writelines("%s\n" % line for line in writeLines)
