from argparse import ArgumentParser
from Bio import Entrez
from ete3 import NCBITaxa
ncbi = NCBITaxa()
Entrez.email = "A.N.Other@example.com"
parser = ArgumentParser()
parser.add_argument('-file', help="txt file to be parsed")
parser.add_argument('-dbType', help="type of database")
parser.add_argument('-filenum', help="output name")
args = parser.parse_args()

filename = args.file
dbType = args.dbType
#filename = 'C:/Users/Andrew.Hwang/Desktop/one_click/output/blastout.txt'
#dbType = 'nt'

memory = {}
writeLines = []
with open(filename, 'r') as f:
    for line in f:
        temp = line.split(",")
        ID=temp[1].split(';')[0]
        pos=int(round(100*((int(temp[2])+int(temp[3]))/2.0) / int(temp[4])))
        if dbType == 'viruses':
            if not memory.has_key(ID):
                handle = Entrez.efetch(db='nucleotide', id=ID, retmode="xml")
                records = Entrez.read(handle)
                lineage = records[0]['GBSeq_taxonomy']
                memory[ID] = lineage
                writeLines.append( lineage + "::" + str(pos))
            else:
                writeLines.append( memory[ID] + "::" + str(pos))
        elif dbType == 'nt':
            lineage = ncbi.get_lineage(ID)
            names = ncbi.get_taxid_translator(lineage)
            strList = [str(names[taxid]) for taxid in lineage]
            writeLines.append( '; '.join(strList) + "::" + str(pos))
			
with open('output/lineage.'+filenum+'.txt', 'w') as the_file:
    for line in writeLines:
        the_file.write(line)
        