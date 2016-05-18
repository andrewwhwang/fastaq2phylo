from Bio import SeqIO
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-blastout', help="blastout file")
parser.add_argument('-fa', help="fasta file")
args = parser.parse_args()

blastout = args.blastout
fa = args.fa
human = []
blast = []
fasta = []
with open(blastout, 'r') as f:
    for line in f:
        items = line.split(",")
        if items[1] == "9606":
            human.append(items[0])
        else:
            blast.append(line)
            
with open(fa) as f:
    record_iter = SeqIO.parse(f,"fasta")
    for read in record_iter:
        if read.id not in human:
            fasta.append(">"+read.id)
            fasta.append(read.seq)
               
with open('output/filtered_blast.txt', "w") as file:
    for item in blast:
        file.write("%s" % item)
        
with open('output/filtered_fasta.fasta', "w") as file:
    for item in fasta:
        file.write("%s\n" % item)