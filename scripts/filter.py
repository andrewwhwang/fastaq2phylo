from Bio import SeqIO
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-blastout', help="blastout file")
parser.add_argument('-fa', help="fasta file")
parser.add_argument('-num', help="number for naming")
args = parser.parse_args()

blastout = args.blastout
fa = args.fa
num = args.num

blastnum = []
fasta = []
with open(blastout, 'r') as f:
    blastnum = [line[:-1] for line in f]

with open(fa) as f:
    record_iter = SeqIO.parse(f,"fasta")
    for i, read in enumerate(record_iter):
        if str(i+1) in blastnum:
            fasta.append(">"+read.id)
            fasta.append(read.seq)

# with open('output/filtered_blast.txt', "w") as file:
#     for item in blast:
#         file.write("%s" % item)
#
with open('output/filtered_fasta.'+num+'.fasta', "w") as file:
    for item in fasta:
        file.write("%s\n" % item)
