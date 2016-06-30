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

queryID = []
fasta = []
with open(blastout, 'r') as f:
    queryID = [line[:-1] for line in f]

with open(fa) as f:
    record_iter = SeqIO.parse(f,"fasta")
    for record in record_iter:
        if str(record.id) in queryID:
            fasta.append(">"+record.id)
            fasta.append(record.seq)

# with open('output/filtered_blast.txt', "w") as file:
#     for item in blast:
#         file.write("%s" % item)
#
with open('output/filtered_fasta.'+num+'.fasta', "w") as file:
    for item in fasta:
        file.write("%s\n" % item)
