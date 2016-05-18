from argparse import ArgumentParser
from Bio import SeqIO
parser = ArgumentParser()
parser.add_argument('-file', help="txt file to be parsed")
args = parser.parse_args()
#makes fasta such that description and sequence are one line each
filename = args.file
with open(filename) as f:
    record_iter = SeqIO.parse(f,"fasta")
    for i in record_iter:
        print ">"+i.id.replace(" ", "")
        print i.seq
