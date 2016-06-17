from argparse import ArgumentParser
from Bio import SeqIO
import random
parser = ArgumentParser()
parser.add_argument('-file', help="txt file to be parsed")
parser.add_argument('-num', help="number of random sequences selected")
parser.add_argument('-total', help="total number of sequences in fasta")
parser.add_argument('-sampleNum', help="sample number for filename")
args = parser.parse_args()

filename = args.file
sampleNum = args.sampleNum
num = int(args.num)
total = int(args.total)


rand = random.sample(range(total), num)
lines = []
with open(filename) as f:
    record_iter = SeqIO.parse(f,"fasta")
    for count, i in enumerate(record_iter):
        if count in rand:
            lines.append(">"+i.id)
            lines.append(i.seq)

with open(filename[:-12]+sampleNum+'.fasta', 'w') as f:
    f.writelines("%s\n" % line for line in lines)
