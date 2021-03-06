from Bio import SeqIO
from argparse import ArgumentParser
from math import ceil
parser = ArgumentParser()
parser.add_argument('-file', help="fasta file to be split")
parser.add_argument('-num', help="number of fasta files to be split into")
parser.add_argument('-total', help="number of sequences")
parser.add_argument('-filenum', help="fasta file's index #")
args = parser.parse_args()

filename = args.file
filenum = args.filenum
num = float(args.num)
total = float(args.total)


def batch_iterator(iterator, batch_size) :
    entry = True
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch



with open(filename) as f:
    seqPer = ceil(total/num)
    record_iter = SeqIO.parse(f,"fasta")
    for i, batch in enumerate(batch_iterator(record_iter, seqPer)):
        fileOut = "output/%s.%i.fasta" % (filenum,i)
        handle = open(fileOut, "w")
        count = SeqIO.write(batch, handle, "fasta")
        handle.close()
        print "Wrote %i reads to %s" % (count, fileOut)
