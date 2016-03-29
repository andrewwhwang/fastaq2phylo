from argparse import ArgumentParser
import random
import sys
parser = ArgumentParser()
parser.add_argument('-file', help="txt file to be parsed")
parser.add_argument('-num', help="number of random sequences selected")
args = parser.parse_args()

filename = args.file
num = int(args.num)

seqID = {}
num_lines = int(sum(.5 for line in open(filename)))
rand = random.sample(range(num_lines), num)
with open(filename, 'r') as f:
    lines=f.readlines()
    for num in rand:
        sys.stdout.write(lines[num*2])
        sys.stdout.write(lines[num*2+1])