from collections import defaultdict
from ete3 import Tree, TreeStyle, NodeStyle, ClusterTree, ProfileFace
from ete3.treeview.faces import add_face_to_node
import re
from argparse import ArgumentParser
import numpy as np
import matplotlib.cm as cm
from math import log

parser = ArgumentParser()
parser.add_argument('-file', help="lineage file to be parsed")
parser.add_argument('-thres', help="threshold of num of reads")
parser.add_argument('-param', help="list of params. Used to make filenames")
parser.add_argument('-samples', help="number of times run")
args = parser.parse_args()

filename = args.file
thres = float(args.thres)
param = args.param
samples = int(args.samples)
#filename = 'C:/Users/Andrew.Hwang/Desktop/one_click/output/lineage.txt'
#thres = 0
def tree(): 
    return defaultdict(tree)

def add(t, path):
    for node in path:
        t = t[node]
    t[path[-1]] = path[-1]

def parseLineage(filename):
    count = {}
    with open(filename, 'r') as f:
        for line in f:
            lineage = line.split('::')[0]
            pos=line.split('::')[1]
            pos = int(pos[:-1]) #gets rid of \n character
            count.setdefault(lineage, []).append(pos)
    return count

def getSuffixandMatrixandNewick(count):
    matrix = '#NAMES\t0\t5\t10\t15\t20\t25\t30\t35\t40\t45\t50\t55\t60\t65\t70\t75\t80\t85\t90\t95\n'
    total = float(sum([len(x) for x in count.values()]))
    suffix = {}
    taxo = tree()
    regex = re.compile('[^a-zA-Z0-9- ]')
    for key in count.keys():
        fullname = key.split('; ')
        fullname = [regex.sub('', x) for x in fullname] 
        species = fullname[-1]
        num = len(count[key])
        
        #getting suffix
        percent = "%.2f" % (100*num/total) + '%'
        suffix[species] = ' %s{%s}' % (num/samples,percent)
        
        #getting matrix
        hist,junk = np.histogram(count[key],range(0,105,5),density=False)
        hist = np.array([.01 if x==0 else x for x in hist])
        hist = species + suffix[species]+ '\t' + '\t'.join([str(x) for x in hist]) + '\n'
        matrix += hist
        
        #getting newick
        if float(percent[:-1]) > thres:
            add(taxo,fullname)
    return suffix, matrix, taxo

def addColors(t):
    pattern = re.compile('[0-9]+\.[0-9]+')
    maxPer = max([float(pattern.findall(node.name)[0]) for node in t.get_leaves()])
    minPer = min([float(pattern.findall(node.name)[0]) for node in t.get_leaves()])
    dif = max(1,maxPer-minPer)
    for n in t.get_leaves():
        nst1 = NodeStyle() 
        nst1["bgcolor"] = getColorStr((float(pattern.findall(n.name)[0])-minPer)/dif)
        n.set_style(nst1)
        
def getColorStr(x):
    cmap = cm.get_cmap('Blues')
    y = lambda x: .15 * (log(x+.018)+4)
    r,g,b,a = cmap(y(x))
    #y = (ln(x+.018) + 4) * .15
    r=int(r*255)
    g=int(g*255)
    b=int(b*255)
    return '#' + "%0.2X" % r + "%0.2X" % g + "%0.2X" % b

def convert(taxo,suffix):
    newickStr = []
    for node, children in taxo.iteritems():
        nodename = str(node)
        if type(children)==defaultdict:
            newickStr.append(convert(children,suffix)+nodename)
        else:
            newickStr.append(nodename + suffix[nodename])
    return '('+','.join(newickStr)+')'

count = parseLineage(filename)
suffix, matrix, taxo = getSuffixandMatrixandNewick(count)
newick = convert(taxo,suffix)
newick += ';'

t = Tree(newick, format=1)
ct = ClusterTree(t.write(), matrix)


addColors(ct)


# nodes are linked to the array table
array = ct.arraytable

# Calculates some stats on the matrix. Needed to establish the color
# gradients.
matrix_dist = [i for r in xrange(len(array.matrix))for i in array.matrix[r] if np.isfinite(i)]
matrix_max = np.max(matrix_dist)
matrix_min = np.min(matrix_dist)
matrix_avg = matrix_min+((matrix_max-matrix_min)/2)

# Creates a profile face that will represent node's profile as a
# heatmap
profileFace  = ProfileFace(matrix_max, matrix_min, matrix_avg, 200, 14, "heatmap",colorscheme=3)
#nameFace = AttrFace("name", fsize=8)
# Creates my own layout function that uses previous faces
def mylayout(node):
    # If node is a leaf
    if node.is_leaf():
        # And a line profile
        add_face_to_node(profileFace, node, 0, aligned=True)
        node.img_style["size"]=2
#        add_face_to_node(nameFace, node, 1, aligned=True)

# Use my layout to visualize the tree
ts = TreeStyle()
ts.layout_fn = mylayout
#ct.show(tree_style=ts)
t.write(format=9, outfile="output/newick/"+param+".nw")
ct.render('output/pngs/'+param+'.png',tree_style=ts)
#print t
