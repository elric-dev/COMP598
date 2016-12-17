'''
@uthor: Haji Mohammad Saleem
Date  : December 9, 2016

Objective 8

Implementation of the Sankoff Algorithm.
Regular: penalties -1 intra group and -2 inter group.
Updated: penalties -1 intra group and -2 inter group and -2 for gaps.
'''

import sys
sys.dont_write_bytecode = True
import os
import paths
from Bio import Phylo
from cStringIO import StringIO
from stockholm_parser import Stockholm
from pprint import pprint
import viennaRNAwrap as vwrap

DATA_PATH = paths.data_path
SEED_PATH = os.path.join(DATA_PATH, 'seed')
TREE_PATH = os.path.join(DATA_PATH, 'tree')


def heading(text):
    print '{:40}{}'.format('', text)
    return


def subheading(text):
    print '{:20}{}'.format('', text)
    return


def hline():
    print '-' * 100
    return


class Sankoff:

    def __init__(self, family, stk):
        '''
        Initialise the Sankoff Class
        '''
        self.family = family
        self.neuclicAcids = ['A', 'C', 'G', 'U']
        self.costMatrix = [[0, 2, 1, 2], [
            2, 0, 2, 1], [1, 2, 0, 2], [2, 1, 2, 0]]
        self.purines = ['A', 'G']
        self.pyrimidines = ['C', 'U']
        self.SS_cons = stk.SS_cons
        self.SS_cons_INFO = None
        self.SQ_align = stk.SQ_align
        self.treeDict = None  # Dict of tree parent and their children
        self.ancestorWeightDict = None  # Dict with the weight matrix
        self.ancestorSeq = None  # Dict with ancestors and their sequences
        self.secStructure = None
        self.consDist = None
        return

    def updateSSCONS(self):
        '''
        Converting the consensus structure to dot bracket notation
        '''
        ss_cons = self.SS_cons
        ss_cons = ss_cons.replace('<', '(')
        ss_cons = ss_cons.replace('>', ')')
        ss_cons = ss_cons.replace('-', '.')
        ss_cons = ss_cons.replace('_', '.')
        ss_cons = ss_cons.replace(':', '.')
        ss_cons = ss_cons.replace(',', '.')
        self.SS_cons = ss_cons
        return

    def includeGaps(self):
        '''
        Update the neuclic acids and cost matrix to include gaps
        '''
        self.neuclicAcids.append('-')
        self.costMatrix = [[0, 2, 1, 2, 2], [2, 0, 2, 1, 2], [
            1, 2, 0, 2, 2], [2, 1, 2, 0, 2], [2, 2, 2, 2, 0]]
        return

    def readTree(self):
        '''
        Read a tree sequence file and convert it into a dictionary
        '''
        tree_file = os.path.join(TREE_PATH, self.family + '.nhx')
        with open(tree_file, 'r') as fin:
            tree_str = fin.read().strip()

        tree_list = []
        node_str = ''
        for char in tree_str:
            if char == ')' or char == '(' or char == ',':
                if node_str:
                    if '_' in node_str:
                        node_str = node_str.split('_')[1]
                        tree_list.append(node_str)
                    node_str = ''
                tree_list.append(char)
            else:
                node_str += char

        # converting the tree requence into a dictionary where parents are key
        # and left and right child are given in a list as value.
        treeDict = {}
        seqStack = []
        ancestor = 0

        for item in tree_list:
            seqStack.append(item)
            if item == ')':
                if seqStack[-5] == '(':
                    ancestor += 1
                    join = seqStack[-5:]
                    seqStack = seqStack[:-5]
                    seqStack.append(ancestor)
                    treeDict[ancestor] = [join[1], join[3]]
                else:
                    if len(seqStack) == 7:
                        node1 = seqStack[1]
                        node2 = seqStack[3]
                        node3 = seqStack[5]
                        ancestor += 1
                        treeDict[ancestor] = [node1, node2]
                        final_ancestor = ancestor + 1
                        treeDict[final_ancestor] = [ancestor, node3]
        self.treeDict = treeDict
        return

    def drawTree(self):
        '''
        Draw the phylogeny tree
        '''
        tree_file = os.path.join(TREE_PATH, self.family + '.nhx')
        with open(tree_file, 'r') as fin:
            tree_str = fin.read().strip()

        tree_list = []
        node_str = ''
        for char in tree_str:
            if char == ')' or char == '(' or char == ',':
                if node_str:
                    if '_' in node_str:
                        node_str = node_str.split('_')[1]
                        tree_list.append(node_str)
                    node_str = ''
                tree_list.append(char)
            else:
                node_str += char

        tree_str = ''.join(tree_list)
        tree = Phylo.read(StringIO(tree_str), "newick")
        tree.rooted = True
        subheading('%s Family Tree' % self.family)
        hline()
        Phylo.draw_ascii(tree)
        hline()
        return

    def familyStats(self):
        seqLength = len(self.SS_cons)
        leafCount = len(self.SQ_align.keys())
        ancsCount = len(self.treeDict.keys())
        subheading(
            'Leaf Nodes: %s, Ancestor Nodes: %s, Sequence Length: %s' %
            (leafCount, ancsCount, seqLength))
        hline()
        return

    def printAncestors(self):
        heading('Ancestors')
        hline()
        print '{:40} {:40} {:40}'.format('Ancestor', 'Left Child', 'Right Child')
        for key in sorted(self.treeDict.keys()):
            print '{:40} {:40} {:40}'.format(str(key), str(self.treeDict[key][0]), str(self.treeDict[key][1]))
        hline()
        return

    def printLeafSeq(self):
        heading('Leaf Nodes')
        hline()
        print "{:30} {}".format('Accession', 'Sequence')
        for key in sorted(self.SQ_align.keys()):
            print "{:30} {}".format(key, self.SQ_align[key])
        hline()
        return

    def printCostMat(self):
        heading('Cost Matrix')
        hline()
        cm = self.costMatrix
        na = self.neuclicAcids
        row1 = [' ']
        row1.extend(na)
        subheading('\t'.join(row1))
        for i in xrange(len(na)):
            cost = [str(x) for x in cm[i]]
            row = [na[i]]
            row.extend(cost)
            subheading('\t'.join(row))
        hline()
        return

    def getAncestorWeight(self):
        '''
        Calculate the ancestoral weight matrix for each possible neucleotide
        '''
        seqLength = len(self.SS_cons)
        ancestorWeightDict = {}

        # initializing the leaf nodes
        n = len(self.neuclicAcids)
        for node in self.SQ_align.iterkeys():  # for all leaf nodes
            nodeCost = []
            for neuclotide in self.SQ_align[
                    node]:  # for all the neucleotides in the sequence
                nucCost = []
                for k in xrange(0, n):  # for each possible neucleotide
                    if k == self.neuclicAcids.index(neuclotide):
                        nucCost.append(0)
                    else:
                        nucCost.append(float("inf"))
                nodeCost.append(nucCost)
            ancestorWeightDict[node] = nodeCost

        # calculating the ancestral sequence weights
        for ancestor in sorted(self.treeDict.keys()):
            leftChild = self.treeDict[ancestor][0]
            rghtChild = self.treeDict[ancestor][1]

            leftWeight = ancestorWeightDict[leftChild]
            rghtWeight = ancestorWeightDict[rghtChild]

            ancsWeight = []
            for i in xrange(
                    seqLength):  # for all the neucleotides in the sequence
                nucCost = []
                for k in xrange(0, n):  # for each possible neucleotide
                    minLeft = float("inf")
                    for l in range(0, n):  # for each possible neucleotide
                        thisCost = self.costMatrix[k][l] + leftWeight[i][l]
                        minLeft = min(minLeft, thisCost)
                    minRght = float("inf")
                    for l in range(0, n):  # for each possible neucleotide
                        thisCost = self.costMatrix[k][l] + rghtWeight[i][l]
                        minRght = min(minRght, thisCost)
                    nucCost.append(minLeft + minRght)
                ancsWeight.append(nucCost)
            ancestorWeightDict[ancestor] = ancsWeight
        self.ancestorWeightDict = ancestorWeightDict
        return

    def getAncestorSeq(self):
        '''
        Calculating the ancestral sequence
        '''
        seqLength = len(self.SS_cons)
        ancestorSeq = {}
        for node in sorted(self.ancestorWeightDict.keys()):
            sequence = ""
            if isinstance(node, int):
                for i in xrange(seqLength):
                    nucCost = self.ancestorWeightDict[node][i]
                    minIdx = nucCost.index(min(nucCost))
                    sequence += self.neuclicAcids[minIdx]
                ancestorSeq[node] = sequence
        self.ancestorSeq = ancestorSeq
        return

    def getSSinfo(self):
        '''
        Retrieve the base pair and bond site information
        '''
        seqLength = len(self.SS_cons)
        stack = []
        index_list = []
        index_dict = {}

        for i in xrange(seqLength):
            char = self.SS_cons[i]
            if char == "(":
                stack.append(i)
                index_list.append(i)
            if char == ")":
                index_dict[i] = stack.pop(-1)
                index_list.append(i)

        self.SS_cons_INFO = {'list': index_list, 'dict': index_dict}
        return

    def getAncestorSeq_extended(self):
        '''
        Calculating the ancestral sequence while trying to preserve the base pairing dependencies
        of the consensus structure by forcing the neucleotide at the closing binding index to be
        something that can form a pair with corresponding opening binding index.

        We have to make sure of two things:
            1. There are no gaps at the indices of binding pair
            2. The closing binding pair is aligned with opening binding pair, i.e., A-U, G-C and G-U
        '''
        allowed_pair = {}
        allowed_pair["A"] = [3]
        allowed_pair["U"] = [0, 2]
        allowed_pair["G"] = [1, 3]
        allowed_pair["C"] = [2]

        BP_list = self.SS_cons_INFO['list']
        BP_dict = self.SS_cons_INFO['dict']

        seqLength = len(self.SS_cons)
        ancestorSeq = {}
        for node in sorted(self.ancestorWeightDict.keys()):
            sequence = ""
            if isinstance(node, int):
                for i in xrange(seqLength):
                    nucCost = self.ancestorWeightDict[node][i]
                    if i in BP_list:  # i is part of a base pair in the consensus structure
                        nucCost.pop(-1)  # removing the gaps from such index
                    if i in BP_dict.iterkeys():  # i is a closing base pair
                        openingidx = BP_dict[i]
                        # neucleotide at opening index
                        openingnuc = sequence[openingidx]

                        if len(allowed_pair[openingnuc]) == 2:
                            a = allowed_pair[openingnuc][0]
                            b = allowed_pair[openingnuc][1]
                            if nucCost[a] > nucCost[b]:
                                sequence += self.neuclicAcids[b]
                            else:
                                sequence += self.neuclicAcids[a]

                        else:
                            sequence += self.neuclicAcids[
                                allowed_pair[openingnuc][0]]
                    else:
                        minIdx = nucCost.index(min(nucCost))
                        sequence += self.neuclicAcids[minIdx]
                ancestorSeq[node] = sequence
        self.ancestorSeq = ancestorSeq
        return

    def printAncsSeq(self):
        heading('Ancestor Nodes')
        hline()
        print "{:30} {}".format('Accession', 'Sequence')
        for key in sorted(self.ancestorSeq.keys()):
            print "{:30} {}".format(str(key), self.ancestorSeq[key])
        hline()
        return

    def getSecondaryStruct(self):
        secStructure = {}
        consDist = {}
        ss_cons = self.SS_cons

        seqLen = len(ss_cons)
        seedid = self.family
        ancesCount = len(self.ancestorSeq.keys())
        distsum = 0
        GCcont = 0
        freqcont = 0.

        for key in sorted(self.ancestorSeq.keys()):
            a, b = vwrap.runRNAFold(self.ancestorSeq[key])
            c = vwrap.runRNADistance(ss_cons, a)
            d = vwrap.mfeRNAFold(self.ancestorSeq[key])
            secStructure[key] = a
            consDist[key] = c
            distsum += int(c)
            freqcont += d
            temp_seq = self.ancestorSeq[key]
            for char in temp_seq:
                if char == "G" or char == "C":
                    GCcont+=1
        print '{:20}{:20}{:20}{:20}{:20}{:20}'.format(seedid, str(ancesCount), str(seqLen), str(distsum), str(GCcont), str(freqcont))
        #print '%s\t%s\t%s\t%s\t%s\t%s'%(seedid, ancesCount, seqLen, distsum, GCcont, freqcont)
        self.secStructure = secStructure
        self.consDist = consDist
        return

    def printSecondaryStruct(self):
        heading('Secondary Structures')
        hline()
        print "{:30} {}".format('Accession', 'Sequence')
        for key in sorted(self.ancestorSeq.keys()):
            print "{:30} {}".format(str(key), self.ancestorSeq[key])
            print "Dist: {:24} {}".format(self.consDist[key], self.secStructure[key])
        hline()

def main(seedFamily):
    # Read and parse the Stockholm file
    stk = Stockholm(seedFamily)
    stk.parse()

    # Parse the phylogeny tree and run Sankoff algorithm
    skf = Sankoff(seedFamily, stk)
    skf.updateSSCONS()
    skf.includeGaps()
    skf.readTree()
    skf.getSSinfo()
    skf.getAncestorWeight()
    skf.getAncestorSeq_extended()
    #skf.getAncestorSeq()
    skf.getSecondaryStruct()
    return


if __name__ == "__main__":
    seeds = os.listdir(SEED_PATH)
    seeds = [x.split('.')[0] for x in seeds]

    
    print '{:20}{:20}{:20}{:20}{:20}{:20}'.format('RNA family', '# of Ancestors', 'Sequence Length', 'RNA Distance', 'GC content', 'MFE freq')
    for seed in seeds:
        main(seed)
