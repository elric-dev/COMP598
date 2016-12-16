'''
@uthor: Haji Mohammad Saleem
Date  : December 9, 2016

Objective 2-5

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
    print '-'*100
    return


class Sankoff:

    def __init__(self, family, stk):
        '''
        Initialise the Sankoff Class
        '''
        heading('Sankoff Analysis for %s Family'%family)
        hline()
        self.family = family 
        self.neuclicAcids = ['A', 'C', 'G', 'U']
        self.costMatrix = [[0,2,1,2],[2,0,2,1],[1,2,0,2],[2,1,2,0]]
        self.purines = ['A', 'G']
        self.pyrimidines = ['C', 'U']
        self.SS_cons = stk.SS_cons
        self.SQ_align = stk.SQ_align
        self.treeDict = None #Dict of tree parent and their children
        self.ancestorSeq = None #Dict with ancestors and their sequences
        self.secStructure = None
        self.consDist = None
        return
   
    def includeGaps(self):
        '''
        Update the neuclic acids and cost matrix to include gaps
        '''
        self.neuclicAcids.append('-')
        self.costMatrix = [[0,2,1,2,2],[2,0,2,1,2],[1,2,0,2,2],[2,1,2,0,2],[2,2,2,2,0]]
        return
    
    def readTree(self):
        '''
        Read a tree sequence file and convert it into a dictionary
        '''
        tree_file = os.path.join(TREE_PATH, self.family+'.nhx')
        with open(tree_file, 'r') as fin:
            tree_str = fin.read().strip()
        
        tree_list = []
        node_str = ''
        for char in tree_str:
            if char == ')'  or char == '(' or char == ',':
                if node_str:
                    if '_' in node_str:
                        node_str = node_str.split('_')[1]
                        tree_list.append(node_str)
                    node_str = ''
                tree_list.append(char)
            else:
                node_str+=char

        #converting the tree requence into a dictionary where parents are key 
        #and left and right child are given in a list as value.
        treeDict = {}
        seqStack = []
        ancestor = 0
        
        for item in tree_list:
            seqStack.append(item)
            if item == ')':
                if  seqStack[-5] == '(':
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
                        final_ancestor = ancestor+1
                        treeDict[final_ancestor] = [ancestor, node3]
        self.treeDict = treeDict
        return
   
    def drawTree(self):
        '''
        Draw the phylogeny tree
        '''
        tree_file = os.path.join(TREE_PATH, self.family+'.nhx')
        with open(tree_file, 'r') as fin:
            tree_str = fin.read().strip()

        tree_list = []
        node_str = ''
        for char in tree_str:
            if char == ')'  or char == '(' or char == ',':
                if node_str:
                    if '_' in node_str:
                        node_str = node_str.split('_')[1]
                        tree_list.append(node_str)
                    node_str = ''
                tree_list.append(char)
            else:
                node_str+=char

        tree_str = ''.join(tree_list)
        tree = Phylo.read(StringIO(tree_str), "newick")
        tree.rooted = True
        subheading('%s Family Tree'%self.family)
        hline()
        Phylo.draw_ascii(tree)
        hline()
	return

    def familyStats(self):
        seqLength = len(self.SS_cons)
        leafCount = len(self.SQ_align.keys())
        ancsCount = len(self.treeDict.keys())
        subheading('Leaf Nodes: %s, Ancestor Nodes: %s, Sequence Length: %s'%(leafCount, ancsCount, seqLength))
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

    def getAncestorSeq(self):
        '''
        Calculate the ancestoral sequence
        '''
        seqLength = len(self.SS_cons)
        ancestorDict = {}
        
        #initializing the leaf nodes
        n = len(self.neuclicAcids)
        for node in self.SQ_align.iterkeys(): #for all leaf nodes
            nodeCost = []
            for neuclotide in self.SQ_align[node]: #for all the neucleotides in the sequence
                nucCost = []
                for k in xrange(0,n): #for each possible neucleotide
                    if k == self.neuclicAcids.index(neuclotide):
                        nucCost.append(0)
                    else:
                        nucCost.append(float("inf"))
                nodeCost.append(nucCost)
            ancestorDict[node] = nodeCost
        
        #calculating the ancestral sequence weights
        for ancestor in sorted(self.treeDict.keys()):
            leftChild = self.treeDict[ancestor][0]
            rghtChild = self.treeDict[ancestor][1]

            leftWeight = ancestorDict[leftChild]
            rghtWeight = ancestorDict[rghtChild]

            ancsWeight = []
            for i in xrange(seqLength): #for all the neucleotides in the sequence
                nucCost = []
                for k in xrange(0,n): #for each possible neucleotide
                    minLeft = float("inf")
                    for l in range(0, n): #for each possible neucleotide
                        thisCost = self.costMatrix[k][l] + leftWeight[i][l]
                        minLeft = min(minLeft, thisCost)
                    minRght = float("inf")
                    for l in range(0, n): #for each possible neucleotide
                        thisCost = self.costMatrix[k][l] + rghtWeight[i][l]
                        minRght = min(minRght, thisCost)
                    nucCost.append(minLeft + minRght)
                ancsWeight.append(nucCost)
            ancestorDict[ancestor] = ancsWeight

        #calculating the ancestral sequence
        ancestorSeq = {}
        for node in sorted(ancestorDict.keys()):
            sequence = ""
            if isinstance( node, int ):
                for i in xrange(seqLength):
                    nucCost = ancestorDict[node][i]
                    minIdx = nucCost.index(min(nucCost))
                    sequence+=self.neuclicAcids[minIdx]
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
        ss_cons = ss_cons.replace('<', '(')
        ss_cons = ss_cons.replace('>', ')')
        ss_cons = ss_cons.replace('-', '.')
        ss_cons = ss_cons.replace('_', '.')
        ss_cons = ss_cons.replace(':', '.')
        ss_cons = ss_cons.replace(',', '.')
        
        for key in sorted(self.ancestorSeq.keys()):
            a,b = vwrap.runRNAFold(self.ancestorSeq[key])
            c = vwrap.runRNADistance(ss_cons, a)
            secStructure[key] = a
            consDist[key] = c
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

if __name__ == "__main__":
    seedFamily = sys.argv[1]

    # Read and parse the Stockholm file
    stk = Stockholm(seedFamily)
    stk.parse()

    # Parse the phylogeny tree and run Sankoff algorithm 
    skf = Sankoff(seedFamily, stk)
    skf.includeGaps()
    skf.readTree()
    #Draw the ascii version of the tree
    skf.drawTree()
    #Get the family stats
    skf.familyStats()
    #Print the cost matrix used in the analysis
    skf.printCostMat()
    #Print the anecstor for each children pair
    skf.printAncestors()
    #Print the sequence for leaf nodes
    skf.printLeafSeq()
    skf.getAncestorSeq()
    #Print the sequence for ancestor nodes
    skf.printAncsSeq()
    skf.getSecondaryStruct()
    #Print the secondary structure
    skf.printSecondaryStruct()
