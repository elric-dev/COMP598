'''
@uthor: Haji Mohammad Saleem
Date  : December 9, 2016

Objective 1
A parser to extract the sequence alignment and the consensus secondary structure from a stockholm file.
'''

import sys
sys.dont_write_bytecode = True
import os
import paths

DATA_PATH = paths.data_path
SEED_PATH = os.path.join(DATA_PATH, 'seed')

class Stockholm:

    def __init__(self, family):
        self.family = family
        self.SS_cons  = None
        self.SQ_align = {}
        return

    def parse(self):
        sq_align = {}
        filepath = os.path.join(SEED_PATH, self.family+'.stockholm')
        with open(filepath, 'r') as fin:
            all_lines = fin.readlines()
        all_lines = [x.strip() for x in all_lines]
        all_lines = [x for x in all_lines if x]
        for line in all_lines:
            if line.startswith("#=GC SS_cons"):
                ss_cons = line.split()[-1]
            if not line.startswith("#") and not line.startswith("/"):
                key = line.split()[0]
                seq = line.split()[1]
                sq_align[key] = seq
        self.SS_cons = ss_cons
        self.SQ_align = sq_align
        return    

if __name__ == "__main__":
    
    rfamFamily = 'RF02375'

    stk = Stockholm(rfamFamily)
    stk.parse()
    
    print "Consensus Secondary Structure"
    print "{:40} :: {}".format('', stk.SS_cons)
    
    print "Sequence Alignment"
    for key in stk.SQ_align.keys():
        print "{:40} :: {}".format(key, stk.SQ_align[key])
