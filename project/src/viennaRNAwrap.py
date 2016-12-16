'''
@uthor: Haji Mohammad Saleem
Date  : December 9, 2016

Objective 5
Python Wrapper for RNAfold and RNADistance
Adapted from Vienna Wrappers for Python by Yann Ponty
'''
import os
import sys
sys.dont_write_bytecode = True

RNAFOLDwrapper = "RNAfold"
RNADISTwrapper = "RNAdistance"

def runRNAFold(seq):
  (tmpFileOut,tmpFileIn,struct,energy) = ("tmpOut.dat","tmpIn.dat",None,None)
  inFile = open(tmpFileIn,"w")
  inFile.write(seq+"\n")
  inFile.close()
  os.system("%s > %s < %s"%(RNAFOLDwrapper,tmpFileOut,tmpFileIn))
  lineno = 0
  for l in open(tmpFileOut,"r"):
    if lineno==1:
      data = l[:-1].split()
      struct = data[0]
      energy = float(" ".join(data[1:])[1:-1])
    lineno += 1
  os.remove("rna.ps")
  os.remove(tmpFileOut)
  os.remove(tmpFileIn)
  return (struct,energy)

def runRNADistance(seq1, seq2):
  (tmpFileOut,tmpFileIn,dist) = ("tmpOut.dat","tmpIn.dat",None)
  inFile = open(tmpFileIn,"w")
  inFile.write(seq1+"\n")
  inFile.write(seq2+"\n")
  inFile.close()
  os.system("%s > %s < %s"%(RNADISTwrapper,tmpFileOut,tmpFileIn))
  with open(tmpFileOut,"r") as f_temp:
    l = f_temp.readline()
    dist = l[3:]
  os.remove(tmpFileOut)
  os.remove(tmpFileIn)
  return dist.strip()
