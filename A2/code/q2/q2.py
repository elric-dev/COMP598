'''
@uthor: Haji Mohamamd Saleem
COMP 599, Assignment 2, Q2
'''

import numpy as np
import pandas as pd
from ggplot import *
import matplotlib.pyplot as plt
from sklearn import metrics
from operator import itemgetter
#------------------------------------------------------------------------------------------

#DATA
alpha_prop = {"A":0., "R":0.21, "N":0.65, "D":0.69, "C":0.68, "E":0.4,\
              "Q":0.39, "G":1., "H":0.61, "I":0.41, "L":0.21, "K":0.26,\
              "M":0.24, "F":0.54, "P":3.16, "S":0.5, "T":0.66, "W":0.49,\
              "Y":0.53, "V":0.61}

#Benchmark 1MBN chain A primary and secondary structure
primary_1MBN=   'VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG'
secondary_1MBN= '   HHHHHHHHHHHHHHTTSHHHHHHHHHHHHHHH HHHHHT HHHHT  SHHHHHH HHHHHHHHHHHHHHHHHHTTTT  HHHHHHHHHHHHHTT   HHHHHHHHHHHHHHHHHH TTTTSHHHHHHHHHHHHHHHHHHHHHHHHHT   '
#------------------------------------------------------------------------------------------

class AlphaPredictor:

    def __init__(self, alpha_prop, primary_seq, secondary_seq, thold):
        self.alpha_prop = alpha_prop
        self.thold = thold
        self.primary = primary_seq
        self.secondary = secondary_seq
        self.predicted = None

    def ah_predictor(self):
        predicted_structure = ''
        for residue in self.primary:
            allele_penalty = self.alpha_prop[residue]
            if allele_penalty < self.thold:
                predicted_structure+='H'
            else:
                predicted_structure+='C'
        self.predicted = predicted_structure
        return 
    
    def get_window(self, seq, index):
        start_index = index - 4
        if start_index < 0:
                start_index = 0
        end_index = index+4
        if end_index >= len(seq):
                end_index = len(seq)-1

        return seq[start_index:end_index+1]

    def ah_predictor_updated(self):
        predicted_structure = ''
	for index in xrange(len(self.primary)):
		seq_window = self.get_window(self.primary, index)
		window_score = sum([self.alpha_prop[residue] for residue in seq_window])
		window_avg = float(window_score)/len(seq_window)
		if window_avg < self.thold:
                	predicted_structure+='H'
            	else:
                	predicted_structure+='C'
	self.predicted = predicted_structure
        return

    def ah_metrics(self):
        TP = 0.
        FP = 0.
        TN = 0.
        FN = 0.
        for actual, predicted in zip(self.secondary, self.predicted):
            if actual == 'H' and predicted == 'H':
                TP+=1.
            if actual == 'H' and predicted != 'H':
                FN+=1.
            if actual != 'H' and predicted == 'H':
                FP+=1.
            if actual != 'H' and predicted != 'H':
                TN+=1.
        TPR = TP/(TP+FN)
        FPR = FP/(FP+TN)
        return TPR, FPR

def get_AUC(TPR_list, FPR_list):
    return np.trapz(TPR_list,FPR_list)

def plotit(tpr, fpr, predictor):
    plt.style.use('ggplot')
    plt.figure()
    plt.plot( fpr, tpr )
    x = np.linspace(0,1)
    plt.plot(x, x, 'b--')
    auc = get_AUC(tpr, fpr)
    print 'AUC %s'%auc
    plt.xlabel( 'False Positive Rate' )
    plt.ylabel( 'True Positive Rate' )
    plt.title("ROC curve with %s predictor, auc:%s"%(predictor, auc))
    plt.savefig(predictor+'.png')
    return


#------------------------------------------------------------------------------------------

#Replacing everything apart from H as C
secondary_1MBN_update = secondary_1MBN.replace('T', 'C').replace('S', 'C').replace(' ', 'C')
#------------------------------------------------------------------------------------------
#Normal predictor
TPR = []
FPR = []
for thold in np.arange(0, 4, 0.005):
    alpha_pred = AlphaPredictor(alpha_prop, primary_1MBN, secondary_1MBN_update, thold)
    alpha_pred.ah_predictor()
    tpr, fpr = alpha_pred.ah_metrics()
    TPR.append(tpr)
    FPR.append(fpr)

plotit(TPR, FPR, 'normal')

#------------------------------------------------------------------------------------------
#Updated predictor
TPR = []
FPR = []
for thold in np.arange(0, 4, 0.005):
    alpha_pred = AlphaPredictor(alpha_prop, primary_1MBN, secondary_1MBN_update, thold)
    alpha_pred.ah_predictor_updated()
    tpr, fpr = alpha_pred.ah_metrics()
    TPR.append(tpr)
    FPR.append(fpr)

plotit(TPR, FPR, 'updated')
