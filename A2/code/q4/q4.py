import operator
import matplotlib.pyplot as plt
import networkx as nx
import itertools
import numpy as np
from scipy import linalg as LA
import sys
import os
import itertools
from pprint import pprint
import pickle
#---------------------------------------------------------------------------------
def read_edge_list(filename):
    G=nx.Graph()
    with open(filename, 'r') as fin:
        all_edges = fin.readlines()
    edge_list = [tuple(x.split()[:2]) for x in all_edges]
    G.add_edges_from(edge_list)
    all_nodes = sorted(nx.nodes(G))

    #print the Graph
    #nx.draw_spring(G, with_labels = True, node_color='b', alpha=0.5, node_size=100)
    #plt.show()
    return G, all_nodes 

def get_index(lst, item):
    return lst.index(item)

def get_neighbours(network, node):
    return network.neighbors(node)

def compare_networks(file1, file2):
    N1, nodes_N1 = read_edge_list(file1)
    N2, nodes_N2 = read_edge_list(file2)
    #print 'N1: %s nodes'%len(nodes_N1)
    #print 'N2: %s nodes'%len(nodes_N2)
    A = np.zeros((len(nodes_N1)*len(nodes_N2), len(nodes_N1)*len(nodes_N2)), dtype=np.float64)

    node_combos = [nodes_N1, nodes_N2]
    node_combos = list(itertools.product(*node_combos))
    
    for ij in xrange(len(node_combos)):
        node_i = node_combos[ij][0]
        node_j = node_combos[ij][1]
        neighbours_i = get_neighbours(N1, node_i)
        neighbours_j = get_neighbours(N2, node_j)
        neighbours_combos = [neighbours_i, neighbours_j]
        neighbours_combos = list(itertools.product(*neighbours_combos))
        for uv in xrange(len(neighbours_combos)):
                node_u = neighbours_combos[uv][0]
                node_v = neighbours_combos[uv][1]
                a_uv = 1.0/(len(get_neighbours(N1, node_u))*len(get_neighbours(N2, node_v))) 
                A[ij, get_index(node_combos, neighbours_combos[uv])] =  a_uv       

    #print 'A generated'
    return A, node_combos

def remove_nodes(node_list, node_pair, eigen):
    node_i = node_pair[0]
    node_j = node_pair[1]
    idx = []
    for i in xrange(len(node_list)):
        if node_list[i][0] == node_i or node_list[i][1] == node_j:
            idx.append(i)
    for i in sorted(idx, reverse=True):
            del node_list[i]
            del eigen[i]
    return node_list, eigen

def find_max(node_list, eigen):
    index, value = max(enumerate(eigen), key=operator.itemgetter(1))
    return value, node_list[index]

def allign_networks(node_list, eigen):
    assert len(node_list) == len(eigen)
    thold = 0.01
    while max(eigen) > thold:
        max_eigen, node_pair = find_max(node_list, eigen)
        print '{:20} {}'.format(node_pair[0], node_pair[1])
        node_list, eigen = remove_nodes(node_list, node_pair, eigen)
        if not node_list:
            break
    return

#--------------------------------------


if __name__ == "__main__":
    network_file1 = sys.argv[1]
    network_file2 = sys.argv[2]

    A, node_combos = compare_networks(network_file1, network_file2)
    
    pklfile = "%s_%s.pkl"%(network_file1.split('.')[0], network_file2.split('.')[0])
    if not os.path.isfile(pklfile):
        e_vals, e_vecs = LA.eig(A)
        eigen = e_vecs[0].real.tolist()
        eigen = [abs(x) for x in eigen]
        with open(pklfile, 'wb') as fout:
            pickle.dump(eigen, fout)
    with open(pklfile, 'rb') as fin:
        eigen = pickle.load(fin)
        #print 'Eigen Vectors Loaded'

    allign_networks(node_combos, eigen)

    #for i in xrange(A.shape[0]):
    #    for j in xrange(A.shape[1]):
    #        if A[i, j] != 0.0:
    #            print i, j, A[i, j]

