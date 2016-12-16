import vienna_wrap as vw
import random
import math 
import numpy as np
import matplotlib.pyplot as plt

#define a class CELL for every element in the population
class Cell:
    def __init__(self, seq, struc):
        self.sequence = seq
        self.structure = struc
        self.fitness = None
        self.id = None
        self.parent = None

#setup a random population
def populate(target, pop_size):
    population = []
    seq_len = len(target)
    #pop_size = N

    for i in range(pop_size):
        sequence = "".join([random.choice("AUCG") for idx in range(seq_len)])
        structure = vw.runRNAFold(sequence)[0]
        new_cell = Cell(sequence, structure)
        new_cell.id = i
        new_cell.parent = i
        population.append(new_cell)
            
    return population

def compute_fitness(population, target):
    beta = -2./len(target)
    tot = []
    
    for cell in population:
        d = int(vw.runRNADistance(cell.structure, target).strip())
        fitness = math.exp(beta * d)
        cell.fitness = fitness
        cell.bp_distance = d
        tot.append(cell.fitness)

    norm = sum(tot)
    for cell in population:
        cell.fitness = float(cell.fitness) / norm

    return None

def mutate(sequence, mutation_rate):
    new_sequence = ""
    mutated = False
    for bp in sequence:
        r = random.random()
        if r < mutation_rate:
            new_sequence = new_sequence + random.choice("AUCG")
            mutated = True
        else:
            new_sequence = new_sequence + bp
            
    return (new_sequence, mutated)

def selection(population, target, mutation_rate):
    beta = -2./len(target)

    parents = np.random.choice(population, len(population), p=[rna.fitness for rna in population], replace=True)

    next_generation = []    
    for i, p in enumerate(parents):
        new_cell = Cell(p.sequence, p.structure)
        new_cell.id = i
        new_cell.parent = p.id
        
        next_generation.append(new_cell)

    for rna in next_generation:      
        mutated_sequence, mutated = mutate(rna.sequence, mutation_rate)
        
        if mutated:
            rna.sequence = mutated_sequence
            rna.structure = vw.runRNAFold(mutated_sequence)[0]
        else:
            continue

    compute_fitness(next_generation, target)

    return next_generation

def evolve(target, generations, pop_size, mutation_rate):
    beta = -2./len(target)
    
    populations = []
    population_stats = []
    initial_population = populate(target, pop_size=pop_size)
    compute_fitness(initial_population, target)

    current_generation = initial_population

    for i in range(generations):
        basepair_dis = []
        for cell in current_generation:
            basepair_dis.append(cell.bp_distance)
        population_stats.append(np.mean(basepair_dis))

        populations.append(current_generation)
        new_gen = selection(current_generation, target, mutation_rate=mutation_rate)
        current_generation = new_gen 
    
    return (populations, population_stats)
#------------------------------------------------------------------
plt.style.use('seaborn-paper')
#------------------------------------------------------------------
targets = [ '((((((((....))))))))', '((((..(((....)))))))', '(((....)))(((....)))']
mut_rate = [ 0.01, 0.02, 0.05, 0.1 ]
N = 100
generations = 500

def plotit(target):
    print target
    plt.figure()
    colors = ['crimson', 'teal', 'steelblue', 'purple']
    i = 0
    for mu in mut_rate:
    	print mu
    	pops, pops_stats = evolve(target, generations=generations, pop_size=N, mutation_rate=mu)
        plt.scatter( range(0,generations), pops_stats, color = colors[i], label = 'mu: '+str(mu))
    	i+=1

    filename = 'Target%s.png'%target
    plt.legend(loc='upper right')
    plt.xlabel('Generations')
    plt.ylabel('Average Base Pair Distance')
    plt.title( 'Target : %s'%target)
    plt.savefig(filename, dpi=700)
    return

for target in targets:
    plotit(target)
