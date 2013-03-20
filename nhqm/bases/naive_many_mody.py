from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations, permutations

# THIS IS FOR FERMIONS

def gen_states(num_sp_states, num_particles=1):
    if num_sp_states < num_particles:
        raise ValueError(
        "There cannot be more particles than states"
        )
    return map(set, combinations(range(num_sp_states), num_particles))

def get_naive_combinations(num_states, num_particles=1):

    states=gen_states(num_states, num_particles)
    states2=gen_states(num_states, num_particles)
    single_state=[]
    for i in xrange(num_states):
        single_state.append( set([i]) )
        
    combinations = []        

    for i, a in enumerate(states):
        for j, b in enumerate(states2):
                #testa alfa och beta och se om det blir nagot
                for beta in b:
                    for alpha in single_state:
                        if (a == (b - set([beta]) | alpha )):
                            combinations.append( (a, b) )
                        
    return combinations                        
            
