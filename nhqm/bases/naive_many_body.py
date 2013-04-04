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

def get_single_particle_combinations(num_states, num_particles=1):

    states=gen_states(num_states, num_particles)
    states2=gen_states(num_states, num_particles)
    single_state=gen_states(num_states, num_particles=1)
        
    combs = []        

    for i, a in enumerate(states):
        for j, b in enumerate(states2):
                #testa alfa och beta och se om det blir nagot
                for beta in b:
                    for alpha in single_state:
                        if (a == (b - set([beta]) | alpha )):
                            combs.append( \
                            (a, b, alpha, set([beta])) )
                        
    return combs                       
            
def get_two_particle_combinations(num_states, num_particles=2):

    states=gen_states(num_states, num_particles)
    states2=gen_states(num_states, num_particles)
    double_state=gen_states(num_states, num_particles=2)
        
    combs = []        

    for i, a in enumerate(states):
        for j, b in enumerate(states2):
                #testa alfa1,2 och beta1,2 och se om det blir nagot
                gamma_delta_permut = combinations(b,2)
                for gamma_delta in gamma_delta_permut:
                    for alpha_beta in double_state:
                        if (a == (b - set(gamma_delta) | alpha_beta )):
                            combs.append( \
                            (a, b, alpha_beta, set(gamma_delta) ) )
                        
    return combs                        
            
