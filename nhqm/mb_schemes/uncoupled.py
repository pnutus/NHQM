from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import product, combinations
from collections import namedtuple, Iterable
from nhqm.mb_schemes.fermion_state import FermionState
from nhqm.QM_helpers import clebsch_gordan, matrix_from_function

name = "Uncoupled"

def hamiltonian(Q, eigvals, eigvecs, interaction, num_particles=2):
    # Generate many-body states
    mb_states = gen_uncoupled_states(Q, num_particles)
    
    # Only take states with m1 + m2 = M
    def total_M(mb_state):
        return Q.M == sum(state.m for state in mb_state)
    mb_states = filter(total_M, mb_states)
    
    # Generate hamiltonian
    matrix_size = len(mb_states)
    def H_func(i, j):
        return H_elem(mb_states[i], mb_states[j], eigvals, interaction)
    return matrix_from_function(H_func, matrix_size, symmetric=True)
            
def H_elem(bra, ket, eigvals, interaction):
    # H_1, one-body interaction
    one_body = 0
    if bra == ket:
        one_body = sum(eigvals[sp.E] for sp in bra)
    
    # H_2, two-body interaction
    a, b = bra
    c, d = ket
    two_body = 0
    if a.m == c.m and b.m == d.m:
        two_body += interaction(a.E, b.E, c.E, d.E)
    if a.m == d.m and b.m == c.m:
        two_body -= interaction(a.E, b.E, d.E, c.E)
        
    return one_body + two_body

def gen_uncoupled_states(Q, num_particles=2):
    q_nums = [(name, value) for name, value in zip(Q._fields, Q) 
                            if isinstance(value, Iterable)]
    names, values = zip(*q_nums)
    SP = namedtuple('sp', names)
    sp_states = [SP(*tup) for tup in product(*values)]
    return map(FermionState, combinations(sp_states, num_particles))
    
if __name__ == '__main__':
    # example
    basis_size = 5
    QNums = namedtuple('qnums', 'l J M E j m')
    Q = QNums(l=1, j=[0.5, 1.5], J=2, M=0, 
              m=[-1.5, -0.5, 0.5, 1.5], 
              E=range(basis_size))
    
           
    mb_states = gen_uncoupled_states(Q, 2)
    for mb_state in mb_states:
        print mb_state
    print "Number of many-body states:", len(mb_states)