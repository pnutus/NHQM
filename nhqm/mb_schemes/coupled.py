from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import product, combinations_with_replacement
from collections import namedtuple, Iterable
from nhqm.QM_helpers import clebsch_gordan, matrix_from_function
 
name = "Coupled"

def hamiltonian(Q, eigvals, eigvecs, interaction, num_particles=2):
    
    # tuples of sp states (sp(E=1, j=1.5), sp(E=7, j=1.5))
    mb_states = gen_coupled_states(Q, num_particles)
    order = len(mb_states)
    
    def H_func(i, j):
        return coupled_H_elem(mb_states[i], mb_states[j], 
                              eigvals, interaction, Q)
    return matrix_from_function(H_func, order)

def coupled_H_elem(bra, ket, eigvals, interaction, Q):
    # H_1, one-body interaction
    one_body = 0
    a, b = ket
    c, d = bra
    if bra == ket:
        one_body = sum(eigvals[sp.E] for sp in bra) * (1 + (a == b)*(-1)**Q.J)
    
    # H_2, two-body interaction
    two_body = (interaction(a.E, b.E, c.E, d.E) 
              + interaction(a.E, b.E, d.E, c.E) * (-1)**Q.J )
    
    return (one_body  + two_body) * N(a, b, Q) * N(c, d, Q)

def N(a, b, Q):
    return sp.sqrt(1 + (a == b)*(-1)**Q.J)/(1 + (a == b))


def gen_coupled_states(Q, num_particles=2):
    q_nums = [(name, value) for name, value in zip(Q._fields, Q) 
                            if isinstance(value, Iterable) and name is not "m"]
    names, values = zip(*q_nums)
    SP = namedtuple('sp', names)
    sp_states = [SP(*tup) for tup in product(*values)]
    return list(combinations_with_replacement(sp_states, num_particles))
    
if __name__ == '__main__':
    # example
    basis_size = 5
    QNums = namedtuple('qnums', 'l J M E j m')
    Q = QNums(l=1, j=[0.5, 1.5], J=2, M=0, 
              m=[-1.5, -0.5, 0.5, 1.5], 
              E=range(basis_size))
    
    mb_states = gen_coupled_states(Q, 2)
    for mb_state in mb_states:
        print mb_state
    print "Number of many-body states:", len(mb_states)
