from itertools import combinations, product
from collections import namedtuple

class FermionState:
    """
    Represents a fermionic many-body state as a list of
    tuples, each tuple containing the quantum numbers for a
    single particle. For example: sp(k=3, j=1.5, m=-0.5)
    """
    def __init__(self, states = []):
        self.states = []
        self.sign = 1
        for state in reversed(states):
            self._create(state)
            
    def _create(self, new_state):
        try:
            i = state_index(self.states, new_state)
            self.states.insert(i, new_state)
            self.sign *= (-1)**i
        except ValueError: # The state already existed
            self.states = []
            self.sign = 0
    
    def _annihilate(self, kill_state):
        try:
            i = self.states.index(kill_state)
            self.states.pop(i)
            self.sign *= (-1)**i
        except ValueError: # The state didn't exist
            self.states = []
            self.sign = 0
    
    def create(self, states):
        """Acts with creation operator on the fermionic state."""
        new_fermion = self.copy()
        for state in states:
            new_fermion._create(state)
        return new_fermion
    
    def annihilate(self, states):
        """Acts with annihilation operator on the fermionic state."""
        new_fermion = self.copy()
        for state in states:
            new_fermion._annihilate(state)
        return new_fermion
            
    def copy(self):
        new_fermion = FermionState()
        new_fermion.states = self.states[:]
        new_fermion.sign = self.sign
        return new_fermion
        
    def __iter__(self):
        return iter(self.states)
        
    def __contains__(self, value):
        return value in self.states
        
    def __str__(self):
        if self.sign == 0:
            return "0"
        sign = "-" if self.sign < 0 else "+" 
        return sign + "| {} >".format(", ".join(str(x) for x in self.states))
    
    def __repr__(self):
        return str(self)

    def __eq__(self,other):
        return self.states == other.states
        
    def __len__(self):
        return len(self.states)
        
    def __getitem__(self, key):
        return self.states[key]
    
def state_index(states, new_state):
    for i, state in enumerate(states):
        if state == new_state:
            raise ValueError
        if state > new_state:
            return i
    return len(states)

def gen_mb_states(quantum_numbers, num_particles=2):
    """
    Generates many-body states given a list of 
    quantum numbers and their possible values.
    Like so [('m', [-1.5, -0.5, 0.5, 1.5]),
             ('k', [0, 1, 2, 3, 4, 5, 6])]
    """
    names, values = zip(*quantum_numbers)
    SP = namedtuple('sp', names)
    sp_states = [SP(*tup) for tup in product(*values)]
    return map(FermionState, combinations(sp_states, num_particles))

if __name__ == '__main__':
    # example
    order = 10
    q_nums = [('m', [-1.5, -0.5, 0.5, 1.5]),
              ('k', range(order))]
    mb_states = gen_mb_states(q_nums, 2)
    for mb_state in mb_states:
        print mb_state
    print "Number of many-body states:", len(mb_states)


# TESTS?