from itertools import combinations

class FermionState:
    """
    Represents a fermionic many-body state as a list of
    tuples, each tuple containing the quantum numbers for a
    single particle. For example: (k-index, j, m) = (3, 1.5, -.5)
    """
    def __init__(self, states = []):
        self.states = []
        self.sign = 1
        for state in states:
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
            
    def _operator(self, operator, states):
        new_fermion = self.copy()
        for state in states:
            new_fermion.operator(state)
        return new_fermion

    
    def create(self, states):
        """Acts with creation operator on the fermionic state."""
        return self._operator(_create, states)
    
    def annihilate(self, states):
        """Acts with annihilation operator on the fermionic state."""
        return self._operator(_annihilate, states)
            
    def copy(self):
        new_fermion = FermionState()
        new_fermion.states = self.states[:]
        new_fermion.sign = self.sign
        return new_fermion
        
    def __iter__(self):
        return iter(self.states)
        
    def __str__(self):
        if self.sign == 0:
            return "0"
        sign = "-" if self.sign < 0 else "+" 
        return sign + "| {} >".format(", ".join(str(x) for x in self.states))
    
    def __repr__(self):
        return str(self)

    def __eq__(self,other):
        return self.sign == other.sign and self.states == other.states
        
    def __len__(self):
        return len(self.states)
    
def state_index(states, new_state):
    for i, state in enumerate(states):
        if new_state == state:
            raise ValueError
        if new_state < state:
            return i
    return len(states)

# TESTS?

print FermionState([(1,-.5), (1, .5), (7, 1.5), (0, -1.5)])