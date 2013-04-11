from itertools import combinations

class FermionState:
    """Represents a fermionic many-body state."""    
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
        except ValueError:
            self.states = []
            self.sign = 0
    
    def _annihilate(self, kill_state):
        try:
            i = self.states.index(kill_state)
            self.states.pop(i)
            self.sign *= (-1)**i
        except ValueError:
            self.states = []
            self.sign = 0
            
    def copy(self):
        new_fermion = FermionState()
        new_fermion.states = self.states[:]
        new_fermion.sign = self.sign
        return new_fermion
    
    def create(self, new_states):
        """Acts with creation operator on the fermionic state."""
        if isinstance(new_states, int):
            new_states = [new_states]
        new_fermion = self.copy()
        for state in new_states:
            new_fermion._create(state)
        return new_fermion
    
    def annihilate(self, kill_states):
        """Acts with annihilation operator on the fermionic state."""
        if isinstance(kill_states, int):
            kill_states = [kill_states]
        new_fermion = self.copy()
        for state in kill_states:
            new_fermion._annihilate(state)
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