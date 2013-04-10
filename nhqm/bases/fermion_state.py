from itertools import combinations

class FermionState:
    """Represents a fermionic many-body state."""    
    def __init__(self, states = [], sign = 1):
        self.states = states
        self.sign = sign
    
    def create(self, new_states):
        """Acts with creation operator on the fermionic state."""
        result_states = self.states
        result_sign = self.sign
        if isinstance(new_states, int):
            new_states = [new_states]
        for new_state in new_states:
            try:
                i = state_index(result_states, new_state)
                result_states.insert(i, new_state)
                result_sign *= (-1)**i            # 
            except IndexError:
                return FermionState(sign = 0)
        return FermionState(result_states, result_sign)
    
    def annihilate(self, kill_states):
        """Acts with annihilation operator on the fermionic state."""
        result_states = self.states
        result_sign = self.sign
        if isinstance(kill_states, int):
            kill_states = [kill_states]
        for kill_state in kill_states:
            try:
                i = result_states.index(kill_state)
                result_states.pop(i)
                result_sign *= (-1)**i
            except ValueError:
                return FermionState(sign = 0)
        return FermionState(result_states, result_sign)
    
    def __iter__(self):
        return iter(self.states)
        
    def __str__(self):
        if self.sign == 0:
            return "0"
        sign = "-" if self.sign < 0 else "" 
        return sign + "| {} >".format(", ".join(str(x) for x in self.states))
    
    def __repr__(self):
        return str(self)
    
def state_index(states, new_state):
    for i, state in enumerate(states):
        if new_state == state:
            raise IndexError
        if new_state < state:
            return i
    return len(states)

# TESTS?