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
        for state in reversed(states):
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

def two_body_indexes(bra, ket):
    result = []
    if len(set(bra) - set(ket)) > 2:
        return []
    for annihilated in combinations(ket, 2):
        for created in combinations(bra, 2):
            new_ket = (ket
                        .annihilate(annihilated)
                        .create(created))
            if new_ket.states == bra.states:
                sign = new_ket.sign
                a,b = created
                c,d = annihilated
                result.append( (a,b,c,d,sign) )
    return result


"""
TESTS
"""
    
import unittest
class RedTests(unittest.TestCase):
    
    def setUp(self):
        self.bra32 = FermionState([1,2,3]) 
        self.ket32 = FermionState([1,3,4])
        self.res32 = [(1, 2, 1, 4, -1), (2, 3, 3, 4, 1)]
        
        self.bra31 = FermionState([1,2,3]) 
        self.ket31 = FermionState([3,5,6])
        self.res31 = [(1,2, 5,6, +1)]

        self.bra33 = FermionState([1,2,3]) 
        self.ket33 = FermionState([1,2,3])
        self.res33 = [(1,2, 1,2, +1), (1,3,1,3,1), (2,3,2,3,1)]
        
        
    def test32(self):
        res = two_body_indexes(self.bra32, self.ket32)
        
        self.assertEquals(res, self.res32 )
        
    def test31(self):
        res = two_body_indexes(self.bra31, self.ket31)
        
        self.assertEquals(res, self.res31 ) 
        
    def test33(self):
        res = two_body_indexes(self.bra33, self.ket33)
        
        self.assertEquals(res, self.res33 ) 