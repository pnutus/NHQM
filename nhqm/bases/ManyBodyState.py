from itertools import combinations

class ManyBodyState:
    """Represents a many body state"""
    minus = False
    states = []
    
    def __init__(self, states_in_list):
        self.state_list = []
    
    def create(self, new_states):
        for new_state in new_states:
            for i, state in enumerate(states):
                if new_state < state:
                    states.insert()
            
            self.state_list.append(state)
            
        self.state_list.sort()
        neg = 0
        for state in states:
            neg += len(self.state_list) - self.state_list.index(state)
        minus = minus  bool(neg % 2)    
        #wtf?
        minus = (minus and not flip) and (not minus and flip)
    