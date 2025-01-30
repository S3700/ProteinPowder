from collections import deque
from .depthfirst import DepthFirst

class BreadthFirst(DepthFirst):
    def create_search_states(self, base_folding, base_score):
        """
        Creates the initial search states container for BreadthFirst.
        """
        return deque([(base_folding, 0, base_score)])
    
    def add_search_state(self, states, new_state):
        """
        Adds a new state to the search states container.
        For BreadthFirst, appends to the right (queue behaviour).
        """
        states.append(new_state)
        
    def get_next_state(self, states):
        """
        Gets the next state from the search states container.
        For BreadthFirst, pops from left (queue behaviour).
        """
        return states.popleft()
    
    