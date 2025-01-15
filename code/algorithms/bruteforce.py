import itertools
from classes.graph import Graph

class Bruteforce:
    def __init__(self, protein_sequence, dimension=3):
        """
        Initialize the Bruteforce algorithm for finding the best folding.
        """
        self.protein_sequence = protein_sequence
        self.dimension = dimension
        self.graph = Graph(protein_sequence, dimension)

    def find_best_solution(self):
        """
        Try all possible foldings and return the one with the best score.
        """
        # Step options based on the dimension (3D or 2D)
        steps = [1, -1, 2, -2, 3, -3] if self.dimension == 3 else [1, -1, 2, -2]
        
        best_score = float('inf')
        best_folding = None

        # Iterate through all possible folding combinations
        for folding in itertools.product(steps, repeat=len(self.protein_sequence) - 1):
            # Update the graph with the current folding configuration
            if not self.graph.apply_folding(folding):  # Skip if crossing detected
                continue
            
            # Calculate the score for the current folding
            score = self.graph.calculate_score()

            # If the score is better (lower), store folding
            if score < best_score:
                best_score = score
                best_folding = folding

        return best_folding, best_score

