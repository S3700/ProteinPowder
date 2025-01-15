import random
from code.classes.graph import Graph

class RandomSolution:
    def __init__(self, protein_sequence, dimension=3, num_random_folds=1):
        """
        Initialize the RandomSolution algorithm for finding the best folding.
        """
        self.protein_sequence = protein_sequence
        self.dimension = dimension
        self.num_random_folds = num_random_folds
        self.graph = Graph(protein_sequence, dimension)  # Create a Graph object with the protein sequence

    def generate_random_folding(self):
        """
        Generate a random folding for the protein sequence.
        """
        # Possible folding steps: move along X, Y, Z axes, or negative directions
        folding = [random.choice([1, -1, 2, -2, 3, -3]) for _ in range(len(self.protein_sequence) - 1)]
        
        return folding


    def find_best_solution(self):
        """
        Perform random search to find the best folding for the protein sequence.
        """
        best_folding = None
        best_score = float('inf')
    
        # Try a number of random foldings and select the best one based on score
        for _ in range(self.num_random_folds):
            folding = self.generate_random_folding()  # Generate a random folding

            # Apply the folding and check if it's valid (no crossings)
            if not self.graph.apply_folding(folding):  # Skip if the folding is invalid (crossing detected)
                continue

            # Calculate the score for the current folding
            score = self.graph.calculate_score()

            # If the score is better (lower), store the folding
            if score < best_score:
                best_score = score
                best_folding = folding

        return best_folding, best_score
