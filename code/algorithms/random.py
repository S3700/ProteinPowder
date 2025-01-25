import random
from code.classes.graph import Graph

class RandomSolution:
    def __init__(self, protein_sequence, dimension=3, num_valid_folds=1):
        """
        Initialize the RandomSolution algorithm for finding the best folding.
        """
        self.protein_sequence = protein_sequence
        self.dimension = dimension
        self.num_valid_folds = num_valid_folds
        self.graph = Graph(protein_sequence, dimension)
        self.all_scores = []  # List to store all valid folding scores

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
        Returns the best folding, its score, and all valid scores encountered.
        """
        best_folding = None
        best_score = float('inf')
        valid_attempts = 0
        
        while valid_attempts < self.num_valid_folds:
            folding = self.generate_random_folding()

            # Apply the folding and check if it's valid (no crossings)
            if self.graph.apply_folding(folding):  # Only proceed if folding is valid
                score = self.graph.calculate_score()
                self.all_scores.append(score)  # Store the score
                valid_attempts += 1

                # If the score is better (lower), store the folding
                if score < best_score:
                    best_score = score
                    best_folding = folding

        return best_folding, best_score, self.all_scores
