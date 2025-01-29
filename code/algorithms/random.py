import random
from code.classes.graph import Graph

class RandomSolution:
    def __init__(self, protein_sequence, num_valid_folds=1):
        """
        Initialize the RandomSolution algorithm for finding the best folding.
        """
        self.protein_sequence = protein_sequence
        self.num_valid_folds = num_valid_folds
        self.graph = Graph(protein_sequence)
        self.all_scores = []

    def generate_random_folding(self):
        """
        Generate a random folding for the protein sequence in 3d space. 
        """
        return [random.choice([1, -1, 2, -2, 3, -3]) for _ in range(len(self.protein_sequence) - 1)]

    def get_valid_folding(self):
        """
        Generate a valid random folding and its score.
        """
        folding = self.generate_random_folding()
        if self.graph.apply_folding(folding):
            score = self.graph.calculate_score()
            return folding, score
        return None, None

    def find_solutions(self):
        """
        Perform random search to find foldings for the protein sequence.
        Returns the best folding, its score, and all valid scores encountered.
        """
        best_folding = None
        best_score = float('inf')
        valid_attempts = 0
        
        while valid_attempts < self.num_valid_folds:
            folding, score = self.get_valid_folding()
            if folding is not None:
                self.all_scores.append(score)
                valid_attempts += 1
                
                if score < best_score:
                    best_score = score
                    best_folding = folding

        return best_folding, best_score, self.all_scores
