import random
from code.classes.graph import Graph

class HillClimber:
    def __init__(self, protein_sequence, dimension=3, max_iterations=10000, num_valid_folds=1):
        """
        Initialize the HillClimber algorithm for finding the best folding.
        """
        self.protein_sequence = protein_sequence
        self.dimension = dimension
        self.max_iterations = max_iterations
        self.num_valid_folds = num_valid_folds
        self.graph = Graph(protein_sequence, dimension)
        self.best_score = float('inf')
        self.best_folding = None
        self.all_scores = []

    def generate_random_folding(self):
        """
        Generate a random initial folding for the protein sequence.
        Returns both the folding and its score if valid, otherwise returns None, None.
        """
        folding = [random.choice([1, -1, 2, -2, 3, -3]) for _ in range(len(self.protein_sequence) - 1)]
        
        # Check if the folding is valid
        if self.graph.apply_folding(folding):
            score = self.graph.calculate_score()
            return folding, score
        return None, None

    def get_initial_state(self):
        """
        Get a valid initial state for the algorithm.
        Keeps trying until a valid folding is found.
        """
        folding, score = None, None
        while folding is None:
            folding, score = self.generate_random_folding()
        return list(folding), score

    def generate_neighbors(self, current_folding):
        """
        Generate valid neighboring states by making small modifications to the current folding.
        Returns a list of (folding, score) tuples for valid neighbors.
        """
        neighbors = []
        possible_directions = [1, -1, 2, -2, 3, -3]
        
        # Try changing each position in the folding
        for i in range(len(current_folding)):
            for direction in possible_directions:
                if direction != current_folding[i]:  # Only try different directions
                    neighbor = current_folding.copy()
                    neighbor[i] = direction
                    
                    # Check if the neighbor folding is valid
                    if self.graph.apply_folding(neighbor):
                        score = self.graph.calculate_score()
                        neighbors.append((neighbor, score))
                    
                    # Reapply the current folding to reset the graph state
                    self.graph.apply_folding(current_folding)
        return neighbors

    def find_best_solution(self):
        """
        Perform Hill Climbing to find the best folding for the protein sequence.
        """
        valid_attempts = 0
        best_score = float('inf')
        best_folding = None
        self.all_scores = []  # Reset scores list

        while valid_attempts < self.num_valid_folds:
            # Get random initial state for each valid fold
            current_folding, current_score = self.get_initial_state()
            initial_score = current_score  # Store initial score
            iterations = 0

            # Perform hill climbing for this fold
            while iterations < self.max_iterations:
                neighbors = self.generate_neighbors(current_folding)
                
                if not neighbors:
                    break

                best_neighbor = min(neighbors, key=lambda x: x[1])
                neighbor_folding, neighbor_score = best_neighbor

                if neighbor_score >= current_score:
                    break

                current_folding = neighbor_folding
                current_score = neighbor_score
                iterations += 1

            # After hill climbing is complete for this fold
            self.all_scores.append(current_score)  # Store only the final score for this fold
            
            # Update best overall solution if current is better
            if current_score < best_score:
                best_score = current_score
                best_folding = current_folding.copy()

            valid_attempts += 1

        # After all folds are complete, apply the best folding to the graph
        self.graph.apply_folding(best_folding)
        
        return best_folding, best_score, self.all_scores
    
    def __repr__(self):
        return f"HillClimber(sequence={self.protein_sequence}, dimension={self.dimension})"