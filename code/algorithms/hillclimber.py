from code.algorithms.random import RandomSolution

class HillClimber(RandomSolution):
    def __init__(self, protein_sequence, max_iterations=10000, num_valid_folds=1):
        """
        Initialize the HillClimber algorithm for finding the best folding, inheriting from RandomSolution.
        """
        super().__init__(protein_sequence, num_valid_folds)
        self.max_iterations = max_iterations

    def generate_neighbors(self, current_folding):
        """
        Generate valid neighboring states by making small modifications to the current folding.
        Returns a list of (folding, score) tuples for valid neighbors.
        """
        neighbors = []
        possible_directions = [1, -1, 2, -2, 3, -3]
        
        for i in range(len(current_folding)):
            for direction in possible_directions:
                if direction != current_folding[i]:
                    neighbor = current_folding.copy()
                    neighbor[i] = direction
                    
                    if self.graph.apply_folding(neighbor):
                        score = self.graph.calculate_score()
                        neighbors.append((neighbor, score))
                    
                    self.graph.apply_folding(current_folding)
        return neighbors

    def find_solutions(self):
        """
        Perform Hill Climbing to find the best folding for the protein sequence.
        """
        valid_attempts = 0
        best_score = float('inf')
        best_folding = None
        self.all_scores = []

        while valid_attempts < self.num_valid_folds:
            # Get a valid initial state
            current_folding, current_score = self.get_valid_folding()
            if current_folding is None:
                continue
   
            iterations = 0
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

        return best_folding, best_score, self.all_scores