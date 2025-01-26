import random
import numpy as np
from code.classes.graph import Graph

class HillClimberv4:
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
        """
        folding = [random.choice([1, -1, 2, -2, 3, -3]) for _ in range(len(self.protein_sequence) - 1)]
        
        # Check if the folding is valid
        if self.graph.apply_folding(folding):
            return folding
        return None

    def calculate_centralization_bonus(self, amino_acids):
        """
        Calculate a bonus for centralizing C amino acids.
        """
        # Find all C amino acids
        c_aminos = [aa for aa in amino_acids if aa.aa_type == 'C']
        
        if not c_aminos:
            return 0
        
        # Calculate the centroid of all amino acids
        all_coords = np.array([aa.position for aa in amino_acids])
        centroid = np.mean(all_coords, axis=0)
        
        # Calculate distances of C amino acids from the centroid
        c_coords = np.array([aa.position for aa in c_aminos])
        distances = np.linalg.norm(c_coords - centroid, axis=1)
        
        # Reward shorter distances (lower is better)
        centralization_bonus = -10 * (1 / (np.mean(distances) + 1))
        
        return centralization_bonus

    def calculate_temp_score(self):
        """
        Calculate the temporary score for a folding including the centralization bonus.
        """
        # Get the base score (original score)
        base_score = self.graph.calculate_score()
        
        # Add the centralization bonus
        centralization_bonus = self.calculate_centralization_bonus(self.graph.amino_acids)
        
        # Return the combined temporary score
        temp_score = base_score + centralization_bonus
        return temp_score

    def generate_neighbors(self, current_folding):
        """
        Generate valid neighboring states by making small modifications to the current folding.
        Returns a list of (folding, temporary_score) tuples for valid neighbors.
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
                        temporary_score = self.calculate_temp_score()
                        neighbors.append((neighbor, temporary_score))
                    
                    # Reapply the current folding to reset the graph state
                    self.graph.apply_folding(current_folding)
        return neighbors

    def find_best_solution(self):
        """
        Perform Hill Climbing to find the best folding for the protein sequence.
        """
        valid_attempts = 0
        
        while valid_attempts < self.num_valid_folds:
            # Generate initial folding
            current_folding = None
            while current_folding is None:
                current_folding = self.generate_random_folding()
            
            # Initial temporary score
            current_temporary_score = self.calculate_temp_score()
            iterations = 0

            # Perform hill climbing for this fold
            while iterations < self.max_iterations:
                # Generate neighbors
                neighbors = self.generate_neighbors(current_folding)
                
                if not neighbors:
                    break

                # Find best neighbor
                best_neighbor = min(neighbors, key=lambda x: x[1])
                neighbor_folding, neighbor_temporary_score = best_neighbor

                # If no improvement, we've reached local optimum
                if neighbor_temporary_score >= current_temporary_score:
                    break

                # Move to the better neighbor
                current_folding = neighbor_folding
                current_temporary_score = neighbor_temporary_score
                iterations += 1

            # Calculate the final true score after optimization
            final_score = self.graph.calculate_score()
            self.all_scores.append(final_score)
            
            # Update best solution if needed
            if final_score < self.best_score:
                self.best_score = final_score
                self.best_folding = current_folding.copy()

            valid_attempts += 1

        # Apply the best folding to the graph
        if self.best_folding:
            self.graph.apply_folding(self.best_folding)
        
        return self.best_folding, self.best_score, self.all_scores
    
    def __repr__(self):
        return f"HillClimber(sequence={self.protein_sequence}, dimension={self.dimension})"

