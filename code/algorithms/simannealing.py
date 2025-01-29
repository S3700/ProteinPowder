import random
import math
from code.classes.graph import Graph

class SimulatedAnnealing:
    def __init__(self, protein_sequence, dimension=3, 
                 max_iterations=10000, 
                 initial_temperature=10, 
                 cooling_rate=0.0012, 
                 num_valid_folds=1):
        """
        Initialize the Simulated Annealing algorithm for finding the best protein folding.
        """
        self.protein_sequence = protein_sequence
        self.dimension = dimension
        self.max_iterations = max_iterations
        self.initial_temperature = initial_temperature
        self.cooling_rate = cooling_rate
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
            score = self.graph.calculate_score()
            return folding, score
        return None, None

    def get_initial_state(self):
        """
        Get a valid initial state for the algorithm.
        """
        folding, score = None, None
        while folding is None:
            folding, score = self.generate_random_folding()
        return list(folding), score

    def generate_neighbor(self, current_folding):
        """
        Generate a neighboring state by making a small modification to the current folding.
        """
        possible_directions = [1, -1, 2, -2, 3, -3]
        
        # Create a copy of the current folding
        neighbor = current_folding.copy()
        
        # Randomly choose a position to modify
        index = random.randint(0, len(neighbor) - 1)
        
        # Choose a new direction different from the current one
        current_direction = neighbor[index]
        possible_new_directions = [d for d in possible_directions if d != current_direction]
        neighbor[index] = random.choice(possible_new_directions)
        
        return neighbor

    def find_best_solution(self):
        valid_attempts = 0
        best_overall_score = float('inf')
        best_overall_folding = None
        self.all_scores = []

        while valid_attempts < self.num_valid_folds:
            # Get random initial state for each valid fold
            current_folding, current_score = self.get_initial_state()
            
            # Initialize temperature
            temperature = self.initial_temperature
            
            # Keep track of the best solution for this fold
            best_folding = current_folding.copy()
            best_score = current_score
            
            # Simulated Annealing iterations
            for _ in range(self.max_iterations):
                # Generate a neighbor
                neighbor = self.generate_neighbor(current_folding)
                
                # Check if neighbor is a valid folding
                if not self.graph.apply_folding(neighbor):
                    continue
                
                # Calculate neighbor's score
                neighbor_score = self.graph.calculate_score()
                
                # Calculate the energy difference
                delta_energy = neighbor_score - current_score
                
                # Decide whether to accept the neighbor
                if delta_energy < 0 or random.random() < math.exp(-delta_energy / temperature):
                    current_folding = neighbor
                    current_score = neighbor_score
                
                # Update best solution for this fold if needed
                if current_score < best_score:
                    best_folding = current_folding.copy()
                    best_score = current_score
                
                # Cool down the temperature
                temperature *= (1 - self.cooling_rate)
            
            # Store the best score for this fold
            self.all_scores.append(best_score)
            
            # Update best overall solution if needed
            if best_score < best_overall_score:
                best_overall_score = best_score
                best_overall_folding = best_folding.copy()
            
            valid_attempts += 1

        # Apply the best folding to the graph
        self.graph.apply_folding(best_overall_folding)
        
        return best_overall_folding, best_overall_score, self.all_scores

    def __repr__(self):
        return f"SimulatedAnnealing(sequence={self.protein_sequence}, dimension={self.dimension})"