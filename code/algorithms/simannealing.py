from code.algorithms.hillclimber import HillClimber
import random
import math

class SimulatedAnnealing(HillClimber):
    def __init__(self, protein_sequence, 
                 max_iterations=10000, 
                 initial_temperature=10, 
                 cooling_rate=0.0012, 
                 num_valid_folds=1):
        """
        Initialize the Simulated Annealing algorithm, inheriting from HillClimber.
        """
        super().__init__(protein_sequence, max_iterations, num_valid_folds)
        self.initial_temperature = initial_temperature
        self.cooling_rate = cooling_rate

    def generate_single_neighbor(self, current_folding):
        """
        Generate a single random neighbor state.
        """
        possible_directions = [1, -1, 2, -2, 3, -3]
        neighbor = current_folding.copy()
        index = random.randint(0, len(neighbor) - 1)
        current_direction = neighbor[index]
        possible_new_directions = [d for d in possible_directions if d != current_direction]
        neighbor[index] = random.choice(possible_new_directions)
        
        return neighbor

    def find_solutions(self):
        """
        Perform Simulated Annealing to find the best folding. 
        """
        valid_attempts = 0
        best_overall_score = float('inf')
        best_overall_folding = None
        self.all_scores = []

        while valid_attempts < self.num_valid_folds:
            current_folding, current_score = self.get_valid_folding()
            if current_folding is None:
                continue

            temperature = self.initial_temperature
            best_folding = current_folding.copy()
            best_score = current_score

            for _ in range(self.max_iterations):
                # Generate a neighboring state
                neighbor = self.generate_single_neighbor(current_folding)

                if not self.graph.apply_folding(neighbor):
                    continue

                neighbor_score = self.graph.calculate_score()
                delta_energy = neighbor_score - current_score

                if delta_energy < 0 or random.random() < math.exp(-delta_energy / temperature):
                    current_folding = neighbor
                    current_score = neighbor_score
                
                if current_score < best_score:
                    best_folding = current_folding.copy()
                    best_score = current_score
                
                temperature *= (1 - self.cooling_rate)
            
            self.all_scores.append(best_score)
            
            if best_score < best_overall_score:
                best_overall_score = best_score
                best_overall_folding = best_folding.copy()
            
            valid_attempts += 1

        return best_overall_folding, best_overall_score, self.all_scores