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
        self.c_indices = [i for i, aa in enumerate(protein_sequence) if aa == 'C']

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
        Includes a heuristic to prefer placing C-amino acids next to each other.
        """
        neighbors = []
        possible_directions = [1, -1, 2, -2, 3, -3]
        
        # Directions for moving around C-amino acids
        c_proximity_directions = {
            1: (1, 0, 0), -1: (-1, 0, 0),
            2: (0, 1, 0), -2: (0, -1, 0),
            3: (0, 0, 1), -3: (0, 0, -1)
        }
        
        # Temporary graph to test foldings
        temp_graph = Graph(self.protein_sequence, self.dimension)
        
        # Try changing each position in the folding
        for i in range(len(current_folding)):
            for direction in possible_directions:
                if direction != current_folding[i]:  # Only try different directions
                    neighbor = current_folding.copy()
                    neighbor[i] = direction
                    
                    # Check if the neighbor folding is valid
                    if temp_graph.apply_folding(neighbor):
                        score = temp_graph.calculate_score()
                        
                        # Heuristic: Analyze C-amino proximity
                        # Check how many C-amino acids are close to each other
                        c_proximity_bonus = self.calculate_c_proximity(temp_graph)
                        
                        neighbors.append((neighbor, score, c_proximity_bonus))
                    
                    # Reset the graph state
                    temp_graph.apply_folding(current_folding)
        
        # Sort neighbors: primary sort by score, secondary by C-proximity (descending)
        # This ensures that if two configurations have the same score, 
        # the one with more C-amino proximity is preferred
        neighbors.sort(key=lambda x: (x[1], -x[2]))
        
        # Return just the folding and score, dropping the proximity bonus
        return [(n[0], n[1]) for n in neighbors]

    def calculate_c_proximity(self, graph):
        """
        Calculate how close C-amino acids are to each other.
        This does NOT modify the score, just provides a proximity measure.
        """
        c_proximity = 0
        c_positions = [graph.amino_acids[i].position for i in self.c_indices]
        
        # Check proximity between C-amino acids
        for i in range(len(c_positions)):
            for j in range(i+1, len(c_positions)):
                # Calculate Manhattan distance
                dist = sum(abs(a-b) for a,b in zip(c_positions[i], c_positions[j]))
                
                # Bonus for being adjacent (distance of 1)
                if dist == 1:
                    c_proximity += 1
        
        return c_proximity

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