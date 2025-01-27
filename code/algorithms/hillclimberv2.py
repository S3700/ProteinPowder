import random
from code.classes.graph import Graph

class HillClimberv2:
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
        self.h_indices = [i for i, aa in enumerate(protein_sequence) if aa == 'H']

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
                        proximity_bonus = self.calculate_c_and_h_proximity(temp_graph)
                        
                        neighbors.append((neighbor, score, proximity_bonus))
                    
                    # Reset the graph state
                    temp_graph.apply_folding(current_folding)
        
        # Sort neighbors: primary sort by score, secondary by C-proximity (descending)
        # This ensures that if two configurations have the same score, 
        # the one with more C-amino proximity is preferred
        neighbors.sort(key=lambda x: (x[1], -x[2]))
        
        # Return just the folding and score, dropping the proximity bonus
        return [(n[0], n[1]) for n in neighbors]

    def calculate_c_and_h_proximity(self, graph):
        """
        Calculate proximity and surrounding of C-amino acids by H-amino acids.
        This does NOT modify the score, only provides a proximity measure.
        """
        proximity_score = 0
        
        # Get positions of C and H amino acids
        c_positions = [graph.amino_acids[i].position for i in self.c_indices]
        h_positions = [graph.amino_acids[i].position for i in self.h_indices]
        
        # Directions for checking adjacency (6 possible directions in 3D)
        directions = [
            (1, 0, 0), (-1, 0, 0),
            (0, 1, 0), (0, -1, 0),
            (0, 0, 1), (0, 0, -1)
        ]
        
        # Check proximity between C-amino acids
        for i in range(len(c_positions)):
            for j in range(i+1, len(c_positions)):
                # Calculate Manhattan distance
                dist = sum(abs(a-b) for a,b in zip(c_positions[i], c_positions[j]))
                
                # Bonus for being adjacent (distance of 1)
                if dist == 1:
                    proximity_score += 5  # More weight for direct C-C proximity
        
        # Check H-amino surrounding for each C-amino
        for c_pos in c_positions:
            h_neighbors = 0
            
            # Check each direction for H-amino neighbors
            for direction in directions:
                neighbor_pos = tuple(c + d for c, d in zip(c_pos, direction))
                
                # Check if the neighboring position is occupied by an H-amino
                if neighbor_pos in h_positions:
                    h_neighbors += 3
            
            # More bonus for more H-amino neighbors
            proximity_score += h_neighbors * 3
        
        return proximity_score

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