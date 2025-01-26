import random
from code.classes.graph import Graph

class HillClimberv3:
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
            return folding
        return None

    def calculate_temporary_score(self):
        """
        Calculate a temporary score that rewards indirect bonds for C and H amino acids.
        Aims to maximize the number of indirect bonds up to a maximum of 3 per amino acid.
        """
        # First, calculate the base score
        base_score = self.graph.calculate_score()
        temporary_score = base_score

        # Track indirect bond counts for C and H amino acids
        c_indices = [i for i, aa in enumerate(self.graph.amino_acids) if aa.aa_type == 'C']
        h_indices = [i for i, aa in enumerate(self.graph.amino_acids) if aa.aa_type == 'H']

        # Penalty system for indirect C bonds
        for c_index in c_indices:
            c_amino = self.graph.amino_acids[c_index]
            indirect_c_bonds = 0
            indirect_h_bonds = 0

            for j, other_aa in enumerate(self.graph.amino_acids):
                # Skip direct neighbors and amino acid itself
                if abs(c_index - j) <= 1:
                    continue

                # Check if the other amino acid is adjacent in space
                if self.graph.is_adjacent_in_space(c_amino, other_aa):
                    # Count C-C bonds
                    if other_aa.aa_type == 'C':
                        indirect_c_bonds += 1
                    # Count C-H bonds
                    if other_aa.aa_type == 'H':
                        indirect_h_bonds += 1

            # Penalty based on number of C bonds
            if indirect_c_bonds == 0:
                temporary_score += 20  # Heavy penalty for no bonds
            elif indirect_c_bonds < 3:
                # Linear bonus for C-C bonds
                temporary_score -= indirect_c_bonds * 6
            
            # Penalty based on number of H bonds
            if indirect_h_bonds == 0:
                temporary_score += 10
            elif indirect_h_bonds < 3:
                # Linear bonus for C-H bonds
                temporary_score -= indirect_h_bonds * 2

        # Similar logic for H amino acids
        for h_index in h_indices:
            h_amino = self.graph.amino_acids[h_index]
            indirect_h_bonds = 0

            for j, other_aa in enumerate(self.graph.amino_acids):
                # Skip direct neighbors and amino acid itself
                if abs(h_index - j) <= 1:
                    continue

                # Check if other amino acid is adjacent in space
                if self.graph.is_adjacent_in_space(h_amino, other_aa):
                    # Count H-H and H-C bonds
                    if other_aa.aa_type == 'H' or other_aa.aa_type == 'C':
                        indirect_h_bonds += 1

            # Reward/penalty based on number of H bonds
            if indirect_h_bonds == 0:
                temporary_score += 10  # Penalty for no bonds
            elif indirect_h_bonds < 3:
                # Linear bonus for H bonds, lower penalty for fewer bonds
                temporary_score -= indirect_h_bonds * 3

        return temporary_score

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
                        temporary_score = self.calculate_temporary_score()
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
            current_temporary_score = self.calculate_temporary_score()
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