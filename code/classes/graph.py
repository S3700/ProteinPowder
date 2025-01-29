from .aminoacid import Aminoacid

class Graph:
    def __init__(self, protein_sequence):
        """
        Initialize the Graph to store amino acids and manage folding.
        """
        self.protein_sequence = protein_sequence
        self.amino_acids = self.create_amino_acids(protein_sequence)
        self.score = 0
    
    def create_amino_acids(self, protein_sequence):
        """
        Creates a list of Aminoacid objects based on the protein sequence.
        """
        amino_acids = []
        position = (0, 0, 0)  # Starting position (origin)
    
        # Placing the first amino at the origin
        amino_acid = Aminoacid(protein_sequence[0], position)
        amino_acids.append(amino_acid)
    
        # For remaining amino's, set positions based on folding directions
        for aa_type in protein_sequence[1:]:
            direction = (1, 0, 0)  # Default direction from origin
            position = self.get_next_position(position, direction)
        
            amino_acid = Aminoacid(aa_type, position)
            amino_acids.append(amino_acid)
    
        # Link neighbours for each amino acid
        for i in range(1, len(amino_acids)):
            amino_acids[i-1].add_neighbour(amino_acids[i])
            amino_acids[i].add_neighbour(amino_acids[i-1])
    
        return amino_acids

    
    def get_next_position(self, current_position, direction):
        """
        Get the next position based on the direction and current position.
        """
        x, y, z = current_position
        dx, dy, dz = direction
        return (x + dx, y + dy, z + dz)
    
    def apply_folding(self, folding):
        """
        Apply the folding configuration to the graph.
        """
        if folding is None:
            print("Error: Folding is None.")
            return False
        
        directions = {
            1: (1, 0, 0),
            -1: (-1, 0, 0),
            2: (0, 1, 0),
            -2: (0, -1, 0),
            3: (0, 0, 1),
            -3: (0, 0, -1)
        }

        coords = [(0, 0, 0)]  # Starting at origin (0, 0, 0)
        self.visited_positions = set(coords)  # Track the visited positions, start with origin

        for i, step in enumerate(folding):
            last_coord = coords[-1]  # Get the last position
            direction = directions[step]  # Get the direction for the current fold step

            next_coord = tuple(sum(x) for x in zip(last_coord, direction))  # Calculate the next position

            # Check for crossing (previously visited position)
            if next_coord in self.visited_positions:
                return False  # Invalid folding due to crossing

            self.amino_acids[i + 1].position = next_coord  # Update position of next amino acid
            self.visited_positions.add(next_coord)  # Mark the position as visited

            coords.append(next_coord)  # Add the new coordinate to the list of coordinates
    
        return True  # Folding is valid (no crossings)

    
    def calculate_score(self):
        """
        Calculate the total score for the folding by considering all amino acid pairs.
        """
        total_score = 0
        
        # Check direct and indirect neighbours
        for i, aa in enumerate(self.amino_acids):            
            for j in range(i + 2, len(self.amino_acids)):  # Checking from i+2 to avoid direct neighbours
                indirect_neighbour = self.amino_acids[j]
                # Check if amino acids are adjacent in 3D space
                if self.is_adjacent_in_space(aa, indirect_neighbour):
                    total_score += aa.calculate_score_with_neighbour(indirect_neighbour)
        
        self.score = total_score
        return self.score
    
    def is_adjacent_in_space(self, aa1, aa2):
        """
        Check if two amino acids are adjacent in the 3D space.
        Only returns True if they are adjacent in exactly one dimension.
        """
        x1, y1, z1 = aa1.position
        x2, y2, z2 = aa2.position
    
        # Check if the two amino acids' positions
        diff_x = abs(x1 - x2)
        diff_y = abs(y1 - y2)
        diff_z = abs(z1 - z2)
    
        # Adjacent if they differ by 1 in one direction and 0 in the other two directions
        return  (diff_x == 1 and diff_y == 0 and diff_z == 0) or \
                (diff_x == 0 and diff_y == 1 and diff_z == 0) or \
                (diff_x == 0 and diff_y == 0 and diff_z == 1)

    def __repr__(self):
        return f"Graph({self.protein_sequence}, Score: {self.score})"