class Aminoacid:
    def __init__(self, aa_type, position):
        """
        Aminoacid class representing an amino acid in the folding.
        """
        self.aa_type = aa_type  # Type of the amino acid (H, P, or C)
        self.position = position  # Position of the amino acid in the 3D space (x, y, z)
        self.neighbours = []  # List to store neighbouring amino acids
        self.folding_options = [(0, 0, 1), (0, 0, -1), (0, 1, 0), (0, -1, 0), (1, 0, 0), (-1, 0, 0)]  # Possible folding directions

    def add_neighbour(self, neighbour):
        """
        Adds a neighbour to current amino acid.
        """
        self.neighbours.append(neighbour)

    def calculate_score_with_neighbour(self, neighbour):
        """
        Calculates the score for two neighbouring amino acids based on their types.
        """
        if self.aa_type == 'H' and neighbour.aa_type == 'H':
            return -1  # H-H bond
        elif self.aa_type == 'C' and neighbour.aa_type == 'C':
            return -5  # C-C bond
        elif (self.aa_type == 'H' and neighbour.aa_type == 'C') or (self.aa_type == 'C' and neighbour.aa_type == 'H'):
            return -1  # H-C or C-H bond
        return 0  # No bonds
    
    def __repr__(self):
        return f"Aminoacid({self.aa_type}, {self.position})"