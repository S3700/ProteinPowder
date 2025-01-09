import itertools
from classes.folding import generate_folding, is_valid_folding
from code.classes.scoring import calculate_score

def brute_force_folding(protein, dimension=3):
    """
    Try all possible foldings and return the one with the best score.
    
    :param protein: Protein sequence as a string.
    :param dimension: Dimension for folding (default is 3).
    :return: Best folding and score.
    """
    steps = [1, -1, 2, -2, 3, -3] if dimension == 3 else [1, -1, 2, -2]
    best_score = float('inf')
    best_folding = None
    
    for folding in itertools.product(steps, repeat=len(protein) - 1):
        if is_valid_folding(folding, dimension):
            score = calculate_score(protein, folding, dimension)
            if score < best_score:
                best_score = score
                best_folding = folding
    
    return best_folding, best_score