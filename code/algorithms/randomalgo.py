import random
from code.classes.scoring import calculate_score
from code.classes.folding import is_valid_folding
from code.visualisation.visualise import print_visual
from code.classes.outputformat import output_format

def generate_random_folding(protein, dimension=3):
    """
    Generate a random folding for a protein sequence. The folding consists of random directions (1, -1, 2, -2, 3, -3).
    
    :param protein: Protein sequence as a string.
    :param dimension: The dimension for folding (default is 3).
    :return: Random folding as a list of steps.
    """
    # Randomly generate a list of folding steps (directions)
    folding = [random.choice([1, -1, 2, -2, 3, -3]) for _ in range(len(protein) - 1)]
    
    # Check if the folding is valid, if not, recursively generate a new one
    while not is_valid_folding(folding, dimension):
        folding = [random.choice([1, -1, 2, -2, 3, -3]) for _ in range(len(protein) - 1)]
    
    return folding

def find_best_random_folding(protein, dimension=3, num_random_folds=100):
    """
    Generate `num_random_folds` random foldings for a protein, calculate their scores,
    and return the folding with the highest score.
    
    :param protein: Protein sequence as a string.
    :param dimension: The dimension for folding (default is 3).
    :param num_random_folds: The number of random foldings to generate and evaluate.
    :return: Best folding and its score.
    """
    best_folding = None
    best_score = float('inf')  # Start with a very bad score (high score is preferred)

    # Generate random foldings and check the score
    for _ in range(num_random_folds):
        folding = generate_random_folding(protein, dimension)
        score = calculate_score(protein, folding, dimension)

        # Update best folding if the current one is better (lower score)
        if score < best_score:
            best_score = score
            best_folding = folding
    
    return best_folding, best_score