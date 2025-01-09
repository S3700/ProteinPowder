import random

def is_valid_folding(folding, dimension=3):
    """
    Check if a folding is valid (no overlap).
    
    :param folding: Folding as a list of steps.
    :param dimension: Dimension (default 3).
    :return: True if valid, otherwise False.
    """
    directions = {
        1: (1, 0, 0),
        -1: (-1, 0, 0),
        2: (0, 1, 0),
        -2: (0, -1, 0),
        3: (0, 0, 1),
        -3: (0, 0, -1)
    }

    coords = [(0, 0, 0)]
    for step in folding:
        last_coord = coords[-1]
        direction = directions[step]
        next_coord = tuple(sum(x) for x in zip(last_coord, direction))
        if next_coord in coords:
            return False
        coords.append(next_coord)

    return len(coords) == len(set(coords))

def generate_folding(protein, dimension=3):
    """
    Generate a random folding for a protein.

    :param protein: Protein sequence as a string.
    :param dimension: Dimension for folding (default is 3).
    :return: A randomly generated folding.
    """
    folding = [random.choice([1, -1, 2, -2, 3, -3]) for _ in range(len(protein) - 1)]
    while not is_valid_folding(folding, dimension):
        folding = [random.choice([1, -1, 2, -2, 3, -3]) for _ in range(len(protein) - 1)]
    return folding
