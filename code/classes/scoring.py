from code.classes.utilities import is_adjacent

def calculate_score(protein, folding, dimension=3):
    """
    Calculate the score of a protein folding in 3D.

    :param protein: Protein sequence as a string (e.g., "HHPHPH").
    :param folding: Folding as a list of steps (e.g., [1, 2, -1, -2, 1]).
    :param dimension: Number of dimensions (default 3).
    :return: Total score.
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
        coords.append(next_coord)

    score = 0
    for i, coord1 in enumerate(coords):
        for j, coord2 in enumerate(coords):
            if abs(i - j) > 1 and coord1 == coord2:
                aa1, aa2 = protein[i], protein[j]
                if aa1 == 'H' and aa2 == 'H':
                    score -= 1
                elif aa1 == 'C' and aa2 == 'C':
                    score -= 5
                elif (aa1 == 'H' and aa2 == 'C') or (aa1 == 'C' and aa2 == 'H'):
                    score -= 1

    for i, coord1 in enumerate(coords):
        for j, coord2 in enumerate(coords):
            if abs(i - j) > 1 and coord1 != coord2 and is_adjacent(coord1, coord2):
                aa1, aa2 = protein[i], protein[j]
                if aa1 == 'H' and aa2 == 'H':
                    score -= 1

    return score