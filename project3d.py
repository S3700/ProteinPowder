import itertools
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def calculate_score(protein, folding, dimension=3):
    """
    Calculate the score of a protein folding in 3D.

    :param protein: Protein sequence as a string (e.g., "HHPHPH").
    :param folding: Folding as a list of steps (e.g., [1, 2, -1, -2, 1]).
    :param dimension: Number of dimensions (default 3).
    :return: Total score.
    """
    # Directions for 3D movement
    directions = {
        1: (1, 0, 0),
        -1: (-1, 0, 0),
        2: (0, 1, 0),
        -2: (0, -1, 0),
        3: (0, 0, 1),
        -3: (0, 0, -1)
    }

    # Calculate coordinates of amino acids
    coords = [(0, 0, 0)]
    for step in folding:
        last_coord = coords[-1]
        direction = directions[step]
        next_coord = tuple(sum(x) for x in zip(last_coord, direction))
        coords.append(next_coord)

    # Calculate score
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

    # Bonus: Count H-H interactions for grid neighbors not connected via chain
    for i, coord1 in enumerate(coords):
        for j, coord2 in enumerate(coords):
            if abs(i - j) > 1 and coord1 != coord2 and is_adjacent(coord1, coord2):
                aa1, aa2 = protein[i], protein[j]
                if aa1 == 'H' and aa2 == 'H':
                    score -= 1

    return score

def is_adjacent(coord1, coord2):
    """
    Check if two coordinates are adjacent on a grid.
    
    :param coord1: First coordinate (x, y, z).
    :param coord2: Second coordinate (x, y, z).
    :return: True if adjacent, otherwise False.
    """
    return sum(abs(a - b) for a, b in zip(coord1, coord2)) == 1

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
    Generate an optimal folding for a protein.

    :param protein: Protein sequence as a string.
    :param dimension: Dimension (default 3).
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

def output_format(protein, folding, score):
    """
    Generate the output in the correct format.

    :param protein: Protein sequence.
    :param folding: Folding as a list of steps.
    :param score: Score of the folding.
    :return: Formatted output.
    """
    result = ["amino,fold"]
    for amino, step in zip(protein, list(folding) + [0]):
        result.append(f"{amino},{step}")
    result.append(f"score,{score}")
    return "\n".join(result)

def print_visual(protein, folding):
    """
    Create a visual representation of the protein folding in 3D.
    Includes red dashed lines for H-H bonds.

    :param protein: Protein sequence.
    :param folding: Folding as a list of steps.
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

    h_bonds = []
    for i, coord1 in enumerate(coords):
        for j, coord2 in enumerate(coords):
            if abs(i - j) > 1 and protein[i] == "H" and protein[j] == "H" and is_adjacent(coord1, coord2):
                h_bonds.append((coord1, coord2))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    xs, ys, zs = zip(*coords)
    ax.plot(xs, ys, zs, '-o', color='gray')

    for i, (x, y, z) in enumerate(coords):
        ax.text(x, y, z, protein[i], fontsize=10, ha='center', va='center',
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='circle'))

    for (x1, y1, z1), (x2, y2, z2) in h_bonds:
        ax.plot([x1, x2], [y1, y2], [z1, z2], 'r--', linewidth=1)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

# Test with example sequences
proteins = ["HHPPHPH"]
for protein in proteins:
    folding, score = generate_folding(protein, dimension=3)
    print(output_format(protein, folding, score))
    print("Visual representation:")
    print_visual(protein, folding)
    print("\n")
