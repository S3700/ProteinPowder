import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from code.classes.utilities import is_adjacent

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
