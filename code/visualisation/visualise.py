import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_folding(ax, coords, protein):
    """
    Plot the protein folding as a 3D line with amino acid labels.
    """
    xs, ys, zs = zip(*coords)
    ax.plot(xs, ys, zs, '-o', color='gray', linewidth=2, alpha=0.5)
    
    # Add labels for amino acids
    for i, (x, y, z) in enumerate(coords):
        label_color = 'red' if protein[i] == 'H' else ('yellow' if protein[i] == 'C' else 'blue')
        ax.text(x, y, z, protein[i], fontsize=12, ha='center', va='center',
                bbox=dict(facecolor=label_color, edgecolor='black', alpha=0.7, boxstyle='circle'))

def calculate_coordinates(protein, folding):
    """
    Calculate all coordinates based on the folding sequence.
    Returns both the coordinates and a mapping of positions to indices.
    """
    directions = {
        1: (1, 0, 0), -1: (-1, 0, 0),
        2: (0, 1, 0), -2: (0, -1, 0),
        3: (0, 0, 1), -3: (0, 0, -1)
    }
    
    coords = [(0, 0, 0)]
    position_to_index = {(0, 0, 0): 0}
    
    for i, step in enumerate(folding):
        last_coord = coords[-1]
        direction = directions[step]
        next_coord = tuple(sum(x) for x in zip(last_coord, direction))
        coords.append(next_coord)
        position_to_index[next_coord] = i + 1
        
    return coords, position_to_index

def plot_bonds(ax, coords, protein, position_to_index):
    """
    Plot indirect bonds between amino acids that are adjacent in 3D space.
    """
    # Define bond styles based on amino acid types
    bond_styles = {
        ('H', 'H'): ('red', '--', 2),
        ('C', 'C'): ('green', '--', 2),
        ('H', 'C'): ('blue', '--', 2),
        ('C', 'H'): ('blue', '--', 2)
    }
    
    # Check each pair of coordinates for potential bonds
    for i, coord1 in enumerate(coords):
        x1, y1, z1 = coord1
        for j, coord2 in enumerate(coords):
            # Skip if same amino acid or direct neighbors
            if abs(i - j) <= 1:
                continue
                
            x2, y2, z2 = coord2
            # Check if amino acids are adjacent in 3D space
            if (abs(x1 - x2) + abs(y1 - y2) + abs(z1 - z2)) == 1:
                aa1_type = protein[i]
                aa2_type = protein[j]
                
                # Get bond style if there is a bond
                bond_style = bond_styles.get((aa1_type, aa2_type))
                if bond_style:
                    color, style, width = bond_style
                    ax.plot([x1, x2], [y1, y2], [z1, z2], 
                           color=color, linestyle=style, linewidth=width, alpha=0.7)

def print_visual(protein, folding, graph):
    """
    Create visual representation of the protein folding in 3D.
    """
    # Calculate coordinates and create position mapping
    coords, position_to_index = calculate_coordinates(protein, folding)
    
    # Create figure size and view angle
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot main folding structure
    plot_folding(ax, coords, protein)
    
    # Plot indirect bonds
    plot_bonds(ax, coords, protein, position_to_index)
    
    # Set labels and viewing angle
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_zlabel('Z', fontsize=12)
    
    ax.view_init(elev=20, azim=45)
    
    plt.tight_layout()
    plt.show()
