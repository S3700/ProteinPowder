def is_adjacent(coord1, coord2):
    """
    Check if two coordinates are adjacent on a grid.
    
    :param coord1: First coordinate (x, y, z).
    :param coord2: Second coordinate (x, y, z).
    :return: True if adjacent, otherwise False.
    """
    return sum(abs(a - b) for a, b in zip(coord1, coord2)) == 1