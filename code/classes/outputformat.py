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