from code.classes.graph import Graph

class DepthFirst:
    def __init__(self, protein_sequence, dimension=3, max_folds=10000):
        """
        Initialize the DepthFirst algorithm for finding the best folding.
        """
        self.protein_sequence = protein_sequence
        self.dimension = dimension
        self.max_folds = max_folds
        self.graph = Graph(protein_sequence, dimension)
        self.best_score = float('inf')
        self.best_folding = None
        self.all_foldings = []
        self.foldings_generated = 0

    def dfs(self, folding, depth, direction = None):
        """
        Perform depth-first search to explore all possible foldings.
        """
        if self.foldings_generated >= self.max_folds:
            return

        if depth == len(self.protein_sequence) - 1:
            if self.graph.apply_folding(folding):
                score = self.graph.calculate_score()
                self.all_foldings.append((folding.copy(), score))
                self.foldings_generated += 1
                if score < self.best_score:
                    self.best_score = score
                    self.best_folding = folding.copy()
            return

        possible_directions = ([1, -1, 2, -2, 3, -3])
        if direction:
            opposite_direction = tuple(-x for x in direction)
            possible_directions = [d for d in possible_directions if d != opposite_direction]

        for direction in possible_directions:
            if self.foldings_generated >= self.max_folds:
                break
            folding.append(direction)
            self.dfs(folding, depth + 1)
            folding.pop()

    def find_best_solution(self):
        """
        Find the best folding for the protein sequence using DFS.
        """
        self.dfs([], 0)
        self.graph.apply_folding(self.best_folding)
        all_scores = [score for _, score in self.all_foldings]
        return self.best_folding, self.best_score, all_scores

    def __repr__(self):
        return f"DepthFirst(sequence={self.protein_sequence}, dimension={self.dimension})"