from code.classes.graph import Graph
import random
import math

class BreadthFirstSearch:
    def __init__(self, protein_sequence, frame_size=5, max_folds=500, keep=0.005):
        self.protein_sequence = protein_sequence
        self.frame_size = frame_size
        self.max_folds = max_folds
        self.keep = keep
        self.graph = Graph(protein_sequence, 3) 

    def generate_frames(self):
        frame_sequences = {}
        sequence_length = len(self.protein_sequence)

        for i in range(1, sequence_length, self.frame_size):
            end = min(i + self.frame_size, sequence_length)
            frame_sequences[i // self.frame_size] = self.protein_sequence[i:end]
        
        return frame_sequences

    def bfs_r_frames(self, current_frame, num_r_frames):
        """
        Generate random folding sequences for the current frame.
        """
        random_frames = []
        possible_directions = [1, -1, 2, -2, 3, -3]

        for _ in range(num_r_frames):
            folding = []
            for _ in current_frame:
                direction = random.choice(possible_directions)
                folding.append(direction)
            random_frames.append(folding)

        return random_frames

    def evaluate_folding(self, folding):
        """
        Evaluate a folding and return its score.
        A lower score represents a better folding.
        """
        self.graph.apply_folding(folding)
        score = self.graph.calculate_score()  
        return score

    def bfs_cull(self, foldings):
        """
        Cull all sequences except the top 5% based on their scores.
        """
        if not foldings:
            print("No foldings to cull.")
            return []

        scores = [self.evaluate_folding(folding) for folding in foldings]
        num_to_keep = max(1, int(len(scores) * self.keep))  
        sorted_foldings = sorted(zip(foldings, scores), key=lambda x: x[1])
        culled_foldings = [folding for folding, score in sorted_foldings[:num_to_keep]]

        print(f"Culled {len(foldings) - len(culled_foldings)} foldings below top {self.keep * 100}% scores.")

        return culled_foldings

    def bfs_recursive(self, frames, current_index, current_foldings):
        """
        Recursive breadth-first search to explore all possible foldings.
        """
        if current_index == len(frames):
            return current_foldings

        foldings = []
        random_foldings = self.bfs_r_frames(frames[current_index], self.max_folds)
        culled_foldings = self.bfs_cull(random_foldings)

        possible_directions = [1, -1, 2, -2, 3, -3]

        for folding in culled_foldings:
            for prev_folding in current_foldings:
                best_new_folding = None
                best_score = float('inf')
                for direction in possible_directions:
                    new_folding = prev_folding + [direction] + folding[1:]
                    score = self.evaluate_folding(new_folding)
                    if score < best_score:
                        best_score = score
                        best_new_folding = new_folding
                foldings.append(best_new_folding)

        return self.bfs_recursive(frames, current_index + 1, foldings)

    def bfs(self):
        """
        Perform breadth-first search to explore all possible foldings.
        """
        frames = self.generate_frames()
        initial_foldings = self.bfs_r_frames(frames[0], self.max_folds)
        culled_foldings = self.bfs_cull(initial_foldings)
        all_foldings = self.bfs_recursive(frames, 1, culled_foldings)
        return all_foldings

    def find_best_solution(self):
        """
        Find the best folding for the protein sequence using BFS.
        """
        all_foldings = self.bfs()
        best_folding = None
        best_score = float('inf')
        all_scores = []

        for folding in all_foldings:
            score = self.evaluate_folding(folding)
            all_scores.append(score)
            if score < best_score:
                best_score = score
                best_folding = folding

        self.best_folding = best_folding
        self.best_score = best_score

        self.graph.apply_folding(self.best_folding)

        print(f"Best selected folding: {self.best_folding}")
        print(f"Protein sequence: {self.protein_sequence}")
        print(f"Length of folding: {len(self.best_folding)}")
        print(f"Length of protein: {len(self.protein_sequence)}")

        return self.best_folding, self.best_score, all_scores

    def __repr__(self):
        return f"BreadthFirst(sequence={self.protein_sequence}, dimension={self.dimension})"