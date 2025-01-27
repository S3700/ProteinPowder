from code.classes.graph import Graph
import random
import math

class DepthFirstSearch:
    def __init__(self, protein_sequence, frame_size=5, max_folds=120, keep=0.025):
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
        
        print(f"Generated {len(frame_sequences)} windows: {frame_sequences}")
        return frame_sequences

    def dfs_r_frames(self, current_frame, num_r_frames):
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

    def dfs_cull(self, foldings, current_index):
        """
        Cull all sequences except the top 10% based on their scores.
        """
        if not foldings:
            print("No foldings to cull.")
            return []

        print(f"Culling at window {current_index + 1}")

        scores = [self.evaluate_folding(folding) for folding in foldings]
        num_to_keep = max(1, int(len(scores) * self.keep)) 
        sorted_foldings = sorted(zip(foldings, scores), key=lambda x: x[1])
        culled_foldings = [folding for folding, score in sorted_foldings[:num_to_keep]]

        print(f"Culled {len(foldings) - len(culled_foldings)} foldings below top {self.keep * 100}% scores.")

        return culled_foldings

    def dfs_recursive(self, frames, current_index, current_foldings):
        """
        Recursive depth-first search to explore all possible foldings.
        """
        if current_index == len(frames):
            return current_foldings

        foldings = []
        random_foldings = self.dfs_r_frames(frames[current_index], self.max_folds)
        culled_foldings = self.dfs_cull(random_foldings, current_index)

        for folding in culled_foldings:
            for prev_folding in current_foldings:
                new_folding = prev_folding.copy()

                if len(new_folding) > 0 and len(folding) > 0:
                    new_folding[-1] = folding[0]

                new_folding.extend(folding[1:])
                foldings.append(new_folding)

        return self.dfs_recursive(frames, current_index + 1, foldings)

    def dfs(self):
        """
        Perform depth-first search to explore all possible foldings.
        """
        frames = self.generate_frames()
        initial_foldings = self.dfs_r_frames(frames[0], self.max_folds)
        culled_foldings = self.dfs_cull(initial_foldings, 0)
        all_foldings = self.dfs_recursive(frames, 1, culled_foldings)
        return all_foldings

    def find_best_solution(self):
        """
        Find the best folding for the protein sequence using DFS.
        """
        all_foldings = self.dfs()
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
        return f"DepthFirst(sequence={self.protein_sequence}, dimension={self.dimension})"
