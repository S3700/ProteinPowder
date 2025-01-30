import time
from code.classes.graph import Graph

class DepthFirst:
    def __init__(self, protein_sequence, num_valid_folds=1):
        """
        Initialize the DepthFirstSearch algorithm for finding the best protein folding.
        """
        self.protein_sequence = protein_sequence
        self.num_valid_folds = num_valid_folds
        self.graph = Graph(protein_sequence)
        self.best_folding = None
        self.best_score = float('inf')
        self.all_scores = []
        self.chunk_size = 5
        
        # Initialize timing
        self.start_time = None
        self.max_runtime = float('inf')
        
        # Store best foldings for each chunk
        self.best_chunk_foldings = [[] for _ in range((len(protein_sequence) - 1 + self.chunk_size - 1) // self.chunk_size)]

    def is_time_exceeded(self):
        """
        Checks if the maximum runtime has been exceeded.
        """
        return time.time() - self.start_time >= self.max_runtime

    def get_possible_directions(self, position_in_chunk):
        """
        Heuristic adds constraints to the fifth amino in a chunk.
        Returns possible directions based on position in chunk.
        """
        if position_in_chunk < 4:
            return [1, -1, 2, -2, 3, -3]  # All directions allowed except 'up' for first four positions
        else:
            return [3]  # Only 'up' for fifth position

    def create_search_states(self, base_folding, base_score):
        """
        Creates the initial search states container for DepthFirst.
        """
        return [(base_folding, 0, base_score)]

    def add_search_state(self, states, new_state):
        """
        Adds a new state to the search states container.
        For DepthFirst, appends to end (stack behavior).
        """
        states.append(new_state)

    def get_next_state(self, states):
        """
        Gets the next state from the search states container.
        For DepthFirst, pops from end (stack behavior).
        """
        return states.pop()

    def process_valid_folding(self, best_foldings, best_score, current_folding, current_score):
        """
        Processes a valid folding and updates best_foldings if necessary.
        """
        if current_score <= best_score:
            if current_score < best_score:
                best_foldings = []
                best_score = current_score
            best_foldings.append((current_folding.copy(), current_score))
            
            # Keep only top 10 best foldings
            if len(best_foldings) > 10:
                best_foldings.sort(key=lambda x: x[1])
                best_foldings = best_foldings[:10]
                
        return best_foldings, best_score

    def explore_chunk(self, chunk_index, base_foldings):
        """
        Explores all possibilities for a chunk using the search algorithm,
        building upon the best foldings from the previous chunk.
        """
        start = chunk_index * self.chunk_size
        end = min(start + self.chunk_size, len(self.protein_sequence) - 1)
        chunk_length = end - start
        
        best_foldings = []
        best_score = float('inf')
        
        for base_folding, base_score in base_foldings:
            states = self.create_search_states(base_folding, base_score)
            
            while states and not self.is_time_exceeded():
                current_folding, pos_in_chunk, current_score = self.get_next_state(states)
                
                # If we've completed this chunk
                if pos_in_chunk == chunk_length:
                    # Only add valid folds to best_foldings
                    if self.graph.apply_folding(current_folding):
                        score = self.graph.calculate_score()
                        best_foldings, best_score = self.process_valid_folding(
                            best_foldings, best_score, current_folding, score
                        )
                    continue
                
                # Try each possible direction
                for direction in self.get_possible_directions(pos_in_chunk):
                    new_folding = current_folding + [direction]
                    if self.graph.apply_folding(new_folding):
                        score = self.graph.calculate_score()
                        self.add_search_state(states, (new_folding, pos_in_chunk + 1, score))
        
        return best_foldings

    def find_solutions(self):
        """
        Runs the search algorithm chunk by chunk.
        """
        # Initialize timing
        self.start_time = time.time()
        
        # Get max_runtime from experiment.py
        import inspect
        for frame in inspect.stack():
            if 'experiment.py' in frame.filename:
                locals_dict = frame.frame.f_locals
                if 'self' in locals_dict and hasattr(locals_dict['self'], 'runtime'):
                    self.max_runtime = locals_dict['self'].runtime
                    break
        
        # Start with empty folding
        current_foldings = [([], 0)]
        
        # Process each chunk
        num_chunks = (len(self.protein_sequence) - 1 + self.chunk_size - 1) // self.chunk_size
        
        for chunk_index in range(num_chunks):
            if self.is_time_exceeded():
                break
            
            # Explore current chunk based on best previous foldings
            best_chunk_foldings = self.explore_chunk(chunk_index, current_foldings)
            
            if not best_chunk_foldings:
                print(f"No valid foldings found for chunk {chunk_index + 1}")
                break
            
            # Update current foldings for next chunk
            current_foldings = best_chunk_foldings
            
            # Update best overall solution if we're at the last chunk
            if chunk_index == num_chunks - 1:
                self.best_folding = best_chunk_foldings[0][0]
                self.best_score = best_chunk_foldings[0][1]
                self.all_scores.append(self.best_score)
        
        # If no solution was found, return empty results
        if not self.best_folding:
            self.best_folding = []
            self.best_score = 0
        
        return self.best_folding, self.best_score, self.all_scores


