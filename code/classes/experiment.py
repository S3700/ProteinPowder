import time
import csv
import os
from datetime import datetime
from functools import wraps

class TimedExperiment:
    def __init__(self, algorithm, max_runtime):
        """
        Initialize TimedExperiment with specific algorithm and runtime.
        Runs the algorithm for the given runtime and stores results.
        """
        self.algorithm = algorithm
        self.runtime = max_runtime
        self.output_dir = "experiment_results"
        self.start_time = None
        
        # Create output directory if it doesn't exist
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def is_time_exceeded(self):
        """
        Check if the maximum runtime has been exceeded.
        """
        return time.time() - self.start_time >= self.runtime

    def time_remaining(self):
        """
        Get remaining runtime in seconds.
        """
        return max(0, self.runtime - (time.time() - self.start_time))

    def run_with_timeout(self, func):
        """
        Decorator to run a function with timeout checking.
        """
        @wraps(func)
        def wrapper(*args, **kwargs):
            if self.is_time_exceeded():
                return None
            return func(*args, **kwargs)
        return wrapper

    def run(self, protein):
        """
        Run the experiment for the specified runtime and collect scores.
        Only save results that were completed within the runtime limit.
        """
        self.start_time = time.time()
        solver = self.algorithm(protein)
        
        # Inject timing methods into the solver
        solver.is_time_exceeded = self.is_time_exceeded
        solver.time_remaining = self.time_remaining
        
        final_runtime = 0 
        best_folding = None
        best_score = float('inf')
        all_scores = []
        last_update_time = self.start_time
        update_interval = 10

        while not self.is_time_exceeded():
            current_time = time.time()

            # Update progress every interval
            if current_time - last_update_time >= update_interval:
                progress_percentage = ((current_time - self.start_time) / self.runtime) * 100
                print(f"Progress: {progress_percentage:.2f}% complete", end='\r')
                last_update_time = current_time

            # Run the algorithm and get results
            folding, score, scores_from_run = solver.find_solutions()
            
            # Break if no valid solution was found
            if folding is None:
                break

            # Calculate elapsed time
            elapsed_time = time.time() - self.start_time

            # Only collect results if within time limit
            if elapsed_time < self.runtime:
                if scores_from_run:
                    all_scores.extend(scores_from_run)
                final_runtime = elapsed_time

                if score is not None and score < best_score:
                    best_score = score
                    best_folding = folding
            else:
                break

        # Save results if any valid solutions were found
        if all_scores:
            self.save_results(all_scores, best_score, len(all_scores), final_runtime, protein)
        else:
            print("\nNo valid solutions found within the time limit.")
            return None, None, []

        return best_folding, best_score, all_scores

    def save_results(self, all_scores, best_score, num_solutions, final_runtime, protein):
        """
        Save experiment results to a CSV file.
        """
        if not all_scores:
            return
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{self.algorithm.__name__}_{timestamp}.csv"
        filepath = os.path.join(self.output_dir, filename)

        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Write experiment metadata
            writer.writerow(['Experiment Info', 'Values'])
            writer.writerow(['Algorithm', self.algorithm.__name__])
            writer.writerow(['Protein', protein])
            writer.writerow(['Final Runtime (s)', f"{final_runtime:.2f}"])
            writer.writerow(['Total Solutions', num_solutions])
            writer.writerow(['Best Score', best_score])
            writer.writerow(['Average Score', f"{sum(all_scores) / len(all_scores):.2f}"])
            writer.writerow([])
            writer.writerow(['Score Index', 'Score'])

            # Write all recorded scores
            for i, score in enumerate(all_scores):
                writer.writerow([i, score])