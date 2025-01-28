import time
import csv
import os
from datetime import datetime

class TimedExperiment:
    def __init__(self, algorithm, runtime=100):
        """
        Initialize TimedExperiment with specific algorithm and runtime.
        Runs the algorithm for the given runtime and stores results.
        """
        self.algorithm = algorithm
        self.runtime = runtime
        self.output_dir = "experiment_results"
        
        # Create output directory if it doesn't exist
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def run(self, protein):
        """
        Run the experiment for the specified runtime and collect scores.
        """
        start_time = time.time()
        solver = self.algorithm(protein)  # Initialize algorithm
    
        best_folding = None
        best_score = float('inf')
        all_scores = []

        last_print_time = 0  # Track time for periodic updates

        while (time.time() - start_time) < self.runtime:
            folding, score, scores_from_run = solver.find_best_solution()
        
            all_scores.extend(scores_from_run)
        
            if score < best_score:
                best_score = score
                best_folding = folding
        
            # Print progress every 10% of total runtime
            elapsed_time = time.time() - start_time
            if elapsed_time - last_print_time > self.runtime * 0.1:
                last_print_time = elapsed_time
                print(f"{elapsed_time:.1f}/{self.runtime}s elapsed... Current folding score: {score}")

        # Save results
        self.save_results(all_scores, best_score, len(all_scores), elapsed_time, protein)
    
        print(f"{self.algorithm.__name__} finished. Best score: {best_score}")

        return best_folding, best_score, all_scores


    def save_results(self, all_scores, best_score, num_solutions, runtime, protein):
        """
        Save experiment results (periodically) to a CSV file.
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{self.algorithm.__name__}_{timestamp}.csv"
        filepath = os.path.join(self.output_dir, filename)

        # Append data instead to save progress periodically
        with open(filepath, 'a', newline='') as f:
            writer = csv.writer(f)
            if f.tell() == 0:
                writer.writerow(['Experiment Info', 'Values'])
                writer.writerow(['Algorithm', self.algorithm.__name__])
                writer.writerow(['Protein', protein])
                writer.writerow(['Runtime (s)', f"{runtime:.2f}"])
                writer.writerow(['Total Solutions', num_solutions])
                writer.writerow(['Best Score', best_score])
                writer.writerow(['Average Score', f"{sum(all_scores) / len(all_scores):.2f}"])
                writer.writerow([])
                writer.writerow(['Score Index', 'Score'])
        
            # Append only new results
            for i, score in enumerate(all_scores[-5:]):  # Save last 5 results periodically
                writer.writerow([i, score])
