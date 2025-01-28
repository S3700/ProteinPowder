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

        while True:
            current_time = time.time()
            elapsed_time = current_time - start_time

            # Stop immediately if runtime has exceeded
            if elapsed_time >= self.runtime:
                break  

            # Run the algorithm
            folding, score, scores_from_run = solver.find_best_solution()

            # Ensure all scores are collected
            all_scores.extend(scores_from_run)

            # Update best score if applicable
            if score < best_score:
                best_score = score
                best_folding = folding

        # Save results to CSV
        self.save_results(all_scores, best_score, len(all_scores), elapsed_time, protein)

        return best_folding, best_score, all_scores


    def save_results(self, all_scores, best_score, num_solutions, runtime, protein):
        """
        Save experiment results to a CSV file.
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{self.algorithm.__name__}_{timestamp}.csv"
        filepath = os.path.join(self.output_dir, filename)

        # Write full results to ensure no data loss
        with open(filepath, 'w', newline='') as f:  # 'w' to overwrite file each time
            writer = csv.writer(f)

            # Write experiment metadata
            writer.writerow(['Experiment Info', 'Values'])
            writer.writerow(['Algorithm', self.algorithm.__name__])
            writer.writerow(['Protein', protein])
            writer.writerow(['Runtime (s)', f"{runtime:.2f}"])
            writer.writerow(['Total Solutions', num_solutions])
            writer.writerow(['Best Score', best_score])
            writer.writerow(['Average Score', f"{sum(all_scores) / len(all_scores):.2f}"])
            writer.writerow([])
            writer.writerow(['Score Index', 'Score'])

            # Write all recorded scores
            for i, score in enumerate(all_scores):  # Save **all** scores
                writer.writerow([i, score])

