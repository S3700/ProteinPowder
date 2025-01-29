import time
import csv
import os
from datetime import datetime

class TimedExperiment:
    def __init__(self, algorithm, max_runtime):
        """
        Initialize TimedExperiment with specific algorithm and runtime.
        Runs the algorithm for the given runtime and stores results.
        """
        self.algorithm = algorithm
        self.runtime = max_runtime  # Total allowed runtime in seconds
        self.output_dir = "experiment_results"
        
        # Create output directory if it doesn't exist
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def run(self, protein):
        """
        Run the experiment for the specified runtime and collect scores.
        Only save results that were completed within the runtime limit.
        """
        start_time = time.time()
        solver = self.algorithm(protein)  # Initialize algorithm
        final_runtime = 0 
        best_folding = None
        best_score = float('inf')
        all_scores = []
        last_update_time = start_time
        update_interval = 10

        while True:
            current_time = time.time()
            elapsed_time = current_time - start_time

            # Stop the experiment immediately if runtime has exceeded
            if elapsed_time >= self.runtime:
                break  

            # Update progress every folding or interval in seconds
            if current_time - last_update_time >= update_interval:
                progress_percentage = (elapsed_time / self.runtime) * 100
                print(f"Progress: {progress_percentage:.2f}% complete", end='\r')
                last_update_time = current_time  # Reset the last update time

            # Run the algorithm and get folding, score, and scores from the current run
            folding, score, scores_from_run = solver.find_best_solution()

            # Calculate elapsed time after each run
            elapsed_time = time.time() - start_time

            # Check if the solution was found within the runtime limit
            if elapsed_time < self.runtime:
                # Only collect scores within the time limit
                all_scores.extend(scores_from_run)
                final_runtime = elapsed_time

                # Update best score if applicable and within time limit
                if score < best_score:
                    best_score = score
                    best_folding = folding
            else:
                break  # Exit the loop if algorithm exceeds the runtime

        # Save results to CSV, ensuring only valid scores are saved
        self.save_results(all_scores, best_score, len(all_scores), final_runtime, protein)

        return best_folding, best_score, all_scores

    def save_results(self, all_scores, best_score, num_solutions, final_runtime, protein):
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
            writer.writerow(['Final Runtime (s)', f"{final_runtime:.2f}"])
            writer.writerow(['Total Solutions', num_solutions])
            writer.writerow(['Best Score', best_score])
            writer.writerow(['Average Score', f"{sum(all_scores) / len(all_scores):.2f}"])
            writer.writerow([])
            writer.writerow(['Score Index', 'Score'])

            # Write all recorded scores (that fit within the runtime)
            for i, score in enumerate(all_scores):
                writer.writerow([i, score])

