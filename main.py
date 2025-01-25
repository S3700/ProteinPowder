from code.algorithms.hillclimber import HillClimber  # Change this based on the algorithm you want to use
from code.visualisation.visualise import print_visual
import matplotlib.pyplot as plt
import time

def main():
    start_time = time.time()
    proteins = ["CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC"] # Input protein sequence
    
    for protein in proteins:
        print(f"\nProcessing protein: {protein}")
        random_solver = HillClimber(protein, dimension=3, num_valid_folds=10)
        
        # Find the best folding and get all scores
        best_folding, best_score, all_scores = random_solver.find_best_solution()
        
        # Output the best folding and its score
        print(f"Best folding found:")
        print(f"Folding: {best_folding}")
        print(f"Score: {best_score}")
        
        # Visualize the best folding
        print_visual(protein, best_folding, random_solver.graph)

        # Create histogram of all scores
        plt.figure(figsize=(10, 6))
        plt.hist(all_scores, bins=50)
        plt.title(f'Distribution of Scores for {protein}\nBest Score: {best_score}')
        plt.xlabel('Score')
        plt.ylabel('Frequency')
        plt.grid(True)
        plt.show()

        # Print some statistics
        print(f"\nStatistics:")
        print(f"Number of valid foldings found: {len(all_scores)}")
        print(f"Average score: {sum(all_scores) / len(all_scores):.2f}")
        print(f"Best score: {best_score}")
        print(f"Worst score: {max(all_scores)}")
        print(f"Runtime in seconds: {(time.time() - start_time)}")

if __name__ == "__main__":
    main()

