from code.algorithms.bruteforce import Bruteforce
from code.algorithms.random import RandomSolution
from code.algorithms.hillclimber import HillClimber
from code.algorithms.simannealing import SimulatedAnnealing
from code.algorithms.breadthfirst import BreadthFirstSearch
from code.algorithms.depthfirst2 import DepthFirst
from code.visualisation.visualise import print_visual
from code.classes.experiment import TimedExperiment
import matplotlib.pyplot as plt


def main():
    proteins = ["HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"]  # Input protein sequence
    algorithm =  SimulatedAnnealing # Input algorithm
    runtime = 10 # Input runtime

    # Initialize experiment handler
    experiment = TimedExperiment(algorithm, runtime)

    for protein in proteins:
        print(f"\nProcessing protein: {protein}")
        print(f"Using Algorithm: {algorithm.__name__}")

        # Run the experiment
        best_folding, best_score, all_scores = experiment.run(protein)

        # Output the best folding and its score
        print(f"Best folding found:\nFolding: {best_folding}\nScore: {best_score}")

        # Use the existing solver instead of creating a new instance
        experiment_solver = experiment.algorithm(protein)
        experiment_solver.graph.apply_folding(best_folding)
        print_visual(protein, best_folding, experiment_solver.graph)

        # Create histogram of all scores
        plt.figure(figsize=(10, 6))
        plt.hist(all_scores, bins=min(50, len(set(all_scores)) or 1))
        plt.title(f'Distribution of Scores for {protein}\nBest Score: {best_score}')
        plt.xlabel('Score')
        plt.ylabel('Frequency')
        plt.grid(True)
        plt.show()

if __name__ == "__main__":
    main()

