from code.algorithms.hillclimber import HillClimber  # Change this based on the algorithm you want to use
from code.visualisation.visualise import print_visual

def main():
    proteins = ["HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"]  # Protein sequences to fold HCPHHPCH
    
    for protein in proteins:
        random_solver = HillClimber(protein, dimension=3)
        
        # Find the best folding and its corresponding score
        best_folding, best_score = random_solver.find_best_solution()
        
        # Check if we got a valid folding
        if best_folding is None:
            print(f"No valid folding found for {protein}.")
            continue
        
        # Output the best folding and its score
        print(f"Best folding for {protein}:")
        print(f"Folding: {best_folding}")
        print(f"Score: {best_score}")
        
        # Visualize the result using print_visual
        print_visual(protein, best_folding, random_solver.graph)

        print("\n")

if __name__ == "__main__":
    main()

