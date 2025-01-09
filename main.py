from code.algorithms.bruteforce import brute_force_folding
from code.visualisation.visualise import print_visual

def main():
    proteins = ["HHPPHPP"]
    num_random_folds = 100
    
    for protein in proteins:
        best_folding, best_score = brute_force_folding(protein, dimension=3)
        print(f"Best folding for {protein}:")
        print(f"Folding: {best_folding}")
        print(f"Score: {best_score}")
        print_visual(protein, best_folding)
        print("\n")

if __name__ == "__main__":
    main()
