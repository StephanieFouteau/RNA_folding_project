#!/usr/bin/python3

import os
import argparse
import matplotlib.pyplot as plt

def plot_scoring_profile(scores_file):
    """
    Plot the scoring profile from the pairwise score file.
    
    Args:
        scores_file (str): file containing pairwise scores (one score per line for each distance bin).
    """
    # Step 1: Read the scores from the file
    scores = []
    with open('scores/' + scores_file, 'r') as file:
        for line in file:
            # Find the score after the colon and strip whitespace
            if ':' in line:
                score = line.split(':')[-1].strip()  # Get text after ':'
                scores.append(float(score))  # Add score to the list
    
    # Step 2: Define the distance bins (0 to 20 Å, 1 Å intervals)
    distance_bins = list(range(len(scores)))  # Each score corresponds to a bin (0-19 Å, 20 bins)
    
    # Step 3: Plot the scores as a function of distance
    plt.figure(figsize=(10, 6))  # Set figure size
    plt.plot(distance_bins, scores, marker='o', linestyle='-', color='darkblue', label='estimated Gibbs free energy')
    
    # Step 4: Add labels, title, and legend
    plt.xlabel("Distance (Å)", fontsize=14)
    plt.ylabel("Score", fontsize=14)
    plt.title("Scoring Profile", fontsize=16)
    plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)  # Add a reference line at score = 0
    plt.legend(fontsize=12)
    
    # Step 5: Save the plot to a file
    plt.savefig('plots/' + str(scores_file) + '.png')
    plt.close()  # Close the figure to avoid overlap in subsequent plots

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot the scoring profile from the pairwise score file')
    parser.add_argument('-scores', required=True, help="Give the path of the folder containing the energy file for the 10 base pairs.")
    
    args = parser.parse_args()
    
    # Validate arguments
    scores = args.scores
    
    # Create the 'plot' folder if it does not exist
    try:
        os.mkdir('plots')
    except OSError:
        pass

    for scores_file in os.listdir(scores):
        plot_scoring_profile(scores_file)
