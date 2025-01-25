#!/usr/bin/python3

import os
import argparse
import math
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from collections import defaultdict
from training import parse_pdb, compute_distances

def linear_interpolation(pair, distance, scores_file):
    """
    Perform linear interpolation to compute the score for a specific base pair at a given distance.
    Each line is taken as the average of the interval. (e.g. line 1 corresponding to the interval [0, 1] is
    considered to be the value 0.5 in the interpolation.)
    Only values between 0.5 and 19.5 are taken into account.
    Args:
        scores_file:  Path to the file containing scores for the pair.
        pair: one base pair (e.g. 'GU')
        distance: Distance between the two bases
        
    Returns:
        Energy (score) associated with base pair and distance
        
    """
    energy_values = []
    #file_path = os.path.join(scores_file, f"{pair}.txt")
    #with open(file_path, 'r') as energy_file:
    with open(scores_file, 'r') as energy_file:
        content = energy_file.readlines()

    # Parse energy values from file
    for line in content:
    # Find the score after the colon and strip whitespace
        if ':' in line:
            score = line.split(':')[-1].strip()  # Get text after ':'
            energy_values.append(float(score))  # Add score to the list
    #energy_values = [float(line.split()[-1]) for line in content]
    
    # Ensure the distance is within valid range
    if not (0.5 <= distance <= 19.5):
        raise ValueError(f"Distance {distance} is out of interpolation range (0.5-19.5).")

    # Perform linear interpolation
    lower_index = math.floor(distance) - 1  # Indices are 0-based
    upper_index = lower_index + 1

    energy_before = energy_values[lower_index]
    energy_after = energy_values[upper_index]
    interpolated_energy = energy_before + (distance - (lower_index + 0.5)) * (energy_after - energy_before)
    return interpolated_energy

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform linear interpolation to compute the score for a specific base pair at a given distance')
    parser.add_argument('-pdb', required=True, help='Path to the PDB file')
    parser.add_argument('-scores', required=True, help='Path to the scores folder')
    
    args = parser.parse_args()

    # Validate arguments
    pdb_file = args.pdb
    scores = args.scores

    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]

    bins = [i for i in range(21)]  # Default: 20 bins from 0 to 20 Å
    min_residue_separation = 3 # Minimum sequence separation

    # Parse PDB file
    atoms = parse_pdb(pdb_file)

    # Compute distances
    pair_distances, all_distances = compute_distances(atoms, min_residue_separation)

    # Perform interpolation for each base pair and distance
    total_gibbs_energy = 0
    for pair, distances in pair_distances.items():
        scores_file = os.path.join(scores, f"{pdb_name}_{pair}.txt")
        for distance in distances:
            if 0.5 <= distance <= 19.5:
                interpolated_energy = linear_interpolation(pair, distance, scores_file)
                print(f"Interpolated energy for pair {pair} at distance {distance:.2f}: {interpolated_energy:.6f}")
            else:
                print(f"Skipping distance {distance:.2f} for pair {pair} (out of range).")
                total_gibbs_energy += interpolated_energy
    print(f"Gibbs Free Energy for structure {pdb_name} is: {total_gibbs_energy:.3f}")
