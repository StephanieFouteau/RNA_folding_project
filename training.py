#!/usr/bin/python3

import os
import argparse
import math
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from collections import defaultdict

def parse_pdb (pdb_file):
    """
    Parse a PDB file and extract atomic coordinates.
    Args:
        pdb_file: path to the PDB file.
    
    Returns:
        atoms: list of tuples (residue_name, chain, residue_index, x, y, z).
    """
    atoms = [] # Initialize empty list to store atoms coordinates
    
    pdb = open(pdb_file, 'r')
    for line in pdb :
        column1 = line[0:6].strip()
        atom_name = line[12:16].strip()
        if column1 == 'ATOM' and atom_name == "C3'":
            residue_name = line[17:20].strip()
            chain = line[21:22].strip()
            residue_index = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        
            atoms.append((residue_name, chain, residue_index, x, y, z))
                
    return atoms

def compute_distances(atoms, min_residue_separation):
    """
    Compute pairwise distances between C3 atoms within the same chain,
    considering only residues separated by at least min_residue_separation.

    Args:
        atoms (list of tuples): List of atoms as tuples (atom_id, x, y, z).
        min_residue_separation (int): Minimum separation between residue indices.
        
    Returns:
        pair_distances (dict): Dictionary with base pair keys and list of distances.
        all_distances (list): List of distances between all atom pairs.
    """
    pair_distances = defaultdict(list)
    all_distances = []
    
    for (atom1, atom2) in combinations(atoms, 2):
        #(residue_name, chain, residue_index, x, y, z)
        residue_name1 = atom1[0]
        chain1 = atom1[1]
        residue_index1 = atom1[2]
        x1 = atom1[3]
        y1 = atom1[4]
        z1 = atom1[5]
        residue_name2 = atom2[0]
        chain2 = atom2[1]
        residue_index2 = atom2[2]
        x2 = atom2[3]
        y2 = atom2[4]
        z2 = atom2[5]
       
        if chain1 == chain2 and abs(residue_index1 - residue_index2) >= min_residue_separation:
            base_pair = ''.join(sorted([residue_name1, residue_name2]))
            distance = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
            pair_distances[base_pair].append(distance)
            all_distances.append(distance)
            
    return pair_distances, all_distances

def plot_distances_distribution(pair_distances, all_distances):
    """
    Plot histograms of distance distributions for each base pair.

    Args:
        pair_distances (dict): Dictionary with base pair keys and list of distances.
        all_distances (list): List of distances between all atom pairs.
    """
    base_pairs = pair_distances.keys()
    plt.figure(figsize=(15, 10))
    
    for i, bp in enumerate(base_pairs):
        plt.subplot(4, 3, i + 1)
        plt.hist(pair_distances[bp], bins=20, color="plum", edgecolor="black")
        plt.title(f"{bp} Distances Distribution")
        plt.xlabel("Distance (Å)")
        
    plt.subplot(4, 3, 11)
    plt.hist(all_distances, bins=20, color="plum", edgecolor="black")
    plt.title("All Distances Distribution")
    
    plt.tight_layout()
    plt.show()
    
def bin_distances(distances, bins):
    """
    Count distances within predefined intervals (or "bins") specified by bins.
    
    Args:
        distances (list of float): List of distances.
        bins (list of float): predefined intervals
    Returns:
        counts (list of int): distances counts for each bin.
    """
    counts = [0] * (len(bins) - 1)
    for distance in distances:
        for i in range(1, len(bins)):
            if bins[i - 1] <= distance < bins[i]:
                counts[i - 1] += 1
                break
    return counts

def compute_frequencies_and_scores(pair_distances, all_distances, bins):
    """
    Compute observed frequencies, reference frequencies, and log-ratio scores.
    Args:
        pair_distances (dict): Dictionary with base pair keys and list of distances.
        all_distances (list): List of distances between all atom pairs.
        bins (list of float): predefined intervals.

    Returns:
        scores (dict): Dictionary with base pair keys and list of scores (estimated Gibbs free energy).

    """
    ref_counts = bin_distances(all_distances, bins)
    total_ref_count = sum(ref_counts)
    ref_freqs = [count / total_ref_count for count in ref_counts]

    scores = {}
    for pair, distances in pair_distances.items():
        pair_counts = bin_distances(distances, bins)
        total_pair_count = sum(pair_counts)
        if total_pair_count == 0:
            pair_freqs = [0] * len(ref_freqs)
        else:
            pair_freqs = [count / total_pair_count for count in pair_counts]

        pair_scores = []
        for obs, ref in zip(pair_freqs, ref_freqs):
            if ref > 0:  # Avoid division by zero
                score = -math.log((obs + 1e-10) / (ref + 1e-10))  # Add small value to avoid log(0)
                pair_scores.append(min(score, 10))  # Cap the score at 10
            else:
                pair_scores.append(10)  # Maximum score if reference frequency is zero
        scores[pair] = pair_scores

    return scores

def save_scores_to_files(scores, bins, output_prefix="scores"):
    """
    Save the computed scores to text files, one for each base pair.
    """
    for pair, scores_list in scores.items():
        with open(f"{output_prefix}_{pair}.txt", "w") as f:
            for i, score in enumerate(scores_list):
                f.write(f"Bin {i + 1} ({bins[i]}-{bins[i + 1]} Å): {score:.6f}\n")

def main(pdb_dir, min_residue_separation=3, bins=None, output_dir="scores"):
    """
    Main function to process multiple PDB files and compute scores.
    """
    if bins is None:
        bins = [i for i in range(21)]  # Default: 20 bins from 0 to 20 Å

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over all PDB files in the input directory
    pdb_files = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if f.endswith(".pdb")]

    for pdb_file in pdb_files:
        try:
            # Parse PDB file
            print(f"Processing PDB file: {pdb_file}")
            atoms = parse_pdb(pdb_file)

            # Compute distances
            pair_distances, all_distances = compute_distances(atoms, min_residue_separation)

            # Plot distances distribution
            plot_distances_distribution(pair_distances, all_distances)

            # Compute frequencies and log-ratio scores
            scores = compute_frequencies_and_scores(pair_distances, all_distances, bins)

            # Save scores to files
            file_name = os.path.splitext(os.path.basename(pdb_file))[0]
            output_prefix = os.path.join(output_dir, file_name)
            save_scores_to_files(scores, bins, output_prefix)

            print(f"Scores for {pdb_file} saved in {output_dir}")


        except Exception as e:
            print(f"Error processing {pdb_file}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute distances, frequencies and estimated Gibbs free energy for a dataset of known 3D structures')
    parser.add_argument('-pdb', required=True, help='Path to the directory containing PDB files')
    parser.add_argument('-out', type=str, default='scores', help='Output directory name')

    args = parser.parse_args()

    # Validate arguments
    pdb_dir = args.pdb
    output_dir = args.out
    
    if os.path.exists(args.out):
        print(f"Error: Output directory {args.out} already exists. Please choose another name.")
        exit(1)
    else:
        os.mkdir(args.out)

    bins = [i for i in range(21)]  # Default: 20 bins from 0 to 20 Å
    min_residue_separation = 3 # Minimum sequence separation
    
    # Iterate over all PDB files in the input directory
    pdb_files = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if f.endswith(".pdb")]

    for pdb_file in pdb_files:
        try:
            # Parse PDB file
            print(f"Processing PDB file: {pdb_file}")
            atoms = parse_pdb(pdb_file)

            # Compute distances
            pair_distances, all_distances = compute_distances(atoms, min_residue_separation)

            # Plot distances distribution
            plot_distances_distribution(pair_distances, all_distances)

            # Compute frequencies and log-ratio scores
            scores = compute_frequencies_and_scores(pair_distances, all_distances, bins)

            # Save scores to files
            file_name = os.path.splitext(os.path.basename(pdb_file))[0]
            output_prefix = os.path.join(output_dir, file_name)
            save_scores_to_files(scores, bins, output_prefix)

            print(f"Scores for {pdb_file} saved in {output_dir}")


        except Exception as e:
            print(f"Error processing {pdb_file}: {e}")















