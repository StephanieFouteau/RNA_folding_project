*Project developped for the GENIOMHE master degree.*  

This project contains scripts for train and use an objective function to estimate\
the Gibbs free energy which can be used to evaluate a predicted RNA structures.   

## Description

1. `training.py`: Calculate the estimated Gibbs free energy for several RNA structures in PDB format.\
-Parse a PDB file and extract atomic coordinates\
-Compute pairwise distances between C3 atoms within the same chain\
-Bins distances into intervals from 0 to 20 Å\
-Compute observed frequencies, reference frequencies, and log-ratio scores    
2. `plot_scores.py`: Plot the scoring profiles, the score (or estimated Gibbs free energy) as a function of the interatomic distance.  
3. `scoring.py`: Perform linear interpolation to compute the score for a specific base pair at a given distance and calculate the estimated Gibbs free energy of the evaluated RNA conformation.  

`PDB_Files`: Folder containing the datas.\
`scores`: Example of `training.py` output.\
`plots`: Example of `plot_scores.py` output.  


---

## Installation

#### 1. Clone the repository from Github:
```
git clone https://github.com/StephanieFouteau/RNA_folding_project.git
```

#### 2. Navigate to the directory where the repository was cloned:
```
cd RNA_folding_project
```
---

## Usage

### training.py

Compute distances, frequencies and estimated Gibbs free energy for a dataset of known 3D structures.

```
training.py [-h] -pdb PDB [-out OUT]
```
Options\
  `-h`, `--help`  show this help message and exit\
  `-pdb` PDB    Path to the directory containing PDB files\
  `-out` OUT    Output directory name

### plot_scores.py

Plot the scoring profile from the pairwise score file.

```
plot_scores.py [-h] -scores SCORES
```
options:\
  `-h`, `--help`      show this help message and exit\
  `-scores` SCORES  Give the path of the folder containing the energy file for the 10 base pairs

### scoring.py

```
scoring.py [-h] -pdb PDB -scores SCORES
```
options:\
  `-h`, `--help`      show this help message and exit\
  `-pdb` PDB        Path to the PDB file\
  `-scores` SCORES  Path to the scores folder

## Outputs

`scores`: Scores generated by training.py.\
`plots`: Plots generated by plot_scores.py.
