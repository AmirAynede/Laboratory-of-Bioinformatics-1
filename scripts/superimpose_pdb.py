# Here’s a Python script that uses SVDSuperimposer from BioPython to perform structure superimposition. The script will:
#	1.	Download the PDB files for Human Cytochrome C (3ZCF:A) and Equine Cytochrome C (3020:A).
#	3.	Perform superimposition using Singular Value Decomposition (SVD).
#	4.	Calculate and print the RMSD (Root Mean Square Deviation).

# This script will download the PDB files, extract the Cα atoms, perform superimposition, 
# and output the RMSD, rotation matrix, and translation vector.

from Bio import PDB
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np
import requests

# Function to download PDB file
def download_pdb(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    with open(f"{pdb_id}.pdb", "wb") as file:
        file.write(response.content)

# Function to extract C-alpha atom coordinates
def get_ca_atoms(pdb_id, chain_id):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, f"{pdb_id}.pdb")

    # Check available chains
    available_chains = [chain.id for chain in structure[0]]
    print(f"PDB {pdb_id} available chains: {available_chains}")

    if chain_id not in available_chains:
        raise ValueError(f"Chain {chain_id} not found in PDB {pdb_id}. Available: {available_chains}")

    return np.array([res["CA"].get_coord() for res in structure[0][chain_id] if "CA" in res])

# Correct PDB IDs
human_pdb, human_chain = "3ZCF", "A"   # Human Cytochrome C
equine_pdb, equine_chain = "3O20", "A" # Equine Cytochrome C (corrected)

# Download PDB files
download_pdb(human_pdb)
download_pdb(equine_pdb)

# Extract C-alpha coordinates
human_atoms = get_ca_atoms(human_pdb, human_chain)
equine_atoms = get_ca_atoms(equine_pdb, equine_chain)

# Align based on the smallest set of atoms
min_atoms = min(len(human_atoms), len(equine_atoms))
human_atoms, equine_atoms = human_atoms[:min_atoms], equine_atoms[:min_atoms]

# Perform superimposition
sup = SVDSuperimposer()
sup.set(human_atoms, equine_atoms)
sup.run()

# Output results
print(f"RMSD: {sup.get_rms():.3f}")
print(f"Rotation Matrix:\n{sup.get_rotran()[0]}")
print(f"Translation Vector:\n{sup.get_rotran()[1]}")



# python3 superimpose file1 file2 (run in terminal)
# python3 superimpose.py 3ZCF 3O20