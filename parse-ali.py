import sys
import numpy as np

def parse_alignment(alignment_file, id_position=1, separator='|'):
    """
    Parses a FASTA multiple sequence alignment file.
    Extracts sequences and stores them in a dictionary with a unique ID.
    """
    alignment_dict = {}  # Dictionary to store sequences with unique identifiers
    sequence_id = ''  # Variable to hold the current sequence ID
    
    with open(alignment_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                # Extract sequence identifier
                sequence_id = line.split(separator)[id_position].strip()
                alignment_dict[sequence_id] = ''
            else:
                # Append sequence data to the corresponding ID
                if sequence_id:
                    alignment_dict[sequence_id] += line.strip()
    
    return alignment_dict

def calculate_frequency_profile(alignment_dict, amino_acids='ACDEFGHILMNPQRSTVWY'):
    """
    Computes the frequency profile for each column in the alignment.
    Also calculates the entropy of each column.
    """
    num_amino_acids = len(amino_acids)  # Number of amino acid types considered
    sequences = list(alignment_dict.values())  # Extract sequences from dictionary
    if not sequences:
        return
    
    num_columns = len(sequences[0])  # Assume all sequences have the same length
    profile = np.zeros((num_amino_acids, num_columns))  # Initialize frequency matrix
    
    # Populate the frequency matrix
    for sequence in sequences:
        for column_index in range(num_columns):
            amino_acid = sequence[column_index]  # Get amino acid at current position
            aa_index = amino_acids.find(amino_acid)  # Find index in reference list
            if aa_index != -1:
                profile[aa_index, column_index] += 1.0  # Increment frequency count
    
    # Normalize the profile and compute entropy
    for column_index in range(num_columns):
        column_sum = profile[:, column_index].sum()  # Total occurrences in column
        if column_sum > 0:
            profile[:, column_index] /= column_sum  # Convert counts to probabilities
        entropy = calculate_entropy(profile[:, column_index])  # Compute entropy
        
        # Print column statistics: index, total frequency, entropy, reference residue, and frequencies
        frequency_vector = '\t'.join([str(freq) for freq in profile[:, column_index]])
        print(column_index, column_sum, entropy, sequences[0][column_index], frequency_vector)
    
    return profile

def calculate_entropy(frequencies):
    """
    Computes the Shannon entropy of a given probability distribution.
    """
    entropy = 0.0  # Initialize entropy value
    for frequency in frequencies:
        if frequency > 0.0:
            entropy -= frequency * np.log2(frequency)  # Shannon entropy formula
    return entropy

if __name__ == '__main__':
    alignment_file = sys.argv[1]  # Get input alignment file from command line argument
    alignment_data = parse_alignment(alignment_file)  # Parse alignment file
    print(alignment_data)  # Print parsed sequences
    calculate_frequency_profile(alignment_data)  # Compute and display frequency profile
