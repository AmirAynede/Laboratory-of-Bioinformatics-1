# Laboratory-of-Bioinformatics-1

Multiple Sequence Alignment (MSA) Entropy Analysis

This script analyzes `multiple sequence alignments (MSA)` to compute amino acid frequency profiles and Shannon entropy for each alignment column. It helps assess sequence conservation by identifying the most conserved residues and their variability.

How the Script Works

### 1. Parsing the Alignment File

The script first reads the MSA file (in FASTA format) and extracts sequences. It stores them in a dictionary where:

- The `key` is the unique sequence identifier.
- The `value` is the corresponding amino acid sequence.

### 2. Building the Frequency Profile

A frequency matrix is created for the alignment:

- Each column in the alignment corresponds to an amino acid position.
- The matrix stores the frequency of each amino acid type (ACDEFGHILMNPQRSTVWY) at every column.
- This helps determine which residues are most conserved.

### 3. Calculating Shannon Entropy

For each column:

-The probability distribution of amino acids is computed.
-The Shannon entropy formula is applied:
 
The Shannon entropy is defined as $H = -\sum p_i \log_2(p_i)$.

where p_i is the probability of an amino acid occurring at that position.

- A `low entropy` value means `high conservation`, while `high entropy` suggests `variability`.

### 4. Determining the Most Conserved Residue

The script identifies:

- The most frequent amino acid at each column.
- Its relative frequency in the column.

Input & Output:

Input:

- A `FASTA-formatted MSA file` (e.g., alignment.fasta).

Output:

```
Column  Total_Count  Entropy  Consensus_Residue  Frequency_Vector

0       5           0.92     A                  0.4   0.2   0.2   0.2
1       5           0.56     G                  0.6   0.1   0.1   0.2
```

where:

- Column: The alignment position.
- Total_Count: Number of residues (excluding gaps).
- Entropy: Shannon entropy value.
- Consensus_Residue: Most frequent amino acid.
- Frequency_Vector: Relative frequencies of all amino acids.

 Usage:

Run the script from the terminal:

``` bash
python3 parse-ali.py alignment.fasta
```
Replace `alignment.fasta` with your `MSA file`.

Applications:

- Identifying highly conserved regions in proteins.
- Detecting variable sites in protein families.
- Assessing sequence conservation in evolutionary studies.
