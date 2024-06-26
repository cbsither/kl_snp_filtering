# Python Script for KL Divergence Computation from Fasta Files

## Description
This script processes a directory of FASTA files to compute the Kullback-Leibler (KL) divergence between two Dirichlet distributions for each site in a sequence alignment. It merges the alignment with reference data from a CSV file, performs computations, and saves the results to a specified output file.

## Dependencies
To run this script, the following Python libraries are required:
- `pandas`: For data manipulation and analysis.
- `numpy`: For numerical operations.
- `Bio`: Part of Biopython, used for parsing FASTA files.
- `scipy`: Used for special functions like gamma and psi functions.
- `argparse`: For parsing command-line options and arguments.
- `os`: For interacting with the operating system, like listing files in a directory.

Make sure you have Python installed on your system (Python 3.x is recommended). You can install these dependencies using pip. For example:
```bash
pip install pandas numpy biopython scipy argparse
```

## Usage
To use this script, you need to have a directory containing FASTA files, a reference CSV file with labcodes and prediction y values, and specify an output file path and a target column.

Run the script from the command line as follows:
```bash
python [script_name].py <fasta_dir> <reference_csv> <output_file> <fasta_code> <target>
```
Replace `[script_name]` with the name of the script file.

### Arguments:
- `fasta_dir`: Directory containing FASTA files.
- `reference_csv`: CSV file with labcodes and prediction y values.
- `output_file`: Path where the output CSV file will be saved.
- `fasta_code`: Link column containing the name of the fasta sequence in the .fasta files
- `target`: Target column in the reference CSV file.

### Reference CSV Layout

The reference CSV should contain at minimum the following:

A 'fasta_code' column:
  -This column is verbatim the name of the fasta sequence in the .fasta files

A 'target' column:
  -This column contains 1 (target) and 0 (non-target) labels identifying whether the fasta sequence is your target or non-target

### Example
```bash
python script.py /path/to/fasta/dir reference_data.csv output_results.csv fasta_link_column target_column
```

## Output
The script will output a CSV file containing the alignment, site, KL score, and frequency counts of nucleotides for target and non-target sequences.
