import argparse
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from scipy.special import gammaln, psi
import warnings
from datetime import datetime

# Ignore FutureWarnings
warnings.filterwarnings('ignore', category=FutureWarning)


def kl_divergence_dirichlet(alpha, beta):
    """
    Compute the KL divergence between two Dirichlet distributions parameterized by alpha and beta.
    """
    sum_alpha = np.sum(alpha)
    sum_beta = np.sum(beta)
    kl = gammaln(sum_alpha) - gammaln(sum_beta) - np.sum(gammaln(alpha)) + np.sum(gammaln(beta))
    kl += np.sum((alpha - beta) * (psi(alpha) - psi(sum_alpha)))
    return kl

def main(fasta_dir, reference_csv, output_file, fasta_code, target):
    # Read reference data
    code_ref = pd.read_csv(reference_csv)

    # Set target species
    target_sp = target
    code_ref['y_'] = code_ref[target_sp].copy()

    # Get fasta files
    fasta_files = [f.split('.')[0] for f in os.listdir(fasta_dir) if f.endswith('.fasta')]
    fasta_files = list(set(fasta_files))

    nuc_key = {'A': 0, 'T': 1, 'G': 2, 'C': 3, '-': 4, 'N': 5}

    # Initialize results DataFrame
    synapo_results = pd.DataFrame(columns=['alignment', 'site', 'KL_Score',
                                           'target_A', 'target_T', 'target_G', 'target_C', 'target_Gap', 'target_N',
                                           'non_target_A', 'non_target_T', 'non_target_G', 'non_target_C', 'non_target_Gap', 'non_target_N'])

    # Process each fasta file
    for gg, i in enumerate(fasta_files):

        print(f'Begin processing: {i} @ {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

        try:
            ex_fasta = SeqIO.parse(os.path.join(fasta_dir, f'{i}.ntogap.nt.fasta'), 'fasta')

            # Create a dictionary to hold all sequences
            data_dict = {seq.id: list(seq.seq) for seq in ex_fasta}

            # Convert the dictionary to a dataframe
            ex_df = pd.DataFrame.from_dict(data_dict, orient='index')
            merged_df = ex_df.merge(code_ref[[fasta_code, 'y_']], left_index=True, right_on=fasta_code, how='left')
            merged_df = merged_df[~pd.isna(merged_df['y_'])].drop_duplicates(fasta_code)
            ex_df = merged_df.rename(columns={'y_': 'y_pred'})

            kl_div = np.zeros(len(ex_df.columns) - 2)

            # Compute KL divergence for each site
            for dd, j in enumerate(ex_df.columns[:-2]):
                non_target = np.ones(6)
                target = np.ones(6)

                for ii, base in ex_df[j].items():
                    y_val = int(ex_df.at[ii, 'y_pred'])
                    index = nuc_key[base]
                    if y_val == 1:
                        target[index] += 1
                    else:
                        non_target[index] += 1

                kl_div[dd] = kl_divergence_dirichlet(alpha=non_target[:4], beta=target[:4])

                if np.sum(target[:4] > 1) == 1 and non_target[np.argmax(target[:4])] == 1:
                    data_list = [[i, dd, kl_div[dd]] + target.tolist() + non_target.tolist()]
                    add_df = pd.DataFrame(data_list, columns=synapo_results.columns)
                    synapo_results = pd.concat([synapo_results, add_df], ignore_index=True)

            print(f'Finished processing: {i} @ {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

        except Exception as e:
            print(f"Error processing file {i}: {e}")

    # Save results
    synapo_results.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process fasta files and compute KL divergence.')
    parser.add_argument('fasta_dir', type=str, help='Directory containing fasta files')
    parser.add_argument('reference_csv', type=str, help='CSV file with labcodes and prediction y values')
    parser.add_argument('output_file', type=str, help='Output file path')
    parser.add_argument('fasta_code', type=str, help='fasta file link column')
    parser.add_argument('target', type=str, help='target column')

    args = parser.parse_args()
    main(args.fasta_dir, args.reference_csv, args.output_file, args.fasta_code, args.target)
