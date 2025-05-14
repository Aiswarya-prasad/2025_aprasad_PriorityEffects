
#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
from glob import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools
import subprocess
import argparse

'''
Accept arguments:
    -s: path to the ASV table - comes from dada2 seq table after chimera removal
    -o: path to the output file (also used as pre-fix for other accompanying output files)
    -l: path to the log file
    -a: path to files with amplicon sequences from isolate genomes to compare against e.g. 03_Analysis/SeqeuncesDB/OnlyIsolateSeqs_amplicons_updated_presence.csv
    -c: path to the file with asv seq vs number of copies expected per strain it contains uids to identify amplicons and the strain names on the column names e.g. 03_Analysis/SeqeuncesDB/OnlyIsolateSeqs_updated_copies.csv
    -i: path to file containing uids, strain name associated and amplicon sequence e.g. 03_Analysis/SeqeuncesDB/OnlyIsolateSeqs_updated.csv

get ASVs from the seqtable provided and infer the species table!

check ASVs against eg. 03_Analysis/SeqeuncesDB/OnlyIsolateSeqs_amplicons_updated_presence.csv if asked for and log results

get C from eg. 03_Analysis/SeqeuncesDB/OnlyIsolateSeqs_updated_copies_amplicons.csv
    
    A    =   C   .   S
(nx1) = (nxm) . (mx1)

    S    =   C^-1    .    A
(mx1)    =  (mXn)    .   (nx1)

A: For each sample, a matrix A can be obtained from the seqtable
where A is all ASVs counts for that sample. It is a column matrix
with ASV counts in the rows. n rows if there are n ASVs. nx1 matrix

C: A matrix C is obtained from the output of 03_Analysis/summarise_genomes_16S.py
it contains the number of copies of each ASV in each strain. It needs to have no
indistinguishable strains with identical ASVs (underdetermined) as many ASVs
as in the seqtable. It is a rectangular matrix with n rows and m columns where
m is the number of strains. n x m matrix

S: It is the matrix we want to solve for. It is a column matrix with the count of strains
for the sample in question. It is a column matrix with m rows and 1 column. Where
m is the number of strains. m x 1 matrix

Now, if C were a square matrix, we could just invert it and multiply it with A to get S.
But C is not a square matrix. So we need to use the least squares method to solve for S.
We can use the numpy function np.linalg.lstsq to solve for S.

What it does is, it finds the S that minimizes the sum of the squares of the differences
between the observed values in A and the values predicted by the model C*S.

Unless we have exactly 1 unique ASV in each strain, C will not be a square matrix. So, we
cannot invert it in the usual sense. Instead, we can use the least squares method
to find the best approximation:

S = argmin_S ||CS - A||^2
where ||.||^2 is the squared L2 norm.
The solution S will be the one that minimizes the sum of the squares of the differences
between the observed values in A and the values predicted by the model C*S.

We need to do this because different strains have different number of ASVs often more than one.
If it have 4 identical copies we trivially divide it by 4 and get the count of the strain. However,
if it has 3 identical and 1 other, we expect that we divide the one that is present in 3 copies
by 3 and this will be equal to the number of copies of the other ASV and also the strain. However,
in practice, this is not the case. The multi-copy ASV is not always present in the expected ratio exactly..
Also consider cases where there are 4 unique ASVs in the strain but their counts are not equal. So, we can instead
use the least squares method to find the best approximation of the strain counts that explains the observed ASV counts.


'''
args = argparse.ArgumentParser(description='Get species table from ASV table and barrnap files')
args.add_argument('-s', type=str, help='Path to the ASV table')
args.add_argument('-o', type=str, help='Path to the output file')
# optional argument for logfile path
args.add_argument('-l', type=str, help='Path to the log file')
args.add_argument('-a', type=str, help='Path to the amplicon sequences information file from sequencing isolate genome amplicons')
args.add_argument('-c', type=str, help='Path to the file with expected copies per strain')
args.add_argument('-i', type=str, help='Path to the file with uid, strain name and sequence for copies file')
args = args.parse_args()

# seqtable_original_path = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/IsolateAmplicons/seqtable_nochim.csv'
# output_path = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/IsolateAmplicons/species_table_updated.csv'
seqtable_original_path = args.s
seqtable_original = pd.read_csv(seqtable_original_path, header = 0, index_col = 0)
output_path = args.o
if args.l:
    log_path = args.l
else:
    log_path = output_path.replace('.csv', '.log')
asv_table_path = output_path.replace('.csv', '_all_asvs.csv')
my_fasta_path = output_path.replace('.csv', '_unmatched_asv_sequences.fasta')

if args.a:
    amplicon_path = args.a
else:
    amplicon_path = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/SeqeuncesDB/OnlyIsolateSeqs_amplicons_updated_presence.csv'
if args.c:
    copies_path = args.c
else:
    copies_path = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/SeqeuncesDB/OnlyIsolateSeqs_updated_copies.csv'
if args.i:
    uid_strain_path = args.i
else:
    uid_strain_path = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/SeqeuncesDB/OnlyIsolateSeqs_updated.csv'


def seq_identical(seq1, seq2):
    if seq1 in seq2 or seq1 in seq2.reverse_complement():
        return True
    if seq2 in seq1 or seq2 in seq1.reverse_complement():
        return True
    else:
        return False


with open(log_path, 'w') as log:
    success = log.write(f'Log file created on {pd.Timestamp.now()}\n')
    success = log.write(f'Input file: {seqtable_original_path}\n')
    success = log.write(f'          contains rows {seqtable_original.index}\n')
    if os.path.exists(amplicon_path):
        success = log.write(f'          Using amplicon sequences from {amplicon_path}\n')
    else:
        success = log.write(f'          No amplicon sequences provided, using only ASV table\n')
    if os.path.exists(uid_strain_path):
        success = log.write(f'Using uid strain information from {uid_strain_path}\n')
    else:
        success = log.write(f'No uid strain information provided, I don\'t know what strains ASV uids correspond to\n')
        sys.exit(1)
    if os.path.exists(copies_path):
        success = log.write(f'Using copies information from {copies_path}\n')
    else:
        success = log.write(f'No copies information provided, using only ASV table\n')
        sys.exit(1)
    success = log.write(f'Output file: {output_path}\n')
    success = log.write(f'Also writing:\n')
    success = log.write(f'          Fasta sequences unique identity for this samples in header: {my_fasta_path}\n')
    success = log.write(f'          ASV table with all ASVs: {asv_table_path}\n')
    success = log.write(f'\n')
    success = log.write(f'\n')
    success = log.write(f'starting to process the ASV table ...\n')
    success = log.write(f'\n')
    
    # get martrix C
    C = pd.read_csv(copies_path)
    C = C.dropna(axis=1, how='all')
    C.set_index('uid', inplace=True)
    # info about known sequences
    df_seq_info = pd.read_csv(uid_strain_path).set_index('uid')
    uid_seq_dict = df_seq_info['seq'].to_dict()
    uid_seq_dict = {int(k): Seq(v) for k, v in uid_seq_dict.items()}
    uid_strain_dict = df_seq_info['strain'].to_dict()

    # clean seqtable - gives As
    seqtable = seqtable_original.copy()
    # seqtable.index = [x.split('_filt_reads.fastq.gz')[0] for x in seqtable_original.index]
    # seqtable.index = [f'{x}_isolate' for x in seqtable.index]
    
    # now, make a version of seqtable with only the ASVs that are in the C matrix index
    # this is the one that will be used to get the species table

    unmatched_asvs = {}
    unk = 0
    asv_matched = False
    seqtable_mod = pd.DataFrame(0, index=seqtable.index, columns=C.index)
    asv_table = pd.DataFrame(0, index=seqtable.index, columns=C.index)
    for column in seqtable.columns:
        asv_seq_table = Seq(column)
        # each asv should only match one uid as uid sequeces are unique
        for uid in uid_seq_dict.keys():
            if seq_identical(asv_seq_table, uid_seq_dict[uid]):
                # get the uid of the ASV
                seqtable_mod[uid] = seqtable[column]
                asv_table[uid] = seqtable[column]
                asv_matched = True
                break
        if not asv_matched:
            unmatched_asvs[f'unk_{unk}'] = asv_seq_table
            asv_table[f'unk_{unk}'] = seqtable[column]
            unk += 1

        asv_matched = False

    seqtable_mod['unknown'] = seqtable.sum(axis=1) - seqtable_mod.sum(axis=1)
    asv_table['total'] = seqtable.sum(axis=1)
    asv_table.columns = [f'{uid_strain_dict[x]}_{x}' if x in uid_strain_dict else x for x in asv_table.columns]
    
    # remove columns with all zeros
    asv_table = asv_table.loc[:, (asv_table != 0).any(axis=0)]
    asv_table.to_csv(asv_table_path, index=True)

    success = log.write(f'ASV table written {asv_table_path}\n')

    # final output is sepcies table
    species_table_columns = [x for x in C.columns]
    species_table_columns.append('unknown')
    species_table =  pd.DataFrame(0, index=seqtable.index, columns=species_table_columns)

    for sample in seqtable.index:
        A = seqtable_mod.loc[sample].values[:-1]
        S = np.linalg.lstsq(C, A, rcond=None)[0]
        S = np.round(S)
        S[S < 0] = 0
        S = np.append(S, seqtable_mod.loc[sample]['unknown'])
        species_table.loc[sample] = S
    
    species_table.to_csv(output_path, index=True)

    success = log.write(f'Species table written {output_path}\n')
    
    # screen unmatched asvs against the amplicon sequences and report!
    log.write(f'Unmatched ASVs: {len(unmatched_asvs)} / {len(seqtable.columns)}\n')
    
    df_amplicon_data = pd.read_csv(amplicon_path)
    for asv_id, asv_seq in unmatched_asvs.items():
        asv_found = False
        for amplicon_seq in df_amplicon_data['seq']:
            if seq_identical(asv_seq, Seq(amplicon_seq)):
                asv_found = True
                associated_strain = df_amplicon_data.loc[df_amplicon_data['seq'] == amplicon_seq, 'strain'].values[0]
                count = df_amplicon_data.loc[df_amplicon_data['seq'] == amplicon_seq, 'count'].values[0]
                total = df_amplicon_data.loc[df_amplicon_data['seq'] == amplicon_seq, 'total'].values[0]
                perc = round(count/total* 100, 2)
                success = log.write(f'{asv_id} found in amplicon sequences in {associated_strain} at {perc}% of {total} reads\n')

    with open(my_fasta_path, 'w') as fasta_file:
        for asv_id, asv_seq in unmatched_asvs.items():
            # full sequence in one line
            record = SeqRecord(asv_seq, id=asv_id, description='')
            success = fasta_file.write(f'>{record.id}\n')
            success = fasta_file.write(f'{record.seq}\n')
        success = log.write(f'Unmatched ASVs written to {my_fasta_path}\n')

    success = log.write(f'Finished writing unmatched ASVs\n')
# not writing fasta anymore..