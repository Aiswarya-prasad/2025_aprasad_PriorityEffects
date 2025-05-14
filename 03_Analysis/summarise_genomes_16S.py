import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from glob import glob
from Bio import SeqIO
from Bio.Seq import Seq
import itertools


def seq_identical(seq1, seq2):
    if seq1 in seq2 or seq1 in seq2.reverse_complement():
        return True
    if seq2 in seq1 or seq2 in seq1.reverse_complement():
        return True
    else:
        return False

def checkmatches(a_seq, seq_list):
    '''
    checks for at lest one identical sequence in the list
    and returns True if found
    '''
    for b_seq in seq_list:
        if seq_identical(a_seq, b_seq):
            return True
    return False


'''
First, collects barrnap sequences from genomes
    00_StrainSelection/Extract_compare_16S/
also collect ASVs from Amplicon sequncing of Isolate genomes
    03_Analysis/IsolateAmplicons/asv_sequences.fasta - writted out from seqtable of the sequencing run into this location by 03_Analysis.394.misannotated/visualize_sequencing_result.py
'''


'''
Genomes from 3 sources Nanopore (by Malick), PacBio (for this experiment) and Old Genomes (Illumina from NCBI)
are found in 00_StrainSelection/Extract_compare_16S/Genomes
barrnap was run on them using 00_StrainSelection/Extract_compare_16S/identify_ribo_seqs.sh

Outputs are collected into a dictionary
barrnap_seqs = {
                Strain: [Seq(...), Seq(...)],
                ...
            }

when the same genome copies several identical copies of the same sequence, they are
all included in the list.
'''

barrnap_files = {os.path.basename(x).split('_')[0]: x for x in glob('00_StrainSelection/Extract_compare_16S/barrnap/*_all_ribo_seqs.fasta')}
barrnap_seqs = {}
for strain in barrnap_files.keys():
    recs = list(SeqIO.parse(open(barrnap_files[strain], 'r'), 'fasta'))
    # only keep those with "16S_rRNA" in the description
    recs = [rec for rec in recs if '16S_rRNA' in rec.description]
    barrnap_seqs[strain] = [rec.seq for rec in recs]

barrnap_uniq_seqs = {}
for strain_iter in barrnap_seqs.keys():
    all_seqs = barrnap_seqs[strain_iter]
    barrnap_uniq_seqs[strain_iter] = []
    all_seqs = barrnap_seqs[strain_iter]
    for i in range(len(all_seqs)):
        if checkmatches(all_seqs[i], all_seqs[:i]):
            continue
        else:
            barrnap_uniq_seqs[strain_iter].append(all_seqs[i])

# strains = barrnap_seqs.keys()
# the above works but
# we want to mantain a fixed order of the strains so,
strains_list = ['ESL0825','ESL0820',
'ESL0822','ESL0170',
'ESL0824','ESL0198',
'ESL0827','ESL0199',
'ESL0200','ESL0819',
'ESL0197',
'ESL0294','ESL0394','ESL0295',
'ESL1028',
'ESL0353','ESL0263','ESL0185',
'ESL0183','ESL0262','ESL0835','ESL0354',
'ESL0393','ESL0259','ESL0260','ESL0184','ESL0350',
'ESL0186','ESL0261','ESL0351']



'''
For each of the ASVs in the IsolateAmplicons experiment, check the strains/samples
in which they are detected if there are >10 in that sample, consider the ASV detected
    Read seqtable from that sequencing run (remove '_filt_reads.fastq.gz' from sample names) 
    The uids for ASVs can be obtained from 03_Analysis/IsolateAmplicons/asv_sequences.fasta
'''

seqtable_isolate_amplicons = pd.read_csv('03_Analysis/IsolateAmplicons/seqtable_nochim.csv', 
                                        sep=',', header=0, index_col=0)
seqtable_isolate_amplicons.index = [x.split('_filt_reads.fastq.gz')[0] for x in seqtable_isolate_amplicons.index]

# low depth if total reads < 100

low_depth_strains = []
for strain in seqtable_isolate_amplicons.index:
    if sum(seqtable_isolate_amplicons.loc[strain]) < 100:
        low_depth_strains.append(strain)

amplicon_id_dict = {}
with open('03_Analysis/IsolateAmplicons/asv_sequences.fasta', 'r') as f:
    for rec in SeqIO.parse(f, 'fasta'):
        amplicon_id_dict[rec.seq] = rec.id

# # check uniqueness of ASVs
# for i, asv in enumerate(asv_seqs):
#     for j in range(i + 1, len(asv_seqs)):
#         if seq_identical(asv, asv_seqs[j]):
#             print(f'ASV {i} and ASV {j} are identical')
#             print(f'{asv} and {asv_seqs[j]}')

amplicon_seqs = {}
for a_strain in seqtable_isolate_amplicons.index:
    amplicon_seqs[a_strain] = []
    detected_asvs = seqtable_isolate_amplicons.loc[a_strain][seqtable_isolate_amplicons.loc[a_strain] > 10].index
    for asv in detected_asvs:
        amplicon_seqs[a_strain].append(Seq(asv))

'''
Compare the detected amplicon sequences to the barrnap sequences
print those that match and make a dictionary of the ones
that do not. Just print for each strain how many amplicon sequences
were detected and their count and whether they matched the barrnap
sequences or not.
'''

amplicon_not_in_barrnap = {}
for strain_iter in amplicon_seqs.keys():
    if strain_iter not in barrnap_seqs:
        print(f'{strain_iter} not in barrnap_seqs')
        continue
    amps_found = 0
    for amplicon_seq in amplicon_seqs[strain_iter]:
        found = False
        for barrnap_seq in barrnap_seqs[strain_iter]:
            if seq_identical(amplicon_seq, barrnap_seq):
                found = True
                amps_found += 1
                # print(f'{strain_iter} {amplicon_id_dict[amplicon_seq]} found in barrnap')
                break
        if not found:
            if strain_iter not in amplicon_not_in_barrnap:
                amplicon_not_in_barrnap[strain_iter] = []
            amplicon_not_in_barrnap[strain_iter].append(amplicon_seq)
            # print(f'{strain_iter} an amplicon_seq NOT found in barrnap')
    uniq_copies_strain_iter = len(set(barrnap_seqs[strain_iter]))
    # if you do not trust barrnap that it will have them be identical,
    # use the seq_identical function to check for identical sequences !
    if strain_iter in low_depth_strains:
        print(f'{strain_iter} copies detected: {amps_found} out of {uniq_copies_strain_iter} - (due to low depth)')
    else:
        print(f'{strain_iter} copies detected: {amps_found} out of {uniq_copies_strain_iter} depth: {sum(seqtable_isolate_amplicons.loc[strain_iter])}')
        # print(f'{strain_iter} copies detected: {amps_found} out of {uniq_copies_strain_iter} depth: {sum(seqtable_isolate_amplicons.loc[strain_iter])} and {len(amplicon_seqs[strain_iter])} amplicon sequences detected')




'''
typically we only want to consider barrnap sequences but as we learnt through this experiment
there can be errors in the genome particularly with Nanopore genomes so, also keep track of the 
amplicon sequences if and when available for the same genomes. This approach also works if there is no
barrnap sequence available for a genome but only an amplicon sequence and vice versa. Another advantage of using
amplicon sequences is that for new strains even if there is no genome sequence available, the amplicon
sequences can be used to identify the strain. - although this is a problem when copy number is unkown - this is usually conserved among species though!
'''

'''
The next goal is to make a database of sequences and the strains they match
we want to only do this based on barrnap but use the amplicon sequences when there
is something matching but not covered by barrnap to point towards a manual investigation
and modification of the database

we need for the 03_Analysis/infer_species_counts.py script,
    the script will receive a seqtable with samples as rows and amplicon sequences as columns
    we need to make for it a database of strain matches for its sequences of interest
     and a strain-sequence-copy matrix that it can use to infer strain counts from ASV counts
'''

os.makedirs('03_Analysis/SeqeuncesDB', exist_ok=True)

# there are duplicate sequences within strains
# but also between strains but for infer_species_counts.py
# we need to give it a copy number table with at least one asv unique to each strain
# each strain it ok for the matrix to be overdetermined but not underdetermined
# i.e. if there are two strains with the same sequence, we need to have a third sequence
# that is unique to one of the strains

'''
make a set of unique sequences from the whole database (barrnap_uniq_seqs) might still contain
duplicates between strains!

so while writing, keep track!
'''

uid_seq_master = {}

unique_seqs = set()

# make database only using unique barrnap sequences
uid = 0
for strain_iter in strains_list:
    for seq in barrnap_uniq_seqs[strain_iter]:
        if seq in unique_seqs:
            print(f'{strain_iter}\'s seq already in unique_seqs!')
            continue
        else:
            unique_seqs.add(seq)
            uid_seq_master[uid] = seq
            uid += 1

# write info table with uid, seq and strain association
with open('03_Analysis/SeqeuncesDB/OnlyIsolateSeqs.csv', 'w') as f:
    success = f.write('uid,strain,seq\n')
    for uid in uid_seq_master.keys():
        strains_associated = []
        uid_seq_iter = uid_seq_master[uid]
        for strain_iter in strains_list:
            if uid_seq_iter in barrnap_uniq_seqs[strain_iter]:
                strains_associated.append(strain_iter)
        if len(strains_associated) > 1:
            print(f'uid {uid} is associated with more than one strain: {strains_associated}')
        strain_print = ';'.join(strains_associated)
        success = f.write(f'{uid},{strain_print},{uid_seq_master[uid]}\n')

# write copies table including all strains for visualization
with open('03_Analysis/IsolateGenomes/strain_ASV_copies_table.csv', 'w') as f:
    # column names strains, row names uids and content number of copies
    success = f.write('uid,')
    for strain_iter in strains_list:
        pass
        success = f.write(f'{strain_iter},')
    success = f.write('\n')
    for uid in uid_seq_master.keys():
        success = f.write(f'{uid},')
        for strain_iter in strains_list:
            copies = 0
            uid_seq_iter = uid_seq_master[uid]
            for barrnap_seq in barrnap_seqs[strain_iter]:
                if seq_identical(uid_seq_iter, barrnap_seq):
                    copies += 1
            success = f.write(f'{copies},')
        success = f.write('\n')


'''
After this was used to identify strains, ESL0820 and ESL0827 were not found
even though they were added to these samples. Upon investigation, it was found
unfortunately due to low depth, ESL0820 and ESL0827 were not detected in the amplicon
sequencing run. However, they were found among the ASVs of the experiment in exactly the
samples they were added to. So, we need to add them to the database manually.

Upon further investigation, this is because of a nanopore genome assembly error
and the barrnap sequence is actually a chimera of two different strains.
adding further confidence, they were found both in the pilot experiment as
well as the complete experiment *2022 and 2024) respectively.
'''

# # special sequences - to consider - not a priori
# close to both uid:2 and uid:3 (not identical) - contains the differentiating SNV of both. Might be chimeric..!
# so considering its copy number to be 2
ESL0820_spl_seq = 'GATGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGGATCCGCCAAGCTTGCTTGGCGGTGAGAGTGGCGAACGGGTGAGTAATGCGTGACCAACCTGCCCTATGCTTCGGAATAGCTCCTGGAAACGGGTGGTAATGCCGGATGCTCCGCGCTGTCGCATGATGGTGCGGGAAAGGGTTTACCGGCATGGGATGGGGTCGCGTCCTATCAGCTTGTTGGCGGGGTGATGGCCTGCCAAGGCTTCGACGGGTAGCCGGCCTGAGAGGGCGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCGACGCCGCGTGCGGGATGACGGCCTTCGGGTTGTAAACCGCTTTTGATTGGGAGCAAGCGAGAGTGAGTGTACCTTTCGAATAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTATCCGGATTTATTGGGCGTAAAGAGCTCGTAGGCGGTTCGTCGCGTCTGGTGTGAAAGTCCATCGCTTAACGGTGGATCGGCGCCGGGTACGGGCGGACTGGAGTGCGGTAGGGGAGACTGGAATTCCCGGTGTAACGGTGGAATGTGTAGATATCGGGAAGAACACCGATGGCGAAGGCAGGTCTCTGGGCCGTCACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGGTGGATGCTGGATGTGGGGCCCGTTCCACGGGTTCCGTGTCGGAGCTAACGCGTTAAGCATCCCGCCTGGGGAGTACGGCCGCAAGGCTAAAACTCAAAGAAATTGACGGGGGCCCGCACAAGCGGCGGAGCATGCGGATTAATTCGATGCAACGCGAAGAACCTTACCTGGGCTTGACATGTGCCGGACGGCCGCGGAGACGCGGCTTCCCTTCGGGGCCGGTTCACAGGTGGTGCATGGTCGTCGTCAGCTCGTGTCGTGAGATGTTGGGTCAAGTCCCGCAACGAGCGCAACCCTCGCCTCGTGTTGCCAGCACGTTATGGTGGGAACTCACGGGGGACCGCCGGGGTTAACCCGGAGGAAGGTGGGGATGACGTCAGATCATCATGCCCCTTACGTCCAGGGCTTCACGCATGCTACAATGGCCGGTACAACGGGATGCGACATGGTGACATGGAGCGGATCCCTGAAAACCGGTCTCAGTTCGGATCGGAGCCTGCAACCCGGCTCCGTGAAGGCGGAGTCGCTAGTAATCGCGGATCAGCAACGCCGCGGTGAATGCGTTCCCGGGCCTTGTACACACCGCCCGTCAAGTCATGAAAGTGGGCAGCACCCGAAGCCGGTGGCCCAACCCGTTTGGGGGGGAGCCGTCTAAGGTGAGGTCCGCGATTGGGACT'
# called unk_1 before - matches to uid:11 which was not otherwise detected (uid:10 also from ESL0827 on the other hand was detected)
# so considering its copy number to be 1
ESL0827_spl_seq = 'GATGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGGATCCAGGCAGCTTGCTGTCTGGTGAGAGTGGCGAACGGGTGAGTAATGCGTGACCAACCTGCCCCATGCTTCGGAATAGCTCCTGGAAACGGGTGGTAATGCCGGATGCTCCGCATCATCGCATGATGGTGTGGGAAAGGGTTTACCGGCATGGGATGGGGTCGCGTCCTATCAGCTTGTTGGCGGGGTGATGGCCTGCCAAGGCTTCGACGGGTAGCCGGCCTGAGAGGGCGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGCGGGATGACGGCCTTCGGGTTGTAAACCGCTTTTGATTGGGAGCAAGCGAGAGTGAGTGTACCTTTCGAATAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTATCCGGATTTATTGGGCGTAAAGAGCTCGTAGGCGGTTCGTCGCGTCTGGTGTGAAAGTCCATCGCTTAACGGTGGATCGGCGCCGGGTACGGGCGGACTGGAGTGCGGTAGGGGAGACTGGAATTCCCGGTGTAACGGTGGAATGTGTAGATATCGGGAAGAACACCGATGGCGAAGGCAGGTCTCTGGGCCGTCACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGGTGGATGCTGGATGTGGGGCCCGTTCCACGGGTTCCGTGTCGGAGCTAACGCGTTAAGCATCCCGCCTGGGGAGTACGGCCGCAAGGCTAAAACTCAAAGAAATTGACGGGGGCCCGCACAAGCGGCGGAGCATGCGGATTAATTCGATGCAACGCGAAGAACCTTACCTGGGCTTGACATGTGCCGGACGGCCGCGGAGACGTGGCTTCCCTTCGGGGCCGGTTCACAGGTGGTGCATGGTCGTCGTCAGCTCGTGTCGTGAGATGTTGGGTCAAGTCCCGCAACGAGCGCAACCCTCGCCTCGTGTTGCCAGCACGTTATGGTGGGAACTCACGGGGGACCGCCGGGGTTAACCCGGAGGAAGGTGGGGATGACGTCAGATCATCATGCCCCTTACGTCCAGGGCTTCACGCATGCTACAATGGCCGGTACAGCGGGATGCGACATGGTGACATGGAGCGGATCCCTGAAAACCGGTCTCAGTTCGGATCGGAGCCTGCAACCCGGCTCCGTGAAGGCGGAGTCGCTAGTAATCGCGGATCAGCAACGCCGCGGTGAATGCGTTCCCGGGCCTTGTACACACCGCCCGTCAAGTCATGAAAGTGGGCAGCACCCGAAGCCGGTGGCCCAACCCGCGAGGGGGGGAGCCGTCTAAGGTGAGGTCCGCGATTGGGACT'

uid_seq_master_updated = {}

unique_seqs_updated = set()

# add the special sequences to the database
uid = 0
for strain_iter in strains_list:
    for seq in barrnap_uniq_seqs[strain_iter]:
        if seq in unique_seqs_updated:
            print(f'{strain_iter}\'s seq already in unique_seqs!')
            continue
        else:
            unique_seqs_updated.add(seq)
            uid_seq_master_updated[uid] = seq
            uid += 1
# add the special sequences to the database
uid_seq_master_updated[uid] = Seq(ESL0820_spl_seq)
uid += 1
uid_seq_master_updated[uid] = Seq(ESL0827_spl_seq)

# write info table with uid, seq and strain association
with open('03_Analysis/SeqeuncesDB/OnlyIsolateSeqs_updated.csv', 'w') as f:
    success = f.write('uid,strain,seq\n')
    for uid in uid_seq_master_updated.keys():
        strains_associated = []
        uid_seq_iter = uid_seq_master_updated[uid]
        for strain_iter in strains_list:
            if uid_seq_iter in barrnap_uniq_seqs[strain_iter]:
                strains_associated.append(strain_iter)
        if len(strains_associated) > 1:
            print(f'uid {uid} is associated with more than one strain: {strains_associated}')
        strain_print = ';'.join(strains_associated)
        if uid_seq_iter == ESL0820_spl_seq:
            strain_print = 'ESL0820'
        if uid_seq_iter == ESL0827_spl_seq:
            strain_print = 'ESL0827'
        success = f.write(f'{uid},{strain_print},{uid_seq_master_updated[uid]}\n')


uid_seq_master_amplicons = {}

unique_seqs_amplicons = set()

with open('03_Analysis/SeqeuncesDB/OnlyIsolateSeqs_amplicons_updated.csv', 'w') as f:
    success = f.write('uid,strain,seq\n')
    uid = 0
    for strain_iter in strains_list:
        for seq in barrnap_uniq_seqs[strain_iter]:
            if seq in unique_seqs_amplicons:
                print(f'{strain_iter}\'s seq already in unique_seqs!')
                continue
            else:
                unique_seqs_amplicons.add(seq)
                uid_seq_master_amplicons[uid] = seq
                uid += 1
    # add the special sequences to the database
    uid_seq_master_amplicons[uid] = Seq(ESL0820_spl_seq)
    uid += 1
    uid_seq_master_amplicons[uid] = Seq(ESL0827_spl_seq)
    uid += 1
    # add the amplicon sequences to the database
    for strain_iter in strains_list:
        if strain_iter not in amplicon_seqs:
            print(f'{strain_iter} not in amplicon_seqs')
            continue
        for seq in amplicon_seqs[strain_iter]:
            if seq in unique_seqs_amplicons:
                print(f'{strain_iter}\'s seq already in unique_seqs!')
                continue
            else:
                unique_seqs_amplicons.add(seq)
                uid_seq_master_amplicons[uid] = seq
                uid += 1

# write info table with uid, seq and strain association
with open('03_Analysis/SeqeuncesDB/OnlyIsolateSeqs_amplicons_updated.csv', 'w') as f:
    success = f.write('uid,strain,seq\n')
    for uid in uid_seq_master_amplicons.keys():
        strains_associated = []
        uid_seq_iter = uid_seq_master_amplicons[uid]
        for strain_iter in strains_list:
            if uid_seq_iter in barrnap_uniq_seqs[strain_iter]:
                # g denotes it matches a sequence found in an isolate genome (by barrnap)
                strains_associated.append(f'g_{strain_iter}')
        for strain_iter in strains_list:
            if strain_iter in low_depth_strains:
                continue
            if strain_iter not in amplicon_seqs:
                continue
            if uid_seq_iter in amplicon_seqs[strain_iter]:
                strains_associated.append(strain_iter)
        if len(strains_associated) > 1:
            print(f'uid {uid} is associated with more than one strain: {strains_associated}')
        strain_print = ';'.join(strains_associated)
        if uid_seq_iter == ESL0820_spl_seq:
            strain_print = 'ESL0820'
        if uid_seq_iter == ESL0827_spl_seq:
            strain_print = 'ESL0827'
        success = f.write(f'{uid},{strain_print},{uid_seq_master_amplicons[uid]}\n')


'''
Now we make the copies table for the 03_Analysis/infer_species_counts.py
using only the barrnap sequences (updated)
This needs to be done for each project as needed because strains that have
identical / shared ASVs need to be "handled"

Look at the copy numbers heatmap and avoid including strains that we know
are identical and strains that were not in the experiment. This was we avoid
attributing reads to another strain that shares an ASV with a strain we know
is not present in the sample! This is only possible for synthetic community
experiments where we know the strains that were added!
'''

# heatmap made using 03_Analysis/Visualize_results.Rmd
# saved in 03_Analysis/IsolateGenomes/strain_ASV_copies_heatmap.pdf


'''
uid 26 is only sequence - 4 copies in both ESL0353 and ESL0185 - only ESL0185 is in the experiment
uid 29-32 are present in both ESL0183 and ESL0262 - only ESL0183 is in the experiment (ESL0262 is its pair only in the pilot but cannot distinguish these two anyway)
uid 57 is shared by ESL0351 (3 copies) and ESL0186 (2 copies) - only ESL0186 is in the experiment
** If ESL0351 were included it would have been fine as it does have an additional 4th seq in 1 copy! uid: 60**

so, exclude ESL0351, ESL0353, ESL0262 from the copies table
'''

exclude_strains = ['ESL0351', 'ESL0353', 'ESL0262']

# write copies table
# be careful of sequences not in the barrnap_seqs

all_barrnap_seqs = set(list(itertools.chain.from_iterable(barrnap_seqs.values())))

with open('03_Analysis/SeqeuncesDB/OnlyIsolateSeqs_updated_copies.csv', 'w') as f:
    # column names strains, row names uids and content number of copies
    success = f.write('uid,')
    for strain_iter in strains_list:
        if strain_iter in exclude_strains:
            continue
        success = f.write(f'{strain_iter},')
    success = f.write('\n')
    for uid in uid_seq_master_updated.keys():
        # skip the replaced uids
        if uid == 2 or uid == 3 or uid == 11:
            continue
        success = f.write(f'{uid},')
        for strain_iter in strains_list:
            if strain_iter in exclude_strains:
                continue
            copies = 0
            uid_seq_iter = uid_seq_master_updated[uid]
            if uid_seq_iter not in all_barrnap_seqs:
                print(f'checking {strain_iter}: {uid} not found because seq not in barrnap_seqs!')
            else:
                for barrnap_seq_iter in barrnap_seqs[strain_iter]:
                    if seq_identical(uid_seq_iter, barrnap_seq_iter):
                        copies += 1
            # since we added special sequences, replace 
            # corresponding predicted ones with 0 for copies
            if uid_seq_iter == ESL0820_spl_seq and strain_iter == 'ESL0820':
                copies = 2
            if uid_seq_iter == ESL0827_spl_seq and strain_iter == 'ESL0827':
                copies = 1
            success = f.write(f'{copies},')
        success = f.write('\n')


# write out the amplicon discovered sequences so infer_species_counts.py
# can check if there are any that match!
# more information on where the amplicon is detected for context!
with open('03_Analysis/SeqeuncesDB/OnlyIsolateSeqs_amplicons_updated_presence.csv', 'w') as f:
    success = f.write('uid,strain,count,total,seq\n')
    for uid in uid_seq_master_amplicons.keys():
        if uid in uid_seq_master_updated:
            continue
        else:
            found_in = []
            for strain_iter in seqtable_isolate_amplicons.index:
                if strain_iter in low_depth_strains:
                    continue
                counts = 0
                uid_seq_iter = uid_seq_master_amplicons[uid]
                counts = seqtable_isolate_amplicons.loc[strain_iter][uid_seq_iter]
                if counts > 10:
                    found_in.append(strain_iter)
                    success = f.write(f'{uid},{strain_iter},{counts},{sum(seqtable_isolate_amplicons.loc[strain_iter])},{uid_seq_master_amplicons[uid]}\n')