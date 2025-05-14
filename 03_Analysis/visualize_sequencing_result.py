import os
import pandas as pd
import argparse
from statistics import mean
from statistics import median
import math
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob

'''
This script is written to judge the success of the nirmalization approach
Plate1_10 contains the numbers used for normalization but this is not correct
because the dilution factor was considered to be 10 instead of 100
knowing this we compare predicted and actual values for number of reads
If the dilution factor were correct, the concentration of the sample would be
10 times more than the predicted value so we should also have 10 times more reads
'''

plate = 'Plate6'
plate_dir = f'/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/PriorityEffectsExperiment/RawData/Plate6/'
sample_sheet = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/PriorityEffectsExperiment/RawData/Plate6/user_biosamples.csv'
predictions_sheet = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/01_LabNotes/05-PCRAmpliPrep/Plate6/plate6_all_norm_data.csv'
# plate = 'Plate5'
# plate_dir = f'/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/PriorityEffectsExperiment/RawData/Plate5/'
# sample_sheet = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/PriorityEffectsExperiment/RawData/Plate5/user_biosamples.csv'
# predictions_sheet = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/01_LabNotes/05-PCRAmpliPrep/Plate5/plate5_all_norm_data.csv'
# plate = 'Plate4'
# plate_dir = f'/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/PriorityEffectsExperiment/RawData/Plate4/'
# sample_sheet = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/PriorityEffectsExperiment/RawData/Plate4/user_biosamples.csv'
# predictions_sheet = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/01_LabNotes/05-PCRAmpliPrep/Plate4/plate4_all_norm_data.csv'
# plate = 'Plate3'
# plate_dir = f'/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/PriorityEffectsExperiment/RawData/Plate3/'
# sample_sheet = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/PriorityEffectsExperiment/RawData/Plate3/user_biosamples.csv'
# predictions_sheet = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/01_LabNotes/05-PCRAmpliPrep/Plate3/plate3_all_norm_data.csv'
# plate = 'Plate2'
# plate_dir = f'/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/RawData/Plate2'
# sample_sheet = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/RawData/Plate2/AP_16S_Pl2_Circular_Consensus_Sequencing_Reads/SAMPLESHEET.csv'
# predictions_sheet = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/01_LabNotes/05-PCRAmpliPrep/Plate2/Plate2_norm_all_norm_data.csv'
# plate = 'Plate1'
# plate_dir = f'/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/RawData/Plate1'
# sample_sheet = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/RawData/Plate1/AP_16S_Pl1_Circular_Consensus_Sequencing_Reads/SAMPLESHEET.csv'
# predictions_sheet = '/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/01_LabNotes/05-PCRAmpliPrep/Plate1/Plate1_norm_all_norm_data.csv'

sample_sheet_df = pd.read_csv(sample_sheet)
raw_file_dict = sample_sheet_df.set_index('Biosample').to_dict()['Barcode']
# raw_file_dict = sample_sheet_df.set_index('Bio Sample').to_dict()['Barcode']


# os.chdir('/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/20230531_PacBio_Sequencing')

# collect a df of number of reads per sample from 
# SequencingData-Plate2/00-raw_reads_renamed
# df_summary = pd.DataFrame(columns=['Sample Name', 'HiFi Reads'])
# raw_files = glob('SequencingData-Plate2/00-raw_reads_renamed/Apis_bee_guts/*.fastq') + \
#             glob('SequencingData-Plate2/00-raw_reads_renamed/AQ_comms/*.fastq') + \
#             glob('SequencingData-Plate2/00-raw_reads_renamed/ESL_strains/*.fastq')
# raw_paths_dir = {os.path.basename(x).split('.fastq')[0]: x for x in raw_files}
# for i, sample in enumerate(raw_paths_dir.keys()):
#     reads = sum(1 for line in open(raw_paths_dir[sample])) / 4
#     df_summary = df_summary._append({'Sample Name': sample, 'HiFi Reads': reads}, ignore_index=True)
#     print(f'{i+1}/{len(raw_paths_dir.keys())}', end='\r')
barcode_summary_file = glob(f'{plate_dir}/*barcodes_summary.csv')
# barcode_summary_file = glob(f'{plate_dir}/*/*barcodes_summary.csv') 
df_summary = pd.read_csv(barcode_summary_file[0])
total_lib_size = sum(df_summary['HiFi Reads'])
# exclude the row called No Name
df_summary = df_summary[df_summary['Sample Name'] != 'No Name']
df_summary = df_summary[['Sample Name', 'HiFi Reads']]
df_predictions = pd.read_csv(predictions_sheet)
df_predictions = df_predictions.sort_values(by=['sample'])
df_predictions['vol_dil_added'] = [f'{x},{y}' for x,y in zip(df_predictions['volume_to_be_added'], df_predictions['dilution'])]
df_predictions['dilution'] = 1/df_predictions['dilution']
df_predictions['reads sequenced'] = df_predictions['sample'].map(df_summary.set_index('Sample Name')['HiFi Reads'])
df_predictions['reads predicted'] = df_predictions['ratio_in_lib'] * total_lib_size

os.makedirs(f'{plate_dir}/visualize_seq_result', exist_ok=True)

# plot a bar plot of df_summary for each sample
# x-axis: sample name
# y-axis: number of reads
# outline in black but not thick
plt.figure(figsize=(12, 6))
plt.axhline(y=10000, color='black', linestyle='dotted')
plt.axhline(y=5000, color='black', linestyle='solid')
# color bars red if < 5000 else blue
for i in range(len(df_summary['HiFi Reads'])):
    if df_summary['HiFi Reads'][i] > 10000:
        plt.bar(df_summary['Sample Name'][i], df_summary['HiFi Reads'][i], color='#33a02c', edgecolor='black', linewidth=0.5)
    if df_summary['HiFi Reads'][i] < 10000 and df_summary['HiFi Reads'][i] > 1000:
        plt.bar(df_summary['Sample Name'][i], df_summary['HiFi Reads'][i], color='#b2df8a', edgecolor='black', linewidth=0.5)
    if df_summary['HiFi Reads'][i] < 1000:
        plt.bar(df_summary['Sample Name'][i], df_summary['HiFi Reads'][i], color='#e31a1c', edgecolor='black', linewidth=0.5)
plt.xticks(rotation=90, size = 6)
plt.yscale('log')
plt.ylabel('Number of Reads')
plt.xlabel('Sample')
plt.title('Number of Reads per Sample')
# plt.show()
plt.savefig(f'{plate_dir}/visualize_seq_result/reads_per_sample.png', dpi=300)


# plot dilution and reads sequenced for each sample
plt.figure(figsize=(12, 6))
plt.axhline(y=10000, color='black', linestyle='dotted')
plt.axhline(y=5000, color='black', linestyle='solid')
# plot sample vs reads sequenced but color by dilution
sns.barplot(x='sample', y='reads sequenced', hue='dilution', data=df_predictions)
plt.xticks(rotation=90, size = 6)
plt.yscale('log')
plt.ylabel('Number of Reads')
plt.xlabel('Sample')
plt.title('Number of Reads per Sample')
# plt.show()
plt.savefig(f'{plate_dir}/visualize_seq_result/reads_per_sample_dilution.png', dpi=300)

# make a box plot of predicted and sequenced reads with line joining them
plt.figure(figsize=(12, 6))
df_plot = df_predictions.melt(id_vars=['sample', 'dilution'], value_vars=['reads predicted', 'reads sequenced'])
ax = sns.lineplot(x='variable', y='value', data=df_plot, hue='sample', style='dilution', dashes=True, markers=True)
ax = sns.boxplot(x='variable', y='value', data=df_plot, hue='variable')
hand, labl = ax.get_legend_handles_labels()
handout=[]
lablout=[]
for h,l in zip(hand,labl):
    if l in ['0.01', '0.1', '1.0', 'reads sequenced', 'reads predicted']:
        lablout.append(l)
        handout.append(h)
ax.legend(handout, lablout)
plt.xticks(rotation=0, size = 10)
plt.yscale('log')
plt.ylabel('Number of Reads')
plt.xlabel('Sample')
plt.title('Number of Reads per Sample')
# do not show hue legend for lineplot
# plt.show()
plt.savefig(f'{plate_dir}/visualize_seq_result/reads_per_sample_boxplot.png', dpi=300)

# plot predicted and sequenced reads as a scatter plot
plt.figure(figsize=(12, 6))
sns.scatterplot(x='reads predicted', y='reads sequenced', data=df_predictions, hue = 'volume_to_be_added', style='dilution')
plt.xlabel('Predicted Reads')
plt.ylabel('Sequenced Reads')
# plt.axhline(y=10000, color='black', linestyle='dashed')
# plt.axhline(y=5000, color='black', linestyle='dotted')
plt.axhline(y=1000, color='black', linestyle='dotted')
# add an x = y line
plt.plot([0, 1000000], [0, 1000000], color='black', linestyle='dashed')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.legend(loc='best', ncol = 2)
# plt.show()
plt.savefig(f'{plate_dir}/visualize_seq_result/predicted_vs_sequenced_reads.png', dpi=300)

# plot predicted and sequenced reads as a scatter plot
plt.figure(figsize=(12, 6))
sns.scatterplot(x='copies', y='reads sequenced', data=df_predictions, hue = 'volume_to_be_added', style='dilution')
plt.xlabel('qPCR copies')
plt.ylabel('Sequenced Reads')
plt.axhline(y=10000, color='black', linestyle='dashed')
plt.axhline(y=5000, color='black', linestyle='dotted')
plt.axhline(y=1000, color='black', linestyle='dotted')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.legend(loc='best', ncol = 3)
# plt.show()
plt.savefig(f'{plate_dir}/visualize_seq_result/copies_vs_sequenced_reads.png', dpi=300)

# plot ng_added to library and reads sequenced
plt.figure(figsize=(12, 6))
sns.scatterplot(x='ng_added', y='reads sequenced', data=df_predictions, hue = 'volume_to_be_added', style='dilution', s=60)
axvlines = [0.1, 1, 5]
for i in axvlines:
    plt.axvline(x=i, color='black', linestyle='dotted')
plt.xlabel('ng added to library')
plt.ylabel('Sequenced Reads')
plt.axhline(y=10000, color='black', linestyle='dashed')
plt.axhline(y=5000, color='black', linestyle='dotted')
plt.axhline(y=1000, color='black', linestyle='dotted')
plt.yscale('log')
plt.xscale('log')
plt.tight_layout()
plt.legend(loc='best', ncol = 3)
# plt.show()
plt.savefig(f'{plate_dir}/visualize_seq_result/ng_added_vs_sequenced_reads.png', dpi=300)

# plot conc to library and reads sequenced
df_predictions_plot = df_predictions.copy()
df_predictions_plot['conc_final'] = df_predictions_plot['conc'] * df_predictions_plot['dilution']
plt.figure(figsize=(12, 6))
sns.scatterplot(x='conc_final', y='reads sequenced', data=df_predictions_plot, hue = 'volume_to_be_added', style='dilution', s=60)
for i in axvlines:
    plt.axvline(x=i, color='black', linestyle='dotted')
plt.xlabel('Concentration of PCR product (after dilution)')
plt.ylabel('Sequenced Reads')
plt.axhline(y=10000, color='black', linestyle='dashed')
plt.axhline(y=5000, color='black', linestyle='dotted')
plt.axhline(y=1000, color='black', linestyle='dotted')
plt.yscale('log')
plt.xscale('log')
plt.tight_layout()
plt.legend(loc='best', ncol = 3)
# plt.show()
plt.savefig(f'{plate_dir}/visualize_seq_result/conc_final_vs_sequenced_reads.png', dpi=300)


plt.figure(figsize=(12, 6))
sns.scatterplot(x='conc', y='reads sequenced', data=df_predictions, hue = 'volume_to_be_added', style='dilution', s=60)
for i in axvlines:
    plt.axvline(x=i, color='black', linestyle='dotted')
plt.xlabel('Concentration of PCR product (before dilution)')
plt.ylabel('Sequenced Reads')
plt.axhline(y=10000, color='black', linestyle='dashed')
plt.axhline(y=5000, color='black', linestyle='dotted')
plt.axhline(y=1000, color='black', linestyle='dotted')
plt.yscale('log')
plt.xscale('log')
plt.tight_layout()
plt.legend(loc='best', ncol = 3)
# plt.show()
plt.savefig(f'{plate_dir}/visualize_seq_result/conc_vs_sequenced_reads.png', dpi=300)

# plot ct to library and reads sequenced
plt.figure(figsize=(12, 6))
sns.scatterplot(x='ct', y='reads sequenced', data=df_predictions, hue = 'volume_to_be_added', style='dilution', s=60)
for i in [30]:
    plt.axvline(x=i, color='black', linestyle='dotted')
plt.xlabel('Ct')
plt.ylabel('Sequenced Reads')
plt.axhline(y=10000, color='black', linestyle='dashed')
plt.axhline(y=5000, color='black', linestyle='dotted')
plt.axhline(y=1000, color='black', linestyle='dotted')
plt.yscale('log')
plt.tight_layout()
plt.legend(loc='best', ncol = 3)
# plt.show()
plt.savefig(f'{plate_dir}/visualize_seq_result/ct_vs_sequenced_reads.png', dpi=300)

# are ct value of orginal sample correlated with ct value of pcr prooduct?
df_original_ct = pd.read_csv('../qPCR_values_of_Original_samples.csv')
df_predictions['original_ct'] = df_predictions['sample'].map(df_original_ct.set_index('Sample.Name')['Ct_mean'])
plt.figure(figsize=(12, 6))
sns.scatterplot(x='original_ct', y='ct', data=df_predictions, hue = 'volume_to_be_added', style='dilution', s=60)
plt.xlabel('Ct of original sample')
plt.ylabel('Ct of PCR product')
plt.axhline(y=30, color='black', linestyle='dotted')
plt.axvline(x=30, color='black', linestyle='dotted')
plt.tight_layout()
plt.legend(loc='best', ncol = 3)
# plt.show()
plt.savefig(f'{plate_dir}/visualize_seq_result/original_ct_vs_ct.png', dpi=300)

# are ct value of orginal sample correlated with ct value of pcr prooduct and seq reads?
plt.figure(figsize=(12, 6))
sns.scatterplot(x='original_ct', y='ct', data=df_predictions, hue = 'reads sequenced', style='dilution', s=60)
plt.xlabel('Ct of original sample')
plt.ylabel('Ct of PCR product')
plt.axhline(y=30, color='black', linestyle='dotted')
plt.axvline(x=30, color='black', linestyle='dotted')
plt.tight_layout()
plt.legend(loc='best', ncol = 3)
# plt.show()
plt.savefig(f'{plate_dir}/visualize_seq_result/original_ct_vs_ct_vs_reads.png', dpi=300)

# original ct vs reads
plt.figure(figsize=(12, 6))
sns.scatterplot(x='original_ct', y='reads sequenced', data=df_predictions, hue = 'volume_to_be_added', style='dilution', s=60)
plt.xlabel('Ct of original sample')
plt.ylabel('Sequenced Reads')
plt.yscale('log')
plt.axhline(y=10000, color='black', linestyle='dashed')
plt.axhline(y=5000, color='black', linestyle='dotted')
plt.axhline(y=1000, color='black', linestyle='dotted')
plt.tight_layout()
plt.legend(loc='best', ncol = 3)
# plt.show()
plt.savefig(f'{plate_dir}/visualize_seq_result/original_ct_vs_reads.png', dpi=300)