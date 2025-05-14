import os
import sys
import pandas as pd
import numpy as np

thresholds = ['999', '998', '995', '99', '97', '95', '90']

for threshold in thresholds:
    cdhit_file = f'03_Analysis/FindNatural_16S/cd-hit_output/asv_ext_{threshold}_clusters.clstr'
    if not os.path.exists(cdhit_file):
        print(f"File {cdhit_file} does not exist.")
        sys.exit(1)

    out_summary = f'03_Analysis/FindNatural_16S/cd-hit_output/asv_ext_mod_{threshold}_summary.csv'
    out_file = f'03_Analysis/FindNatural_16S/cd-hit_output/asv_ext_{threshold}_clusters_mod.clstr'

    strain_name_spec_dict = {'ESL0825': 'Bifidobacterium apousia',
                            'ESL0820': 'Bifidobacterium apousia',
                            'ESL0822': 'Bifidobacterium asteroides',
                            'ESL0170': 'Bifidobacterium asteroides',
                            'ESL0824': 'Bifidobacterium polysaccharolyticum',
                            'ESL0198': 'Bifidobacterium polysaccharolyticum',
                            'ESL0827': 'Bifidobacterium sp1.',
                            'ESL0199': 'Bifidobacterium sp1.',
                            'ESL0200': 'Bifidobacterium sp2.',
                            'ESL0819': 'Bifidobacterium sp2.',
                            'ESL0197': 'Bifidobacterium coryneforme',
                            'ESL0294': 'Bombilactobacillus mellis',
                            'ESL0394': 'Bombilactobacillus mellis',
                            'ESL0295': 'Bombilactobacillus mellis',
                            'ESL1028': 'Bombilactobacillus mellifer',
                            'ESL0353': 'Lactobacillus apis',
                            'ESL0263': 'Lactobacillus apis',
                            'ESL0185': 'Lactobacillus apis',
                            'ESL0185_ESL0353': 'Lactobacillus apis',
                            'ESL0183': 'Lactobacillus helsingborgensis',
                            'ESL0262': 'Lactobacillus helsingborgensis',
                            'ESL0262_ESL0183': 'Lactobacillus helsingborgensis',
                            'ESL0835': 'Lactobacillus helsingborgensis',
                            'ESL0354': 'Lactobacillus helsingborgensis',
                            'ESL0393': 'Lactobacillus melliventris',
                            'ESL0259': 'Lactobacillus melliventris',
                            'ESL0260': 'Lactobacillus melliventris',
                            'ESL0184': 'Lactobacillus melliventris',
                            'ESL0350': 'Lactobacillus melliventris',
                            'ESL0186': 'Lactobacillus kullabergensis',
                            'ESL0351_ESL0186': 'Lactobacillus kullabergensis',
                            'ESL0261': 'Lactobacillus kullabergensis',
                            'ESL0351': 'Lactobacillus kullabergensis'}

    # goal is to identify threshold that puts only strain together that are
    # from the same species but also does not separate strains of the same species into
    # many different cluster

    # lines not defining the cluster look like '4\t1007aa, >ASV2510... *\n'

    # Read the CD-HIT output file
    cluster_members_dict = {}
    with open(cdhit_file, 'r') as f:
        for line in f:
            if line.startswith('>Cluster'):
                cluster_id = line.split(' ')[1].strip()
                cluster_members_dict[cluster_id] = []
            else:
                # Extract the strain name from the line
                strain_name = line.split('>')[1].split('...')[0].strip()
                # Add the strain name to the current cluster
                if cluster_id in cluster_members_dict:
                    cluster_members_dict[cluster_id].append(strain_name)

    # Create a DataFrame to store the cluster information
    clusters = []
    for cluster_id, members in cluster_members_dict.items():    
        # Create a dictionary for each cluster
        species = []
        for member in members:
            if 'ESL' in member:
                a_species = strain_name_spec_dict.get(member.split('-')[1])
                if a_species not in species:
                    species.append(a_species)
        cluster_info = {
            'Cluster ID': cluster_id,
            # 'Members': ', '.join(members),
            'Species': ';'.join(species),
            'Num Members': len(members),
            'Num ESL Members': len([m for m in members if 'ESL' in m]),
            'ESL Members': ';'.join([m for m in members if 'ESL' in m]),
            'Num Species': len(set(species))
        }
        clusters.append(cluster_info)
    clusters_df = pd.DataFrame(clusters)
    clusters_df.to_csv(out_summary, index=False)
