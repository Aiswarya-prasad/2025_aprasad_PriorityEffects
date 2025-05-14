import os
import sys
import logging
import pandas as pd
import numpy as np
import argparse


# counts_table_file = '03_Analysis/Results/tables/priority_effect_strength/temp/counts_all_full_combos.csv'
# low_detection_samples_file = '03_Analysis/PriorityEffectsExperiment/samples_low_detection.csv'
# output_file = '03_Analysis/Results/tables/df_priority_effect_strength_simple.csv'
# output_file_sim = '03_Analysis/Results/tables/df_priority_effect_strength_simple_by_sim.csv'

argparser = argparse.ArgumentParser(description='Calculate priority effect strength')
argparser.add_argument('--counts_table_file', type=str, help='Path to the counts table file', required=True)
argparser.add_argument('--low_detection_samples_file', type=str, help='Path to the low detection samples file', required=True)
argparser.add_argument('--output_file', type=str, help='Path to the output file', required=True)
args = argparser.parse_args()
counts_table_file = args.counts_table_file
low_detection_samples_file = args.low_detection_samples_file

args = argparser.parse_args()

counts_table_file = args.counts_table_file
low_detection_samples_file = args.low_detection_samples_file
output_file = args.output_file
output_file_sim = output_file.replace('.csv', '_by_sim.csv')


df_counts_all = pd.read_csv(counts_table_file)
df_counts = df_counts_all[['ID', 'Strain', 'CommunityCombo', 'spec_name', 'abs_counts', 'arrival', 'limit_passed', 'geq_passed', 'sample_geq', 'Type']]

community_combos = ['1', '2', '3', '4', '5', '6']
# community_combos = ['1', '2', '3', '4', '5', '6', '1_rep']

low_detection_ids = [x for x in pd.read_csv(low_detection_samples_file)['x'].values]

strains_to_skip = [
    'ESL0197', # no partner
    'ESL1028', # no partner
    'ESL0294', # no partner
    'unknown', # no partner
    'ESL0353', # not in experiment
    'ESL0351' # not in experiment
]

with open (output_file, 'w+') as f_simple:
    with open(output_file_sim, 'w+') as f_simple_sim:
        success = f_simple.write('Strain,CommunityCombo,priority_effect_strength_median,num_samples_first,num_samples_second,spec_name\n')
        success = f_simple_sim.write('Strain,CommunityCombo,priority_effect_strength_median,null_f,null_s,num_samples_first,num_samples_second,spec_name\n')
        # success = f_simple_sim.write('Strain,CommunityCombo,priority_effect_strength_median,confint_bot,confint_top,num_samples_first,num_samples_second,spec_name\n')
        for strain in df_counts['Strain'].unique():
            spec_name_strain = df_counts[df_counts['Strain'] == strain]['spec_name'].values[0]
            df_strain = df_counts[df_counts['Strain'] == strain]
            df_strain = df_strain[(df_strain['arrival'] == 'first') | (df_strain['arrival'] == 'second')]
            # empty df
            if strain in strains_to_skip or df_strain.empty:
                print('no data for strain: ', strain)
                continue
            else:
                 print(strain)
            for community_combo in community_combos:
                    df_community_strain = df_strain[df_strain['CommunityCombo'] == community_combo]
                    # median of those that say first
                    first_values = df_community_strain[df_community_strain['arrival'] == 'first']['abs_counts'].values
                    second_values = df_community_strain[df_community_strain['arrival'] == 'second']['abs_counts'].values
                    num_first = len(first_values)
                    num_second = len(second_values)
                    # print(f'Community: {community_combo}')
                    # print(f'log10(first): {np.log10(np.median(first_values))} (bees - {num_first})')
                    # print(f'log10(second): {np.log10(np.median(second_values))} (bees - {num_second})')
                    # print(f'priority_effect_strength: {np.log10(np.median(first_values)/np.median(second_values))}')
                    if np.median(second_values) != 0 and np.median(first_values) != 0:
                        priority_eff_ratio = np.log10(np.median(first_values)/np.median(second_values))
                    else:
                        if np.median(second_values) == 0 and np.median(first_values) != 0:
                            priority_eff_ratio = np.log10(np.median(first_values)/1)
                        elif np.median(second_values) != 0 and np.median(first_values) == 0:
                            priority_eff_ratio = -1
                        else:
                            priority_eff_ratio = 0
                            
                    success = f_simple.write(f'{strain},{community_combo},{priority_eff_ratio},{num_first},{num_second},{spec_name_strain}\n')
                    # only up to here relevant for Fig. 5
                    
                    # simulate n=100 times?
                    n = 1000
                    priority_effect_strength = []
                    for i in range(n):
                        first_a = np.random.choice(first_values)
                        second_a = np.random.choice(second_values)
                        if second_a != 0 and first_a != 0:
                            strength = np.log10((first_a)/(second_a))
                        else:
                            if second_a == 0 and first_a != 0:
                                strength = np.log10((first_a)/(1))
                            elif second_a != 0 and first_a == 0:
                                strength = -1
                            else:
                                strength = 0
                        priority_effect_strength.append(strength)
                        # if strength < 0:
                        #     print(f'strength ({strain}): {strength}')
                    priority_effect_strength_null_f = []
                    for i in range(n):
                        first_a = np.random.choice(first_values)
                        second_a = np.random.choice(first_values)
                        if second_a != 0 and first_a != 0:
                            strength = np.log10((first_a)/(second_a))
                        else:
                            if second_a == 0 and first_a != 0:
                                strength = np.log10((first_a)/(1))
                            elif second_a != 0 and first_a == 0:
                                strength = -1
                            else:
                                strength = 0
                        priority_effect_strength_null_f.append(strength)
                    priority_effect_strength_null_s = []
                    for i in range(n):
                        first_a = np.random.choice(second_values)
                        second_a = np.random.choice(second_values)
                        if second_a != 0 and first_a != 0:
                            strength = np.log10((first_a)/(second_a))
                        else:
                            if second_a == 0 and first_a != 0:
                                strength = np.log10((first_a)/(1))
                            elif second_a != 0 and first_a == 0:
                                strength = -1
                            else:
                                strength = 0
                    priority_effect_strength_null = priority_effect_strength_null_f + priority_effect_strength_null_s
                    # print(f'priority_effect_strength(n=100): {np.median(priority_effect_strength)}')
                    # confidences = np.percentile(priority_effect_strength, [5, 95])
                    # print(f'confidence interval: {confidences[0]} - {confidences[1]}')
                    success = f_simple_sim.write(f'{strain},{community_combo},{np.median(priority_effect_strength)},{np.median(priority_effect_strength_null_f)},{np.median(priority_effect_strength_null_s)},{num_first},{num_second},{spec_name_strain}\n')

"""
obsolete!
"""

'''
with open('03_Analysis/Results/tables/df_priority_effect_strength_null_and_actual_simulations.csv', 'w+') as f_sim_null:
    success = f_sim_null.write('Strain,priority_effect_strength,CommunityCombo,comparison,arrival,spec_name\n')
    for strain in df_counts['Strain'].unique():
        spec_name_strain = df_counts[df_counts['Strain'] == strain]['spec_name'].values[0]
        if strain in strains_to_skip:
            continue
        df_strain = df_counts[df_counts['Strain'] == strain]
        df_strain = df_strain[(df_strain['arrival'] == 'first') | (df_strain['arrival'] == 'second')]
            # dict of dicts of dicts
            # draw from first and second, same community combo,
            # draw from first and seconf from different community combos
            # draw from first and first from same community combo
            # draw from second and second from same community combo
            # have a column to track if it is first and second, first and first or second and second and
            # same or different community combo
        dict_community_strain = df_strain[['CommunityCombo', 'abs_counts', 'arrival']].groupby(['arrival', 'CommunityCombo'])['abs_counts'].apply(list).to_dict()
    # do first vs second same community combo
        for community_combo in community_combos:
            first_values = dict_community_strain[('first', community_combo)]
            second_values = dict_community_strain[('second', community_combo)]
            num_first = len(first_values)
            num_second = len(second_values)
            if num_first == 0 or num_second == 0:
                continue
            priority_effect_strength = []
            for i in range(100):
                first_a = np.random.choice(first_values)
                second_a = np.random.choice(second_values)
                strength = np.log10((first_a)/(second_a))
                priority_effect_strength.append(strength)
                success = f_sim_null.write(f'{strain},{strength},{community_combo},same,actual,{spec_name_strain}\n')

            # do first vs second different community combo
        for community_combo in community_combos:
            for community_combo_2 in community_combos:
                if community_combo == community_combo_2:
                    continue
                first_values = dict_community_strain[('first', community_combo)]
                second_values = dict_community_strain[('second', community_combo_2)]
                num_first = len(first_values)
                num_second = len(second_values)
                if num_first == 0 or num_second == 0:
                    continue
                priority_effect_strength = []
                for i in range(100):
                    first_a = np.random.choice(first_values)
                    second_a = np.random.choice(second_values)
                    strength = np.log10((first_a)/(second_a))
                    priority_effect_strength.append(strength)
                    success = f_sim_null.write(f'{strain},{strength},{community_combo},different,actual,{spec_name_strain}\n')

            # do first vs first same community combo
        for community_combo in community_combos:
            first_values = dict_community_strain[('first', community_combo)]
            num_first = len(first_values)
            if num_first < 2:
                continue
            priority_effect_strength = []
            for i in range(100):
                first_a = np.random.choice(first_values)
                second_a = np.random.choice(first_values)
                strength = np.log10((first_a)/(second_a))
                priority_effect_strength.append(strength)
                success = f_sim_null.write(f'{strain},{strength},{community_combo},same,null_first,{spec_name_strain}\n')

            # do second vs second same community combo
        for community_combo in community_combos:
            second_values = dict_community_strain[('second', community_combo)]
            num_second = len(second_values)
            if num_second < 2:
                continue
            priority_effect_strength = []
            for i in range(100):
                first_a = np.random.choice(second_values)
                second_a = np.random.choice(second_values)
                strength = np.log10((first_a)/(second_a))
                priority_effect_strength.append(strength)
                success = f_sim_null.write(f'{strain},{strength},{community_combo},same,null_second,{spec_name_strain}\n')

            # do first vs first different community combo
        for community_combo in community_combos:
            for community_combo_2 in community_combos:
                if community_combo == community_combo_2:
                    continue
                first_values = dict_community_strain[('first', community_combo)]
                num_first = len(first_values)
                if num_first < 2:
                    continue
                priority_effect_strength = []
                for i in range(100):
                    first_a = np.random.choice(first_values)
                    second_a = np.random.choice(first_values)
                    strength = np.log10((first_a)/(second_a))
                    priority_effect_strength.append(strength)
                    success = f_sim_null.write(f'{strain},{strength},{community_combo},different,null_first,{spec_name_strain}\n')

            # do second vs second different community combo
        for community_combo in community_combos:
            for community_combo_2 in community_combos:
                if community_combo == community_combo_2:
                    continue
                second_values = dict_community_strain[('second', community_combo)]
                num_second = len(second_values)
                if num_second < 2:
                    continue
                priority_effect_strength = []
                for i in range(100):
                    first_a = np.random.choice(second_values)
                    second_a = np.random.choice(second_values)
                    strength = np.log10((first_a)/(second_a))
                    priority_effect_strength.append(strength)
                    success = f_sim_null.write(f'{strain},{strength},{community_combo},different,null_second,{spec_name_strain}\n')





# the dataframe contains the following columns all the counts of all the
# treatments across all the samples
# I want priority effect strength to be calculated for each
# strain in each community combo and I also want to know two values one, considering
# only limi_passed samples and another considering all the samples
# Also want the number of samples from which the priority effect strength was calculated
# instead of summarizing the counts as median etc, I want to take a random sample from a list
# of all for the numerator and denominator several times and then calculate the priority effect strength
# reporting the median, mean and standard deviation of the priority effect strength (repeat n=100 times)
# to calculate the priority effect strength, I will do log ratio of Ci_j and Cj_i where i is the focal strain
# and it is first arriver and j is the other strain and it is the second arriver
# alternative measure is, ln( C_ji / C_0i ) - ln( C_ij / C_i0 ) where C_0i and C_i0 are cases where i was 
# the only strain fed and considered equivalent but this is for implementing later..

# I will first create a dictionary of all the counts for each strain in each community combo
# I will then use this dictionary to calculate the priority
# effect strength for each strain in each community combo

strains = df_counts['Strain'].unique()
# for strain in strains:
#     if '_' in strain:
#         strains = np.append(strains, strain.split('_'))
community_combos = ['1', '2', '3', '4', '5', '6', '1_rep']
community_combos_all = ['1', '2', '3', '4', '5', '6', 'm1', 'm4', 'm6', 'm7', '1_rep']
# community_combos = df_counts['CommunityCombo'].unique()
arrivals = df_counts['arrival'].unique()

# create a dictionary of all the counts for each strain in each community combo
dict_counts_all = {}
dict_geq_passed_all = {}
dict_sample_geq_all = {}
for strain in strains:
    dict_counts_all[strain] = {}
    dict_geq_passed_all[strain] = {}
    dict_sample_geq_all[strain] = {}
    for community_combo in community_combos_all:
        dict_counts_all[strain][community_combo] = {}
        dict_geq_passed_all[strain][community_combo] = {}
        dict_sample_geq_all[strain][community_combo] = {}
        for arrival in ['first', 'second', 'only', 'dropout']:
            dict_counts_all[strain][community_combo][arrival] = {}
            dict_geq_passed_all[strain][community_combo][arrival] = {}
            dict_sample_geq_all[strain][community_combo][arrival] = {}
            df_temp = df_counts[(df_counts['Strain'] == strain) & (df_counts['CommunityCombo'] == community_combo) & (df_counts['arrival'] == arrival)]
            dict_counts_all[strain][community_combo][arrival] = df_temp['abs_counts'].values
            dict_geq_passed_all[strain][community_combo][arrival] = df_temp['geq_passed'].values
            dict_sample_geq_all[strain][community_combo][arrival] = df_temp['sample_geq'].values

# make df from dict of dicts
df_dict_counts_all = pd.DataFrame(columns = ['Strain', 'CommunityCombo', 'arrival', 'counts', 'geq_passed', 'sample_geq'])
for strain in strains:
    for community_combo in community_combos_all:
        for arrival in ['first', 'second', 'only', 'dropout']:
            for the_count, geq_passed, sample_geq in zip(dict_counts_all[strain][community_combo][arrival], dict_geq_passed_all[strain][community_combo][arrival], dict_sample_geq_all[strain][community_combo][arrival]):
                df_dict_counts_all = df_dict_counts_all._append({'Strain': strain, 'CommunityCombo': community_combo, 'arrival': arrival, 'counts': the_count, 'geq_passed': geq_passed, 'sample_geq': sample_geq}, ignore_index = True)

df_dict_counts_all.to_csv('03_Analysis/Results/tables/df_dict_counts_all.csv', index = False)

# create a dictionary of all the counts for each strain in each community combo
dict_counts = {}
dict_geq_passed = {}
dict_sample_geq = {}
for strain in strains:
    dict_counts[strain] = {}
    dict_geq_passed[strain] = {}
    dict_sample_geq[strain] = {}
    for community_combo in community_combos:
        dict_counts[strain][community_combo] = {}
        dict_geq_passed[strain][community_combo] = {}
        dict_sample_geq[strain][community_combo] = {}
        for arrival in ['first', 'second', 'only', 'dropout']:
            dict_counts[strain][community_combo][arrival] = {}
            dict_geq_passed[strain][community_combo][arrival] = {}
            dict_sample_geq[strain][community_combo][arrival] = {}
            df_temp = df_counts[(df_counts['Strain'] == strain) & (df_counts['CommunityCombo'] == community_combo) & (df_counts['arrival'] == arrival)]
            dict_counts[strain][community_combo][arrival] = df_temp['abs_counts'].values
            dict_geq_passed[strain][community_combo][arrival] = df_temp['geq_passed'].values
            dict_sample_geq[strain][community_combo][arrival] = df_temp['sample_geq'].values

# make df from dict of dicts
df_dict_counts = pd.DataFrame(columns = ['Strain', 'CommunityCombo', 'arrival', 'counts', 'geq_passed', 'sample_geq'])
for strain in strains:
    for community_combo in community_combos:
        for arrival in ['first', 'second', 'only', 'dropout']:
            for the_count, geq_passed, sample_geq in zip(dict_counts[strain][community_combo][arrival], dict_geq_passed[strain][community_combo][arrival], dict_sample_geq[strain][community_combo][arrival]):
                df_dict_counts = df_dict_counts._append({'Strain': strain, 'CommunityCombo': community_combo, 'arrival': arrival, 'counts': the_count, 'geq_passed': geq_passed, 'sample_geq': sample_geq}, ignore_index = True)

df_dict_counts.to_csv('03_Analysis/Results/tables/df_dict_counts.csv', index = False)

# now I will calculate the priority effect strength for each strain in each community combo
# I will take a random sample from the list of counts for the numerator and denominator
# and calculate the priority effect strength
# I will repeat this 100 times and report the median, mean and standard deviation of the priority effect strength
# for each strain in each community combo

def calculate_priority_effect_strength(by_simulation = False, write_bootstraps = False, n = 100, list_first = [], list_second = [], list_only = [], verbose = False, strain = '', community_combo = ''):
    list_first_filt = [x for x in list_first if str(x) != 'nan']
    list_second_filt = [x for x in list_second if str(x) != 'nan']
    list_only_filt = [x for x in list_only if str(x) != 'nan']
    if by_simulation:
        if write_bootstraps:
            with open('03_Analysis/Results/tables/bootstraps_priority_effect_strength.csv', 'w+') as f:
                f.write('Strain,CommunityCombo,arrival,second_arriver,priority_effect_strength\n')
        if verbose:
            print('Running by simulation')
        priority_effect_strength = []
        if list_only_filt == []:
            if len(list_first_filt) == 0 or len(list_second_filt) == 0:
                return {'median': np.nan, 'mean': np.nan, 'std': np.nan}
        else:
            if len(list_first_filt) == 0 or len(list_second_filt) == 0 or len(list_only_filt) == 0:
                return {'median': np.nan, 'mean': np.nan, 'std': np.nan}
        for i in range(n):
            first_a = np.random.choice(list_first_filt)
            second_a = np.random.choice(list_second_filt)
            if len(list_only) == 0:
                strength = np.log10((first_a)/(second_a))
                priority_effect_strength.append(strength)
                if write_bootstraps:
                    with open('03_Analysis/Results/tables/bootstraps_priority_effect_strength.csv', 'a') as f:
                        f.write(str(first_a) + ',' + str(second_a) + ',' + str(strength) + '\n')
                if verbose:
                    print('first_a: ', first_a, 'second_a: ', second_a, 'strength: ', strength)
                continue
            else:
                only_a_1 = np.random.choice(list_only_filt)
                only_a_2 = np.random.choice(list_only_filt)
                strength = np.log10((first_a)/(only_a_1)) - np.log10((second_a)/(only_a_2))
                priority_effect_strength.append(strength)
                if verbose:
                    print('first_a: ', first_a, 'second_a: ', second_a, 'only_a_1: ', only_a_1, 'only_a_2: ', only_a_2, 'strength: ', strength)

        return {'median': np.median(priority_effect_strength), 'mean': np.mean(priority_effect_strength), 'std': np.std(priority_effect_strength)}
    else:
        if verbose:
            print('Running without simulation')
        if list_only_filt == []:
            if len(list_first_filt) == 0 or len(list_second_filt) == 0:
                return {'median': np.nan, 'mean': np.nan}
        else:
            if len(list_first_filt) == 0 or len(list_second_filt) == 0 or len(list_only_filt) == 0:
                return {'median': np.nan, 'mean': np.nan}
        
        if len(list_only) == 0:
            mean_first = np.mean(list_first_filt)
            mean_second = np.mean(list_second_filt)
            priority_effect_strength_mean = np.log10((mean_first)/(mean_second))
            median_first = np.median(list_first_filt)
            median_second = np.median(list_second_filt)
            priority_effect_strength_median = np.log10((median_first)/(median_second))

            if verbose:
                print('mean')
                print('first: ', mean_first, 'second: ', mean_second)
                print('median')
                print('first: ', median_first, 'second: ', median_second)
                print('median: ', priority_effect_strength_median, 'mean: ', priority_effect_strength_mean)

        else:
            mean_first = np.mean(list_first_filt)
            mean_second = np.mean(list_second_filt)
            mean_only = np.mean(list_only_filt)
            priority_effect_strength_mean = np.log10((mean_first)/(mean_only)) - np.log10((mean_second)/(mean_only))
            median_first = np.median(list_first_filt)
            median_second = np.median(list_second_filt)
            median_only = np.median(list_only_filt)
            priority_effect_strength_median = np.log10((median_first)/(median_only)) - np.log10((median_second)/(median_only))
            
            if verbose:
                print('mean')
                print('first: ', mean_first, 'second: ', mean_second, 'only: ', mean_only)
                print('median')
                print('first: ', median_first, 'second: ', median_second, 'only: ', median_only)
                print('median: ', priority_effect_strength_median, 'mean: ', priority_effect_strength_mean)
        
        return {'median': priority_effect_strength_median, 'mean': priority_effect_strength_mean}

strains_to_skip = [
    # 'ESL0197', # no partner
    # 'ESL1028', # no partner
    'unknown', # no partner
    'ESL0353', # not in experiment
    'ESL0351' # not in experiment
]

dict_priority_effect_strength = {}
for strain in strains:
    if strain in strains_to_skip:
        continue
    dict_priority_effect_strength[strain] = {}
    for community_combo in community_combos:
        if community_combo == 'SW':
            continue
        dict_priority_effect_strength[strain][community_combo] = {}
        first_arriver = dict_counts[strain][community_combo]['first']
        second_arriver = dict_counts[strain][community_combo]['second']
        only_arriver = dict_counts[strain][community_combo]['only']
        len_first = len(first_arriver)
        len_second = len(second_arriver)
        len_only = len(only_arriver)
        dict_priority_effect_strength[strain][community_combo] = calculate_priority_effect_strength(list_first = first_arriver, list_second = second_arriver)
        # dict_priority_effect_strength[strain][community_combo] = calculate_priority_effect_strength(list_first = first_arriver, list_second = second_arriver, list_only = only_arriver)
        dict_priority_effect_strength[strain][community_combo]['num_samples_first'] = len_first
        dict_priority_effect_strength[strain][community_combo]['num_samples_second'] = len_second
        dict_priority_effect_strength[strain][community_combo]['num_samples_only'] = len_only

# write as a table
df_priority_effect_strength = pd.DataFrame(columns = ['Strain', 'CommunityCombo', 'priority_effect_strength_median', 'priority_effect_strength_mean', 'num_samples_first', 'num_samples_second'])
for strain in strains:
    if strain in strains_to_skip:
        continue
    for community_combo in community_combos:
        if community_combo == 'SW':
            continue
        df_priority_effect_strength = df_priority_effect_strength._append({'Strain': strain, 'CommunityCombo': community_combo, 'priority_effect_strength_median': dict_priority_effect_strength[strain][community_combo]['median'], 'priority_effect_strength_mean': dict_priority_effect_strength[strain][community_combo]['mean'], 'num_samples_first': dict_priority_effect_strength[strain][community_combo]['num_samples_first'], 'num_samples_second': dict_priority_effect_strength[strain][community_combo]['num_samples_second']}, ignore_index = True)

df_priority_effect_strength.to_csv('03_Analysis/Results/tables/df_priority_effect_strength_simple.csv', index = False)
# df_priority_effect_strength.to_csv('03_Analysis/Results/tables/df_priority_effect_strength_adv.csv', index = False)


dict_priority_effect_strength = {}
for strain in strains:
    if strain in strains_to_skip:
        continue
    dict_priority_effect_strength[strain] = {}
    for community_combo in community_combos:
        if community_combo == 'SW':
            continue
        dict_priority_effect_strength[strain][community_combo] = {}
        first_arriver = dict_counts[strain][community_combo]['first']
        second_arriver = dict_counts[strain][community_combo]['second']
        only_arriver = dict_counts[strain][community_combo]['only']
        len_first = len(first_arriver)
        len_second = len(second_arriver)
        len_only = len(only_arriver)
        dict_priority_effect_strength[strain][community_combo] = calculate_priority_effect_strength(list_first = first_arriver, list_second = second_arriver, by_simulation = True)
        # dict_priority_effect_strength[strain][community_combo] = calculate_priority_effect_strength(list_first = first_arriver, list_second = second_arriver, list_only = only_arriver, by_simulation = True)
        dict_priority_effect_strength[strain][community_combo]['num_samples_first'] = len_first
        dict_priority_effect_strength[strain][community_combo]['num_samples_second'] = len_second
        dict_priority_effect_strength[strain][community_combo]['num_samples_only'] = len_only

# write as a table
df_priority_effect_strength = pd.DataFrame(columns = ['Strain', 'CommunityCombo', 'priority_effect_strength_median', 'priority_effect_strength_mean', 'std', 'num_samples_first', 'num_samples_second'])
for strain in strains:
    if strain in strains_to_skip:
        continue
    for community_combo in community_combos:
        if community_combo == 'SW':
            continue
        df_priority_effect_strength = df_priority_effect_strength._append({'Strain': strain,
                                                                           'CommunityCombo': community_combo,
                                                                           'priority_effect_strength_median': dict_priority_effect_strength[strain][community_combo]['median'],
                                                                           'priority_effect_strength_mean': dict_priority_effect_strength[strain][community_combo]['mean'],
                                                                           'std': dict_priority_effect_strength[strain][community_combo]['std'],
                                                                           'num_samples_first': dict_priority_effect_strength[strain][community_combo]['num_samples_first'],
                                                                           'num_samples_second': dict_priority_effect_strength[strain][community_combo]['num_samples_second']}, ignore_index = True)

df_priority_effect_strength.to_csv('03_Analysis/Results/tables/df_priority_effect_strength_simple_by_simulation.csv', index = False)
# df_priority_effect_strength.to_csv('03_Analysis/Results/tables/df_priority_effect_strength_adv_by_simulation.csv', index = False)'
'''