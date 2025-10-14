# Study Title

PRIORITY EFFECTS AT THE STRAIN LEVEL IN HONEYBEE
GUT MICROBIOTA

# Abstract

Gut microbial communities differ at the strain level among individual hosts, but the mechanisms driving this variation remain poorly understood. One potential factor is priority effects, a process in which differences in the timing and order of microbial colonization influence subsequent community assembly ("first come, first served" dynamics). We hypothesize that these stochastic/neutral processes act at the strain level within species due to the niche overlap of closely related bacteria, and may drive community divergence even under similar environmental conditions. To test these predictions, we sequentially colonized microbiota-depleted honeybees with two distinct microbial communities composed of the same twelve species but different strains, ensuring that individuals shared species-level composition but differed at the strain level. We found that firstcomer strains consistently dominated the resulting communities, suggesting strong priority effects. Dropout experiments in which the firstcomer strain of a species was removed led to only partial increase in the colonization success of the conspecific latecomer, suggesting that inter-species interactions also play a role. Our results underscore the importance of priority effects for gut microbial community assembly at the strain level and in shaping the specialized gut microbiota of bees.

This is the code base accompanying the study described above.

# Priority Effects Analysis

This repository contains code and scripts for analyzing priority effects in microbial communities. All analysis is contained within the 03_Analysis folder, which includes Python and R scripts for data parsing, statistical analysis, and figure generation.

03_Analysis/ \
├── infer_species_counts.py \
├── Make_Figures.Rmd \
├── parse_cd-hit_clusters.py \
├── priority_effect_strength.py \
├── qpcr_data_parse.R \
├── summarise_genomes_16S.py \
├── Utilities.R \
├── Visualize_results.Rmd \
├── visualize_sequencing_result.py \

## Dependencies
Python 3.x
Common packages: pandas, numpy, matplotlib, seaborn (check each script for specifics)
R (>= 4.0)
Common packages: tidyverse, ggplot2, knitr, rmarkdown

## Project Structure
*Data parsing*: parse_cd-hit_clusters.py, qpcr_data_parse.R, summarise_genomes_16S.py
*Analysis*: infer_species_counts.py, priority_effect_strength.py
*Visualization*: visualize_sequencing_result.py, Make_Figures.Rmd, Visualize_results.Rmd, Utilities.R

The repository is indended for documentation purposes. If you would like to use any scripts or their parts, they are found within the subdirectories. To run the entire pipeline ensure that all mentioned packes are installed and the clone this repository and run the respective scripts after modifying file paths as needed. For more information do not hesitate to reach out!