# Study Title

PRIORITY EFFECTS AT THE STRAIN LEVEL IN HONEYBEE
GUT MICROBIOTA

# Abstract

Gut microbial communities often differ at the strain level even among closely related
individuals, but the ecological mechanisms driving this variation remain poorly understood.
One potential driver is priority effects, differences in the timing and order of microbial
colonization, which can lead to the assembly of distinct communities, even under similar
environmental conditions. Priority effects may specifically play an important role in shaping
microbial communities at the strain level, given that strains of the same species typically
occupy similar ecological niches. To test this, we examined gut microbiota assembly in
honeybees, in which age-matched nestmates are known to host similar microbial
communities at the species level but vary in strain composition. We sequentially colonized
microbiota-depleted honeybees with two distinct microbial communities, each composed of
the same twelve species but different strains. We found that firstcomer strains consistently
dominated the resulting communities, though the strength of these priority effects varied
among closely related strains and species. Dropping out individual strains from the firstcomer
community only partially improved the colonization success of latecomer conspecifics,
suggesting that priority effects also act across species boundaries. Our results underscore the
importance of priority effects for gut microbial community assembly at the strain level and in
shaping the specialized gut microbiota of bees.

This is the code base accompanying the study described above.

# Priority Effects Analysis

This repository contains code and scripts for analyzing priority effects in microbial communities. All analysis is contained within the 03_Analysis folder, which includes Python and R scripts for data parsing, statistical analysis, and figure generation.

03_Analysis/
├── infer_species_counts.py
├── Make_Figures.Rmd
├── parse_cd-hit_clusters.py
├── priority_effect_strength.py
├── qpcr_data_parse.R
├── summarise_genomes_16S.py
├── Utilities.R
├── Visualize_results.Rmd
├── visualize_sequencing_result.py

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