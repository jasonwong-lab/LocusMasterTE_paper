# lasTEq_paper
Scripts for figures and results in lasTEq paper


### Figures_Rscripts
Rscripts used to generate main figures.


### Data_lists
Containing supplementary information about data used.

1. Used_SG-Nex_accession_list.xlsx: Excel file contains information about used samples in SG-NEx study with sample ID.
2. Processed_long_read.xlsx: Links directed to "processed long-read TPM counts" readily used as long-read input in lasTEq.\
6 samples with 5 different cell lines. Two HCT116 samples exist generated from our lab and SG-NEx study.


### preprocessing_commands.sh
Important commands use to preprocess samples. Largely contains three major parts.
1. Processing long-read: Minimap2 and FeatureCounts
2. Aligning short-read: STAR
3. Simulation commands: Spanki
