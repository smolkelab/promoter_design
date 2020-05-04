# Promoter Design
Code and analysis scripts developed for "Model-driven generation of artificial yeast promoters" (Kotopka BJ, Smolke CD. Model-driven generation of artificial yeast promoters. Nat Commun. 2020;11(1):2113.).

## Commit history
Commit a8dc340 corresponds to the BioRxiv preprint (https://www.biorxiv.org/content/10.1101/748616v). 

## Installation/Setup
Python scripts were run using Python 2.7 on Ubuntu 16.04 (xenial); R scripts (in "Figures/RStudio" and "Code") were run within RStudio 1.1.463, using R 3.5.1, on Windows 10. For Python requirements see 'requirements.txt' in the main directory; external R dependencies are data.table, ggplot2, MASS, viridis, UpSetR, gplots, and beanplot. tensorflow-gpu highly recommended for model training.

## What does this do?
Subdirectories 'GPD' and 'ZEV' extract sequences and promoter activity estimates for the GPD and ZEV libraries from FASTQ files generated in FACS-seq experiments (main text Fig. 1); additionally, 'GPD/models' and 'ZEV/models' train neural network models to predict promoter activity. 'joined' combines data from these experiments and trains a model ensemble on the pooled data (Fig. 2). 'seq_design' then uses these models to design sets of novel promoters to fulfill user-specified objectives, and 'oligo_design' takes these designs as input and generates a FASTA file containing oligos which can be ordered from Twist Bioscience to synthesize the designed promoters, along with PCR primers to selectively amplify a particular promoter set. Activities of these promoters were measured by FACS-seq; data analysis for this experiment is in 'design_testing' (Fig. 3). 'model_interpretation' runs in silico mutagenesis experiments to identify trends in sequences which the model identifies as high-performing (Fig. 5). Subdirectories 'mean_extraction' and 'models' contain code invoked elsewhere to analyze FACS-seq results or train models, respectively. Lastly, 'Figures/RStudio' contains R scripts for generating figures and data analysis (including all analysis for Fig. 4 - individual-sequence flow cytometry validation of FACS-seq results), and scripts in 'Code' contain functions used in 'Figures/RStudio'. These scripts were also used to perform statistical analyses whose results are described in the main text.

### Practical notes for re-running figure scripts
To generate figures, obtain the Supplementary Data from Zenodo (https://zenodo.org/record/3735426). You will need to store it in the top level of your D: drive, or edit paths in scripts as appropriate. Create directories 'Figures/PNGs' and 'Figures/Final PNGs'; then source the corresponding R scripts in 'Figures/RStudio'. Outputs will appear in 'PNGs' or 'Final PNGs'; for figures whose output appears in PNGs, I used Abode Illustrator to arrange subpanels and add some axis labels. Illustrator was also used to make schematics, such as Figure 1A; Figure S23 was generated using a Benchling plasmid map.

## Raw data deposition
NGS data deposited to NCBI GEO as Series GSE135464.
Other data available from Zenodo (https://zenodo.org/record/3735426).

## Acknowledgements
Code in mean_extraction/hgbrian_GPU_NW.py is adapted from work by Brian Naughton (https://github.com/hgbrian/nw_align, commit 54221ee). 'bdfacs_setpolys.m' (used to set gates for FACS-seq experiments) by Brent Townshend (https://github.com/btownshend/TwoColor, eddf69e). All code used with permission of creators.
