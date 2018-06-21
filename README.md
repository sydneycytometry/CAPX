# Cytometry analysis pipeline for large and complex datasets (CAPX).


The Cytometry Analysis Pipeline for large and compleX datasets (CAPX) is a set of scripts that brings together existing clustering and data visualisation tools into a single pipeline. This is used as a standard analysis pipeline at the Sydney Cytometry facility for high-dimensional flow and mass cytometry data. 

This analysis pipleine uses three main scripts, described below. 


### Three main scripts used in the CAPX pipeline


**CAPX-1-preprocess**

Pre-processing of data to prepare it for clustering and dimensionality reduction. Output is a 'clustering ready' .csv file (and .fcs file).



**CAPX-2-cluster**

Data (including datasets of tens of millions of cells) is clustered using FlowSOM (https://www.ncbi.nlm.nih.gov/pubmed/25573116), subsampled (with differential downsampling options), and visualisated using tSNE (https://lvdmaaten.github.io/publications/papers/JMLR_2014.pdf, https://www.ncbi.nlm.nih.gov/pubmed/23685480) via the rtsne package (https://cran.r-project.org/web/packages/Rtsne/index.html). Subsequently, this script will also use the code from 'tSNEplots' (https://github.com/sydneycytometry/tSNEplots) to generate coloured tSNE plot images for each marker and each sample. Other scripts can be used on the output data files to give an identity to each cluster (ClusterPlots, SumTables, HeatMap_MFI).



**CAPX-3-annotate**

Following the explorationg marker expression on the cells and clusters generated, this script gives the ability to add a population name to each cluster, and gives the ability to 'merge' redundant clusters if desired.



### Extra scripts for additional functions


**tSNEplots** (https://github.com/sydneycytometry/tSNEplots) can be used to automatically create a coloured tSNE plot for every marker and sample (and group).


**AutoGraph** (https://github.com/sydneycytometry/AutoGraph) can be used to automatically plot dot plots to compare measurements (cells per tissue, median fluorescence intensity (MFI) etc) of each cluster/population between groups.


### The following scripts will be included in the next release of CAPX:


**ClusterPlots** - automated generation of coloured tSNE plots showing clusters.  


**FlexiPlots** - a script with adjustable parameters for visualising plots

**SumTables** - generates a table summarising the analysed dataset: samples vs clusters -- number of cells per cluster per sample, MFI of each marker on each cluster in each sample, etc.


**HeatMaps_MFI** - generates a heatmap for measuring MFI for each marker on each cluster, per sample. Includes 'fold-change' visualisation options.


**HeatMaps_CellNum** - generates a heatmap for measuring the number of cells in each cluster in each sample. Includes 'fold-change' visualisation options.



