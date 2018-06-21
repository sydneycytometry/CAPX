# Cytometry analysis pipeline for large and complex datasets (CAPX).


The Cytometry Analysis Pipeline for large and compleX datasets (CAPX) is a set of scripts that brings together existing clustering and data visualisation tools into a single pipeline. This is used as a standard analysis pipeline at the Sydney Cytometry facility for high-dimensional flow and mass cytometry data. 

This analysis pipleine uses three main scripts, described below. 


## Three main scripts used in the CAPX pipeline

### CAPX-1-preprocess

Pre-processing of data to prepare it for clustering and dimensionality reduction. Output is a 'clustering ready' .csv file (and .fcs file).



### CAPX-2-cluster

Data(including datasets of tens of millions of cells) using FlowSOM (https://www.ncbi.nlm.nih.gov/pubmed/25573116) and visualisated using tSNE (https://lvdmaaten.github.io/publications/papers/JMLR_2014.pdf) via the rtsne package (https://cran.r-project.org/web/packages/Rtsne/index.html), with differential downsampling options. This will also use the code from 'tSNEplots' (https://github.com/sydneycytometry/tSNEplots) to generate coloured tSNE plot images for each marker and each sample. Other scripts can be used on the output data files to give an identity to each cluster (ClusterPlots, SumTables, HeatMap_MFI).



### CAPX-3-annotate

Following the explorationg marker expression on the cells and clusters generated, this script gives the ability to add a population name to each cluster, and gives the ability to 'merge' redundant clusters if desired.



## Extra scripts for additional functions

tSNEplots (https://github.com/sydneycytometry/tSNEplots)

ClusterPlots

FlexiPlots

SumTables

HeatMaps_MFI

HeatMaps_CellNum

AutoGraph (https://github.com/sydneycytometry/AutoGraph)

