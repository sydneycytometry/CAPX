# CAPX
Cytometry analysis pipeline for large and complex datasets (CAPX).




# Five scripts used in the CAPX pipeline

### CAPX-1-preprocess

Pre-processing of data to prepare it for clustering and dimensionality reduction.



### CAPX-2-cluster

Clustering of data using FlowSOM and visualisation using tSNE, with differential downsampling options. Will also use the code from 'tSNEplots' to generate coloured tSNE plot images for each marker and each sample. Other scripts can be used on the output data files to give an identity to each cluster (ClusterPlots, SumTables, HeatMap_MFI).



### CAPX-3-annotate

Explorationg of expression of each marker on each cluster, addition of a population name.





# Extra scripts for additional functions

tSNEplots (https://github.com/sydneycytometry/tSNEplots)

ClusterPlots

FlexiPlots

SumTables

HeatMaps_MFI

HeatMaps_CellNum

AutoGraph (https://github.com/sydneycytometry/AutoGraph)

