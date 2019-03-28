# Cytometry analysis pipeline for large and complex datasets (CAPX).


The Cytometry Analysis Pipeline for large and compleX datasets (CAPX) is workflow for discovery analysis of high-dimensional cytometry data. Specifically the CAPX workflow is provided here as a set of scripts that bring together existing clustering and data visualisation tools into a single pipeline. This is used as an analysis pipeline at the Sydney Cytometry facility for high-dimensional flow and mass cytometry data. This same approach can be undertaken largely in FlowJo or cytofkit (see https://sydneycytometry.org.au/capx for more information).


Data (including datasets of tens of millions of cells) is clustered using FlowSOM (https://www.ncbi.nlm.nih.gov/pubmed/25573116), subsampled (with differential downsampling options), and visualisated using tSNE (https://lvdmaaten.github.io/publications/papers/JMLR_2014.pdf, https://www.ncbi.nlm.nih.gov/pubmed/23685480) via the rtsne package (https://cran.r-project.org/web/packages/Rtsne/index.html). Subsequently, this script will also use the code from 'tSNEplots' (https://github.com/sydneycytometry/tSNEplots) to generate coloured tSNE plot images for each marker and each sample. Other scripts can be used on the output data files to give an identity to each cluster (ClusterPlots, SumTables, HeatMap_MFI).


## Getting started ##


**Download CAPX script**


Go to 'releases' above (https://github.com/sydneycytometry/CAPX/releases) and download source code for the latest version.


**Download supporting scripts**


CytoTools (https://github.com/sydneycytometry/CytoTools) provides a number of supporting scripts, including: SumTables (generates a table summarising the analysed dataset: samples vs clusters -- number of cells per cluster per sample, MFI of each marker on each cluster in each sample, etc), HeatMaps (generates a heatmap for measuring the number of cells in each cluster in each sample or the MFI for each marker on each cluster, per sample. Includes 'fold-change' visualisation options), ClusterPlots (automated generation of coloured tSNE plots showing clusters), and FlexiPlots (a script with adjustable parameters for visualising plots).

AutoGraph (https://github.com/sydneycytometry/AutoGraph) can be used to automatically plot dot plots to compare measurements (cells per tissue, median fluorescence intensity (MFI) etc) of each cluster/population between groups.



**Protocols**


For usage instructions, please see https://github.com/sydneycytometry/CAPX/wiki/CAPX-usage-protocol.


## Citation ##


If you use these scripts in your work, please cite this github using the information below. You can cite the specific version that you used in your work.


Ashhurst, T. M. (2018). Cytometry Analysis Pipeline for large and compleX dataests v2.5. GitHub repository. DOI: TBC, repository: https://github.com/sydneycytometry/CAPX.


**Version history**


See https://github.com/sydneycytometry/CAPX/releases.


## References ##


**Packages used**


'plyr'

'data.table'

'rstudioapi'

'devtools'

'FlowSOM'

'Rtsne'

'ggplot2'

'colorRamps'

'ggthemes'

'scales'

'flowCore'

'Biobase'

'flowViz'


**Specific references**


Sofie Van Gassen, Britt Callebaut and Yvan Saeys (2017). FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data. http://www.r-project.org, http://dambi.ugent.be.

L.J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research 9(Nov):2579-2605, 2008.

L.J.P. van der Maaten. Accelerating t-SNE using Tree-Based Algorithms. Journal of Machine Learning Research 15(Oct):3221-3245, 2014.

Jesse H. Krijthe (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using a Barnes-Hut Implementation, URL: https://github.com/jkrijthe/Rtsne

