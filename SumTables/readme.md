SumTables contains two scripts: SumTables_CellNo and SumTables_MFI. 

SumTables_CellNo will read a .csv file containing an entire experiment (i.e. all samples in one file). Assuming the file contains keywords so that each cell can be attributed to a sample, the script will generate a table showing cluster or population vs sample. Each value will represent the number of cells from each cluster or population in each sample.

SumTables_MFI will read a .csv file containing a single sample (or a mixture of samples). The script will produce a table showing cluster or population vs marker, and each value will represent the MFI of each marker for each population or cluster.
