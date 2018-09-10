### Cytometry Analysis Pipeline for Complex Datasets (CAPX) part 2 - cluster
    
    # Thomas Ashhurst
    # 2018-05-30
    # thomas.ashhurst@sydney.edu.au
    # www.github.com/SydneyCytometry/CAPX

### Summary
    
    # 1 INSTALL AND LOAD PACKAGES
    # 2 USER INPUT -- DATA PREPARATION
    
    ### END USER INPUT ###
     
    # 3 FlowSOM
    # 4 tSNE
    # 5 tSNEplots

###################################################### 1. INSTALL AND-OR LOAD PACKAGES ###################################################### 

    ### 1.1 - Install packages (if not already installed)
        if(!require('flowCore')) {install.packages('flowCore')}
        if(!require('plyr')) {install.packages('plyr')}
        if(!require('Biobase')) {install.packages("Biobase")}
        if(!require('Rtsne')) {install.packages("Rtsne")}
        if(!require('devtools')){install.packages("devtools")}
        if(!require('Rphenograph')){devtools::install_github("JinmiaoChenLab/Rphenograph")}
        if(!require('flowViz')) {source("https://bioconductor.org/biocLite.R") 
          biocLite('flowViz')}
        if(!require('FlowSOM')) {install.packages('FlowSOM')} 
        if(!require('data.table')) {install.packages('data.table')}
        if(!require('rstudioapi')) {install.packages('rstudioapi')}
        if (!require("ggplot2")){install.packages("ggplot2")} # for plotting tSNE graphs
        if (!require("colorRamps")){install.packages("colorRamps")} # for colour scheme management
        if (!require("ggthemes")){install.packages("ggthemes")} # for plot themes
        if (!require("scales")){install.packages("scales")} # for re-scaling if necessary
    
    ### 1.2 Load packages       
        library('flowCore')
        library('plyr')
        library('Biobase')
        library('Rtsne')
        library('devtools')
        library("Rphenograph") # not actually needed
        library('flowViz')    
        library('FlowSOM')
        library('data.table')
        library('rstudioapi')
        library('ggplot2')
        library('colorRamps') # for colour scheme management
        library('ggthemes') # for plot themes
        library('scales') # for re-scaling, only if necessary
    
    ### 1.3 - Set working directory and assign as 'PrimaryDirectory'
    
        ## In order for this to work, a) rstudioapi must be installed and b) the location of this .r script must be in your desired working directory
            dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
            setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
            getwd()
            PrimaryDirectory <- getwd()
            PrimaryDirectory
        
        ## Use this to manually set the working directory
            #setwd("/Users/Tom/Desktop/Experiment")                          # Set your working directory here (e.g. "/Users/Tom/Desktop/") -- press tab when selected after the '/' to see options
            #getwd()                                                         # Check your working directory has changed correctly
            #PrimaryDirectory <- getwd()                                     # Assign the working directory as 'PrimaryDirectory'
            #PrimaryDirectory
        

###################################################### 2. USER INPUT - DATA PREPARATION ###################################################### 
       
    ### 2.1 - Specify options
        
        ## Use to list the .fcs or .csv files in the working directory
            list.files(path=PrimaryDirectory, pattern = ".fcs")     # see a list of FCS files
            list.files(path=PrimaryDirectory, pattern = ".csv")     # see a list of CSV files
        
        ## Specify general options
            File_Type             <- ".csv"           # Specify the type of input file. Use ".fcs" OR ".csv", in lower case      # fcs files not current supported
            FileNos               <- "single"         # Can be "single" for one file, or "multi" for multiple files. If "single" the only the file desired for analysis should be in the working directory, otherwise the script will use the first file it finds.
        
        ## New name for data
            Data_name             <- "demo_data"      # If one single file, what is the name of the file (not required if reading 'multi' files from a directory)
            Run_Groups            <- 1                # An option to write .csv and .fcs files for each group. Turn OFF (0) if you have no group distinguishing keywords
        
        ## FlowSOM options
            Run_FlowSOM           <- 1                # Option to run FlowSOM. Enter 1 for yes, or 0 for no
            FlowSOM_k             <- 20               # Number of metaclusters to derive from FlowSOM (can also be selected automatically, but this often does not perform well)
            FlowSOM_seed          <- 42               # Seed for running FlowSOM
            FlowSOM_meta_seed     <- 42               # Seed for metaclustering
            
            FlowSOM_clus_name     <- "FlowSOM_original_cluster" # Label for metacluster parameter name
            FlowSOM_meta_name     <- "FlowSOM_meta_cluster" # Label for metacluster parameter name
        
        ## tSNE options
            Run_tSNE              <- 1                # Option to run tSNE. Enter a 1 to run tSNE, or 0 to skip 
            tSNE_dwnsmpl          <- "specific"       # can be "uniform" or "specific". 'uniform' means the entire dataset will be evenly downsampled to a target total number of cells. 'specific' means each sample is downsampled to a specific value.
    
            tSNE_nsub_seed        <- 42               # Seed for downsampling for tSNE
            tSNE_seed             <- 42               # Seed for running tSNE
            tSNE_X_name           <- "tSNE1"          # Label for tSNE1 parameter
            tSNE_Y_name           <- "tSNE2"          # Label for tSNE2 parameter

            # IF tSNE_dwnsmpl = "uniform"
            tSNE_nsub             <- 70000            # ONLY used if tSNE_dwnsmpl = "uniform". Specify the total downsample size for the data run in tSNE.
            
        ## tSNEplot options
            Run_tSNEplots         <- 1              # 1 to run tSNEplots and generate coloured tSNE maps
    
        ## ClusterPlot options -- can ONLY be run if both FlowSOM and tSNE have been run
            Run_ClusterPlots      <- 1              # 1 to generate coloured cluster plots and number plots
            cluspl.n.sub          <- 1000           # Define subsampling target and seed for NUMBER PLOT ONLY (default = 1000 per plot, and a set of 42)
            cluspl.n.sub.seed     <- 42 
            
        
    ### 2.2 - Read files into a list of dataframes (run all of 2.2 at once, no user input)
        
        ## Loop for single files
            if(FileNos == "single"){
              if(File_Type == ".csv"){
                file.single <- list.files(path=PrimaryDirectory, pattern = ".csv")
                data <- read.csv(file.single[1], check.names = FALSE)
              }
              if(File_Type == ".fcs"){
                file.single <- list.files(path=PrimaryDirectory, pattern = ".fcs")
                data <- exprs(read.FCS(FileNames[1], transformation = FALSE)) 
                data <- parameters(read.FCS(File)) 
              }
            }
      
        ## Loop for multiple files
            if(FileNos == "multi"){
              DataList=list() # Creates and empty list to start 
              if (File_Type == ".csv"){
                FileNames <- list.files(path=PrimaryDirectory, pattern = ".csv")
                for (File in FileNames) { # Loop to read files into the list
                  tempdata <- read.csv(File, check.names = FALSE)
                  File <- gsub(".csv", "", File)
                  DataList[[File]] <- tempdata
                }
              }
              if (File_Type == ".fcs"){
                FileNames <- list.files(path=PrimaryDirectory, pattern = ".fcs")
                ParaList=list()
                for (File in FileNames) { # Loop to read files into the list
                  tempdata <- exprs(read.FCS(FileNames[1], transformation = FALSE)) 
                  File <- gsub(".fcs", "", File)
                  DataList[[File]] <- tempdata
                  ParaList[[File]] <- parameters(read.FCS(File)) 
                  ParaList[[File]]
                }
              }
              data <- rbindlist(DataList) ## Concatenate into one large data frame -- what if they have a column conflict??
              data
            }
         
        ## Review data
            dim(data)
            head(data)
       
        ## Choose the column that defines the sample NAMES
            as.matrix(names(data))
            samp.col <- "SampleName"
            head(data[samp.col]) # check that this is actually the sample names
            
        ## Choose the column that defines the GROUPS NAMES
            as.matrix(names(data))
            grp.col <- "GroupName"
            head(data[grp.col]) # check that this is actually the group names
            
        
    ### 2.3 - Choose markers for clustering
        
        ## Review column names
            names(data)                            # show data with headings
            ColumnNames <- unname(colnames(data)) # assign reporter and marker names (column names) to 'ColumnNames'
            as.matrix(ColumnNames) # view the column 'number' for each parameter
            
        ## Specify column numbers to be used for clustering
            ClusteringColNos <- c(1,3:8,10,12:14,16:24,26)
            ClusteringColNos
            
            ClusteringCols <- ColumnNames[ClusteringColNos] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10
           
            ClusteringCols  # check that the column names that appear are the ones you want to analyse
            ColumnNames[-ClusteringColNos] # Check which columns are being EXCLUDED!   
            
    ### 2.4 - Downsample targets for tSNE analysis
        
        ## Create a list of the number of rows 'cells' in each sample 
            SampleNames <- unique(data[[samp.col]]) # unique sample names
            
            samp.name <- data[samp.col]
            grp.row.list = list()
            
            for(i in SampleNames){
              grp.row.list[i] <- nrow(subset(data, samp.name == i))
            }
            
        ## Check the number of cells in each file
            as.matrix(grp.row.list) 
        
        ## IF USING tSNE_dwnsmpl = "uniform" -- then the following lines are NOT relevant and do not need to be run
            # --> SKIP to 3.1
            
        ## IF USING tSNE_dwnsmpl = "specific" -- each sample has a unique downsample target
            # --> Set downsample targets -- can use the spreadsheet available from www.github.com/sydneycytometry/capx
            DownSampleTargets <- c(1000,
                                   1000,
                                   1000,
                                   1000,
                                   1000,
                                   1000,
                                   1000,
                                   1000,
                                   1000,
                                   1000,
                                   1000,
                                   1000)
            
            as.matrix(DownSampleTargets) # value of each target
            sum(DownSampleTargets) # sum of all targets -- the total number of samples to appear in tSNE analysis
            
        ## Review you have the right number of DownSampleTargets entires
            length(DownSampleTargets)
            length(DownSampleTargets) == length(1:length(SampleNames)) # Should return TRUE if the length is the same
            

############################################################################################################################
###################################################### END USER INPUT ######################################################         
############################################################################################################################        
                            
                   
###################################################### 3. DATA PREP ######################################################         

    ### 3.1. Save a copy of 'data' and a list of columns for clustering
        
        ## 
        AllSampleNames <- unique(data[[samp.col]])
        AllGroupNames <- unique(data[[grp.col]])
        
        ## Directory creation                    
        setwd(PrimaryDirectory)
        dir.create("Output_cluster", showWarnings = FALSE)
        setwd("Output_cluster")
        OutputDirectory <- getwd()
        OutputDirectory
        
        ## Export data in CSV format
        data_csv <- paste(paste0(Data_name), "_used_for_clustering", ".csv", sep = "")
        write.csv(x = data, file = data_csv, row.names=FALSE)
        
        ## Export list of columns for clustering
        ClusteringCols_csv <- paste(paste0("Columns_for_Clustering"), ".csv", sep = "")
        write.csv(x = ClusteringCols, file = ClusteringCols_csv)

                
###################################################### 4. FlowSOM ######################################################         
                
      if(Run_FlowSOM == 1) {
        
        ### 4.1 - Create output directories
        
        setwd(OutputDirectory)
        dir.create("Output_FlowSOM", showWarnings = FALSE)
        
        setwd("Output_FlowSOM")
        dir.create("Output_FlowSOM_info", showWarnings = FALSE)
        dir.create("Output_FlowSOM_data", showWarnings = FALSE)
        
        ### 4.2 - Prep data for FlowSOM
        
        ## Check data and data column names
        head(data)
        dimnames(data)[[2]]
        
        ## Create FCS file metadata - column names with descriptions
        metadata <- data.frame(name=dimnames(data)[[2]], desc=paste('column',dimnames(data)[[2]],'from dataset'))
        
        ## Create FCS file metadata - ranges, min, and max settings -- by default, they are commented out (adjust ranges manually in FlowJo)
        #metadata$range <- apply(apply(data,2,range),2,diff) # throws an error because of word entry -- hopefully is ok
        metadata$minRange <- apply(data,2,min)
        metadata$maxRange <- apply(data,2,max)
        
        ## Create flowframe with data
        data.ff <- new("flowFrame",
                       exprs=as.matrix(data), # in order to create a flow frame, data needs to be read as matrix
                       parameters=AnnotatedDataFrame(metadata))
        
        head(exprs(data.ff))
        
        data_FlowSOM <- data.ff
        
        # choose markers for FlowSOM analysis
        FlowSOM_cols <- ClusteringCols
        
        ### 4.3. - Run FlowSOM              
        
        ## set seed for reproducibility
        set.seed(FlowSOM_seed)
        
        ## run FlowSOM (initial steps prior to meta-clustering)
        FlowSOM_out <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
        FlowSOM_out <- FlowSOM::BuildSOM(FlowSOM_out, colsToUse = FlowSOM_cols)
        FlowSOM_out <- FlowSOM::BuildMST(FlowSOM_out)
        
        ### some warnings will be returned because of the 'SampleName' and 'GroupName' entries
        
        ## Optional visualization
        
        #FlowSOM::PlotStars(FlowSOM_out) # won't plot names if arial font problems still exists, # if 'sample name' is in there, can't plot
        
        #set.seed(42)
        #FlowSOM::PlotStars(FlowSOM_out,view="tSNE")
        
        #print(colnames(FlowSOM_out$map$medianValues))
        #FlowSOM::PlotMarker(FlowSOM_out,"BUV395.CD11b")
        #FlowSOM::PlotNumbers(UpdateNodeSize(FlowSOM_out,reset=TRUE))
        
        ## extract cluster labels (pre meta-clustering) from output object
        labels_pre <- FlowSOM_out$map$mapping[, 1]
        labels_pre
        length(labels_pre)
        res.original <- labels_pre
        
        ## run meta-clustering
        FlowSOM_out_meta <- FlowSOM::metaClustering_consensus(FlowSOM_out$map$codes, k = FlowSOM_k, seed = FlowSOM_meta_seed)
        
        # note: In the PREVIOUS version of FlowSOM, the meta-clustering function 
        # FlowSOM::metaClustering_consensus() does not pass along the seed argument 
        # correctly, so results are not reproducible. We use the internal function 
        # ConsensusClusterPlus::ConsensusClusterPlus() to get around this. 
        
        # seed <- 1234
        # out <- ConsensusClusterPlus::ConsensusClusterPlus(t(out$map$codes), maxK = FlowSOM_kvalue, seed = seed)
        # out <- out[[FlowSOM_kvalue]]$consensusClass
        
        # However, this
        # HAS BEEN fixed in the next update of FlowSOM (version 1.5); then the following 
        # (simpler) code can be used instead:
        
        ## Optional visualisation
        # FlowSOM::PlotMarker(FlowSOM_out,"BUV395.CD11b", backgroundValues = as.factor(FlowSOM_out_meta))
        
        
        ## extract META (?) cluster labels from output object
        labels <- FlowSOM_out_meta[labels_pre]
        
        ## summary of cluster sizes and number of clusters
        table(labels)
        length(table(labels))
        
        ## save ORIGINAL cluster labels
        res.original <- data.frame("labels_pre" = labels_pre)
        colnames(res.original)[grepl('labels_pre',colnames(res.original))] <- FlowSOM_clus_name
        
        dim(data)
        dim(res.original)
        head(res.original)
        
        ## save META cluster labels
        res.meta <- data.frame("labels" = labels)
        colnames(res.meta)[grepl('labels',colnames(res.meta))] <- FlowSOM_meta_name
        
        dim(data)
        dim(res.meta)
        head(res.meta)
        
        
        ### 4.4. - Add FlowSOM cluster numbers to data, then export data to .csv and .fcs 
        
        ## Add FlowSOM original cluster number to data
        output_data <- cbind(data, res.original)
        data <- output_data
        head(data)
        dim(data)            
        
        ## Add FlowSOM meta cluster number to data
        output_data <- cbind(data, res.meta) 
        data <- output_data
        head(data)
        dim(data)
        
        ### testing: could make a 'list' with data and res, and use rbindlist(list) to merge them much faster than cbind
        #   test <- rbindlist(list(data, res))
        #   dim(test)
        #   not designed for adding columns, but for adding rows -- could use fill = TRUE to add the missing cols with NA's, then remove them, but seems messy
        
        ## Save cluster labels to CSV
        setwd(OutputDirectory)
        setwd("Output_FlowSOM/Output_FlowSOM_info")
        getwd()
        
        write.table(res.original, file = "original_cluster_labels_FlowSOM.txt", row.names = FALSE, quote = FALSE, sep = "\t")
        write.table(res.meta, file = "meta_cluster_labels_FlowSOM.txt", row.names = FALSE, quote = FALSE, sep = "\t")
        
        ## Save FlowSOM_out as a .csv file
        capture.output(FlowSOM_out, file = "FlowSOM_output.txt")
        
        # These don't work because FlowSOM_out can't be coerced into a dataframe  
        # write.csv(as.list(FlowSOM_out), file = paste(paste0(Data_name), "_FlowSOM_details.csv", sep = ""))
        # lapply(as.list(FlowSOM_out), function(x) write.table( data.frame(x), 'test.csv'  , append= T, sep=',' ))
        # capture.output(summary(FlowSOM_out), file = "FlowSOM_output.txt")
        # could try: write.list(FlowSOM_out, file = "test")
        
        ## Write to .csv
        setwd(OutputDirectory)
        setwd("Output_FlowSOM/Output_FlowSOM_data")
        getwd()
        
        write.csv(data, file = paste(paste0(Data_name), "_with_FlowSOM.csv", sep = ""), row.names=FALSE)
        
        ## Write to .fcs
        metadata <- data.frame(name=dimnames(data)[[2]],desc=paste('column',dimnames(data)[[2]],'from dataset'))
        metadata ## Double check metadata
        
        ## Create FCS file metadata - ranges, min, and max settings
        #metadata$range <- apply(apply(data,2,range),2,diff)
        metadata$minRange <- apply(data,2,min)
        metadata$maxRange <- apply(data,2,max)
        
        data.ff <- new("flowFrame",exprs=as.matrix(data), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
        head(data.ff)
        write.FCS(data.ff, paste0(Data_name, "_with_FlowSOM", ".fcs")) ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
        
        head(exprs(data.ff))                  
        
        ### 4.5. - Export each group as a separate .csv and .fcs file             
        
        if(Run_Groups == 1){
            for(a in AllGroupNames){
              data_subset <- subset(data, data[[grp.col]] == a)
              dim(data_subset)
              
              ## write .csv
              write.csv(data_subset, file = paste(Data_name, "_", a, "_with_FlowSOM", ".csv"), row.names=FALSE)
              
              ## write .fcs
              metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
              metadata
              
              ## Create FCS file metadata - ranges, min, and max settings
              #metadata$range <- apply(apply(data_subset,2,range),2,diff)
              metadata$minRange <- apply(data_subset,2,min)
              metadata$maxRange <- apply(data_subset,2,max)
              
              data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
              head(data_subset.ff)
              write.FCS(data_subset.ff, paste0(Data_name, "_", a, "_with_FlowSOM", ".fcs"))
            }
        }
        
        ### 4.6. - Export each sample as a separate .csv and .fcs file             
        
        for (a in AllSampleNames) {
          
          data_subset <- subset(data, data[[samp.col]] == a)
          dim(data_subset)
          
          ## write .csv
          write.csv(data_subset, file = paste(Data_name, "_", a, "_with_FlowSOM", ".csv"), row.names=FALSE)
          
          ## write .fcs
          metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
          metadata
          
          ## Create FCS file metadata - ranges, min, and max settings
          #metadata$range <- apply(apply(data_subset,2,range),2,diff)
          metadata$minRange <- apply(data_subset,2,min)
          metadata$maxRange <- apply(data_subset,2,max)
          
          data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
          head(data_subset.ff)
          write.FCS(data_subset.ff, paste0(Data_name, "_", a, "_with_FlowSOM", ".fcs"))
        }
        
        setwd(PrimaryDirectory)
        getwd()
        
      }


###################################################### 5. tSNE ######################################################                       
                
    if(Run_tSNE == 1){
      
      ### 5.1 - Create output directories
      
      setwd(OutputDirectory)
      dir.create("Output_tSNE", showWarnings = FALSE)
      setwd("Output_tSNE")
      dir.create("Output_tSNE_info", showWarnings = FALSE)
      dir.create("Output_tSNE_data", showWarnings = FALSE)
      setwd(OutputDirectory)
      
      ### 5.2 - Subsample full dataset (including FlowSOM cluster no, if FlowSOM was run)
      
      ## take res dataframe 'output_data'
      ## subsample to target no 'output_data_subsample'
      ## reduce to data_specific 'output_data_subsample_specific'       
      
      ## output_data
      head(data)
      dim(data)
      
      ## Subsampling
      
      if(tSNE_dwnsmpl == "uniform"){
        dim(data)
        set.seed(tSNE_nsub_seed)
        data_subsampled <- data[sample(1:nrow(data), tSNE_nsub), ]
        dim(data_subsampled)
        data_subsampled <- as.data.frame(data_subsampled)
      }
      
      if(tSNE_dwnsmpl == "specific"){
        
        data_subsampled <- data.frame()
        for (i in c(1:length(AllSampleNames))) {
          nam <- AllSampleNames[i]
          nsub <- DownSampleTargets[i] 
          data.temp <- subset(data, data[[samp.col]] == nam) # works
          nrow(data.temp)
          set.seed(tSNE_nsub_seed)
          data.temp <- data.temp[sample(1:nrow(data.temp), nsub), ]
          nrow(data.temp)
          data_subsampled <- rbind(data_subsampled, data.temp)
        }
        dim(data_subsampled)
      }

      ### 2.3 - Create set of column (parameter) names
      head(data_subsampled)                            # show data with headings
      
      ColumnNames_for_tSNE <- unname(colnames(data_subsampled))[ClusteringColNos] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10
      ColumnNames_for_tSNE  # check that the column names that appear are the ones you want to analyse
      
      dim(data_subsampled) # Review dimensionality of 'data' -- N events and N parameters
      data_specific <- data_subsampled[, ColumnNames_for_tSNE] # Prepare data for Rtsne --  select columns to use
      dim(data_specific) # Review dimensionality of 'data' -- N events and N parameters
      
      ### 5.3 - Save a copy of 'data' and 'data specific'        
      
      setwd(OutputDirectory)
      setwd("Output_tSNE/Output_tSNE_info")
      
      ## Export subsampled data in CSV format (record of what was run)
      subsampled_data_name <- paste0(Data_name, "_subsampled_data")
      subsampled_data <- paste(subsampled_data_name, ".csv", sep = "")
      write.csv(x = data_subsampled, file = subsampled_data)
      
      ## Export a list of parameters used to run tSNE and Phenograph
      tSNEparametersname <- paste0(Data_name, "_data_and_columns_used_for_tSNE")
      tSNEparameters <- paste(tSNEparametersname, ".csv", sep = "")
      write.csv(x = data_specific, file = tSNEparameters)      
      
      ### 5.4 - Run tSNE ('data' at this point is scaled, transformed) -- RUN ALL OF STEP 5
      
      ## Set working directory to "Output_info"
      setwd(OutputDirectory)
      setwd("Output_tSNE/Output_tSNE_info")
      
      ## Turn sink on (send the verbose text from the tSNE algorithm progress updates to a .txt file)
      verbose_name <- paste0("tSNE_verbose_", Data_name, ".txt")
      sink(file = verbose_name, append=TRUE, split=FALSE, type = c("output", "message"))
      
      ## Run tSNE algorithm
      set.seed(tSNE_seed)                             # default = 42 -- sets seed for reproducibility, delete or comment out for random tSNE
      tsne_out <- Rtsne(as.matrix(data_specific),     # dataset (assigned above) read as a matrix
                        dims = 2,                     # default = 2 -- only input 2 or 3 (2 = plot in 2D, 3 = plot in 3D)
                        initial_dims = 50,            # default = 50 -- number of dimensions retained in initial PCA step
                        perplexity = 30,              # default = 30
                        theta = 0.5,                  # default = 0.5 -- use 0.5 for Barnes-Hut tSNE, 0.0 for exact tSNE (takes longer)
                        check_duplicates = FALSE,      # default = TRUE, can set to FALSE if problematic
                        pca = TRUE,                   # default = TRUE -- performs a PCA prior to actual tSNE run
                        max_iter = 1000,              # default = 1000 -- default total iterations
                        verbose = TRUE,               # default = FALSE -- best to set as TRUE to receive feedback
                        is_distance = FALSE,          # default = FALSE -- experimental, using X as a distance matrix
                        Y_init = NULL,                # default = NULL -- recommend to use NULL for random initialisation
                        stop_lying_iter = 250,        # default = 250 -- number of iterations of early exaggeration
                        mom_switch_iter = 250,        # default = 250 -- number of iterations before increased momentum of spread
                        momentum = 0.5,               # default = 0.5 -- initial momentum of spread
                        final_momentum = 0.8,         # default = 0.8 -- momentum of spread at 'final_momentum'
                        eta = 200,                    # default = 200 -- learning rate
                        exaggeration_factor = 12.0    # default = 12.0 -- factor used during early exaggeration
      )
      
      ## turn sink off
      sink()
      
      ## if tSNE detected duplicates, the following errors should occur: 
      # "Error in Rtsne.default(as.matrix(data_specific), dims = 2, initial_dims = 50,  : Remove duplicates before running TSNE."
      # "Error in inherits(.data, "split") : object 'tsne_out' not found"
      
      
      ### 5.5 - Save info from tSNE run and add tSNE parameters to dataset
      
      ## Save tsne_out output info as .csv (i.e., output settings as raw results)
      tsne_out.df <- ldply(tsne_out, data.frame)
      output_name <- paste0("tSNE_parameter_values_", Data_name, ".csv")
      write.csv(x = tsne_out.df, file = output_name) # pretty blood good -- doesn't give row number for Y, costs, or itercosts -- but easy to figure out
      
      ## Add Y (tSNE coordinate) values to starting data
      tsne_out_Y <- tsne_out$Y
      head(tsne_out_Y)
      
      colnames(tsne_out_Y) <- c(tSNE_X_name, tSNE_Y_name)
      head(tsne_out_Y)
      
      ## save cluster labels
      data_subsampled <- cbind(data_subsampled, tsne_out_Y)
      head(data_subsampled)
      
      ### 3.5 - Write all data (with tSNE and FlowSOM parameters) to .csv and .fcs
          
          ## Set working directory
          setwd(OutputDirectory)
          setwd("Output_tSNE/Output_tSNE_data")
          
          ## Save data (with new tSNE parameters) as CSV
          write.csv(x = data_subsampled, file = paste0(Data_name, "_with_tSNE", ".csv"), row.names=FALSE)
          
          ## Check data and data column names
          head(data_subsampled)
          dimnames(data_subsampled)[[2]]
          
          ## Create FCS file metadata - column names with descriptions
          metadata <- data.frame(name=dimnames(data_subsampled)[[2]],
                                 desc=paste('column',dimnames(data_subsampled)[[2]],'from dataset')
          )
          
          ## Create FCS file metadata - ranges, min, and max settings
          #metadata$range <- apply(apply(data_subsampled,2,range),2,diff)
          metadata$minRange <- apply(data_subsampled,2,min)
          metadata$maxRange <- apply(data_subsampled,2,max)
          
          ## Create flowframe with tSNE data
          data_subsampled.ff <- new("flowFrame",
                                    exprs=as.matrix(data_subsampled), # in order to create a flow frame, data needs to be read as matrix
                                    parameters=AnnotatedDataFrame(metadata)
          )
          
          ## Check flow frame
          data_subsampled.ff
          
          ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
          new_file_name_fcs <- paste0(Data_name, "_with_tSNE", ".fcs")
          write.FCS(data_subsampled.ff, new_file_name_fcs)
          
          ### There is a delay here -- fcs file ends up in primary directory
          
          
      ### 3.7 - Write GROUPED data (with tSNE and FlowSOM parameters) to .csv and .fcs
      if(Run_Groups == 1){
          for(a in AllGroupNames){
            data_subset <- subset(data_subsampled, data_subsampled[[grp.col]] == a)
            dim(data_subsampled)
            
            ## write .csv
            write.csv(data_subset, file = paste(Data_name, "_", a, "_with_tSNE.csv", sep = ""), row.names=FALSE)
            
            ## write .fcs
            metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
            metadata
            
            ## Create FCS file metadata - ranges, min, and max settings
            #metadata$range <- apply(apply(data_subset,2,range),2,diff)
            metadata$minRange <- apply(data_subset,2,min)
            metadata$maxRange <- apply(data_subset,2,max)
            
            data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
            head(data_subset.ff)
            write.FCS(data_subset.ff, paste0(Data_name, "_", a, "_with_tSNE", ".fcs"))
          }
      }
      
      ### 3.8 - Write individual files (with tSNE and FlowSOM parameters) to .csv and .fcs     
      
      for (a in AllSampleNames) {
        data_subset <- subset(data_subsampled, data_subsampled[[samp.col]] == a)

        ## write .csv
        write.csv(data_subset, file = paste(Data_name, "_", a, "_with_tSNE.csv", sep = ""), row.names=FALSE)
        
        ## write .fcs
        metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
        
        ## Create FCS file metadata - ranges, min, and max settings
        #metadata$range <- apply(apply(data_subset,2,range),2,diff)
        metadata$minRange <- apply(data_subset,2,min)
        metadata$maxRange <- apply(data_subset,2,max)
        
        data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
        head(data_subset.ff)
        write.FCS(data_subset.ff, paste0(Data_name, "_", a, "_with_tSNE", ".fcs"))
      }
      
      
      ## Move back to PrimaryDirectory
      setwd(PrimaryDirectory)
      getwd()
      
      # note -- in the final data output, all parameters are included, but only the subsampled and/or transformed cells
      
    }

                
###################################################### 6. tSNEplots ######################################################  
                
    if(Run_tSNEplots == 1){    
      
      ## Create 'jet' colour scheme (not available by default in R)
      jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      
      ## Set your working directory here (e.g. "/Users/Tom/Desktop/")
      setwd(OutputDirectory)
      setwd("Output_tSNE/Output_tSNE_data") 
      
      ## Assign the working directory as 'PrimaryDirectory'
      PlotDirectory <- getwd()
      PlotDirectory
      
      ## Create a list of file names (names of the samples) and check file names
      PlotFiles <- list.files(path=PlotDirectory, pattern = ".csv") # ALTERNATIVE: list.files(pattern = '\\.csv')
      PlotFiles
      
      ## In the output of the previous line, you will see the names for the tSNE parameters -- insert them in between the "" below
      plotXname <- tSNE_X_name
      plotYname <- tSNE_Y_name
      
      
      ##### STEP 2b: ESTABLISH GLOBAL SCALE LIMITS FOR COLOUR and XY #####

      ## Create a 'list' of the data from all CSV files, then combine data into one large dataframe
      tables <- lapply(PlotFiles, read.csv, header = TRUE)
      combined.df <- do.call(rbind , tables)
      
      numeric.only <- sapply(combined.df, is.numeric)
      combined.df <- combined.df[ , numeric.only] # removes any non 'numeric' values
      
      ## Find column names for whole dataset
      names(combined.df)
      
      ## find maximum and minimum tSNE-X value --> define these
      Xmax <- max(combined.df[[plotXname]]) # double check function
      Ymax <- max(combined.df[[plotYname]])
      
      ## find maximum and minimum tSNE-Y value --> define these
      Xmin <- min(combined.df[[plotXname]]) # double check function
      Ymin <- min(combined.df[[plotYname]])

      # Using STEP 2b, the colour scale max and min will be the same for all samples, despite what the individual sample max or min is
      # Also using STEP2b, the X and Y limits will be the same for all samples, despite what the individual sample max or min is
      
      
      ##### STEP 3: Loop with  samples in separate folders #####  
      ## First, change the tSNE parameters (on lines 103 and 104)
      ## Then run all of the script below
      
      ## Set wd
      setwd(PlotDirectory)
      getwd()
      
      ## Create output folder (if a folder called "Output" already exists, nothing will happen)
      dir.create("Output(tSNEplots)", showWarnings = FALSE)
      
      for (File in PlotFiles){ 
        
        ## CSV is read into dataframe 'CurrentSampleCSV'
        CurrentSampleCSV <- read.csv(File)
        CurrentSampleCSV
        
        numeric.only <- sapply(CurrentSampleCSV, is.numeric) 
        CurrentSampleCSV <- CurrentSampleCSV[ , numeric.only] # removes any non 'numeric' values
        
        # Modifications to the name are made here
        File <- gsub(" ", "_", File) # replaces empty spaces in the file name with '_'
        File <- gsub(".csv", "", File) # removes ".csv" from the file name 
        
        ## Defines the 'tSNE' parameters that will be used to set the X and y axis
        plotX <- CurrentSampleCSV[[plotXname]] # defines the tSNE1 (x-axis) parameter name from your file
        plotY <- CurrentSampleCSV[[plotYname]] # defines the tSNE2 (y-axis) parameter name from your file
        
        ## Create subdirectory
        setwd("Output(tSNEplots)") 
        newdir <- paste0(File)
        dir.create(newdir, showWarnings = FALSE)
        setwd(newdir) 
        getwd()
        
        ## Sub-loop to create one image for every parameter
        for (i in names(CurrentSampleCSV)){ 
          
          tSNEplotLoop <- ggplot(
            data = `CurrentSampleCSV`,
            aes(x = plotX, y = plotY)) +
            geom_point(size = 0.5, mapping=aes_string(color=i))+ # 2 for large # 0.5 for small
            scale_colour_gradientn(colours = jet.colors(50),
                                   limits = c(quantile(combined.df[[i]], probs = c(0.01)), #0.03-01 seems to work well
                                              quantile(combined.df[[i]], probs = c(0.995))), #0.97-995 seems to work well
                                   oob=squish) + 
            ggtitle(i) +
            xlim(Xmin, Xmax)+
            ylim(Ymin, Ymax)+
            # xlab("tSNE1") + # use if desired, must also change theme settings below
            # ylab("tSNE2") + # use if desired, must also change theme settings below
            theme( # panel.background = element_rect(fill = "white", colour = "white", size = 0.5), # change 'colour' to black for informative axis
              axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              panel.grid.major = element_blank(),
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank(),
              legend.position = "right",
              legend.text=element_text(size=15), # large = 30 # small = 8
              legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
              legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
              legend.title=element_blank(),
              plot.title = element_text(color="Black", face="bold", size=15, hjust=0) # size 70 for large, # 18 for small
            )
          ggsave(tSNEplotLoop, filename = paste0(File, "-", i,".jpeg"), width = 4, height = 3) # Large size default width = 14.4 height = 12 (11.35 without title), default = PDF, # small w3.6, h3
        }
        
        ## Move back to PrimaryDirectory
        setwd(PlotDirectory)
        getwd()  
      }
    }
     
        
###################################################### 7. ClusterPlots ######################################################  
    
    if(Run_ClusterPlots == 1){
          
      ### Setup
          
          ## WD and colours
          jet.black <- colorRampPalette(c("black", "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000", "purple"))
          colour.scheme <- jet.black
      
          setwd(OutputDirectory)
          setwd("Output_tSNE/Output_tSNE_data")
          dir.create("Output_NumPlots", showWarnings = FALSE)
          
          FileNames <- list.files(path=getwd(), pattern = ".csv")
          FileNames
          
          ## In the output of the previous line, you will see the names for the tSNE parameters -- insert them in between the "" below
          plotXname <- tSNE_X_name
          plotYname <- tSNE_Y_name
          cluster.name <- FlowSOM_meta_name
          
          cluspl.n.sub # already defined
          cluspl.n.sub.seed # already defined
          
      ### Read files 
          ## Read all CSV files into a list
          files  <- list.files(pattern = '\\.csv')
          
          ## Create a 'list' of the data from all CSV files, then combine data into one large dataframe
          tables <- lapply(files, read.csv, header = TRUE)
          combined.df <- do.call(rbind , tables)
          numeric.only <- sapply(combined.df, is.numeric)
          combined.df <- combined.df[ , numeric.only] # removes any non 'numeric' values
          
          ## find maximum and minimum tSNE-X value --> define these
          Xmax <- max(combined.df[[plotXname]]) # double check function
          Ymax <- max(combined.df[[plotYname]])
          
          ## find maximum and minimum tSNE-Y value --> define these
          Xmin <- min(combined.df[[plotXname]]) # double check function
          Ymin <- min(combined.df[[plotYname]])
          
      ### Loop for Number Plot
          
          for (File in FileNames){ 
            setwd(OutputDirectory)
            setwd("Output_tSNE/Output_tSNE_data")
            
            ## CSV is read into dataframe 'CurrentSampleCSV'
            CurrentSampleCSV <- read.csv(File)
            numeric.only <- sapply(CurrentSampleCSV, is.numeric) 
            CurrentSampleCSV <- CurrentSampleCSV[ , numeric.only] # removes any non 'numeric' values
            
            ## Subsample to pre-defined number manageable number of points (only if samples has more cells than the target downsample number)
            if(nrow(CurrentSampleCSV) > cluspl.n.sub){
              set.seed(cluspl.n.sub.seed)
              CurrentSampleCSV <- CurrentSampleCSV[sample(1:nrow(CurrentSampleCSV), cluspl.n.sub), ]
            }

            # Modifications to the name are made here
            File <- gsub(" ", "_", File) # replaces empty spaces in the file name with '_'
            File <- gsub(".csv", "", File) # removes ".csv" from the file name 
            
            ## Defines the 'tSNE' parameters that will be used to set the X and y axis
            plotX <- CurrentSampleCSV[[plotXname]] # defines the tSNE1 (x-axis) parameter name from your file
            plotY <- CurrentSampleCSV[[plotYname]] # defines the tSNE2 (y-axis) parameter name from your file
            clust.no <- CurrentSampleCSV[[cluster.name]]
            
            ## Subdirectory
            setwd(OutputDirectory)
            setwd("Output_tSNE/Output_tSNE_data/Output_NumPlots")
            
            tSNEplotLoop <- ggplot(
              data = `CurrentSampleCSV`,
              aes(x = plotX, y = plotY)) +
              scale_colour_gradientn(colours = colour.scheme(length(unique(clust.no)))) +
              geom_text(aes(label=clust.no, colour=clust.no)) +
              
              #aes(x = plotX, y = plotY)) +
              #geom_point(size = 0.5)+ # 2 for large # 0.5 for small
              #scale_colour_manual(name = cluster.name, values = c(colour.scheme(length(unique(clust.no))))) +
              
              xlim(Xmin, Xmax)+
              ylim(Ymin, Ymax)+
              ggtitle(cluster.name) +
              # xlab("tSNE1") + # use if desired, must also change theme settings below
              # ylab("tSNE2") + # use if desired, must also change theme settings below
              theme( # panel.background = element_rect(fill = "white", colour = "white", size = 0.5), # change 'colour' to black for informative axis
                axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                panel.grid.major = element_blank(),
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank(),
                legend.position = "right",
                legend.text=element_text(size=15), # large = 30 # small = 8
                legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
                legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
                legend.title=element_blank(),
                plot.title = element_text(color="Black", face="bold", size=15, hjust=0) # size 70 for large, # 18 for small
              )
            ggsave(tSNEplotLoop, filename = paste0(File, ".jpeg"), width = 6, height = 4.5) # Large size default width = 14.4 height = 12 (11.35 without title), default = PDF, # small w3.6, h3
            
            ## Move back to PrimaryDirectory
            setwd(PrimaryDirectory)
            getwd()  
          }
          
          
          
      ### CLUSTER COLOUR PLOT
          
          ## Set wd
          setwd(OutputDirectory)
          setwd("Output_tSNE/Output_tSNE_data")
          dir.create("Output_ClusterPlots", showWarnings = FALSE)
          
          for (File in FileNames){ 
            setwd(OutputDirectory)
            setwd("Output_tSNE/Output_tSNE_data")
            
            ## CSV is read into dataframe 'CurrentSampleCSV'
            CurrentSampleCSV <- read.csv(File)
            numeric.only <- sapply(CurrentSampleCSV, is.numeric) 
            CurrentSampleCSV <- CurrentSampleCSV[ , numeric.only] # removes any non 'numeric' values
            
            # Modifications to the name are made here
            File <- gsub(" ", "_", File) # replaces empty spaces in the file name with '_'
            File <- gsub(".csv", "", File) # removes ".csv" from the file name 
            
            ## Defines the 'tSNE' parameters that will be used to set the X and y axis
            plotX <- CurrentSampleCSV[[plotXname]] # defines the tSNE1 (x-axis) parameter name from your file
            plotY <- CurrentSampleCSV[[plotYname]] # defines the tSNE2 (y-axis) parameter name from your file
            clust.no <- CurrentSampleCSV[[cluster.name]]
            
            ## Create subdirectory
            setwd(OutputDirectory)
            setwd("Output_tSNE/Output_tSNE_data/Output_ClusterPlots")
            
            tSNEplotLoop <- ggplot(
              data = `CurrentSampleCSV`,
              aes(x = plotX, y = plotY, colour = as.factor(clust.no))) +
              geom_point(size = 1)+ # 2 for large # 0.5 for small
              #scale_colour_gradientn(colours = colour.scheme(50)) + 
              scale_colour_manual(name = cluster.name, values = c(colour.scheme(length(unique(clust.no))))) +
              ggtitle(cluster.name) +
              xlim(Xmin, Xmax)+
              ylim(Ymin, Ymax)+
              ggtitle(cluster.name) +
              # xlab("tSNE1") + # use if desired, must also change theme settings below
              # ylab("tSNE2") + # use if desired, must also change theme settings below
              theme( # panel.background = element_rect(fill = "white", colour = "white", size = 0.5), # change 'colour' to black for informative axis
                axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                panel.grid.major = element_blank(),
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank(),
                legend.position = "right",
                #legend.text=element_text(size=10), # large = 30 # small = 8
                #legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
                #legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
                #legend.title=element_blank(),
                plot.title = element_text(color="Black", face="bold", size=15, hjust=0) # size 70 for large, # 18 for small
              )
            ggsave(tSNEplotLoop, filename = paste0(File, ".jpeg"), width = 10.5, height = 6.75) # Large size default width = 14.4 height = 12 (11.35 without title), default = PDF, # small w3.6, h3
            
            ## Move back to PrimaryDirectory
            setwd(PrimaryDirectory)
            getwd()  
          }
          
          ## Move back to PrimaryDirectory
          setwd(PrimaryDirectory)
          getwd() 
    }   
        