### Cytometry Analysis Pipeline for Complex Datasets (CAPX) v2.5_FMW

    # Thomas Ashhurst/Felix Marsh-Wakefield
    # 2019-06-07
    # thomas.ashhurst@sydney.edu.au
    # www.github.com/SydneyCytometry/CAPX

### Summary
    
    # 1 INSTALL AND LOAD PACKAGES
        # 1.1 - Install packages
        # 1.2 - Load packages
        # 1.3 - Set working directory
    
    # 2 USER INPUT -- DATA PREPARATION
        # 2.1 - Preferences
    
    # 3 USER INPUT -- LINE BY LINE
        # 3.1 - Remove unhelpful keywords
        # 3.2 - Remove duplicates
        # 3.3 - Arcsinh transformation
        # 3.4 - Add sample identifying keywords
        # 3.5 - Add group keywords
        # 3.6 - Downsample (options)
        # 3.7 - Merge data and remove duplicates (if selected)

    # 4 END USER INPUT - RUN ALL AT ONCE
        # 4.1 - Write .csv and .fcs files (all data in one file) 
        # 4.2 - Write .csv and .fcs files (individual files)
        # 4.3 - Write .csv and .fcs files (grouped files)

###################################################### 1. INSTALL AND-OR LOAD PACKAGES ###################################################### 


    ### 1.1 - Install packages (if not already installed) from CRAN
            
            ## From CRAN (required)
                if(!require('plyr')) {install.packages('plyr')}
                if(!require('data.table')) {install.packages('data.table')}
                if(!require('rstudioapi')) {install.packages('rstudioapi')}
                if(!require('devtools')){install.packages("devtools")}
            
            ## From CRAN (required for clustering and tSNE)
                if(!require('FlowSOM')) {source("https://bioconductor.org/biocLite.R") # for running FlowSOM ### POSSIBLY INSTALL FROM GITHUB
                  biocLite('FlowSOM')} 
                if(!require('Rtsne')) {install.packages("Rtsne")} # for running tSNE

            ## From CRAN (only required for plotting -- if installation unsuccessful, set Run_tSNEplots and Run_ClusterPlots to 0)
                if (!require("ggplot2")){install.packages("ggplot2")} # for plotting tSNE graphs
                if (!require("colorRamps")){install.packages("colorRamps")} # for colour scheme management
                if (!require("ggthemes")){install.packages("ggthemes")} # for plot themes
                if (!require("scales")){install.packages("scales")} # for re-scaling if necessary

            ## From Bioconductor (required)
                if(!require('flowViz')) {source("https://bioconductor.org/biocLite.R") 
                  biocLite('flowViz')}
                if(!require('flowCore')) {source("https://bioconductor.org/biocLite.R") 
                  biocLite('flowCore')}
                if(!require('Biobase')) {source("https://bioconductor.org/biocLite.R") 
                  biocLite('Biobase')}
    
            ## From Github repositories
                #if(!require('Rphenograph')){devtools::install_github("JinmiaoChenLab/Rphenograph")} # not required

    ### 1.4 Load packages       
                library('plyr')
                library('data.table')
                library('rstudioapi')
                library('devtools')
                library('FlowSOM') ###
                library('Rtsne')
                library('ggplot2')
                library('colorRamps')
                library('ggthemes')
                library('scales')
                library('flowCore') ###
                library('Biobase')
                library('flowViz') ### matrixStats

    ### 1.5 - Set working directory and assign as 'PrimaryDirectory'
    
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
            
     

###################################################### 2. DATA INPUT and PREPARATION ###################################################### 
         
    ### 2.1 - Read data into workspace
            
            ## Use to list the .fcs or .csv files in the working directory -- important, the only FCS/CSV file in the directory should be the one desired for analysis. If more than one are found, only the first file will be used
                list.files(path=PrimaryDirectory, pattern = ".fcs")     # see a list of FCS files
                list.files(path=PrimaryDirectory, pattern = ".csv")     # see a list of CSV files
                
            ## File details
                file.type               <- ".csv"         # Specfiy file type (".csv" or ".fcs") -- readings .fcs files not currently functional
                data.name               <- "Exp_name"    # a new name for the data - suggest name is sampletype_date_time (e.g. liver_20180203_1400)
                
            ## Check the list of files
                FileNames <- list.files(path=PrimaryDirectory, pattern = file.type) # Generate list of files
                as.matrix(FileNames) # See file names in a list

            ## Read data from Files into list of data frames
                DataList=list() # Creates and empty list to start 
                Length_check = list() # creates an empty list
                ColName_check = list() 
                nrow_check = list()
                
                if (file.type == ".csv"){
                  for (File in FileNames) { # Loop to read files into the list
                    tempdata <- fread(File, check.names = FALSE)
                    File <- gsub(".csv", "", File)
                    DataList[[File]] <- tempdata
                  }
                  for(i in c(1:(length(DataList)))){Length_check[[i]] <- length(names(DataList[[i]]))} # creates a list of the number of columns in each sample
                  for(i in c(1:(length(DataList)))){ColName_check[[i]] <- names(DataList[[i]])}
                  name.table <- data.frame(matrix(unlist(ColName_check), nrow = length(DataList), byrow = T))
                  for(i in c(1:(length(DataList)))){nrow_check[[i]] <- nrow(DataList[[i]])}
                }
                
                ###
                rm(tempdata)
                ###

                
                #if (file.type == ".fcs"){
                #  ParaList=list()
                #  for (File in FileNames) { # Loop to read files into the list
                #    tempdata <- exprs(read.FCS(FileNames[File], transformation = FALSE)) 
                #    File.g <- gsub(".fcs", "", File)
                #    DataList[[File.g]] <- tempdata
                #    ParaList[[File.g]] <- parameters(read.FCS(File)) 
                #  }
                # for(i in c(1:(length(DataList)))){Length_check[[i]] <- length(names(DataList[[i]]))} # creates a list of the number of columns in each sample
                # for(i in c(1:(length(DataList)))){ColName_check[[i]] <- names(DataList[[i]])}
                # name.table <- data.frame(matrix(unlist(ColName_check), nrow = length(DataList), byrow = T))
                # for(i in c(1:(length(DataList)))){nrow_check[[i]] <- nrow(DataList[[i]])}
                #}
               
                
    ### 2.2 - Data review
            
            ## Review the data from the first sample
                head(DataList[[1]])
            
            ## Review number of columns (parameters, or features), and the number of rows (cells) in each sample
                as.matrix(Length_check)    # number of columns in each sample
                as.matrix(nrow_check)      # number of rows in each sample
                
    ### 2.3 - Review column names, and remove troubleshom columns (if required)
            
            ## Review number of columns, and then all column names
                as.matrix(Length_check)
                name.table

            ## Remove troublesom columns (if required)
                ############### ONLY IF REQUIRED ###############
                ## Remove any troublesome columns (if required)
                #for (i in c(1:(length(DataList)))) {
                #  DataList[[i]]$SampleID <- NULL # after the '$', put column name here
                #}
                ################################################ 
            
            ## Final check -- ensure the number of columns in each file is consistent
                Length_check = list() # creates an empty list
                for(i in c(1:(length(DataList)))){Length_check[[i]] <- length(names(DataList[[i]]))} # creates a list of the number of columns in each sample
                as.matrix(Length_check) # check that the number of columns in each sample is the same length
                
            ## Final check -- check all column names
                for(i in c(1:(length(DataList)))){ColName_check[[i]] <- names(DataList[[i]])}
                name.table <- data.frame(matrix(unlist(ColName_check), nrow = length(DataList), byrow = T))
                name.table

    ### 2.4 - Remove any columns that represent empty channels -- leave in all cellular markers, live/dead, etc
            
            ## Column selection
                #as.matrix(names(DataList[[1]])) # review the list of columns
                #col.rmv <- c() # select the columns to remove
            
            ## Column review    
                #as.matrix(names(DataList[[1]][-c(col.rmv)])) # Check columns to KEEP
                #as.matrix(names(DataList[[1]][c(col.rmv)])) # Check columns to REMOVE
            
            ## Remove columns and check data
                #for (i in c(1:(length(DataList)))) {
                #  DataList[[i]] <- DataList[[i]][-c(col.rmv)]
                #}
                
            ## Review columns remaining
                #as.matrix(names(DataList[[1]]))
            
    ### 2.5 - Add sample identifiers (necessary to merge sample)
            
            ## Create a list of 'SampleName' and SampleNo' entries -- 1:n of samples, these will be matched in order
                AllSampleNames <- c(names(DataList))
                AllSampleNames # Character (words)
                
                AllSampleNos <- c(1:(length(DataList)))       ### <-- changed to 4+ to match with sample
                AllSampleNos # Intiger/numeric (numbers)
                
            ## Add 'SampleNo' and SampleName' to each 
                for (i in AllSampleNos) {DataList[[i]]$SampleNo <- i}
                for (a in AllSampleNames) {DataList[[a]]$SampleName <- a} # Inserting names doesn't stop the rest working, just messes with some auto visual stuff
                
                head(DataList[[1]]) # check to see that the sample name and number keywords have been entered on each row
                
                as.matrix(AllSampleNos)
                as.matrix(AllSampleNames)
                
    ### 2.6 - Add 'GroupNo' and 'GroupName' to each (required if generating 'group' files) (can skip, if no groups present in sample)
                
            ## REQUIRED IF YOU NEED TO SEPARATE DIFFERENT TREATMENT GROUPS -- OTHERWISE, CAN IGNORE
                
            ## Create empty lists
                group.names = list()
                group.nums = list()
                
            ## Setup group names for each sample [[1]] for sample 1, [[2]] for sample 2, etc
                group.names[[1]] <- "Mock_D7"                    # name of group 1
                group.names[[2]] <- "WNV_D7"     # name of group 2 (repeat line with [[3]] for a third group, continue for other groups)

                group.names # check group names
                
            ## Specify which samples belong to each group
                as.matrix(AllSampleNames)
                
                group.nums[[1]] <- c(1:6)       # samples that belong to group 1
                group.nums[[2]] <- c(7:12)      # samples that belong to group 2 (repeat line with [[3]] for a third group, continue for other groups)

                group.nums # check group names
                
            ## Add 'GroupName' and 'GroupNo' keywords to dataset
                num.of.groups <- c(1:length(group.names))
                for(a in num.of.groups){
                  for(i in c(group.nums[[a]])){
                    DataList[[i]]$GroupNo <- a   
                    DataList[[i]]$GroupName <- group.names[[a]]
                  }
                }
                
            ## Check column names
                head(DataList[[1]]) # check the entries for one of the samples in group 1
                head(DataList[[7]]) # check the entries for one of the samples in group 2
                
    
###################################################### 3. USER INPUT - CLUSTERING AND DIMENSIONALITY REDUCTION PREFERENCES ###################################################### 
        
    ### 3.1 - Downsampling          
            
            ## Review cells in sample 
                as.matrix(nrow_check) 
                
            ## Downsampling and duplicates
                downsample.files        <- 1              # Do you wish to downsample each of the samples? No = 0, Yes = 1
                downsample.seed         <- 42             # Seed for downsampling
 
            ## Set downsample targets (must be less than cells in sample)
                as.matrix(nrow_check)                     # Check the number of cells in each sample
                downsample.targets      <- c(rep(9000, each=length(unique(AllSampleNos))))
                
                # Use to set downsample value to be equal to minimum value
                #downsample.targets <- c(rep(nrow_check[[which.min(nrow_check)]], each=length(unique(AllSampleNos))))
                
                as.matrix(nrow_check) # Cells in sample
                as.matrix(downsample.targets) # Downsample targest   
                length(as.matrix(nrow_check)) == length(as.matrix(downsample.targets)) # If you have the same number of targets as samples, then = TRUE. Else, = FALSE
                
            ## Run downsampling
                if(downsample.files == 1){   
                  for (i in AllSampleNos) {
                    nsub <- downsample.targets[i]
                    set.seed(downsample.seed)
                    datalisttemp <- DataList[[i]]
                    datalisttemp <- datalisttemp[sample(1:nrow(datalisttemp), nsub), ]
                    DataList[[i]] <- datalisttemp
                  }
                  Downsample_check = list()
                  for(i in AllSampleNos){
                    Downsample_check[[i]] <- downsample.targets[i] == nrow(DataList[[i]])
                  }
                  as.matrix(Downsample_check) 
                }
                
            ## If downsampling performed, should return TRUE for each sample successfully downsampled
                
            ## Check number of rows (cells) per sample
                nrow_check = list()
                for(i in c(1:(length(DataList)))){nrow_check[[i]] <- nrow(DataList[[i]])}
                
                as.matrix(nrow_check)

                
    ### 3.2 - Merge data and remove duplicates
                
            ## Concatenate into one large data frame
                rbindlist(DataList, use.names = TRUE, #"use.names" groups columns based on name
                          fill = TRUE) #"fill" creates new columns if they're missing, giving them NA 
                data
                
            ## Remove 'DataList' 
                rm(DataList)
                
            ## Remove duplicates
                remove.duplicates       <- 1              # Removing duplicates aids in tSNE analysis? No = 0, Yes = 1
                
                dim(data)
                
                if(remove.duplicates == 1){
                  data <- data[!duplicated(data), ]  # remove rows containing duplicate values within rounding
                  data <- as.data.frame(data)
                  dim(data) # check to see if the number of parameters has reduced (no change means no duplicates were present)
                }

            ## Remove any troublesome columns (if required)
                as.matrix(colnames(data))
                
                {
                  data[["190BCKG"]] <- NULL # between [[" "]], put column name here
                  data[["140Ce"]] <- NULL
                  data[["133Cs"]] <- NULL
                  data[["157Gd"]] <- NULL
                  data[["113In"]] <- NULL
                  data[["208Pb"]] <- NULL
                  data[["120Sn"]] <- NULL
                  data[["131Xe"]] <- NULL
                  data[["SampleID"]] <- NULL
                  data[["Event.."]] <- NULL
                }
                
                as.matrix(colnames(data))
                
    ### 3.3 - If required, Arcsinh transformation     
                
            ## Do you want to perform an ARCSINH transformation?
                arcsinh.transform       <- 0            # No = 0, Yes = 1 # default = 0
                
            ## Specifications
                asinh_scale             <- 15           # ONLY if arcsinh.transform = 1. For CyTOF choose between 5 and 15, for flow you will need between 200 and 2000.
                
                head(data)
                col.names.dl <- names(data)             # show data with headings
                as.matrix(col.names.dl)                 # view the column 'number' for each parameter
                
                col.nos.scale <- c(2:32)                # specify column numbers to be transformed - e.g. c(11, 23, 10)] 
                
                col.names.dl[col.nos.scale]             # Columns to transform
                col.names.dl[-col.nos.scale]            # Columns NOT to transform
                
            ## Perform transform (if selected in preferences, otherwise transformation will not run)
                if(arcsinh.transform == 1){
                  col.names.SCALE <- col.names.dl[col.nos.scale]
                  data[, col.names.SCALE] <- asinh(data[, col.names.SCALE] / asinh_scale)
                  head(data)
                  summary(data)
                  }
            
                
    ### 3.4 - Write .csv and .fcs files of pre-processed samples   
            write.merged.file       <- 0              # Do you want to write one large merged file? No = 0, Yes = 1
            write.sep.files         <- 0              # Do you also want to write indivdual files for each sample? No = 0, Yes = 1
            write.group.files       <- 0              # Do you also want to write one file for each group? (requires use.groups = 1) No = 0, Yes = 1
                
     
###################################################### 4. USER INPUT - CLUSTERING AND DIMENSIONALITY REDUCTION PREFERENCES ###################################################### 
            
    ### 4.1 - Identify columns you would like to use to separate samples and groups
                
            ## Choose the column that defines the sample NAMES
                as.matrix(names(data))
                samp.col <- "SampleName"
                head(data[samp.col]) # check that this is actually the sample names
                
            ## Choose the column that defines the GROUPS NAMES
                as.matrix(names(data))
                grp.col <- "GroupName"
                head(data[grp.col]) # check that this is actually the group names
                
    ### 4.2 - Specify columns to use for clustering
                
            ## Review column names
                names(data)                            # show data with headings
                ColumnNames <- unname(colnames(data)) # assign reporter and marker names (column names) to 'ColumnNames'
                as.matrix(ColumnNames) # view the column 'number' for each parameter
                
            ## Specify column numbers to be used for clustering and first tSNE plots
                ClusteringColNos <- c(5,6,8,9,11,13,17:19,21:29,32)
                ClusteringColNos
                
                ClusteringCols <- ColumnNames[ClusteringColNos] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10
                
                ClusteringCols  # check that the column names that appear are the ones you want to analyse
                ColumnNames[-ClusteringColNos] # Check which columns are being EXCLUDED! 
                
                #If tSNE plots are to use different markers for algorithm, state them here (e.g. including all markers)
                as.matrix(ColumnNames)
                
                ClusteringColNos_tSNE <- c(2,4:6,8:9,11:13,16:19,21:30,32)
                ClusteringColNos_tSNE
                
                ClusteringCols_tSNE <- ColumnNames[ClusteringColNos_tSNE]
                ClusteringCols_tSNE
                
    ### 4.3 - Specify CLUSTERING options   
            
            ## FlowSOM options
                Run_FlowSOM           <- 1                # Option to run FlowSOM. Enter 1 for yes, or 0 for no
                FlowSOM_k             <- 40               # Number of metaclusters to derive from FlowSOM (can also be selected automatically, but this often does not perform well)
                FlowSOM_seedA          <- 8               # Seed for running FlowSOM
                FlowSOM_meta_seed     <- 42               # Seed for metaclustering
                
                # Repeat FlowSOM run with extra seeds
                Run_FlowSOM_repeat     <- 1               # Option to run FlowSOM with extra seeds. Enter 1 for yes, or 0 for no
                FlowSOM_seedB          <- 12              # 2nd seed for running FlowSOM
                FlowSOM_seedC          <- 42              # 3rd seed for running FlowSOM
                
                FlowSOM_clus_name     <- "FlowSOM_original_cluster" # Label for metacluster parameter name
                FlowSOM_meta_name     <- "FlowSOM_meta_cluster" # Label for metacluster parameter name
                
                write.pre.FSOM       <- 0              # Save a csv of the data that was used for clustering (i.e. including only columns desired for clustering)
                write.FSOM.merged     <- 1              # Do you want to write one large merged file? No = 0, Yes = 1
                write.FSOM.sep         <- 1              # Do you also want to write indivdual files for each sample? No = 0, Yes = 1
                write.FSOM.group       <- 1              # Do you also want to write one file for each group? (requires use.groups = 1) No = 0, Yes = 1
                
    ### 4.4 - Specify tSNE options   
                                
            ## tSNE options for first run
                Run_tSNE              <- 1                # Option to run tSNE. Enter a 1 to run tSNE, or 0 to skip 
                tSNE_seed             <- 42               # Seed for running tSNE
                tSNE_X_name           <- "tSNE_X_select"          # Label for tSNE1 parameter
                tSNE_Y_name           <- "tSNE_Y_select"          # Label for tSNE2 parameter
                
                write.pre.tSNE        <- 1              # Save a csv of the data that was used for clustering (i.e. including only columns desired for clustering)
                write.tSNE.merged     <- 1              # Do you want to write one large merged file? No = 0, Yes = 1
                write.tSNE.sep         <- 1              # Do you also want to write indivdual files for each sample? No = 0, Yes = 1
                write.tSNE.group       <- 1              # Do you also want to write one file for each group? (requires use.groups = 1) No = 0, Yes = 1
                
            ## tSNE options for second run (this can be used to have a second tSNE plot using different markers than what was used for FlowSOM)
                Run_tSNE2             <- 1                # Option to run tSNE. Enter a 1 to run tSNE, or 0 to skip 
                tSNE_seed2            <- 42               # Seed for running tSNE
                tSNE_X_name2          <- "tSNE_X_all"          # Label for tSNE1 parameter
                tSNE_Y_name2          <- "tSNE_Y_all"          # Label for tSNE2 parameter
                
                write.pre.tSNE2       <- 1              # Save a csv of the data that was used for clustering (i.e. including only columns desired for clustering)
                write.tSNE.merged2    <- 1              # Do you want to write one large merged file? No = 0, Yes = 1
                write.tSNE.sep2       <- 1              # Do you also want to write indivdual files for each sample? No = 0, Yes = 1
                write.tSNE.group2     <- 1              # Do you also want to write one file for each group? (requires use.groups = 1) No = 0, Yes = 1
                
            ## tSNE downsampling options
                tSNE_dwnsmpl          <- "specific"       # can be "specific" or "uniform". 'specific' means each sample is downsampled to a specific value. 'uniform' means the entire dataset will be evenly downsampled to a target total number of cells. 
                tSNE_nsub_seed        <- 42               # Seed for downsampling for tSNE
                
                as.matrix(nrow_check)             # review the target number of cells per sample for the pre-processing downsampling
                as.matrix(AllSampleNames)
                
                  # IF tSNE_dwnsmpl = "specific" -- set downsample targets for each sample # Must be less than the number of cells in the original samples
                  tSNE.dwnsmp.targets   <- c(100,
                                             100,
                                             100,
                                             100,
                                             100,
                                             100,
                                             100,
                                             100,
                                             100,
                                             100,
                                             100,
                                             100
                                             )
                
                  # IF tSNE_dwnsmpl = "uniform"
                  #tSNE_nsub             <- 70000            # ONLY used if tSNE_dwnsmpl = "uniform". Specify the total downsample size for the data run in tSNE.
        
    ### 4.5 - Specify PLOTTING and other output options  
                
            ## Generate summary tables?
                #run.sumtables         <- 1              # Generate summary tables for quantitiatve analysis # NOT YET ACTIVE
                
            ## tSNEplot options
                Run_tSNEplots         <- 1              # 1 to run tSNEplots and generate coloured tSNE maps
                
            ## ClusterPlot options -- can ONLY be run if both FlowSOM and tSNE have been run
                Run_ClusterPlots      <- 1              # 1 to generate coloured cluster plots and number plots
                cluspl.n.sub          <- 1000           # Define subsampling target and seed for NUMBER PLOT ONLY (default = 1000 per plot, and a set of 42)
                cluspl.n.sub.seed     <- 42 
         
     
                
############################################################################################################################
###################################################### END USER INPUT ###################################################### 
############################################################################################################################ 

    ## Create output directory
        setwd(PrimaryDirectory)
        dir.create("Output_preprocess", showWarnings = FALSE)
        setwd("Output_preprocess")
        OutputDirectory <- getwd()
        OutputDirectory

    ### 4.1 - Export data in single large file -- both .csv and .fcs format
       
        if(write.merged.file == 1){
        
          ## write .csv
          csv.filename <- paste(paste0(data.name), ".csv", sep = "")
          fwrite(x = data, file = csv.filename, row.names=FALSE)
  
          ## write .fcs
            #metadata <- data.frame(name=dimnames(data)[[2]],desc=paste('column',dimnames(data)[[2]],'from dataset'))
          
          ## Create FCS file metadata - ranges, min, and max settings
          #metadata$range <- apply(apply(data,2,range),2,diff)
            #metadata$minRange <- apply(data,2,min)
            #metadata$maxRange <- apply(data,2,max)
          
            #data.ff <- new("flowFrame",exprs=as.matrix(data), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
            #head(data.ff)
            #write.FCS(data.ff, paste0(data.name, ".fcs"))
        }
        
        
    ### 4.2 - Export data in individual files -- both .csv and .fcs format

        if(write.sep.files == 1){
          
          for (a in AllSampleNames) {
            
            data_subset <- subset(data, SampleName == a)
            dim(data_subset)
            
            ## write .csv
            fwrite(data_subset, file = paste(data.name, "_", a,".csv", sep = ""), row.names=FALSE)
            
            ## write .fcs
              #metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
            
            ## Create FCS file metadata - ranges, min, and max settings
            #metadata$range <- apply(apply(data_subset,2,range),2,diff)
              #metadata$minRange <- apply(data_subset,2,min)
              #metadata$maxRange <- apply(data_subset,2,max)
            
              #data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
              #head(data_subset.ff)
              #write.FCS(data_subset.ff, paste0(data.name, "_", a, ".fcs"))
          }
        }
         
    ### 4.3 - Export data as grouped files -- both .csv and .fcs format
                  
        if(write.group.files == 1){
          
          for(a in group.names){
            data_subset <- subset(data, GroupName == a)
            dim(data_subset)
            
            ## write .csv
            fwrite(data_subset, file = paste(data.name, "_", a,".csv", sep = ""), row.names=FALSE)
            
            ## write .fcs
              #metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
              #metadata
            
            ## Create FCS file metadata - ranges, min, and max settings
            #metadata$range <- apply(apply(data_subset,2,range),2,diff)
              #metadata$minRange <- apply(data_subset,2,min)
              #metadata$maxRange <- apply(data_subset,2,max)
            
              #data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
              #head(data_subset.ff)
              #write.FCS(data_subset.ff, paste0(data.name,"_",a,".fcs"))
          }
        }

        
###################################################### 5. DATA PREP FOR CLUSTERING ######################################################         
        
    ### 5.1. Save a copy of 'data' and a list of columns for clustering
        
        AllSampleNames <- unique(data[[samp.col]])
        AllGroupNames <- unique(data[[grp.col]])
        
        ## Directory creation                    
        setwd(PrimaryDirectory)
        dir.create("Output_cluster", showWarnings = FALSE)
        setwd("Output_cluster")
        OutputDirectory <- getwd()
        OutputDirectory
        
        ## Export list of columns for clustering
        ClusteringCols_csv <- paste(paste0("Columns_for_Clustering"), ".csv", sep = "")
        write.csv(x = ClusteringCols, file = ClusteringCols_csv)
        
        ## Export list of columns for second tSNE plot
        tSNE2_Cols_csv <- paste(paste0("ClusteringColNos_tSNE"), ".csv", sep = "")
        write.csv(x = ClusteringCols_tSNE, file = tSNE2_Cols_csv)
        
        ## Export data in CSV format
        if(write.pre.FSOM == 1){
          data_csv <- paste(paste0(data.name), "_used_for_clustering", ".csv", sep = "")
          write.csv(x = data, file = data_csv, row.names=FALSE)
        }

        
###################################################### 6. FlowSOM ######################################################         
        
        if(Run_FlowSOM == 1) {
          
        ### 6.1 - Create output directories
          
          setwd(OutputDirectory)
          dir.create("Output_FlowSOM", showWarnings = FALSE)
          
          setwd("Output_FlowSOM")
          dir.create("Output_FlowSOM_info", showWarnings = FALSE)
          dir.create("Output_FlowSOM_data", showWarnings = FALSE)
          
        ### 6.2 - Prep data for FlowSOM
          
          ## Check data and data column names
          head(data)
          dimnames(data)[[2]]
          
          ## Create FCS file metadata - column names with descriptions
          metadata <- data.frame(name=dimnames(data)[[2]], desc=paste('column',dimnames(data)[[2]],'from dataset'))
          
          ## Create FCS file metadata - ranges, min, and max settings -- by default, they are commented out (adjust ranges manually in FlowJo)
          #metadata$range <- apply(apply(data,2,range),2,diff) # throws an error because of word entry -- hopefully is ok
            #metadata$minRange <- apply(data,2,min)
            #metadata$maxRange <- apply(data,2,max)
          
          ## Create flowframe with data
          data.ff <- new("flowFrame",
                         exprs=as.matrix(data), # in order to create a flow frame, data needs to be read as matrix
                         parameters=AnnotatedDataFrame(metadata))
          
          head(exprs(data.ff))
          
          data_FlowSOM <- data.ff
          
          # choose markers for FlowSOM analysis
          FlowSOM_cols <- ClusteringCols
          
        ### 6.3. - Run FlowSOM              
          
          ## set seed for reproducibility
          set.seed(FlowSOM_seedA)
          
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
          colnames(res.original)[grepl('labels_pre',colnames(res.original))] <- paste0(FlowSOM_clus_name, FlowSOM_seedA)
          
          dim(data)
          dim(res.original)
          head(res.original)
          
          ## save META cluster labels
          res.meta <- data.frame("labels" = labels)
          colnames(res.meta)[grepl('labels',colnames(res.meta))] <- paste0(FlowSOM_meta_name, FlowSOM_seedA)
          
          dim(data)
          dim(res.meta)
          head(res.meta)
          
        ### 6.4. - Add FlowSOM cluster numbers to data
          
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
          
          
          ## Running FlowSOM with extra seeds
          if (Run_FlowSOM_repeat == 1) {
            
           ### 6.4.1 - Prep data for FlowSOM
            
            ## Check data and data column names
            head(data)
            dimnames(data)[[2]]
            
            ## Create FCS file metadata - column names with descriptions
            metadata <- data.frame(name=dimnames(data)[[2]], desc=paste('column',dimnames(data)[[2]],'from dataset'))
            
            ## Create flowframe with data
            data.ff <- new("flowFrame",
                           exprs=as.matrix(data), # in order to create a flow frame, data needs to be read as matrix
                           parameters=AnnotatedDataFrame(metadata))
            
            head(exprs(data.ff))
            
            data_FlowSOM <- data.ff
            
            # choose markers for FlowSOM analysis
            FlowSOM_cols <- ClusteringCols
            
           ### 6.4.2. - Run FlowSOM              
            
            ## set seed for reproducibility
            set.seed(FlowSOM_seedB)
            
            ## run FlowSOM (initial steps prior to meta-clustering)
            FlowSOM_out <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
            FlowSOM_out <- FlowSOM::BuildSOM(FlowSOM_out, colsToUse = FlowSOM_cols)
            FlowSOM_out <- FlowSOM::BuildMST(FlowSOM_out)
            
            ### some warnings will be returned because of the 'SampleName' and 'GroupName' entries

            ## extract cluster labels (pre meta-clustering) from output object
            labels_pre <- FlowSOM_out$map$mapping[, 1]
            labels_pre
            length(labels_pre)
            res.original <- labels_pre
            
            ## run meta-clustering
            FlowSOM_out_meta <- FlowSOM::metaClustering_consensus(FlowSOM_out$map$codes, k = FlowSOM_k, seed = FlowSOM_meta_seed)

            ## extract META (?) cluster labels from output object
            labels <- FlowSOM_out_meta[labels_pre]
            
            ## summary of cluster sizes and number of clusters
            table(labels)
            length(table(labels))
            
            ## save ORIGINAL cluster labels
            res.original <- data.frame("labels_pre" = labels_pre)
            colnames(res.original)[grepl('labels_pre',colnames(res.original))] <- paste0(FlowSOM_clus_name, FlowSOM_seedB)
            
            dim(data)
            dim(res.original)
            head(res.original)
            
            ## save META cluster labels
            res.meta <- data.frame("labels" = labels)
            colnames(res.meta)[grepl('labels',colnames(res.meta))] <- paste0(FlowSOM_meta_name, FlowSOM_seedB)
            
            dim(data)
            dim(res.meta)
            head(res.meta)
            
           ### 6.4.3 - Add FlowSOM cluster numbers to data
            
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
            
           ### 6.4.4 - Prep data for FlowSOM
            
            ## Check data and data column names
            head(data)
            dimnames(data)[[2]]
            
            ## Create FCS file metadata - column names with descriptions
            metadata <- data.frame(name=dimnames(data)[[2]], desc=paste('column',dimnames(data)[[2]],'from dataset'))
            
            ## Create flowframe with data
            data.ff <- new("flowFrame",
                           exprs=as.matrix(data), # in order to create a flow frame, data needs to be read as matrix
                           parameters=AnnotatedDataFrame(metadata))
            
            head(exprs(data.ff))
            
            data_FlowSOM <- data.ff
            
            # choose markers for FlowSOM analysis
            FlowSOM_cols <- ClusteringCols
            
           ### 6.4.5. - Run FlowSOM              
            
            ## set seed for reproducibility
            set.seed(FlowSOM_seedC)
            
            ## run FlowSOM (initial steps prior to meta-clustering)
            FlowSOM_out <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
            FlowSOM_out <- FlowSOM::BuildSOM(FlowSOM_out, colsToUse = FlowSOM_cols)
            FlowSOM_out <- FlowSOM::BuildMST(FlowSOM_out)
            
            ### some warnings will be returned because of the 'SampleName' and 'GroupName' entries
            
            ## extract cluster labels (pre meta-clustering) from output object
            labels_pre <- FlowSOM_out$map$mapping[, 1]
            labels_pre
            length(labels_pre)
            res.original <- labels_pre
            
            ## run meta-clustering
            FlowSOM_out_meta <- FlowSOM::metaClustering_consensus(FlowSOM_out$map$codes, k = FlowSOM_k, seed = FlowSOM_meta_seed)

            ## extract META (?) cluster labels from output object
            labels <- FlowSOM_out_meta[labels_pre]
            
            ## summary of cluster sizes and number of clusters
            table(labels)
            length(table(labels))
            
            ## save ORIGINAL cluster labels
            res.original <- data.frame("labels_pre" = labels_pre)
            colnames(res.original)[grepl('labels_pre',colnames(res.original))] <- paste0(FlowSOM_clus_name, FlowSOM_seedC)
            
            dim(data)
            dim(res.original)
            head(res.original)
            
            ## save META cluster labels
            res.meta <- data.frame("labels" = labels)
            colnames(res.meta)[grepl('labels',colnames(res.meta))] <- paste0(FlowSOM_meta_name, FlowSOM_seedC)
            
            dim(data)
            dim(res.meta)
            head(res.meta)
            
           ### 6.4.6 - Add FlowSOM cluster numbers to data
            
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
            
            
          }
          
          
        ### 6.5. - Export data to .csv and .fcs 
          
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
          # write.csv(as.list(FlowSOM_out), file = paste(paste0(data.name), "_FlowSOM_details.csv", sep = ""))
          # lapply(as.list(FlowSOM_out), function(x) write.table( data.frame(x), 'test.csv'  , append= T, sep=',' ))
          # capture.output(summary(FlowSOM_out), file = "FlowSOM_output.txt")
          # could try: write.list(FlowSOM_out, file = "test")
          
          ## Write to .csv
          setwd(OutputDirectory)
          setwd("Output_FlowSOM/Output_FlowSOM_data")
          getwd()
          

          ### 6.6 - Write merged sample file
          
              if(write.FSOM.merged == 1){
              
              fwrite(data, file = paste(paste0(data.name), "_with_FlowSOM.csv", sep = ""), row.names=FALSE)
                
              ## Write to .fcs
                #metadata <- data.frame(name=dimnames(data)[[2]],desc=paste('column',dimnames(data)[[2]],'from dataset'))
                #metadata ## Double check metadata
              
              ## Create FCS file metadata - ranges, min, and max settings
              #metadata$range <- apply(apply(data,2,range),2,diff)
                #metadata$minRange <- apply(data,2,min)
                #metadata$maxRange <- apply(data,2,max)
              
                #data.ff <- new("flowFrame",exprs=as.matrix(data), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
                #head(data.ff)
                #write.FCS(data.ff, paste0(data.name, "_with_FlowSOM", ".fcs")) ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
              
              #head(exprs(data.ff))
              }
          
          ### 6.7. - Export each group as a separate .csv and .fcs file             
          
              if(write.FSOM.group == 1){
                for(a in AllGroupNames){
                  data_subset <- subset(data, data[[grp.col]] == a)
                  dim(data_subset)
                  
                  ## write .csv
                  fwrite(data_subset, file = paste0(data.name, "_", a, "_with_FlowSOM", ".csv"), row.names=FALSE)
                  
                  ## write .fcs
                    #metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
                    #metadata
                  
                  ## Create FCS file metadata - ranges, min, and max settings
                  #metadata$range <- apply(apply(data_subset,2,range),2,diff)
                    #metadata$minRange <- apply(data_subset,2,min)
                    #metadata$maxRange <- apply(data_subset,2,max)
                  
                    #data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
                    #head(data_subset.ff)
                    #write.FCS(data_subset.ff, paste0(data.name, "_", a, "_with_FlowSOM", ".fcs"))
                }
              }
          
          ### 6.8. - Export each sample as a separate .csv and .fcs file             
              
              if(write.FSOM.sep == 1){
                for (a in AllSampleNames) {
                  
                  data_subset <- subset(data, data[[samp.col]] == a)
                  dim(data_subset)
                  
                  ## write .csv
                  fwrite(data_subset, file = paste0(data.name, "_", a, "_with_FlowSOM", ".csv"), row.names=FALSE)
                  
                  ## write .fcs
                    #metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
                    #metadata
                  
                  ## Create FCS file metadata - ranges, min, and max settings
                  #metadata$range <- apply(apply(data_subset,2,range),2,diff)
                    #metadata$minRange <- apply(data_subset,2,min)
                    #metadata$maxRange <- apply(data_subset,2,max)
                  
                    #data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
                    #head(data_subset.ff)
                    #write.FCS(data_subset.ff, paste0(data.name, "_", a, "_with_FlowSOM", ".fcs"))
                }
              }
              
              setwd(PrimaryDirectory)
              getwd()
          
        }
        

###################################################### 7. tSNE-all markers ######################################################                       
        
        if(Run_tSNE2 == 1){
          
        ### 7.1 - Create output directories
          
          setwd(OutputDirectory)
          dir.create("Output_tSNE_all", showWarnings = FALSE)
          setwd("Output_tSNE_all")
          dir.create("Output_tSNE_info", showWarnings = FALSE)
          dir.create("Output_tSNE_data", showWarnings = FALSE)
          setwd(OutputDirectory)
          
        ### 7.2 - Subsample full dataset (including FlowSOM cluster no, if FlowSOM was run)
          
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
              nsub <- tSNE.dwnsmp.targets[i] 
              data.temp <- subset(data, data[[samp.col]] == nam) # works
              nrow(data.temp)
              set.seed(tSNE_nsub_seed)
              data.temp <- data.temp[sample(1:nrow(data.temp), nsub), ]
              nrow(data.temp)
              data_subsampled <- rbind(data_subsampled, data.temp)
            }
            dim(data_subsampled)
          }
          
        ### 7.3 - Create set of column (parameter) names
          head(data_subsampled)                            # show data with headings
          
          ColumnNames_for_tSNE <- unname(colnames(data_subsampled))[ClusteringColNos_tSNE] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10
          ColumnNames_for_tSNE  # check that the column names that appear are the ones you want to analyse
          
          dim(data_subsampled) # Review dimensionality of 'data' -- N events and N parameters
          data_specific <- data_subsampled[, ColumnNames_for_tSNE] # Prepare data for Rtsne --  select columns to use
          dim(data_specific) # Review dimensionality of 'data' -- N events and N parameters
          
        ### 7.4 - Save a copy of 'data' and 'data specific'        
          
          setwd(OutputDirectory)
          setwd("Output_tSNE_all/Output_tSNE_info")
          
          ## Export a list of parameters used to run tSNE and Phenograph
          tSNEparametersname <- paste0(data.name, "_data_and_columns_used_for_tSNE")
          tSNEparameters <- paste(tSNEparametersname, ".csv", sep = "")
          write.csv(x = data_specific, file = tSNEparameters)  
          
          ## Export subsampled data in CSV format (record of what was run)
          if(write.pre.tSNE2 == 1){
            subsampled_data_name <- paste0(data.name, "_subsampled_data")
            subsampled_data <- paste(subsampled_data_name, ".csv", sep = "")
            write.csv(x = data_subsampled, file = subsampled_data)
            
            
          }
          
        ### 7.5 - Run tSNE ('data' at this point is scaled, transformed) -- RUN ALL OF STEP 5
          
          ## Set working directory to "Output_info"
          setwd(OutputDirectory)
          setwd("Output_tSNE_all/Output_tSNE_info")
          
          ## Turn sink on (send the verbose text from the tSNE algorithm progress updates to a .txt file)
          verbose_name <- paste0("tSNE_verbose_", data.name, ".txt")
          sink(file = verbose_name, append=TRUE, split=FALSE, type = c("output", "message"))
          
          ## Run tSNE algorithm
          set.seed(tSNE_seed2)                             # default = 42 -- sets seed for reproducibility, delete or comment out for random tSNE
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
          
          
        ### 7.6 - Save info from tSNE run and add tSNE parameters to dataset
          
          ## Save tsne_out output info as .csv (i.e., output settings as raw results)
          tsne_out.df <- ldply(tsne_out, data.frame)
          output_name <- paste0("tSNE_parameter_values_", data.name, ".csv")
          fwrite(x = tsne_out.df, file = output_name, row.names = TRUE) # pretty blood good -- doesn't give row number for Y, costs, or itercosts -- but easy to figure out
          
          ## Add Y (tSNE coordinate) values to starting data
          tsne_out_Y <- tsne_out$Y
          head(tsne_out_Y)
          
          colnames(tsne_out_Y) <- c(tSNE_X_name2, tSNE_Y_name2)
          head(tsne_out_Y)
          
          ## save cluster labels
          data_subsampled <- cbind(data_subsampled, tsne_out_Y)
          head(data_subsampled)
          
        ### 7.7 - Write all data (with tSNE and FlowSOM parameters) to .csv and .fcs
          
          ## Set working directory
          setwd(OutputDirectory)
          setwd("Output_tSNE_all/Output_tSNE_data")
              
              if(write.tSNE.merged2 == 1){
    
                  ## Save data (with new tSNE parameters) as CSV
                  write.csv(x = data_subsampled, file = paste0(data.name, "_with_tSNE", ".csv"), row.names=FALSE)
                  
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
                  new_file_name_fcs <- paste0(data.name, "_with_tSNE", ".fcs")
                  write.FCS(data_subsampled.ff, new_file_name_fcs)
                  
                  ### There is a delay here -- fcs file ends up in primary directory
              }
          
        ### 7.8 - Write GROUPED data (with tSNE and FlowSOM parameters) to .csv and .fcs
          if(write.tSNE.group2 == 1){
            for(a in AllGroupNames){
              data_subset <- subset(data_subsampled, data_subsampled[[grp.col]] == a)
              dim(data_subsampled)
              
              ## write .csv
              write.csv(data_subset, file = paste0(data.name, "_", a, "_with_tSNE.csv", sep = ""), row.names=FALSE)
              
              ## write .fcs
              metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
              metadata
              
              ## Create FCS file metadata - ranges, min, and max settings
              #metadata$range <- apply(apply(data_subset,2,range),2,diff)
              metadata$minRange <- apply(data_subset,2,min)
              metadata$maxRange <- apply(data_subset,2,max)
              
              data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
              head(data_subset.ff)
              write.FCS(data_subset.ff, paste0(data.name, "_", a, "_with_tSNE", ".fcs"))
            }
          }
          
        ### 7.9 - Write individual files (with tSNE and FlowSOM parameters) to .csv and .fcs     
          
          if(write.tSNE.sep2 == 1){
    
              for (a in AllSampleNames) {
                data_subset <- subset(data_subsampled, data_subsampled[[samp.col]] == a)
                
                ## write .csv
                write.csv(data_subset, file = paste0(data.name, "_", a, "_with_tSNE.csv", sep = ""), row.names=FALSE)
                
                ## write .fcs
                metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
                
                ## Create FCS file metadata - ranges, min, and max settings
                #metadata$range <- apply(apply(data_subset,2,range),2,diff)
                metadata$minRange <- apply(data_subset,2,min)
                metadata$maxRange <- apply(data_subset,2,max)
                
                data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
                head(data_subset.ff)
                write.FCS(data_subset.ff, paste0(data.name, "_", a, "_with_tSNE", ".fcs"))
              }
              
              ## Move back to PrimaryDirectory
              setwd(PrimaryDirectory)
              getwd()
          }
          
          # note -- in the final data output, all parameters are included, but only the subsampled and/or transformed cells
          
        }
        
        
###################################################### 8. tSNEplots-all markers ######################################################  
        
        if(Run_tSNEplots == 1){    
          
          ## Create 'jet' colour scheme (not available by default in R)
          jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
          
          ## Set your working directory here (e.g. "/Users/Tom/Desktop/")
          setwd(OutputDirectory)
          setwd("Output_tSNE_all/Output_tSNE_data") 
          
          ## Assign the working directory as 'PrimaryDirectory'
          PlotDirectory <- getwd()
          PlotDirectory
          
          ## Create a list of file names (names of the samples) and check file names
          PlotFiles <- list.files(path=PlotDirectory, pattern = ".csv") # ALTERNATIVE: list.files(pattern = '\\.csv')
          PlotFiles
          
          ## In the output of the previous line, you will see the names for the tSNE parameters -- insert them in between the "" below
          plotXname <- tSNE_X_name2
          plotYname <- tSNE_Y_name2
          
          
        ##### 8.1 - ESTABLISH GLOBAL SCALE LIMITS FOR COLOUR and XY #####
          
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
          
          
        ### 8.2 - Loop with  samples in separate folders #####  
          ## First, change the tSNE parameters (on lines 103 and 104)
          ## Then run all of the script below
          
          ## Set wd
          setwd(PlotDirectory)
          getwd()
          
          ## Create output folder (if a folder called "Output" already exists, nothing will happen)
          dir.create("Output(tSNEplots-all)", showWarnings = FALSE)
          
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
            setwd("Output(tSNEplots-all)") 
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
        
        
###################################################### 9. ClusterPlots-all ######################################################  
        
        if(Run_ClusterPlots == 1){
          
        ### 9.1 - Setup
          
          ## WD and colours
          jet.black <- colorRampPalette(c("black", "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000", "purple"))
          colour.scheme <- jet.black
          
          setwd(OutputDirectory)
          setwd("Output_tSNE_all/Output_tSNE_data")
          dir.create("Output_NumPlots", showWarnings = FALSE)
          
          FileNames <- list.files(path=getwd(), pattern = ".csv")
          FileNames
          
          ## In the output of the previous line, you will see the names for the tSNE parameters -- insert them in between the "" below
          plotXname <- tSNE_X_name2
          plotYname <- tSNE_Y_name2
          cluster.name <- FlowSOM_meta_name
          
          cluspl.n.sub # already defined
          cluspl.n.sub.seed # already defined
          
        ### 9.2 - Read files 
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
          
        ### 9.3 - Loop for Number Plot
          
          for (File in FileNames){ 
            setwd(OutputDirectory)
            setwd("Output_tSNE_all/Output_tSNE_data")
            
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
            
            cluster.name <- paste0(FlowSOM_meta_name, FlowSOM_seedA)
            clust.no <- CurrentSampleCSV[[cluster.name]]
            
            ## Subdirectory
            setwd(OutputDirectory)
            setwd("Output_tSNE_all/Output_tSNE_data/Output_NumPlots")
            dir.create(paste0("Output_NumPlots", "_seed", FlowSOM_seedA), showWarnings = FALSE)
            
            setwd(OutputDirectory)
            setwd(paste0("Output_tSNE_all/Output_tSNE_data/Output_NumPlots/", "Output_NumPlots", "_seed", FlowSOM_seedA))
            
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
          
          if (Run_FlowSOM_repeat == 1){
           ### 9.3.1 - Second FlowSOM seed plots
            for (File in FileNames){ 
              setwd(OutputDirectory)
              setwd("Output_tSNE_all/Output_tSNE_data")
              
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
              
              cluster.name <- paste0(FlowSOM_meta_name, FlowSOM_seedB)
              clust.no <- CurrentSampleCSV[[cluster.name]]
              
              ## Subdirectory
              setwd(OutputDirectory)
              setwd("Output_tSNE_all/Output_tSNE_data/Output_NumPlots")
              dir.create(paste0("Output_NumPlots", "_seed", FlowSOM_seedB), showWarnings = FALSE)
              
              setwd(OutputDirectory)
              setwd(paste0("Output_tSNE_all/Output_tSNE_data/Output_NumPlots/", "Output_NumPlots", "_seed", FlowSOM_seedB))
              
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
            
           ### 9.3.2 - Third FlowSOM seed plots
            for (File in FileNames){ 
              setwd(OutputDirectory)
              setwd("Output_tSNE_all/Output_tSNE_data")
              
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
              
              cluster.name <- paste0(FlowSOM_meta_name, FlowSOM_seedC)
              clust.no <- CurrentSampleCSV[[cluster.name]]
              
              ## Subdirectory
              setwd(OutputDirectory)
              setwd("Output_tSNE_all/Output_tSNE_data/Output_NumPlots")
              dir.create(paste0("Output_NumPlots", "_seed", FlowSOM_seedC), showWarnings = FALSE)
              
              setwd(OutputDirectory)
              setwd(paste0("Output_tSNE_all/Output_tSNE_data/Output_NumPlots/", "Output_NumPlots", "_seed", FlowSOM_seedC))
              
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
            
          }
          
          
        ### 9.4 - CLUSTER COLOUR PLOT
          
          ## Set wd
          setwd(OutputDirectory)
          setwd("Output_tSNE_all/Output_tSNE_data")
          dir.create("Output_ClusterPlots", showWarnings = FALSE)
          
          for (File in FileNames){ 
            setwd(OutputDirectory)
            setwd("Output_tSNE_all/Output_tSNE_data")
            
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
            
            cluster.name <- paste0(FlowSOM_meta_name, FlowSOM_seedA)
            clust.no <- CurrentSampleCSV[[cluster.name]]
            
            ## Create subdirectory
            setwd(OutputDirectory)
            setwd("Output_tSNE_all/Output_tSNE_data/Output_ClusterPlots")
            dir.create(paste0("Output_ClusterPlots", "_seed", FlowSOM_seedA), showWarnings = FALSE)
            
            setwd(OutputDirectory)
            setwd(paste0("Output_tSNE_all/Output_tSNE_data/Output_ClusterPlots/", "Output_ClusterPlots", "_seed", FlowSOM_seedA))
            
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
          
          if (Run_FlowSOM_repeat == 1) {
           ### 9.4.1 - Second FlowSOM seed plots
            for (File in FileNames){ 
              setwd(OutputDirectory)
              setwd("Output_tSNE_all/Output_tSNE_data")
              
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
              
              cluster.name <- paste0(FlowSOM_meta_name, FlowSOM_seedB)
              clust.no <- CurrentSampleCSV[[cluster.name]]
              
              ## Create subdirectory
              setwd(OutputDirectory)
              setwd("Output_tSNE_all/Output_tSNE_data/Output_ClusterPlots")
              dir.create(paste0("Output_ClusterPlots", "_seed", FlowSOM_seedB), showWarnings = FALSE)
              
              setwd(OutputDirectory)
              setwd(paste0("Output_tSNE_all/Output_tSNE_data/Output_ClusterPlots/", "Output_ClusterPlots", "_seed", FlowSOM_seedB))
              
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
            
           ### 9.4.2 - Second FlowSOM seed plots
            for (File in FileNames){ 
              setwd(OutputDirectory)
              setwd("Output_tSNE_all/Output_tSNE_data")
              
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
              
              cluster.name <- paste0(FlowSOM_meta_name, FlowSOM_seedC)
              clust.no <- CurrentSampleCSV[[cluster.name]]
              
              ## Create subdirectory
              setwd(OutputDirectory)
              setwd("Output_tSNE_all/Output_tSNE_data/Output_ClusterPlots")
              dir.create(paste0("Output_ClusterPlots", "_seed", FlowSOM_seedC), showWarnings = FALSE)
              
              setwd(OutputDirectory)
              setwd(paste0("Output_tSNE_all/Output_tSNE_data/Output_ClusterPlots/", "Output_ClusterPlots", "_seed", FlowSOM_seedC))
              
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
            
          }
          
          ## Move back to PrimaryDirectory
          setwd(PrimaryDirectory)
          getwd() 
        }   
        
###################################################### 10. tSNE-select markers ######################################################                       
        
        if(Run_tSNE == 1){
          
        ### 10.1 - Create output directories
          
          setwd(OutputDirectory)
          dir.create("Output_tSNE_all+select", showWarnings = FALSE)
          setwd("Output_tSNE_all+select")
          dir.create("Output_tSNE_info", showWarnings = FALSE)
          dir.create("Output_tSNE_data", showWarnings = FALSE)
          setwd(OutputDirectory)
          
        ### 10.2 - Subsample full dataset (including FlowSOM cluster no, if FlowSOM was run)
          
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
              nsub <- tSNE.dwnsmp.targets[i] 
              data.temp <- subset(data, data[[samp.col]] == nam) # works
              nrow(data.temp)
              set.seed(tSNE_nsub_seed)
              data.temp <- data.temp[sample(1:nrow(data.temp), nsub), ]
              nrow(data.temp)
              data_subsampled <- rbind(data_subsampled, data.temp)
            }
            dim(data_subsampled)
          }
          
        ### 10.3 - Create set of column (parameter) names
          head(data_subsampled)                            # show data with headings
          
          ColumnNames_for_tSNE <- unname(colnames(data_subsampled))[ClusteringColNos] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10
          ColumnNames_for_tSNE  # check that the column names that appear are the ones you want to analyse
          
          dim(data_subsampled) # Review dimensionality of 'data' -- N events and N parameters
          data_specific <- data_subsampled[, ColumnNames_for_tSNE] # Prepare data for Rtsne --  select columns to use
          dim(data_specific) # Review dimensionality of 'data' -- N events and N parameters
          
        ### 10.4 - Save a copy of 'data' and 'data specific'        
          
          setwd(OutputDirectory)
          setwd("Output_tSNE_all+select/Output_tSNE_info")
          
          ## Export a list of parameters used to run tSNE and Phenograph
          tSNEparametersname <- paste0(data.name, "_data_and_columns_used_for_tSNE")
          tSNEparameters <- paste(tSNEparametersname, ".csv", sep = "")
          write.csv(x = data_specific, file = tSNEparameters)  
          
          ## Export subsampled data in CSV format (record of what was run)
          if(write.pre.tSNE == 1){
            subsampled_data_name <- paste0(data.name, "_subsampled_data")
            subsampled_data <- paste(subsampled_data_name, ".csv", sep = "")
            write.csv(x = data_subsampled, file = subsampled_data)
            
            
          }
          
        ### 10.5 - Run tSNE ('data' at this point is scaled, transformed) -- RUN ALL OF STEP 5
          
          ## Set working directory to "Output_info"
          setwd(OutputDirectory)
          setwd("Output_tSNE_all+select/Output_tSNE_info")
          
          ## Turn sink on (send the verbose text from the tSNE algorithm progress updates to a .txt file)
          verbose_name <- paste0("tSNE_verbose_", data.name, ".txt")
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
          
          
        ### 10.6 - Save info from tSNE run and add tSNE parameters to dataset
          
          ## Save tsne_out output info as .csv (i.e., output settings as raw results)
          tsne_out.df <- ldply(tsne_out, data.frame)
          output_name <- paste0("tSNE_parameter_values_", data.name, ".csv")
          fwrite(x = tsne_out.df, file = output_name, row.names = TRUE) # pretty blood good -- doesn't give row number for Y, costs, or itercosts -- but easy to figure out
          
          ## Add Y (tSNE coordinate) values to starting data
          tsne_out_Y <- tsne_out$Y
          head(tsne_out_Y)
          
          colnames(tsne_out_Y) <- c(tSNE_X_name, tSNE_Y_name)
          head(tsne_out_Y)
          
          ## save cluster labels
          data_subsampled <- cbind(data_subsampled, tsne_out_Y)
          head(data_subsampled)
          
        ### 10.7 - Write all data (with tSNE and FlowSOM parameters) to .csv and .fcs
          
          ## Set working directory
          setwd(OutputDirectory)
          setwd("Output_tSNE_all+select/Output_tSNE_data")
          
          if(write.tSNE.merged == 1){
            
            ## Save data (with new tSNE parameters) as CSV
            write.csv(x = data_subsampled, file = paste0(data.name, "_with_tSNE", ".csv"), row.names=FALSE)
            
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
            new_file_name_fcs <- paste0(data.name, "_with_tSNE", ".fcs")
            write.FCS(data_subsampled.ff, new_file_name_fcs)
            
            ### There is a delay here -- fcs file ends up in primary directory
          }
          
        ### 10.8 - Write GROUPED data (with tSNE and FlowSOM parameters) to .csv and .fcs
          if(write.tSNE.group == 1){
            for(a in AllGroupNames){
              data_subset <- subset(data_subsampled, data_subsampled[[grp.col]] == a)
              dim(data_subsampled)
              
              ## write .csv
              write.csv(data_subset, file = paste0(data.name, "_", a, "_with_tSNE.csv", sep = ""), row.names=FALSE)
              
              ## write .fcs
              metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
              metadata
              
              ## Create FCS file metadata - ranges, min, and max settings
              #metadata$range <- apply(apply(data_subset,2,range),2,diff)
              metadata$minRange <- apply(data_subset,2,min)
              metadata$maxRange <- apply(data_subset,2,max)
              
              data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
              head(data_subset.ff)
              write.FCS(data_subset.ff, paste0(data.name, "_", a, "_with_tSNE", ".fcs"))
            }
          }
          
        ### 10.9 - Write individual files (with tSNE and FlowSOM parameters) to .csv and .fcs     
          
          if(write.tSNE.sep == 1){
            
            for (a in AllSampleNames) {
              data_subset <- subset(data_subsampled, data_subsampled[[samp.col]] == a)
              
              ## write .csv
              write.csv(data_subset, file = paste0(data.name, "_", a, "_with_tSNE.csv", sep = ""), row.names=FALSE)
              
              ## write .fcs
              metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
              
              ## Create FCS file metadata - ranges, min, and max settings
              #metadata$range <- apply(apply(data_subset,2,range),2,diff)
              metadata$minRange <- apply(data_subset,2,min)
              metadata$maxRange <- apply(data_subset,2,max)
              
              data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
              head(data_subset.ff)
              write.FCS(data_subset.ff, paste0(data.name, "_", a, "_with_tSNE", ".fcs"))
            }
            
            ## Move back to PrimaryDirectory
            setwd(PrimaryDirectory)
            getwd()
          }
          
          # note -- in the final data output, all parameters are included, but only the subsampled and/or transformed cells
          
        }
        
        
###################################################### 11. tSNEplots-all markers ######################################################  
        
        if(Run_tSNEplots == 1){    
          
          ## Create 'jet' colour scheme (not available by default in R)
          jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
          
          ## Set your working directory here (e.g. "/Users/Tom/Desktop/")
          setwd(OutputDirectory)
          setwd("Output_tSNE_all+select/Output_tSNE_data") 
          
          ## Assign the working directory as 'PrimaryDirectory'
          PlotDirectory <- getwd()
          PlotDirectory
          
          ## Create a list of file names (names of the samples) and check file names
          PlotFiles <- list.files(path=PlotDirectory, pattern = ".csv") # ALTERNATIVE: list.files(pattern = '\\.csv')
          PlotFiles
          
          ## In the output of the previous line, you will see the names for the tSNE parameters -- insert them in between the "" below
          plotXname <- tSNE_X_name
          plotYname <- tSNE_Y_name
          
          
        ### 11.1 - ESTABLISH GLOBAL SCALE LIMITS FOR COLOUR and XY #####
          
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
          
          
        ### 11.2 - Loop with  samples in separate folders #####  
          ## First, change the tSNE parameters (on lines 103 and 104)
          ## Then run all of the script below
          
          ## Set wd
          setwd(PlotDirectory)
          getwd()
          
          ## Create output folder (if a folder called "Output" already exists, nothing will happen)
          dir.create("Output(tSNEplots-all+select)", showWarnings = FALSE)
          
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
            setwd("Output(tSNEplots-all+select)") 
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
        
        
###################################################### 12. ClusterPlots-select ######################################################  
        
        if(Run_ClusterPlots == 1){
          
        ### 12.1 - Setup
          
          ## WD and colours
          jet.black <- colorRampPalette(c("black", "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000", "purple"))
          colour.scheme <- jet.black
          
          setwd(OutputDirectory)
          setwd("Output_tSNE_all+select/Output_tSNE_data")
          dir.create("Output_NumPlots", showWarnings = FALSE)
          
          FileNames <- list.files(path=getwd(), pattern = ".csv")
          FileNames
          
          ## In the output of the previous line, you will see the names for the tSNE parameters -- insert them in between the "" below
          plotXname <- tSNE_X_name
          plotYname <- tSNE_Y_name
          cluster.name <- FlowSOM_meta_name
          
          cluspl.n.sub # already defined
          cluspl.n.sub.seed # already defined
          
        ### 12.2 - Read files 
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
          
        ### 12.3 - Loop for Number Plot
          
          for (File in FileNames){ 
            setwd(OutputDirectory)
            setwd("Output_tSNE_all+select/Output_tSNE_data")
            
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
            
            cluster.name <- paste0(FlowSOM_meta_name, FlowSOM_seedA)
            clust.no <- CurrentSampleCSV[[cluster.name]]
            
            ## Subdirectory
            setwd(OutputDirectory)
            setwd("Output_tSNE_all+select/Output_tSNE_data/Output_NumPlots")
            dir.create(paste0("Output_NumPlots", "_seed", FlowSOM_seedA), showWarnings = FALSE)
            
            setwd(OutputDirectory)
            setwd(paste0("Output_tSNE_all+select/Output_tSNE_data/Output_NumPlots/", "Output_NumPlots", "_seed", FlowSOM_seedA))
            
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
          
          if (Run_FlowSOM_repeat == 1){
            ### 12.3.1 - Second FlowSOM seed plots
            for (File in FileNames){ 
              setwd(OutputDirectory)
              setwd("Output_tSNE_all+select/Output_tSNE_data")
              
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
              
              cluster.name <- paste0(FlowSOM_meta_name, FlowSOM_seedB)
              clust.no <- CurrentSampleCSV[[cluster.name]]
              
              ## Subdirectory
              setwd(OutputDirectory)
              setwd("Output_tSNE_all+select/Output_tSNE_data/Output_NumPlots")
              dir.create(paste0("Output_NumPlots", "_seed", FlowSOM_seedB), showWarnings = FALSE)
              
              setwd(OutputDirectory)
              setwd(paste0("Output_tSNE_all+select/Output_tSNE_data/Output_NumPlots/", "Output_NumPlots", "_seed", FlowSOM_seedB))
              
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
            
            ### 12.3.2 - Third FlowSOM seed plots
            for (File in FileNames){ 
              setwd(OutputDirectory)
              setwd("Output_tSNE_all+select/Output_tSNE_data")
              
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
              
              cluster.name <- paste0(FlowSOM_meta_name, FlowSOM_seedC)
              clust.no <- CurrentSampleCSV[[cluster.name]]
              
              ## Subdirectory
              setwd(OutputDirectory)
              setwd("Output_tSNE_all+select/Output_tSNE_data/Output_NumPlots")
              dir.create(paste0("Output_NumPlots", "_seed", FlowSOM_seedC), showWarnings = FALSE)
              
              setwd(OutputDirectory)
              setwd(paste0("Output_tSNE_all+select/Output_tSNE_data/Output_NumPlots/", "Output_NumPlots", "_seed", FlowSOM_seedC))
              
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
            
          }
          
          
        ### 12.4 - CLUSTER COLOUR PLOT
          
          ## Set wd
          setwd(OutputDirectory)
          setwd("Output_tSNE_all+select/Output_tSNE_data")
          dir.create("Output_ClusterPlots", showWarnings = FALSE)
          
          for (File in FileNames){ 
            setwd(OutputDirectory)
            setwd("Output_tSNE_all+select/Output_tSNE_data")
            
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
            
            cluster.name <- paste0(FlowSOM_meta_name, FlowSOM_seedA)
            clust.no <- CurrentSampleCSV[[cluster.name]]
            
            ## Create subdirectory
            setwd(OutputDirectory)
            setwd("Output_tSNE_all+select/Output_tSNE_data/Output_ClusterPlots")
            dir.create(paste0("Output_ClusterPlots", "_seed", FlowSOM_seedA), showWarnings = FALSE)
            
            setwd(OutputDirectory)
            setwd(paste0("Output_tSNE_all+select/Output_tSNE_data/Output_ClusterPlots/", "Output_ClusterPlots", "_seed", FlowSOM_seedA))
            
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
          
          if (Run_FlowSOM_repeat == 1) {
            ### 12.4.1 - Second FlowSOM seed plots
            for (File in FileNames){ 
              setwd(OutputDirectory)
              setwd("Output_tSNE_all+select/Output_tSNE_data")
              
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
              
              cluster.name <- paste0(FlowSOM_meta_name, FlowSOM_seedB)
              clust.no <- CurrentSampleCSV[[cluster.name]]
              
              ## Create subdirectory
              setwd(OutputDirectory)
              setwd("Output_tSNE_all+select/Output_tSNE_data/Output_ClusterPlots")
              dir.create(paste0("Output_ClusterPlots", "_seed", FlowSOM_seedB), showWarnings = FALSE)
              
              setwd(OutputDirectory)
              setwd(paste0("Output_tSNE_all+select/Output_tSNE_data/Output_ClusterPlots/", "Output_ClusterPlots", "_seed", FlowSOM_seedB))
              
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
            
          ### 12.4.2 - Second FlowSOM seed plots
            for (File in FileNames){ 
              setwd(OutputDirectory)
              setwd("Output_tSNE_all+select/Output_tSNE_data")
              
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
              
              cluster.name <- paste0(FlowSOM_meta_name, FlowSOM_seedC)
              clust.no <- CurrentSampleCSV[[cluster.name]]
              
              ## Create subdirectory
              setwd(OutputDirectory)
              setwd("Output_tSNE_all+select/Output_tSNE_data/Output_ClusterPlots")
              dir.create(paste0("Output_ClusterPlots", "_seed", FlowSOM_seedC), showWarnings = FALSE)
              
              setwd(OutputDirectory)
              setwd(paste0("Output_tSNE_all+select/Output_tSNE_data/Output_ClusterPlots/", "Output_ClusterPlots", "_seed", FlowSOM_seedC))
              
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
            
          }
          
          ## Move back to PrimaryDirectory
          setwd(PrimaryDirectory)
          getwd() 
        }   
        
        
