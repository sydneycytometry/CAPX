### Cytometry Analysis Pipeline for Complex Datasets (CAPX) part 3 - annotate

    # Thomas Ashhurst
    # 2018-05-30
    # thomas.ashhurst@sydney.edu.au
    # www.github.com/SydneyCytometry/CAPX

### Summary
  
    # 1 INSTALL AND LOAD PACKAGES
    # 2 READ FILES
    # 3 ANNOTATE
    # 4 WRITE FILES

###################################################### 1. INSTALL AND-OR LOAD PACKAGES ###################################################### 
    
    ### 1.1 - Install packages (if not already installed)
        if(!require('flowCore')) {install.packages('flowCore')}
        if(!require('plyr')) {install.packages('plyr')}
        if(!require('data.table')) {install.packages('data.table')}
        if(!require('rstudioapi')) {install.packages('rstudioapi')}
        if(!require('R.utils')) {install.packages('R.utils')}
        if(!require('ggplot2')) {install.packages('ggplot2')}
        if(!require('gplots')) {install.packages('gplots')}
        if(!require('RColorBrewer')) {install.packages('RColorBrewer')}
        if(!require('tidyr')) {install.packages('tidyr')}
        if(!require('Biobase')) {install.packages('Biobase')}

    ### 1.2 Load packages       
        library('flowCore')
        library('plyr')
        library('data.table')
        library('rstudioapi')
        library('R.utils')
        library('ggplot2')
        library('gplots')
        library('RColorBrewer')
        library('tidyr')
        library('Biobase')

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

    ## Preferences
        Data_name             <- "demo_data"
    
    ## Preferences: write .csv and .fcs files    
        write.merged.file       <- 1              # do you want to write one large merged file?
        write.sep.files         <- 1              # do you also want to write indivdual files for each sample?
        write.group.files       <- 1              # do you also want to write one file for each group? (requires creation of group keywords)
        
    
###################################################### 2. READ FILES ###################################################### 
    
    # DO ON MAIN FlowSOM FILES, preferably not tSNE files -- numbers too small, differential downsampling
    
    ## Use to list the .fcs or .csv files in the working directory
        list.files(path=PrimaryDirectory, pattern = ".fcs")     # see a list of FCS files
        list.files(path=PrimaryDirectory, pattern = ".csv")     # see a list of CSV files
        
    ## Specify general options
        File_Type             <- ".csv"               # Specify the type of input file. Use ".fcs" OR ".csv", in lower case      # fcs files not current supported
        FileNos               <- "single"             # Can be "single" for one file, or "multi" for multiple files. If "single" the only the file desired for analysis should be in the working directory, otherwise the script will use the first file it finds.

    ## Loop for single and loop for multiple files
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
    
    ## Choose the column that defines the sample NAMES
        as.matrix(names(data))
        samp.col <- "SampleName"
        head(data[samp.col]) # check that this is actually the sample names
    
    ## Choose the column that defines the GROUPS NAMES
        as.matrix(names(data))
        grp.col <- "GroupName"
        head(data[grp.col]) # check that this is actually the group names
    
    ## Choose the column that defines the CLUSTERS you want to annotate
        as.matrix(names(data))
        clust.col <- "FlowSOM_meta_cluster"
        head(data[grp.col]) # check that this is actually the group names
    
    ## Choose the a name for your new 'population name' column
        pop.col <- "PopulationName"
    
    # Duplicate 'cluster number' column, and rename to population ID column
        data[[pop.col]] = data[[clust.col]] 
        
    
###################################################### 3. Annotate ###################################################### 

    ## EXPLORE the tSNE plots and HEATMAPS to work out what clusters reflect what populations
        
        max(data$FlowSOM_meta_cluster) # max cluster number
        min(data$FlowSOM_meta_cluster) # min cluster number
        
        pop.names <- c(n1 <- "1", 
                       n2 <- "2", 
                       n3 <- "3",
                       n4 <- "4",
                       n5 <- "5",
                       n6 <- "6",
                       n7 <- "7",
                       n8 <- "8", 
                       n9 <- "9",
                       n10 <- "10",
                       n11 <- "11",
                       n12 <- "12",
                       n13 <- "13",
                       n14 <- "14",
                       n15 <- "15",
                       n16 <- "16",
                       n17 <- "17",
                       n18 <- "18",
                       n19 <- "19",
                       n20 <- "20")
        
        pop.names

        for(i in sort(c(unique(data[[clust.col]])))){ 
          data[[pop.col]][data[[clust.col]] == i] <- pop.names[i]
        }
        
        head(data)
            

############################################################################################################################
###################################################### END USER INPUT ######################################################         
############################################################################################################################         
        
        
###################################################### 4. Write files ###################################################### 
        
    ### Write .csv and .fcs files

        AllSampleNames <- unique(data[[samp.col]])
        AllGroupNames <- unique(data[[grp.col]])
        
        setwd(PrimaryDirectory)
        dir.create("Output_rename", showWarnings = FALSE)
        setwd("Output_rename")
        output.directory <- getwd()
        
    ### Write all data (with tSNE and FlowSOM parameters) to .csv and .fcs
        
        if(write.merged.file == 1){
          
          ## Set working directory
          setwd(output.directory)
          
          ## Save data (with new tSNE parameters) as CSV
          write.csv(x = data, file = paste0(Data_name, ".csv"), row.names=FALSE)
          
          ## Check data and data column names
          head(data)
          dimnames(data)[[2]]
          
          ## Create FCS file metadata - column names with descriptions
          metadata <- data.frame(name=dimnames(data)[[2]],
                                 desc=paste('column',dimnames(data)[[2]],'from dataset')
          )
          
          ## Create FCS file metadata - ranges, min, and max settings
          #metadata$range <- apply(apply(data,2,range),2,diff)
          metadata$minRange <- apply(data,2,min)
          metadata$maxRange <- apply(data,2,max)
          
          ## Create flowframe with tSNE data
          data.ff <- new("flowFrame",
                                    exprs=as.matrix(data), # in order to create a flow frame, data needs to be read as matrix
                                    parameters=AnnotatedDataFrame(metadata)
          )
          
          ## Check flow frame
          data.ff
          
          ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
          new_file_name_fcs <- paste0(Data_name, ".fcs")
          write.FCS(data.ff, new_file_name_fcs)
          
          ### There is a delay here -- fcs file ends up in primary directory
        }
        
    ### Write GROUPED data (with tSNE and FlowSOM parameters) to .csv and .fcs
        if(write.group.files == 1){
          for(a in AllGroupNames){
            data_subset <- subset(data, data[[grp.col]] == a)
            dim(data_subset)
            
            ## write .csv
            write.csv(data_subset, file = paste(Data_name, "_", a, ".csv", sep = ""), row.names=FALSE)
            
            ## write .fcs
            metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
            metadata
            
            ## Create FCS file metadata - ranges, min, and max settings
            #metadata$range <- apply(apply(data_subset,2,range),2,diff)
            metadata$minRange <- apply(data_subset,2,min)
            metadata$maxRange <- apply(data_subset,2,max)
            
            data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
            head(data_subset.ff)
            write.FCS(data_subset.ff, paste0(Data_name, "_", a, ".fcs"))
          }
        }
    
    ### Write individual files (with tSNE and FlowSOM parameters) to .csv and .fcs     
        if(write.sep.files == 1){
          for (a in AllSampleNames) {
            data_subset <- subset(data, data[[samp.col]] == a)
            dim(data_subset)
            
            ## write .csv
            write.csv(data_subset, file = paste(Data_name, "_", a, ".csv", sep = ""), row.names=FALSE)
            
            ## write .fcs
            metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
            
            ## Create FCS file metadata - ranges, min, and max settings
            #metadata$range <- apply(apply(data_subset,2,range),2,diff)
            metadata$minRange <- apply(data_subset,2,min)
            metadata$maxRange <- apply(data_subset,2,max)
            
            data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
            head(data_subset.ff)
            write.FCS(data_subset.ff, paste0(Data_name, "_", a, ".fcs"))
          }
        }
        
    ### Move back to PrimaryDirectory
        setwd(PrimaryDirectory)
        getwd()
        


        
        
