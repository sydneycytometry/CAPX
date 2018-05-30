### Cytometry Analysis Pipeline for Complex Datasets (CAPX) part 1 - data preprocessing
    # Thomas Ashhurst
    # 2018-05-17
    # thomas.ashhurst@sydney.edu.au
    # www.github.com/SydneyCytometry

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

    ### 1.1 - Install packages (if not already installed)
        if(!require('flowCore')) {install.packages('flowCore')}
        if(!require('plyr')) {install.packages('plyr')}
        if(!require('Biobase')) {install.packages('Biobase')}
        if(!require('data.table')) {install.packages('data.table')}
        if(!require('rstudioapi')) {install.packages('rstudioapi')}
        
    ### 1.2 Load packages       
        library('flowCore')
        library('plyr')
        library('Biobase')
        library('data.table')
        library('rstudioapi')

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
        
        ## Use to list the .fcs or .csv files in the working directory -- important, the only FCS/CSV file in the directory should be the one desired for analysis. If more than one are found, only the first file will be used
        list.files(path=PrimaryDirectory, pattern = ".fcs")     # see a list of FCS files
        list.files(path=PrimaryDirectory, pattern = ".csv")     # see a list of CSV files
        
        file.type               <- ".csv" # Specfiy file type (".csv" or ".fcs")
        data.name               <- "BM_20180518_1100"   # suggest name is sampletype_date_time (e.g. liver_20180203_1400)

        arcsinh.transform       <- 0      # 0 = no, 1 = yes # default = 0 
        asinh_scale             <- 15     # for CyTOF choose between 5 and 15, for flow you will need between 200 and 2000.

        remove.duplicates       <- 1      # removing duplicates aids in tSNE analysis
        
        write.merged.file       <- 1      # do you want to write one large merged file?
        write.sep.files         <- 1      # do you also want to write indivdual files for each sample?
        write.group.files       <- 1      # do you also want to write one file for each group? (requires creation of group keywords)
        
        #add.group.keywords     <- 1      # inactive
        #downsample.files       <- 1      # inactive
        
###################################################### 3. USER INPUT - LINE BY LINE ###################################################### 
        
    ### 3.1 - Read files into a list of dataframes (run all of 2.2, no user input)
        
        ## Check the list of files
        FileNames <- list.files(path=PrimaryDirectory, pattern = file.type) # Generate list of files
        as.matrix(FileNames) # See file names in a list
        
        ## Read data from Files into list of data frames
        DataList=list() # Creates and empty list to start 
        
        if (file.type == ".csv"){
          for (File in FileNames) { # Loop to read files into the list
            tempdata <- read.csv(File)
            File <- gsub(".csv", "", File)
            DataList[[File]] <- tempdata
          }
        }
        
        if (file.type == ".fcs"){
          ParaList=list()
          for (File in FileNames) { # Loop to read files into the list
            tempdata <- exprs(read.FCS(FileNames[1], transformation = FALSE)) 
            File <- gsub(".fcs", "", File)
            DataList[[File]] <- tempdata
            ParaList[[File]] <- parameters(read.FCS(File)) 
            ParaList[[File]]
          }
        }
    

    ### 3.2 - OPTIONAL - Remove troubleshom columns (if required) and add new columns for SampleName and SampleNo
        
        ## Check the number of columns in each file is consistent
        Length_check = list() # creates an empty list
        for(i in c(1:(length(DataList)))){Length_check[[i]] <- length(names(DataList[[i]]))} # creates a list of the number of columns in each sample
        
        as.matrix(Length_check) # check that the number of columns in each sample is the same length
        
        ## Check the names of the columns in each file is consistent
        ColName_check = list() 
        for(i in c(1:(length(DataList)))){ColName_check[[i]] <- names(DataList[[i]])}
        
        ColName_check # check that the number of columns in each sample is the same length
        
        ############### ONLY IF REQUIRED ###############
        #    # Remove any troublesome columns (if required)
        #    for (i in c(1:(length(DataList)))) {
        #      DataList[[i]]$SampleID <- NULL # after the '$', put column name here
        #    }
        ################################################ 
        
        ## Final check -- ensure the number of columns in each file is consistent
        Length_check = list() # creates an empty list
        for(i in c(1:(length(DataList)))){Length_check[[i]] <- length(names(DataList[[i]]))} # creates a list of the number of columns in each sample
        as.matrix(Length_check) # check that the number of columns in each sample is the same length
        
        ## Review the kind of numerical entries in the data
        head(DataList[[1]])
        
        
    ### 3.3 - IF REQUIRED - Arcsinh transformation (data may have already been transformed, SKIP IF NOT REQUIRED)        

        ## Review paramters (columns)
        col.names.dl <- names(DataList[[1]]) # show data with headings
        as.matrix(col.names.dl) # view the column 'number' for each parameter

        ## IF required, specify column numbers to be transformed - e.g. c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10. Alternatively, e.g. c(1:32) to include all 32 columns
        col.nos.scale <- c(2:32) 
        col.nos.scale
        
        if(arcsinh.transform == 1){
            col.names.SCALE <- col.names.dl[col.nos.scale]
            var <- col.names.SCALE  # check that the column names that appear are the ones you want to analyse
    
            for (i in c(1:(length(DataList)))) {
              DataList[[i]][,var] <- sapply(DataList[[i]][,var], function(x) x = x/asinh_scale)
              }
            #DataList <- asinh(DataList / asinh_scale)    # transforms all columns! including event number etc
        }
        
        DataList[[1]] # check the data has been converted 
        
            
            
         
    ### 3.4 - ## Add sample identifiers (necessary to merge sample)
        
        #### REMOVE the ".csv" at the end of the file name
         
        ## Create a list of 'SampleName' and SampleNo' entries -- 1:n of samples, these will be matched in order
        AllSampleNames <- c(names(DataList))
        AllSampleNames # Character (words)
        
        AllSampleNos <- c(1:(length(DataList)))       ### <-- changed to 4+ to match with sample
        AllSampleNos # Intiger/numeric (numbers)
        
        ## Add 'SampleNo' and SampleName' to each 
        for (i in AllSampleNos) {DataList[[i]]$SampleNo <- i}
        for (a in AllSampleNames) {DataList[[a]]$SampleName <- a} # Inserting names doesn't stop the rest working, just messes with some auto visual stuff
        head(DataList[[1]])
        ## Alternative: matrix - 1st column use for 'i', matrix[1]$i or something -- then use matrix[2]$i for the name?
        
        head(DataList[[1]]) # check to see that the sample name and number keywords have been entered on each row
        
        as.matrix(AllSampleNos)
        as.matrix(AllSampleNames)
        
        
        
    ### 3.5 - OPTIONAL - Add 'GroupNo' and 'GroupName' to each (required if generating 'group' files)
                
        ## Setup 'GroupNo' and 'GroupName' for samples
        Group_1_name <- "Mock_D7" # enter the name for the first group     
        Group_2_name <- "WNV_D7" 
        
        Group_1_range <- c(1:6) # enter the sample numbers that belong to group 1. E.g. c(1:6), or c(1,4,5) etc.
        Group_2_range <- c(7:12) # enter the sample numbers that belong to group 2. E.g. c(7:12), or c(7,9,12) etc.       
        
        Group_names <- c(Group_1_name, Group_2_name)
        Group_names
        
        ####################################################
        ###### Add any additional groups if required #######
        ####################################################
        
        ## Group 1 loop
        for(i in Group_1_range){
          DataList[[i]]$GroupNo <- 1    
          DataList[[i]]$GroupName <- Group_1_name
        }
        
        ## Group 2 loop
        for(i in Group_2_range){
          DataList[[i]]$GroupNo <- 2    
          DataList[[i]]$GroupName <- Group_2_name
        }
        
        ####################################################
        #### Create new loop for each additional group #####
        ####################################################
        
        ## Check column names
        head(DataList[[1]]) # check the entries for one of the samples in group 1
        head(DataList[[10]]) # check the entries for one of the samples in group 2
        

        
    ### 3.6 - OPTIONAL - Specify downsample targets

        ## Check number of cells present in each sample
        nrow.check = list() # creates an empty list
        for(i in c(1:(length(DataList)))){nrow.check[[i]] <- nrow(DataList[[i]])} # creates a list of the number of rows (cells) in each sample
        
        as.matrix(nrow.check) # shows the number of rows (cells) in each sample
        
        ## Specify downsampling targets (can skip if not downsampling data)
        DownSampleTargets <- c(8000,
                               8000,
                               8000,
                               8000,
                               8000,
                               8000,
                               8000,
                               8000,
                               8000,
                               8000,
                               8000,
                               8000
                               )
        
        ## Review you have the right number of DownSampleTargets entires
        length(DownSampleTargets)
        length(DownSampleTargets) == length(AllSampleNos) # Should return TRUE if you have the same number of targets as samples
        
        ## Review the downsample target values
        dwnsmpl.check = list() 
        for(i in c(1:(length(DownSampleTargets)))){dwnsmpl.check[[i]] <- (DownSampleTargets[[i]])} # creates a list of each of the downsample targets
        
        as.matrix(dwnsmpl.check) # shows each of the downsample targets
        
        ## Perform subsampling (if elected to)
        for (i in AllSampleNos) {
          
          nsub <- DownSampleTargets[i]
          set.seed(123)
          datalisttemp <- DataList[[i]]
          datalisttemp <- datalisttemp[sample(1:nrow(datalisttemp), nsub), ]
          DataList[[i]] <- datalisttemp
        }
        
        Downsample_check = list()
        for(i in AllSampleNos){
          Downsample_check[[i]] <- DownSampleTargets[i] == nrow(DataList[[i]])
        }
        
        ## Should all return 'TRUE' if successfully downsampled
        as.matrix(Downsample_check)
 
        
        
    ### 3.7 - Merge data and remove duplicates
        
        ## Concatenate into one large data frame
        data <- rbindlist(DataList) 
        data
        
        ## Remove duplicates
        dim(data) # note the second entry, which is number of parameters

        if(remove.duplicates == 1){
          data <- data[!duplicated(data), ]  # remove rows containing duplicate values within rounding
        }
        
        dim(data) # check to see if the number of parameters has reduced (no change means no duplicates were present)
        

      
###################################################### 4. END USER INPUT - RUN ALL AT ONCE ###################################################### 


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
          write.csv(x = data, file = csv.filename)
  
          ## write .fcs
          metadata <- data.frame(name=dimnames(data)[[2]],desc=paste('column',dimnames(data)[[2]],'from dataset'))
          metadata
          
          ## Create FCS file metadata - ranges, min, and max settings
          #metadata$range <- apply(apply(data,2,range),2,diff)
          metadata$minRange <- apply(data,2,min)
          metadata$maxRange <- apply(data,2,max)
          
          data.ff <- new("flowFrame",exprs=as.matrix(data), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
          head(data.ff)
          write.FCS(data.ff, paste0(data.name, ".fcs"))
        }
        
        
        
    ### 4.2 - Export data in individual files -- both .csv and .fcs format

        if(write.sep.files == 1){
          
          for (a in AllSampleNames) {
            
            data_subset <- subset(data, SampleName == a)
            dim(data_subset)
            
            ## write .csv
            write.csv(data_subset, file = paste(data.name, "_", a,".csv", sep = ""))
            
            ## write .fcs
            metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
            metadata
            
            ## Create FCS file metadata - ranges, min, and max settings
            #metadata$range <- apply(apply(data_subset,2,range),2,diff)
            metadata$minRange <- apply(data_subset,2,min)
            metadata$maxRange <- apply(data_subset,2,max)
            
            data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
            head(data_subset.ff)
            write.FCS(data_subset.ff, paste0(data.name, "_", a, ".fcs"))
          }
        }
         
        
         
    ### 4.3 - Export data as grouped files -- both .csv and .fcs format
                  
        if(write.group.files == 1){
          
          for(a in Group_names){
            data_subset <- subset(data, GroupName == a)
            dim(data_subset)
            
            ## write .csv
            write.csv(data_subset, file = paste(data.name, "_", a, ".csv", sep =))
            
            ## write .fcs
            metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
            metadata
            
            ## Create FCS file metadata - ranges, min, and max settings
            #metadata$range <- apply(apply(data_subset,2,range),2,diff)
            metadata$minRange <- apply(data_subset,2,min)
            metadata$maxRange <- apply(data_subset,2,max)
            
            data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
            head(data_subset.ff)
            write.FCS(data_subset.ff, paste0(data.name, "_", a, ".fcs"))
          }
        }
        
        
    
