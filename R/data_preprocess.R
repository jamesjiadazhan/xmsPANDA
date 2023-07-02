#' @title data_preprocess
#' @description This function performs data transformation, normalization
#' @usage 
#' data_preprocess(Xmat = NA, Ymat = NA, feature_table_file = NA,
#'                 parentoutput_dir, class_labels_file = NA,
#'                 num_replicates = 3, feat.filt.thresh = NA,
#'                 summarize.replicates = TRUE, summary.method = "mean",
#'                 all.missing.thresh = 0.5, group.missing.thresh = 0.7,
#'                 log2transform = TRUE, medcenter = FALSE,
#'                 znormtransform = FALSE, quantile_norm = TRUE,
#'                 lowess_norm = FALSE, rangescaling = FALSE,
#'                 paretoscaling = FALSE, mstus = FALSE, sva_norm =
#'                 FALSE, TIC_norm = FALSE, eigenms_norm = FALSE,
#'                 madscaling = FALSE, vsn_norm = FALSE, cubicspline_norm
#'                 = FALSE, missing.val = 0, samplermindex = NA,
#'                 rep.max.missing.thresh = 0.5, summary.na.replacement =
#'                 "zeros", featselmethod = NA, pairedanalysis = FALSE,
#'                 normalization.method = "none", input.intensity.scale =
#'                 "raw", create.new.folder = TRUE)
#' @param Xmat R object for feature table. If this is given, it will be read directly. If not, read.table function will be used to read the feature table. 
#' @param Ymat R object for response/class labels matrix. If this is given, then class can be set to NA.
#' @param feature_table_file Feature table that includes the mz, retention time, and measured intensity in each sample 
#' for each analyte. The first 2 columns should be the mz and time. The remaining columns
#' should correspond to the samples in the class labels file with each column including the intensity profile
#' of a sample. Full path required.
#' @param parentoutput_dir Provide full path of the folder where you want the results to be written.
#' @param class_labels_file File with class labels information for each sample. Samples should be in the same order
#' as in the feature table. Please use the same format as in the example folder.
#' @param num_replicates Number of technical replicates
#' @param feat.filt.thresh Percent Intensity Difference or Coefficient of variation threshold; feature filtering
#' Use NA to skip this step.
#' @param summarize.replicates Do the technical replicates per sample need to be averaged or median summarized?
#' @param summary.method Method for summarizing the replicates. Options: "mean" or "median"
#' @param summary.na.replacement How should the missing values be represented? 
#' Options: "zeros", "halffeaturemin", "halfsamplemin","halfdatamin", "none"
#' "zeros": replaces missing values by 0
#' "halfsamplemin": replaces missing value by one-half of the lowest signal intensity in the
#' corresponding sample
#' "halfdatamin": replaces missing value by one-half of the lowest signal intensity in the
#' complete dataset
#' "halffeaturemin": replaces missing value by one-half of the lowest signal intensity for the
#' current feature
#' "none": keeps missing values as NAs
#' Users are recommended to perform imputation prior to performing biomarker discovery.
#' @param missing.val How are the missing values represented in the input data? Options: "0" or "NA"
#' @param samplermindex Column index of any additional or irrelevant columns to be deleted.
#' Options: "NA" or list of column numbers. eg: c(1,3,4) Default=NA
#' @param rep.max.missing.thresh What propotion of replicates are allowed to have missing values during the averaging or 
#' median summarization step of each biological sample? If the number of replicates with
#' missing values is greater than the defined threshold, then the summarized value is 
#' represented by the "missing.val" parameter. If the number of replicates with missing values
#' is less than or equal to the defined threshold, then the summarized value is equal to the 
#' mean or the median of the non-missing values. Default: 0.5
#' @param all.missing.thresh What propotion of total number of samples should have an intensity?
#' Default: 0.5
#' @param group.missing.thresh What propotion of samples in either of the two groups should have an intensity?
#' If at least x% of the samples in either group have a signal, then the feature is retained
#' for further analysis. Default: 0.7
#' @param log2transform Data transformation: Please refer to http://www.biomedcentral.com/1471-2164/7/142
#' Try different combinations; such as log2transform=TRUE, znormtransfrom=FALSE
#' or log2transform=FALSE, znormtransfrom=TRUE
#' @param medcenter Median centering of metabolites
#' @param znormtransform Auto scaling; each metabolite will have a mean of 0 and unit variance
#' @param quantile_norm Performs quantile normalization. Normalization options: Please set only one of the options to be TRUE
#' @param lowess_norm Performs lowess normalization. Normalization options: Please set only one of the options to be TRUE
#' @param madscaling Performs median adjusted scale normalization. Normalization options: Please set only one of the options to be TRUE
#' @param vsn_norm Boolean flag for VSN normalization (default is FALSE)
#' @param cubicspline_norm Boolean flag for cubic spline normalization (default is FALSE)
#' @param missing.val Value to replace missing data with (default is 0)
#' @param samplermindex Sample RM Index (default is NA)
#' @param rep.max.missing.thresh Maximum missing threshold for replicate (default is 0.5)
#' @param summary.na.replacement Method to replace NA in summary (default is "zeros")
#' @param featselmethod Feature selection method (default is NA)
#' @param pairedanalysis Boolean flag for paired analysis (default is FALSE)
#' @param input.intensity.scale Intensity scale of the input (default is "raw")
#' @param create.new.folder Boolean flag to create new folder (default is TRUE)
#' @param log2.transform.constant Constant for log2 transformation.
#' @param alphabetical.order Boolean flag for alphabetical order (default is FALSE)
#' @return Pre-processed data matrix.
#' @author Karan Uppal <kuppal3gt@gmail.com>, Jiada James Zhan <jzha832@emory.edu>

data_preprocess <- function(
  Xmat = NA,
  Ymat = NA,
  feature_table_file = NA,
  parentoutput_dir,
  class_labels_file = NA,
  num_replicates = 3,
  feat.filt.thresh = NA,
  summarize.replicates = TRUE,
  summary.method = "mean",
  all.missing.thresh = 0.1,
  group.missing.thresh = 0.8,
  normalization.method = "none",
  log2transform = FALSE,
  medcenter = FALSE,
  znormtransform = FALSE,
  quantile_norm = FALSE,
  lowess_norm = FALSE,
  rangescaling = FALSE,
  paretoscaling = FALSE,
  mstus = FALSE,
  sva_norm = FALSE,
  TIC_norm = FALSE,
  eigenms_norm = FALSE,
  madscaling = FALSE,
  vsn_norm = FALSE,
  cubicspline_norm = FALSE,
  missing.val = 0,
  samplermindex = NA,
  rep.max.missing.thresh = 0.5,
  summary.na.replacement = "zeros",
  featselmethod = NA,
  pairedanalysis = FALSE,
  input.intensity.scale = "raw",
  create.new.folder = TRUE,
  log2.transform.constant = 1,
  alphabetical.order = FALSE
) {

  # Suppress warnings
  options(warn = -1)

  # Read the data file or use the given data matrix
  data_matrix <- ifelse(
    is.na(Xmat),
    utils::read.table(Xmat,sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE),
    Xmat
  )

  # Set working directory and create necessary folders
  setwd(parentoutput_dir)

  if (create.new.folder) {
    dir.create("Stage1", showWarnings = FALSE)
    setwd("Stage1")
  } else {
    dir.create("Tables", showWarnings = FALSE)
    setwd("Tables")
  }

  # Normalize the normalization.method parameter
  normalization.method=tolower(normalization.method)

  # Simplified the logical checks for setting normalization methods
  log2transform = normalization.method %in% c("log2quantilenorm", "log2transform", "sva_norm", "eigenms_norm")
  quantile_norm = normalization.method == "log2quantilenorm" || normalization.method == "quantile_norm"
  znormtransform = normalization.method == "znormtransform"
  lowess_norm = normalization.method == "lowess_norm"
  rangescaling = normalization.method == "rangescaling"
  paretoscaling = normalization.method == "paretoscaling"
  mstus = normalization.method == "mstus"
  sva_norm = normalization.method == "sva_norm"
  eigenms_norm = normalization.method == "eigenms_norm"
  vsn_norm = normalization.method == "vsn_norm"
  TIC_norm = normalization.method == "tic_norm"
  cubicspline_norm = normalization.method == "cubicspline_norm"
  madscaling = normalization.method == "mad_norm"
  
  cnames<-colnames(data_matrix)
  cnames<-tolower(cnames)
  
  
  # If the input intensity scale is "log2", then we will not perform the log2 transformation
  if (input.intensity.scale == "log2") {
    log2transform = FALSE
  }

  # Check if the data matrix has "name" in its column names
  check_names<-grep(cnames,pattern="^name$")
  
  names_with_mz_time<-NA
  
  X<-data_matrix
  if(length(check_names)>0){
    
    if(check_names==1){
      
      check_names1<-grep(cnames,pattern="^mz$")
      check_names2<-grep(cnames,pattern="^time$")
      
      
      if(length(check_names1)<1 & length(check_names2)<1){
        mz<-seq(1.00001,nrow(X)+1,1)
        time<-seq(1.01,nrow(X)+1,1.00)
        check_ind<-gregexpr(cnames,pattern="^name$")
        check_ind<-which(check_ind>0)
        X<-as.data.frame(X)
        
        
        Name<-as.character(X[,check_ind])
        
        X<-cbind(mz,time,X[,-check_ind])
        names_with_mz_time=cbind(Name,mz,time)
        
        names_with_mz_time<-as.data.frame(names_with_mz_time)
        X<-as.data.frame(X)
        
        
        write.table(names_with_mz_time,file="Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
        
      }else{
        
        if(length(check_names1)>0 & length(check_names2)>0){
          
          check_ind<-gregexpr(cnames,pattern="^name$")
          check_ind<-which(check_ind>0)
          Name<-as.character(X[,check_ind])
          X<-X[,-check_ind]
          names_with_mz_time=cbind(Name,X$mz,X$time)
          colnames(names_with_mz_time)<-c("Name","mz","time")
          names_with_mz_time<-as.data.frame(names_with_mz_time)
          X<-as.data.frame(X)
          write.table(names_with_mz_time,file="Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
        }
      }
      
    }
  }else{
    
    
    
    check_names1<-grep(cnames[1],pattern="^mz$")
    check_names2<-grep(cnames[2],pattern="^time$")
    if(length(check_names1)<1 || length(check_names2)<1){
      stop("Invalid feature table format. The format should be either Name in column A or mz and time in columns A and B. Please check example files.")
    }
  }
  
  data_matrix<-X
  #print(head(data_matrix))
  
  # Set a missing threshold if not specified
  if (is.na(all.missing.thresh)) {
    all.missing.thresh = -1
  }
  
  # Remove columns corresponding to 'samplermindex'
  if (!is.na(samplermindex)) {
    data_matrix <- data_matrix[, -samplermindex]
  }
  
  # Keep only unique records
  data_matrix <- unique(data_matrix)
  
  # Replace missing values with NAs
  if (!is.na(missing.val)) {
    data_matrix <- replace(as.matrix(data_matrix), which(data_matrix == missing.val), NA)
  }

  data_matrix_orig <- data_matrix

  # Get sample names
  snames <- colnames(data_matrix)

  # Create directory for output
  dir.create(parentoutput_dir, showWarnings = FALSE)

  fheader = "transformed_log2fc_threshold_"

  data_m <- as.matrix(data_matrix[, -c(1:2)])
  
  
  # Step 2: Average replicates
  if(summarize.replicates == TRUE) {
    if(num_replicates > 1) {
      # Apply the function getSumreplicates to data_matrix
      # This function is assumed to be from another R package
      # If it is not, please replace package_name with the appropriate package
      data_m <- getSumreplicates(data_matrix, alignment.tool = "apLCMS", numreplicates = num_replicates, numcluster = 10, rep.max.missing.thresh = rep.max.missing.thresh, summary.method = summary.method, summary.na.replacement, missing.val = missing.val)
      
      # Replace NA values with the missing_val
      data_m <- replace(data_m, which(is.na(data_m) == TRUE), missing.val)
      
      # Print messages and generate filenames based on the summary.method
      filename <- switch(
        summary.method,
        "mean" = { print("Replicate averaging done"); "Rawdata_averaged.txt"},
        "median" = { print("Replicate median summarization done"); "Rawdata_median_summarized.txt"},
        NA
      )
      
      # Write the resulting data matrix to a file
      write.table(cbind(data_matrix[, 1:2], data_m), file = filename, sep = "\t", row.names = FALSE)
      
      # Update data_matrix
      data_matrix <- cbind(data_matrix[, 1:2], data_m)
    }
  }
  
  
  # Continue with the original dataset if no averaging was done
  data_matrix <- cbind(data_matrix[, 1:2], data_m)
  data_matrix_orig <- data_matrix
  data_subjects <- data_m
  ordered_labels <- list()
  num_samps_group <- new("list")
  
  
  # If class labels are provided, process the labels
  if (!is.na(class_labels_file)) {
    data_matrix <- list()
  
    # Load class labels from the file or use the provided matrix
    classlabels <- if (is.na(Ymat)) {
      readr::read_table(class_labels_file, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      Ymat
    }
    
    cnames1<-colnames(classlabels)
    cnames1[1]<-c("SampleID")
    cnames1[2]<-c("Class")
    
    colnames(classlabels)<-cnames1 #c("SampleID","Class")
    
    f1<-table(classlabels$SampleID)
    
    
    
    classlabels<-as.data.frame(classlabels)
    classlabels<-classlabels[seq(1,dim(classlabels)[1],num_replicates),]
    #print(classlabels)
    class_labels_levels<-levels(as.factor(classlabels[,2]))
    bad_rows<-which(class_labels_levels=="")
    if(length(bad_rows)>0){
      class_labels_levels<-class_labels_levels[-bad_rows]
    }
    
    # Update the matrix based on the class labels
    for(c in 1:length(class_labels_levels)) {
      classlabels_index <- which(classlabels[, 2] == class_labels_levels[c])
      ordered_labels <- c(ordered_labels, as.character(classlabels[classlabels_index, 2]))
      num_samps_group[[c]] <- length(classlabels_index)
    }
  
    data_matrix <- cbind(data_matrix_orig[, 1:2], data_matrix)
    data_m <- as.matrix(data_matrix[, -c(1:2)])
    
    
    } else {
      if(is.na(Ymat)==TRUE)
      {
        classlabels<-rep("A",dim(data_m)[2])
        classlabels<-as.data.frame(classlabels)
        ordered_labels<-classlabels
        num_samps_group[[1]]<-dim(data_m)[2]
        class_labels_levels<-c("A")
        data_m<-as.matrix(data_matrix[,-c(1:2)])
      } else {
        classlabels<-Ymat
        classlabels<-as.data.frame(classlabels)
        
        if(pairedanalysis==TRUE){
          classlabels<-classlabels[,-c(2)]
        }
      
      cnames1<-colnames(classlabels)
      cnames1[1]<-c("SampleID")
      cnames1[2]<-c("Class")
      
      colnames(classlabels)<-cnames1
      #colnames(classlabels)<-c("SampleID","Class")
      f1<-table(classlabels$SampleID)
      
      
      classlabels<-as.data.frame(classlabels)
      classlabels<-classlabels[seq(1,dim(classlabels)[1],num_replicates),]
      #print(classlabels)
      
      
      if(alphabetical.order==FALSE){
        
        classlabels[,2]<-factor(classlabels[,2],levels=unique(classlabels[,2]))
        class_labels_levels<-levels(as.factor(classlabels[,2]))
        
      }else{
          class_labels_levels<-levels(as.factor(classlabels[,2]))
      }
      
      
      bad_rows<-which(class_labels_levels=="")
      if(length(bad_rows)>0){
        class_labels_levels<-class_labels_levels[-bad_rows]
      }
      
      for(c in 1:length(class_labels_levels))
      {
       
        classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
        ordered_labels<-c(ordered_labels,as.character(classlabels[classlabels_index,2]))
        num_samps_group[[c]]<-length(classlabels_index)
      }
    }
  }
  
    
  # Check if sample names match
  rnames_xmat <- colnames(data_matrix[, -c(1:2)])
  
  if (!is.na(classlabels)) {
    rnames_ymat <- as.character(classlabels[, 1])
  } else {
    rnames_ymat <- rnames_xmat
  }
  
  if(all(rnames_xmat==rnames_ymat)==FALSE){
    stop("Sample names do not match between feature table and classlabels. Please try replacing spaces and - with .")
  }

  
  #Step 3a) Remove features if signal is not detected in at least x% of all samples
  ##################################################################################
  metab_zeros={}
  data_clean<-{}
  clean_metabs<-{}
  total_good_metabs<-{}
  
  total_sigs<-apply(data_m,1,function(x){
    if(is.na(missing.val)==FALSE){return(length(which(x!=missing.val)))
    }else{
      return(length(which(is.na(x)==FALSE)))
    }})
  
  
  useful_metabs_1<-which(total_sigs>=dim(data_m)[2])
  
  useful_metabs_2<-which(total_sigs>=dim(data_m)[2]*0.95)
  
  
  if(is.na(all.missing.thresh)==FALSE)
  {
    
    total_sig_thresh<-dim(data_m)[2]*all.missing.thresh
    
    total_good_metabs<-which(total_sigs>total_sig_thresh)
    
  }
  
  #remove bad features based on all missing values criteria
  if(length(total_good_metabs)>0){
    data_m<-data_m[total_good_metabs,]
    data_matrix<-data_matrix[total_good_metabs,]
    #print(paste("Dimension of data matrix after overall ",all.missing.thresh,"% signal threshold filtering",sep=""))
    #print(paste("Dimension of data matrix after using overall ",100*all.missing.thresh, "% signal criteria for filtering:"),sep="")
    #print(dim(data_matrix))
  }else{
    stop(paste("None of the metabolites have signal in ",all.missing.thresh*100, "% of samples",sep=""))
  }
  
  
  #Step 3b) Find features for which the signal is not detected in at least x% of samples in either of the groups
  
  
  data_m<-data_matrix[,-c(1:2)]
  
  
  
  
  if(is.na(group.missing.thresh)==FALSE)
  {
    
    
    if(length(class_labels_levels)>1){
      
      
      
      
      clean_metabs<-lapply(1:dim(data_matrix)[1],function(metab_num)
      {
        clean_metabs<-NA
        for(c in 1:length(class_labels_levels)){
          
          classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
          templabels<-classlabels[,2]
          if(is.na(missing.val)==FALSE){
            num_cursig<-length(which(data_m[metab_num,which(templabels==class_labels_levels[c])]!=missing.val))
          }else{
            num_cursig<-length(which(is.na(data_m[metab_num,which(templabels==class_labels_levels[c])])==FALSE))
          }
          sig_thresh_cur<-length(which(templabels==class_labels_levels[c]))*group.missing.thresh
          if(num_cursig>=sig_thresh_cur)
          {
            clean_metabs<-metab_num
            break   #for(i in 1:4){if(i==3){break}else{print(i)}}
            
          }
          
        }
        return(clean_metabs)
      })
    }
    else{
      
      
      
      if(length(class_labels_levels)==1){
        num_samps_group[[1]]<-num_samps_group[[1]]
        
        
        sig_thresh_groupA<-group.missing.thresh*num_samps_group[[1]]
        
        
        clean_metabs<-lapply(1:dim(data_matrix)[1],function(metab_num)
        {
          if(is.na(missing.val)==FALSE){
            num_sigsA<-length(which(data_m[metab_num,1:num_samps_group[[1]]]!=missing.val))
            
          }else{
            
            num_sigsA<-length(which(is.na(data_m[metab_num,1:num_samps_group[[1]]])==FALSE))
          }
          
          if((num_sigsA>=sig_thresh_groupA) )
          {
            #clean_metabs<-c(clean_metabs,metab_num)
            
            return(metab_num)
          }else{
            
            return(NA)
          }
          
        })
      }
      
      
    }
    
    clean_metabs<-unlist(clean_metabs)
    clean_metabs<-na.omit(clean_metabs)
    
    
    
  }else{
    
    
    
    clean_metabs<-seq(1,dim(data_matrix)[1])
  }
  ####################################################################################
  
  #Step 4) Replace missing values
  if(summarize.replicates==TRUE)
  {
    
    {
      
      if(is.na(missing.val)==FALSE){
        
        #print("Replacing missing values with NAs.")
        data_m<-replace(as.matrix(data_m),which(data_m==missing.val),NA)
      }
      
      
      if(summary.na.replacement=="zeros"){
        data_m<-replace(data_m,which(is.na(data_m)==TRUE),0)
      }else{
        if(summary.na.replacement=="halfsamplemin"){
          data_m<-apply(data_m,2,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(min(x,na.rm=TRUE)/2)}; return(x)})
        }else{
          
          if(summary.na.replacement=="halfdatamin"){
            
            
            min_val<-min(data_m,na.rm=TRUE)*0.5
            data_m<-replace(data_m,which(is.na(data_m)==TRUE),min_val)
            
            #data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
          }else{
            if(summary.na.replacement=="halffeaturemin"){
              data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(min(x,na.rm=TRUE)/2)}; return(x)})
              data_m<-t(data_m)
            }else{
              
              
              if(summary.na.replacement=="bpca"){
                library(pcaMethods)
                
                
                pc1 <- pcaMethods::pca(t(data_m), method="bpca", nPcs=3,scale="uv")
                
                data_m<-pcaMethods::completeObs(pc1)
                
                try(detach("package:pcaMethods",unload=TRUE),silent=TRUE)
                
                data_m<-t(data_m)
                
              }else{
                
                if(summary.na.replacement=="knn"){
                  suppressMessages(library(impute))
                  
                  data_m<-apply(data_m,2,as.numeric)
                  data_m<-impute.knn(data.matrix(data_m),k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
                  data_m<-data_m$data
                  
                }else{
                  
                  if(summary.na.replacement=="featuremean"){
                    data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(mean(x,na.rm=TRUE))}; return(x)})
                    data_m<-t(data_m)
                  }else{
                    
                    if(summary.na.replacement=="randomforest"){
                      
                      #  #save(data_m,file="data_m.Rda")
                      final_set<-RFimpute(t(data_m))
                      
                      final_set<-t(final_set$ximp)
                      final_set<-as.data.frame(final_set)
                      
                      
                    }else{
                      
                      if(summary.na.replacement=="QRILC"){
                        
                        # final_set<-QRILCimpute(data_m)
                        library(tmvtnorm)
                        final_set<-QRILCimpute(data_m) #rows: features; cols: samples
                        ##save(final_set,file="final_set.Rda")
                        final_set<-ldply(final_set,rbind) 
                        final_set<-t(final_set)               
                        final_set<-as.data.frame(final_set)
                        
                        if(length(which(final_set<0))>0){
                          final_set<-replace(as.matrix(final_set),which(final_set<0),0)
                        }
                        
                      }else{
                        data_m<-data_m
                        
                      }
                      
                      
                    }
                    
                    
                  }
                  
                }
              }
              
              
              
            }
          }
        }
        
        
      }
    }
  }else
  {
    data_m<-data_matrix[,-c(1:2)]
    
    if(is.na(missing.val)==FALSE){
      
     # print("Replacing missing values with NAs.")
      data_m<-replace(as.matrix(data_m),which(data_m==missing.val),NA)
    }
    
    if(summary.na.replacement=="zeros"){
      data_m<-replace(data_m,which(is.na(data_m)==TRUE),0)
    }else{
      if(summary.na.replacement=="halfsamplemin"){
        data_m<-apply(data_m,2,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(min(x,na.rm=TRUE)/2)}; return(x)})
      }else{
        
        if(summary.na.replacement=="halfdatamin"){
          
          
          min_val<-min(data_m,na.rm=TRUE)*0.5
          data_m<-replace(data_m,which(is.na(data_m)==TRUE),min_val)
          
          #data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
        }else{
          if(summary.na.replacement=="halffeaturemin"){
            data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(min(x,na.rm=TRUE)/2)}; return(x)})
            data_m<-t(data_m)
          }else{
            
            if(summary.na.replacement=="bpca"){
              
              suppressMessages(library(pcaMethods))
              pc1 <- pcaMethods::pca(t(data_m), method="bpca", nPcs=3, scale="uv")
              
              data_m<-completeObs(pc1)
              data_m<-t(data_m)
            }else{
              if(summary.na.replacement=="knn"){
                suppressMessages(library(impute))
                data_m<-impute.knn(data.matrix(data_m),k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
                data_m<-data_m$data
              }else{
                
                
                if(summary.na.replacement=="randomforest"){
                  
                  final_set<-RFimpute(t(data_m))
                  
                  final_set<-t(final_set$ximp)
                  final_set<-as.data.frame(final_set)
                  
                  
                }else{
                  
                  if(summary.na.replacement=="QRILC"){
                    
                    # final_set<-QRILCimpute(data_m)
                    
                    library(tmvtnorm)
                    final_set<-QRILCimpute(data_m) #rows: features; cols: samples
                    ##save(final_set,file="final_set.Rda")
                    final_set<-ldply(final_set,rbind) 
                    final_set<-t(final_set)               
                    final_set<-as.data.frame(final_set)
                    
                    if(length(which(final_set<0))>0){
                      final_set<-replace(as.matrix(final_set),which(final_set<0),0)
                    }
                    
                  }
                  
                  
                }
                
                
                
              }
            }
            
          }
        }
      }
      
      
    }
    
  }
  
  
  
  
  #group-wise missing values
  if(length(clean_metabs)>0)
  {
    data_m<-data_m[clean_metabs,]
    data_matrix<-data_matrix[clean_metabs,]
    
   # print(paste("Dimension of data matrix after using group-wise (Factor 1) ",100*group.missing.thresh, "% signal criteria for filtering:"),sep="")
    #print(dim(data_matrix))
    
  }
  
  data_m<-as.data.frame(data_m)
  data_matrix<-as.data.frame(data_matrix)
  
  ##save(data_matrix,file="data_matrix.Rda")
  ##save(data_m,file="data_m.Rda")
  
  
  data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
  write.table(data_matrix,file="pretransformation.txt",sep="\t",row.names=FALSE)
  ####################################################################
  #Step 4) Data transformation and normalization
  
  data_m_prescaling<-data_m
  
  
  
  
  
  if(TIC_norm==TRUE){
    
    ####savedata_m,file="data_m_raw.Rda")
    
    
    data_m<-do_TIC_norm(data_m)
    
  }else{
    
    if(mstus==TRUE){
      
      
      data_m<-do_MSTUS_norm(data_m,missing.val=missing.val)
      
      
      
    }
    
  }
  
  if(cubicspline_norm==TRUE){
    
    data_m<-do_cubicspline_norm(data_m)
    
  }
  if(log2transform==TRUE)
  {
    data_m<-do_log2transform_norm(data_m,log2.transform.constant)
    
    
  }
  
  
  if(eigenms_norm==TRUE){
    
    feature_id_vector<-paste(data_matrix[,c(1)],"_",data_matrix[,c(2)],sep="")
    
    
    
    if(featselmethod=="limma2wayrepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="spls2wayrepeat"){
      analysistype="twowayrepeat"
      pairedanalysis=TRUE
    }else{
      
      if(featselmethod=="limma1wayrepeat" | featselmethod=="lm1wayanovarepeat" | featselmethod=="spls1wayrepeat"){
        analysistype="onewayrepeat"
        pairedanalysis=TRUE
      }else{
        pairedanalysis=FALSE
        
      }
      
    }
    
    data_m<-do_eigenms_norm(data_m,classlabels,featselmethod,feature_id_vector,pairedanalysis=pairedanalysis)
    
  }
  
  if(sva_norm==TRUE){
    
    data_m<-do_sva_norm(data_m,classlabels,featselmethod)
    
  }
  
  
  
  
  if(quantile_norm==TRUE)
  {
    data_m<-do_quantile_norm(data_m)
    
    
  }
  
  if(vsn_norm==TRUE)
  {
    data_m<-do_vsn_norm(data_m)
    
    
    
  }
  
  if(madscaling==TRUE)
  {
    
    colmedians=apply(data_m,2,function(x){median(x,na.rm=TRUE)})
    
    Y=sweep(data_m,2,colmedians)
    mad<-apply(abs(Y),2,function(x){median(x,na.rm=TRUE)})
    const<-1 #prod(mad)^(1/length(mad))
    scale.normalized<-sweep(Y,2,const/mad,"*")+mean(colmedians,na.rm=TRUE)
    
    
    data_m<-scale.normalized
    
  }
  
  if(lowess_norm==TRUE)
  {
    data_m<-do_loess_norm(data_m)
    #print("lowess")
  }
  
  
  
  if(medcenter==TRUE)
  {
    data_m<-do_medcenter_norm(data_m)
    
    
  }
  if(znormtransform==TRUE)
  {
    data_m<-do_znormtransform_norm(data_m)
  }
  
  
  
  if(paretoscaling==TRUE){
    #pareto scaling
    data_m<-do_paretoscaling_norm(data_m)
    
  }
  
  if(rangescaling==TRUE){
    #range scaling
    data_m<-do_rangescaling_norm(data_m)
  }
  
  
  
  
  
  data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
  
  
  data_m<-round(data_m,5)
  data_m<-as.data.frame(data_m)
  
  
  num_rows<-dim(data_m)[1]
  num_columns<-dim(data_m)[2]
  
  #print("num rows is ")
  #print(num_rows)
  #for apLCMS:
  rnames<-paste("mzid_",seq(1,num_rows),sep="")
  rownames(data_m)=rnames
  
  rnames_xmat<-colnames(data_m)
  rnames_ymat<-classlabels[,1]
  
  check_ylabel<-regexpr(rnames_ymat[1],pattern="^[0-9]*",perl=TRUE)
  check_xlabel<-regexpr(rnames_xmat[1],pattern="^X[0-9]*",perl=TRUE)
  
  
  if(length(check_ylabel)>0 && length(check_xlabel)>0){
    if(attr(check_ylabel,"match.length")>0 && attr(check_xlabel,"match.length")>0){
      
      rnames_ymat<-paste("X",rnames_ymat,sep="")
    }
  }
  classlabels<-classlabels[match(rnames_xmat,rnames_ymat),]
  
  
  filename<-paste("ordered_classlabels_file.txt",sep="")
  write.table(classlabels, file=filename,sep="\t",row.names=FALSE)
  
  filename<-paste("Normalized_sigthreshfilt_averaged_data.txt",sep="")
  data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
  write.table(data_matrix, file=filename,sep="\t",row.names=FALSE)
  data_matrix_prescaling<-cbind(data_matrix[,c(1:2)],data_m_prescaling)
  setwd("../")
  return(list(data_matrix_afternorm_scaling=data_matrix,data_matrix_prescaling=data_matrix_prescaling,classlabels=classlabels,names_with_mz_time=names_with_mz_time))
  #return(data_matrix)
}
