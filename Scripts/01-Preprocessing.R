#!/usr/bin/env Rscript

# Preprocessed data is already in the repository.
# However, here are the steps to create it.
# First of all, go to the public dataset repository:
# http://fcon_1000.projects.nitrc.org/indi/retro/southwestuni_qiu_index.html
# And get the followging files:
# "swu_slim_addphenodata.tar.gz"            
# "swu_slim_connmats_subs25629-30757.tar.gz"
# "swu_slim_connmats_subs30758-31057.tar.gz"
# "swu_slim_phenodata_time1.tsv"            
# "swu_slim_phenodata_time2.tsv"            
# "swu_slim_phenodata_time3.tsv"            
# "swu_slim_scan_params_anat.pdf"           
# "swu_slim_scan_params_dti.pdf"            
# "swu_slim_scan_params_rest.pdf"

# Then load the next function that only needs the input directory
# where those files are storaged. It will create the output preprocessed data.

########################################################################
# Preprocessing
########################################################################

pp_func <- function(in_dir){
  
  # Input directory 'in_dir' must be a character vector path
  if(!is.character(in_dir)) stop("The input directory must be a character string!")
  # Input directory 'in_dir' must exist.
  if(!dir.exists(in_dir)) stop("The input directory must exist!")
  # Input directory 'in_dir' must contain the following files
  in_files <- c("swu_slim_addphenodata.tar.gz",           
                "swu_slim_connmats_subs25629-30757.tar.gz",
                "swu_slim_connmats_subs30758-31057.tar.gz",
                "swu_slim_phenodata_time1.tsv",           
                "swu_slim_phenodata_time2.tsv",           
                "swu_slim_phenodata_time3.tsv")
  if(!all(in_files %in% dir(in_dir))) stop("The input directory must contain the SWU-SLIM files!")
  
  # Concatenate phenotypic info
  phen_ls <- list.files(path = in_dir, pattern = "swu_slim_phenodata_time", full.names = T)
  for(ii in 1:length(phen_ls)){
    if(!exists("DF")){
      DF <- read.table(phen_ls[ii], header = T)
      names(DF)[3] <- "Age"
    } else{
      phen <- read.table(phen_ls[ii], header = T)
      phen$Age <- NA
      # Concatenate datasets
      DF <- rbind(DF,phen[,match(names(DF),names(phen))])
    }
  }
  
  # Curate dataset
  # Sort data in order to compute longitudinal ages
  DF <- DF[order(DF$Match_ID,DF$Session),]
  # Set dates in standard format
  Sys.setlocale("LC_TIME","en_US.UTF-8") # careful with the local language in order to set the months
  DF$Scan_Time <- as.Date(DF$Scan_Time, "%d-%b-%Y")
  # Compute ages
  DF$Age <- as.numeric(DF$Age)
  for(ii in 1:nrow(DF)){
    # Find NA values in Age variable
    if(is.na(DF$Age[ii])){
      # If its previous session has age, compute it based on scan date difference
      if(!is.na(DF$Age[ii-1]) & DF$Match_ID[ii]==DF$Match_ID[ii-1]){
        # Compute differences
        diffscan <- as.numeric(difftime(DF$Scan_Time[ii],DF$Scan_Time[ii-1]))/365.25
        DF$Age[ii] <- DF$Age[ii-1] + diffscan
      }
    }
  }
  # Curate Sex variable
  DF$Sex[which(DF$Sex == "male\xa1\xa1")] <- "male"
  DF$Sex <- factor(DF$Sex) # Relevel
  # Curate logical variables
  DF$func[which(DF$func=="n/a")] <- NA
  DF$DTI[which(DF$DTI=="n/a")] <- NA
  
  # Add psychometric variables
  # Create temporary directory
  tmp_dir <- file.path(in_dir,paste0("tmp",round(runif(1,100,999))))
  if(!dir.exists(tmp_dir)) dir.create(path = tmp_dir)
  untar(tarfile = file.path(in_dir,"swu_slim_addphenodata.tar.gz"), exdir = tmp_dir)
  # Read average frame-wise displacement (FD) for each fMRI scan
  phen_ls <- list.files(tmp_dir, "qap_functional_temporal", recursive = T, full.names = T)
  DF$mean_fd <- NA
  for(ii in 1:length(phen_ls)){
    # Read dataset
    phen <- read.table(phen_ls[ii], header = T, sep = "\t")
    # Match sessions
    id_intersect <- intersect(DF$Match_ID[which(DF$Session==ii)],phen$Match_ID)
    df_key <- match(id_intersect,DF$Match_ID[which(DF$Session==ii)])
    phen_key <- match(id_intersect,phen$Match_ID)
    # Check if ID coincide
    if(!identical(DF$Match_ID[which(DF$Session==ii)][df_key],phen$Match_ID[phen_key])){
      stop("ID's do not match!!")
    }
    DF$mean_fd[which(DF$Session==ii)][df_key] <- phen$mean_fd[phen_key]
  }
  # Do the same for psychometric variables
  # BDI
  phen_ls <- list.files(tmp_dir, "BDI", recursive = T, full.names = T)
  DF$BDI_Score <- NA
  for(ii in 1:length(phen_ls)){
    # Read dataset
    phen <- read.table(phen_ls[ii], header = T, sep = "\t")
    # Match sessions
    id_intersect <- intersect(DF$Match_ID[which(DF$Session==ii)],phen$Match_ID)
    df_key <- match(id_intersect,DF$Match_ID[which(DF$Session==ii)])
    phen_key <- match(id_intersect,phen$Match_ID)
    # Check if ID coincide
    if(!identical(DF$Match_ID[which(DF$Session==ii)][df_key],phen$Match_ID[phen_key])){
      stop("ID's do not match!!")
    }
    DF$BDI_Score[which(DF$Session==ii)][df_key] <- as.character(phen$BDI_Score[phen_key])
  }
  # Curate variable
  DF$BDI_Score[which(DF$BDI_Score=="n/a")] <- NA
  DF$BDI_Score <- as.numeric(DF$BDI_Score)
  # STAI-State
  phen_ls <- list.files(tmp_dir, "state_anxiety", recursive = T, full.names = T)
  DF$State_anxiety <- NA
  for(ii in 1:length(phen_ls)){
    # Read dataset
    phen <- read.table(phen_ls[ii], header = T, sep = "\t")
    # Match sessions
    id_intersect <- intersect(DF$Match_ID[which(DF$Session==ii)],phen$Match_ID)
    df_key <- match(id_intersect,DF$Match_ID[which(DF$Session==ii)])
    phen_key <- match(id_intersect,phen$Match_ID)
    # Check if ID coincide
    if(!identical(DF$Match_ID[which(DF$Session==ii)][df_key],phen$Match_ID[phen_key])){
      stop("ID's do not match!!")
    }
    DF$State_anxiety[which(DF$Session==ii)][df_key] <- as.character(phen$State_anxiety[phen_key])
  }
  # Curate variable
  DF$State_anxiety[which(DF$State_anxiety=="n/a")] <- NA
  DF$State_anxiety <- as.numeric(DF$State_anxiety)
  # STAI-Trait
  phen_ls <- list.files(tmp_dir, "Trait_anxiety", recursive = T, full.names = T)
  DF$Trait_anxiety <- NA
  for(ii in 1:length(phen_ls)){
    # Read dataset
    phen <- read.table(phen_ls[ii], header = T, sep = "\t")
    # Match sessions
    id_intersect <- intersect(DF$Match_ID[which(DF$Session==ii)],phen[,1])
    df_key <- match(id_intersect,DF$Match_ID[which(DF$Session==ii)])
    phen_key <- match(id_intersect,phen[,1])
    # Check if ID coincide
    if(!identical(DF$Match_ID[which(DF$Session==ii)][df_key],phen[phen_key,1])){
      stop("ID's do not match!!")
    }
    DF$Trait_anxiety[which(DF$Session==ii)][df_key] <- as.character(phen$Trait_anxiety[phen_key])
  }
  # Curate variable
  DF$Trait_anxiety[which(DF$Trait_anxiety=="n/a")] <- NA
  DF$Trait_anxiety <- as.numeric(DF$Trait_anxiety)
  
  # Save results
  write.csv(DF, "phenotypic.csv", row.names = F)
  # Remove intermediate files
  unlink(tmp_dir, recursive = T)
  
  #########################################################################
  # Read Connectivity Matrices
  
  # Read compressed files
  gz_ls <- list.files(path = in_dir, pattern = "connmats", full.names = T)
  gz_n <- length(gz_ls)
  if(gz_n == 0) stop("No .tar.gz files found!!")
  # Uncompress files
  tmp_dir <- file.path(in_dir,paste0("tmp",round(runif(1,100,999))))
  if(!dir.exists(tmp_dir)) dir.create(path = tmp_dir)
  for(ii in gz_ls) untar(tarfile = ii, exdir = tmp_dir)
  # List connectivity matrices
  #cmx_ls <- list.files(path = tmp_dir, pattern = "ROI_FC.mat", full.names = T, recursive = T)
  gsr_ls <- list.files(path = tmp_dir, pattern = "ROI_FC_Globreg.mat", full.names = T, recursive = T)
  # List only Dosenbach160 connectivity matrices
  #cmx_ls <- cmx_ls[grep("Dosenbach_160", cmx_ls)]
  gsr_ls <- gsr_ls[grep("Dosenbach_160", gsr_ls)]
  # Create subject and session record
  str_len <- length(unlist(strsplit(gsr_ls[1],"/")))
  cmxDF <- as.data.frame(t(sapply(1:length(gsr_ls), function(x) unlist(strsplit(gsr_ls[x],"/"))[str_len-3:2])))
  names(cmxDF) <- c("id","session")
  # Concatenate matrices
  # Load 'R.matlab' package
  if(!require("R.matlab")) install.packages("R.matlab"); library(R.matlab)
  # Store matrices
  gsr3D <- array(0, dim = c(160, 160, length(gsr_ls)))
  cmxDF$mean_fd_d160gsr <- as.numeric(NA)
  cmxDF$no_nan_d160gsr <- 0
  cat("Progress:")
  for(ii in 1:length(gsr_ls)){
    # Print progress
    if(ii%%100 == 0) cat(paste0("...",ii))
    # Read every matrix
    cmx <- readMat(gsr_ls[ii])
    cmxDF$mean_fd_d160gsr[ii] <- mean(cmx$fd[-1])
    gsr3D[,,ii] <- cmx$z
    if(sum(is.nan(cmx$z))<1) cmxDF$no_nan_d160gsr[ii] <- 1
  }
  cat("\n")
  
  # Curate data
  # Remove NaN matrices
  gsr3D <- gsr3D[,,which(cmxDF$no_nan_d160==1)]
  cmxDF <- cmxDF[which(cmxDF$no_nan_d160==1),]
  cmxDF$session <- as.integer(cmxDF$session)
  
  # Match datasets
  gsr3D_d160 <- array(as.numeric(NA), dim = c(160, 160, nrow(DF)))
  phen_key <- paste0(DF$Match_ID,"_s",DF$Session)
  cmx_key <- paste0(cmxDF$id,"_s",cmxDF$session)
  gsr3D_d160[,,match(cmx_key,phen_key)] <- gsr3D
  
  # Save data
  #saveRDS(gsr3D_d160, file.path(out_dir, "cmx3D_d160gsr.rds"))
  # Update FD
  DF$cmx_fd <- as.numeric(NA)
  DF$cmx_fd[match(cmx_key,phen_key)] <- cmxDF$mean_fd_d160
  write.csv(DF, "phenotypic.csv", row.names = F)
  # Remove intermediate files
  unlink(tmp_dir, recursive = T)
  
  #########################################################################
  # Connectivity matrices grouped by lobe and network levels
  # Create averaging function
  cmx_avg <- function(vol3D, avg_factor){
    nnodes <- nlevels(avg_factor)
    mat3D <- sapply(1:dim(vol3D)[3], function(x){
      # Averaging rows
      aux <- aggregate(x = vol3D[,,x],
                       by = list(avg_factor),
                       FUN = function(y) mean(y[is.finite(y)]))
      # Averaging columns
      aux2 <- aggregate(x = t(aux[,-1]),
                        by = list(avg_factor),
                        FUN = "mean")
      return(as.matrix(aux2[,-1]))
    })
    # Reshape object and save it
    mat3D[is.nan(mat3D)] <- NA
    mat3D <- array(mat3D, c(nnodes, nnodes, dim(vol3D)[3]))
    return(mat3D)
  }
  
  # Load 'brainGraph' package
  if(!require("brainGraph")) install.packages("brainGraph"); suppressMessages(library(brainGraph))
  
  # Lobe it
  lobe3D_d160gsr <- cmx_avg(gsr3D_d160, dosenbach160$lobe)
  
  # Save matrices
  #saveRDS(lobe3D_d160gsr, "lobe3D_d160gsr.rds")
  out_dir <- "Lobe"
  if(!dir.exists(out_dir)) dir.create(out_dir)
  for(ii in 1:length(phen_key)){
    outname <- paste0(phen_key[ii],".txt")
    write.table(lobe3D_d160gsr[,,ii], file.path(out_dir,outname),
                quote = F, col.names = F, row.names = F)
  }
}
