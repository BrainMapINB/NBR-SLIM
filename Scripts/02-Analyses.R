#!/usr/bin/env Rscript

# Clear workspace
rm(list = ls())

# Read preprocessed data
DF <- read.csv("https://raw.githubusercontent.com/BrainMapINB/NBR-SLIM/main/Data/phenotypic.csv")
# List of connectivity matrices and read them
cmx_ls <- paste0("https://raw.githubusercontent.com/BrainMapINB/NBR-SLIM/main/Data/Lobe/",
                 DF$Match_ID,"_s",DF$Session,".txt")
lobe3D <- sapply(cmx_ls, function(x) as.matrix(read.table(x)))
# Reshape
lobe3D <- array(lobe3D, c(8,8,nrow(DF)))

# Remove subjects with only one session
table(table(DF$Match_ID))
tp1_key <- match(names(which(table(DF$Match_ID)==1)),DF$Match_ID)
DF <- DF[-tp1_key,]
lobe3D <- lobe3D[,,-tp1_key]

# Sample descritive statistics
length(table(DF$Match_ID)) # Number of subjects
table(DF$Sex[match(names(table(DF$Match_ID)), DF$Match_ID)]) # Sex frequencies
range(DF$Age, na.rm = T) # Age range
table(table(DF$Match_ID)) # Session frequencies

########################################################################
# Analyses
########################################################################

# All inputs to run NBS are in the repository.
# However, inside this IF there are the commands to generate them
if(0){
  # First of all, remove NA's (NBS does not allow NA's)
  DFNBS <- DF[which(!is.na(DF$State_anxiety)),]
  lobeNBS <- lobe3D[,,which(!is.na(DF$State_anxiety))]
  DFNBS <- DFNBS[which(!is.na(lobeNBS[1,1,])),]
  lobeNBS <- lobeNBS[,,which(!is.na(lobeNBS[1,1,]))]
  
  # Balance sample (remove those with only one session)
  DFNBS$Match_ID <- factor(DFNBS$Match_ID)
  lobeNBS <- lobeNBS[,,-match(names(which(table(DFNBS$Match_ID)==1)),DFNBS$Match_ID)]
  DFNBS <- DFNBS[-match(names(which(table(DFNBS$Match_ID)==1)),DFNBS$Match_ID),]
  
  # Take two samples: with two and three time points
  # Balance sample for 2tp (remove third session in those with three sessions)
  DFNBS$Match_ID <- factor(DFNBS$Match_ID)
  lobeNBS2 <- lobeNBS[,,-(match(names(which(table(DFNBS$Match_ID)==3)),DFNBS$Match_ID)+2)]
  DFNBS2 <- DFNBS[-(match(names(which(table(DFNBS$Match_ID)==3)),DFNBS$Match_ID)+2),]
  # Balance sample for 3tp (remove third session in those with three sessions)
  lobeNBS3 <- lobeNBS[,,-c(match(names(which(table(DFNBS$Match_ID)==2)),DFNBS$Match_ID),
                           match(names(which(table(DFNBS$Match_ID)==2)),DFNBS$Match_ID)+1)]
  DFNBS3 <- DFNBS[-c(match(names(which(table(DFNBS$Match_ID)==2)),DFNBS$Match_ID),
                     match(names(which(table(DFNBS$Match_ID)==2)),DFNBS$Match_ID)+1),]
  
  # Design matrix for 2tp (including individual slopes)
  DFNBS2$Match_ID <- factor(DFNBS2$Match_ID)
  dmat <- matrix(0, nrow = nrow(DFNBS2), ncol = length(table(DFNBS2$Match_ID)))
  for(ii in 1:nrow(DFNBS2)) dmat[ii,as.integer(DFNBS2$Match_ID)[ii]] <- 1 # Add intercepts
  dmat <- cbind(rep(1,nrow(DFNBS2)),DFNBS2$State_anxiety, dmat)
  write.table(dmat, "NBS-2tp_dmat.txt", quote = F, row.names = F, col.names = F)
  
  # Contrast vector for NBS-2tp
  conmat <- rep(0, ncol(dmat))
  conmat[2] <- 1
  write.table(matrix(conmat, nrow = 1, ncol = ncol(dmat)),
              "NBS-2tp_con.txt",
              quote = F, row.names = F, col.names = F)
  
  # Design matrix for 3tp (including individual slopes)
  DFNBS3$Match_ID <- factor(DFNBS3$Match_ID)
  dmat <- matrix(0, nrow = nrow(DFNBS3), ncol = length(table(DFNBS3$Match_ID)))
  for(ii in 1:nrow(DFNBS3)) dmat[ii,as.integer(DFNBS3$Match_ID)[ii]] <- 1 # Add intercepts
  dmat <- cbind(rep(1,nrow(DFNBS3)),DFNBS3$State_anxiety, dmat)
  write.table(dmat, "NBS-3tp_dmat.txt", quote = F, row.names = F, col.names = F)
  
  # Contrast vector
  conmat <- rep(0, ncol(dmat))
  conmat[2] <- 1
  write.table(matrix(conmat, nrow = 1, ncol = ncol(dmat)),
              "../Analyses/NBS-3tp/NBS-3tp_con.txt",
              quote = F, row.names = F, col.names = F)
  
  # Connectivity matrices NBS-2tp
  for(ii in 1:dim(lobeNBS2)[3]){
    dest_file <- "NBS-2tp_mat"
    if(ii < 100) dest_file <- paste0(dest_file,"0") 
    if(ii < 10) dest_file <- paste0(dest_file,"0")
    dest_file <- paste0(dest_file,ii,".txt")
    write.table(lobeNBS2[,,ii], dest_file,
                quote = F, row.names = F, col.names = F)
  }
  # Block
  write.table(matrix(as.integer(DFNBS2$Match_ID), nrow = nrow(DFNBS2), ncol = 1),
              "NBS-2tp_block.txt", quote = F, row.names = F, col.names = F)
  
  # Connectivity matrices NBS-3tp
  for(ii in 1:dim(lobeNBS3)[3]){
    dest_file <- "NBS-3tp_mat"
    if(ii < 100) dest_file <- paste0(dest_file,"0") 
    if(ii < 10) dest_file <- paste0(dest_file,"0")
    dest_file <- paste0(dest_file,ii,".txt")
    write.table(lobeNBS3[,,ii], paste0("../Analyses/NBS-3tp/cmx/",dest_file),
                quote = F, row.names = F, col.names = F)
  }
  # Block
  write.table(matrix(as.integer(DFNBS3$Match_ID), nrow = nrow(DFNBS3), ncol = 1),
              "NBS-3tp_block.txt", quote = F, row.names = F, col.names = F)
}

# Test NBR-LME
set.seed(18900217)
library(NBR)
library(parallel)
DF$Session <- factor(DF$Session)
DF$Match_ID <- factor(DF$Match_ID)
ncores <- detectCores()
tic <- Sys.time()
nbr_result <- nbr_lme_aov(net = lobe3D, diag = F,
                      nnodes = 8, idata = DF,
                      thrP = NULL, thrF = 4,
                      mod = "~ State_anxiety",
                      rdm = "~ 1+State_anxiety|Match_ID",
                      nperm = 1000, cores = ncores,
                      nudist = T, na.action = na.exclude
)
toc <- Sys.time()
show(toc-tic)
# Save results in plain text format
write.table(nbr_result$components$State_anxiety,
            "NBR-LME_State_Anxiety_comp.txt",
            quote = F, row.names = F)
write.table(nbr_result$fwe$State_anxiety,
            "NBR-LME_State_Anxiety_fwe.txt",
            quote = F, row.names = F)
write.table(nbr_result$nudist,
            "NBR-LME_State_Anxiety_nudist.txt",
            quote = F, row.names = F)
