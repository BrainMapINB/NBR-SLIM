#!/usr/bin/env Rscript

# Clear workspace
rm(list = ls())

# Input parameters
#op <- par() # Save original parameters

########################################################################
# Load data
########################################################################

# Read statistical networks
# Load 'R.matlab' package
if(!require("R.matlab")) install.packages("R.matlab"); library(R.matlab)

# Read input stats
# NBS-2tp
aux <- readMat("https://github.com/BrainMapINB/NBR-SLIM/blob/main/Analyses/NBS-2tp/NBS2tp_lobe_stai_F4.mat?raw=true")
NBS2tp_stat <- aux[[1]][[3]][[4]]
NBS2tp_stat[NBS2tp_stat<4] <- 0
NBS2tp_perm <- aux[[1]][[2]][[5]]
# NBS-3tp
aux <- readMat("https://github.com/BrainMapINB/NBR-SLIM/blob/main/Analyses/NBS-3tp/NBS3tp_lobe_stai_F4.mat?raw=true")
NBS3tp_stat <- aux[[1]][[3]][[4]]
NBS3tp_stat[NBS3tp_stat<4] <- 0
NBS3tp_perm <- aux[[1]][[2]][[5]]
# NBR
aux <- read.table("https://raw.githubusercontent.com/BrainMapINB/NBR-SLIM/main/Analyses/NBR-LME/NBR-LME_State_Anxiety_comp.txt", header = T)
NBR_stat <- matrix(0, 8, 8)
for(ii in 1:nrow(aux)) NBR_stat[aux[ii,2],aux[ii,3]] <- aux[ii,5]+4
NBR_stat <- NBR_stat + t(NBR_stat)
aux <- read.table("https://raw.githubusercontent.com/BrainMapINB/NBR-SLIM/main/Analyses/NBR-LME/NBR-LME_State_Anxiety_fwe.txt", header = T)
NBR_perm <- aux[1,4]
aux <- read.table("https://raw.githubusercontent.com/BrainMapINB/NBR-SLIM/main/Analyses/NBR-LME/NBR-LME_State_Anxiety_nudist.txt", header = T)
NBR_perm <- c(NBR_perm, aux[,2])

# Find components in NBS permuted matrices
source("https://raw.githubusercontent.com/BrainMapINB/NBR-SLIM/main/Scripts/Auxiliary/find_components.R")
# Homogenizate data
all_perm <- matrix(0, nrow = 1001, ncol = 3)
tri_pos <- which(upper.tri(matrix(nrow = 8, ncol = 8)), arr.ind = T)
# NBR-2tp
all_perm[,1] <- sapply(1:1001, function(x){ y <- find_components(NBS2tp_perm[x,], 4, tri_pos, 8)
if(is.null(y)) return(0) else return(max(aggregate(y[,5] ~ y[,4], y, sum)[, 2]))})
# NBR-3tp
all_perm[,2] <- sapply(1:1001, function(x){ y <- find_components(NBS3tp_perm[x,], 4, tri_pos, 8)
if(is.null(y)) return(0) else return(max(aggregate(y[,5] ~ y[,4], y, sum)[, 2]))})
# NBR
all_perm[,3] <- NBR_perm
colnames(all_perm) <- c("NBS2tp","NBS3tp","NBR")
all_perm <- as.data.frame(all_perm)

########################################################################
# Visualization
########################################################################

# Load 'corrplot', 'circlize' & 'scales' package
if(!require("corrplot")) install.packages("corrplot"); library(corrplot)
if(!require("circlize")) install.packages("circlize"); library(circlize)
if(!require("scales")) install.packages("scales"); library(scales)

# Corrplots input parameters
mod_labs <- c("NBS2tp","NBS3tp","NBR")
lob_labs <- c("FRT","PAR","TEMP","OCC","INS","CING","SUB","CBL")
colorcito <- colorRampPalette(c("white","yellow","orange","red"))

# Chord diagram input parameters
lob_fac <- as.factor(1:length(lob_labs))
levels(lob_fac) <- lob_labs

# Generate panel
par(mfrow=c(3,3))
for(nn in 1:3){
  
  # Corrplot
  mat_stat <- get(paste0(mod_labs[nn],"_stat"))
  colnames(mat_stat) <- row.names(mat_stat) <- lob_labs
  corrplot(mat_stat, is.corr = F,type = "lower", tl.col = "black",
           method = "square", cl.cex = 1, cl.lim = c(0,12.5),
           diag = F, mar=c(3,1,1,1)+0.1, col = colorcito(100))
  
  # Chord diagram of significant (uncorrected) intra and inter-module effects
  # Initialize the plot
  circos.clear()
  circos.initialize(factors = lob_fac, xlim=c(-1,1))
  # Build the regions of track #1
  circos.trackPlotRegion(factors = lob_fac, ylim = c(-1,1),
                         bg.col = "aliceblue",
                         bg.border = "black")
  # Add labels
  lob_pos <- ux(29.75, "mm")*(1:nlevels(lob_fac))
  circos.text(lob_pos,rep(0,nlevels(lob_fac)),
              labels = lob_labs,
              facing = "bending.inside", cex=1.4)
  # Add a links between a point and another
  # Find significant (NBS corrected) weigths
  up_tri <- which(upper.tri(mat_stat), arr.ind = T)
  sig_tri <- which(mat_stat[upper.tri(mat_stat)]>4)
  # Draw links
  for(ii in 1:length(sig_tri)){
    # Set central position
    rnum <- runif(1); p1 <- c(rnum,rnum-1)
    rnum <- runif(1); p2 <- c(rnum,rnum-1)
    # Set color
    fval <- mat_stat[up_tri[sig_tri[ii],1],up_tri[sig_tri[ii],2]]
    fper <- round(100*fval/12.5)
    fcol <- colorcito(100)[fper]
    # Draw link
    circos.link(lob_fac[up_tri[sig_tri[ii],1]], p1,
                lob_fac[up_tri[sig_tri[ii],2]], p2,
                col = scales::alpha(fcol,.9),
                border = scales::alpha(fcol,.4),
                h.ratio=0.6)
  }
  # Clear cicle parameters
  circos.clear()
  
  # Plot null distributions
  d <- density(all_perm[-1,nn])
  plot(d, xlim = c(-0.5,20), ylim = c(0,1.8),
       xlab = "Strength", axes = F, main = mod_labs[nn], col = "blue")
  polygon(d, col = "cyan", border = "blue")
  axis(1); axis(2, las=1)
  abline(v = all_perm[1,nn], col = "brown", lty=2, lwd = 2)
  print(sum(all_perm[-1,nn]>=all_perm[1,nn])/1000)
}
