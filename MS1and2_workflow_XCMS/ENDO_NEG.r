### Analysis of diatom Mono- and Co-culture MS data
# MS1 data pre-processing and statistical analysis with xcms
# Parameters for peak detection calculated with IPO
# Ms2 data relation to MS1 precursors and preparation for annotation by MAW

###---- library ----
# Load libraries
#library(parallel)               # Detect number of cpu cores
#library(foreach)                # For multicore parallel
#library(doMC)                   # For multicore parallel
library(RColorBrewer)           # For colors
library(MSnbase)                # MS features
library(xcms)                   # Swiss army knife for metabolomics
library(CAMERA)                 # Metabolite Profile Annotation
library(Spectra)                # Spectra package needed for XCMS3
library(vegan)
library(multcomp)               # For Tukey test
library(Hmisc)                  # For correlation test
library(gplots)                 # For fancy heatmaps
library(circlize)               # For sunburst plot
library(plotrix)                # For sunburst plot
library(caret)                  # Swiss-army knife for statistics
library(pROC)                   # Evaluation metrics
library(PRROC)                  # Evaluation metrics
library(multiROC)               # Evaluation metrics
#library(chemodiv)               # Chemodiversity (Petren 2022)
#library(rcdk)                   # CDK
#library(rinchi)                 # Converting SMILES to InchiKey
library(plotly)                 # For creating html plots
library(htmlwidgets)            # For creating html plots
#library(shiny)                  # HTML in R
#library(sunburstR)              # HTML-sunburst plots
#library(heatmaply)              # HTML heatmaps
library(stringr)              
library(randomForest)           # Random forest
library(MetaboAnalystR)         # Random forest
#library(iESTIMATE)
source("https://raw.githubusercontent.com/ipb-halle/iESTIMATE/main/R/_functions.r")

########## set directory and list files ########
start.time <- Sys.time()

# set data directory for MS1 data files
input_dir_MS1 <- paste(getwd(), "/MS1_pos_neg/", sep = "")
input_dir_MS1

#setwd(input_dir_MS1)# set data directory for MS1 data files
#input_dir_MS1 <- paste(getwd(), "/MS1/", sep = "")
# set data directory for MS2 data files
input_dir_MS2 <- paste(getwd(), "/MS2/", sep = "")
input_dir_MS2
# file lists
files_MS1 <- list.files(input_dir_MS1)
files_MS2 <- list.files(input_dir_MS2)
files_MS2

# create plot directory
if (dir.exists(paste(getwd(), "/ endo_neg_1ms2_plots/", sep = ""))){
  print("plots directory already exists")
  start_time <- Sys.time()
}  else{
  dir.create(" endo_neg_1ms2_plots")
  start_time <- Sys.time()
  print("plots folder has been created")
}

raw_data_MS1_ENDO_neg <- list.files(input_dir_MS1, pattern = "neg_ENDO")
raw_data_MS1_ENDO_neg

raw_data_MS1_ENDO_neg <- list.files(input_dir_MS1, pattern = "neg_ENDO")
raw_data_MS1_ENDO_neg_files <- paste(input_dir_MS1, raw_data_MS1_ENDO_neg, sep ="")
raw_data_MS1_ENDO_neg_files

raw_data_MS2_ENDO_neg <- list.files(input_dir_MS2, pattern = "ENDOneg")
raw_data_MS2_ENDO_neg_files <- paste(input_dir_MS2, raw_data_MS2_ENDO_neg, sep ="")
raw_data_MS2_ENDO_neg_files

all_files <- c(raw_data_MS1_ENDO_neg_files, raw_data_MS2_ENDO_neg_files)
all_files

all_files_names <- str_remove(all_files, ".mzML")
all_files_names <- str_remove_all(all_files_names, input_dir_MS1)
all_files_names[26:27]<- str_remove_all(all_files_names[26:27], input_dir_MS2)
all_files_names[26:(25 + length(raw_data_MS2_ENDO_neg))]<- str_remove_all(all_files_names[26:(25 + length(raw_data_MS2_ENDO_neg))], input_dir_MS2)
all_files_names

# create vector with sample classes according to culture information sheet
samp_groups <- c("CoCuPp", "CoCuSm", "CoCuSm", "CoCuPp", "CoCuPp", "CoCuSm",
                 rep(x = "Sm", times = 8), 
                 rep(x = "Pp", times = 8),
                 "CoCuSm", "CoCuPp", "MB", rep(x = "MS2", times = length(raw_data_MS2_ENDO_neg)))

CoCuPp1 <- "Co-culture sample from Prymnesium parvum"
CoCuSm1 <- "Co-culture sample from Skeletonema marinoi"
Sm1 <- "Mono-culture sample from Skeletonema marinoi"
Pp1 <- "Mono-culture sample from Prymnesium parvum"
MB1 <- "Media Blank"
ms2 <- "MS2"

# create vector with sample classes according to culture information sheet
samp_groups_description <- c(CoCuPp1, CoCuSm1, CoCuSm1, CoCuPp1, CoCuPp1, CoCuSm1,
                 rep(x = Sm1, times = 8), 
                 rep(x = Pp1, times = 8),
                 CoCuSm1, CoCuPp1, MB1, rep(x = ms2, times = length(raw_data_MS2_ENDO_neg)))

samp_groups_description

# create vector with colors 
CoCuPp1 <- ("#3E134F")
CoCuSm1 <- rep("#F36E35", 2)
CoCuPp2 <- rep("#3E134F", 2)
CoCuSm2 <- ("#F36E35")
Sm <- rep("#F8B83C", 8)
Pp <- rep("#C53270", 8)
CoCuSm3 <- ("#F36E35")
CoCuPp3 <- ("#3E134F")
MB <- rep("#040404", 1)
ms2<- rep("aquamarine", length(raw_data_MS2_ENDO_neg))
color <- c(CoCuPp1, CoCuSm1, CoCuPp2, CoCuSm2, Sm, Pp, CoCuSm3, CoCuPp3, MB, ms2)



#all_files



# create phenodata based on culture type
pheno_data_ENDO <- data.frame(sample_name = all_files_names, sample_group = samp_groups, samp_groups_description=samp_groups_description)
pheno_col_ENDO <- data.frame(color)


#pheno_data_ENDO

msd <- readMSData(files = all_files,
                  pdata = new("NAnnotatedDataFrame",pheno_data_ENDO),
                  mode = "onDisk",
                  centroided = TRUE)
#msd 

# inspect data 
table(msLevel(msd))
head(fData(msd)[, c("scanWindowLowerLimit", "scanWindowUpperLimit",
                    "originalPeaksCount", "msLevel", 
                    "polarity", "retentionTime")])

# Restrict data to 650
msd <- filterRt(msd, c(0, 650))

# ONLY FOR MS1 DATA
# subset data for msLevel = 1 and save raw data
#msd <- filterMsLevel(msd, msLevel = 1)
#table(msLevel(msd))

# create plot directory
if (dir.exists(paste(getwd(), "/ endo_neg_1ms2_Results/", sep = ""))){
  print("plots directory already exists")
  start_time <- Sys.time()
}  else{
  dir.create(" endo_neg_1ms2_Results")
  start_time <- Sys.time()
  print("results folder has been created")
}
write.csv(fData(msd), file=paste(filename = " endo_neg_1ms2_Results/ENDO_neg_raw_data.csv", sep = ""), row.names=FALSE)

# Get base peak chromatograms
#register(bpstart(SnowParam()))
# setwd(input_dir_MS1_polarity)

chromas_ENDO_neg <- chromatogram(msd, 
                             aggregationFun="max", 
                             #msLevel = 1,
                             )
#chromas_ENDO_neg

# Plot chromatograms based on phenodata groups
#pdf(file="plots/ENDO_chromas.pdf", encoding="ISOLatin1", pointsize=2, width=6, height=4, family="Helvetica")
jpeg(filename = "endo_neg_1ms2_plots/endo_neg_chromas.jpeg", width = 2000, height = 1200, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=2, cex.lab=2, cex.main=2)
plot(chromas_ENDO_neg, main="Raw chromatograms", xlab="retention time [s]", ylab="intensity", col = color)
legend("topleft", bty="n", pt.cex=3, cex=1.5, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(msd@phenoData@data[["sample_group"]]))
dev.off()


# Get TICs
#pdf(file="plots/endo_tics.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "endo_neg_1ms2_plots/endo_neg_tics.jpeg", width = 2000, height = 1200, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=2, cex.lab=2, cex.main=2)
tics_endo_neg <- split(tic(msd), f=fromFile(msd))
boxplot(tics_endo_neg, col=color, ylab="intensity", xlab="sample", main="Total ion current", outline = FALSE)
legend("topleft", bty="n", pt.cex=3, cex=1.5, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(msd@phenoData@data[["sample_group"]]))
dev.off()

### ---- pre-processing ----

# grouping/binning for peak detection, based on similarity of the base peak chromatogram 
chromas_bin_ENDO_neg <- MSnbase::bin(chromas_ENDO_neg, binSize=2)
chromas_bin_ENDO_neg

chromas_bin_cor_ENDO_neg <- cor(log2(do.call(cbind, lapply(chromas_bin_ENDO_neg, intensity)))) # transformation with log
chromas_bin_cor_ENDO_neg

colnames(chromas_bin_cor_ENDO_neg) <- rownames(chromas_bin_cor_ENDO_neg) <- msd$sample_name

chromas_bin_cor_ENDO_neg[is.na(chromas_bin_cor_ENDO_neg)] <- 0


jpeg(filename = " endo_neg_1ms2_plots/heatmap_chromas_bin_ENDO_neg.jpeg", width = 1000, height = 1000, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(7,0,0,7), cex.axis=0.9, cex=0.6)
heatmap(chromas_bin_cor_ENDO_neg)
dev.off()


# Assess retention times and intensities of first file
head(rtime(chromas_ENDO_neg[1, 1]))
head(intensity(chromas_ENDO_neg[1, 1]))

# check for polarity
#head(fData(msd)[, c("polarity", "filterString", "msLevel", "retentionTime")])
#table(polarity(msd))

ms_params_ENDO_neg <- CentWaveParam(ppm=35, mzCenterFun="wMean", peakwidth=c(14, 59), 
                                     prefilter=c(3, 140), mzdiff=0.0155, snthresh=7, noise=0, 
                                     integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
                                     fitgauss=FALSE, roiList=list(), roiScales=numeric())
ms_params_ENDO_neg

# CHOSE THE PARAMETERS 
ms_data_ENDO_neg <- findChromPeaks(msd, param=ms_params_ENDO_neg)
ms_data_ENDO_neg

# check the detected peaks
head(chromPeaks(ms_data_ENDO_neg))
chromPeakData(ms_data_ENDO_neg)

ms_summary_ENDO_neg <- lapply(split.data.frame(chromPeaks(ms_data_ENDO_neg), 
                                            f=chromPeaks(ms_data_ENDO_neg)[, "sample"]), 
                           FUN=function(z) { c(peak_count=nrow(z), rt=quantile(z[, "rtmax"] - z[, "rtmin"])) } )

ms_summary_ENDO_neg <- do.call(rbind, ms_summary_ENDO_neg)

rownames(ms_summary_ENDO_neg) <- basename(fileNames(ms_data_ENDO_neg))
rownames(ms_summary_ENDO_neg)

print(ms_summary_ENDO_neg)

table(msLevel(ms_data_ENDO_neg))



write.csv(as.data.frame(table(msLevel(ms_data_ENDO_neg))), file=" endo_neg_1ms2_Results/ENDO_neg_ms_data.csv", row.names=FALSE)

jpeg(filename = " endo_neg_1ms2_plots/ENDO_neg_ms_data.jpeg", width = 1000, height = 600, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(4,18,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plotChromPeakImage(ms_data_ENDO_neg, main="Frequency of identified peaks per RT", binSize = 20)
#dev.off()

plotChromPeakImage(ms_data_ENDO_neg, main="Frequency of identified peaks per RT", binSize = 20)

## Group peaks
ms_data_ENDO_neg <- groupChromPeaks(ms_data_ENDO_neg, param=PeakDensityParam(
  sampleGroups=ms_data_ENDO_neg$sample_group, minFraction=0.7, bw=2.5))
ms_data_ENDO_neg

## RT correction
ms_data_ENDO_neg <- adjustRtime(ms_data_ENDO_neg, param=PeakGroupsParam(
  minFraction=0.7,smooth="loess",span=0.8,family="gaussian"))

# Plot the difference of raw and adjusted retention times
#pdf(file="plots/ENDO_ms1_raw_adjusted.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
jpeg(filename = " endo_neg_1ms2_plots/ENDO_neg_ms_raw_adjusted.jpeg", width = 500, height = 1000, quality = 100, bg = "white")
par(mfrow=c(2,1), mar=c(4.5,4.2,4,1), cex=0.8)
plot(chromas_ENDO_neg, peakType="none", main="Raw chromatograms")
plotAdjustedRtime(ms_data_ENDO_neg, lwd=2, main="Retention Time correction")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
dev.off()

## Group peaks
ms_data_ENDO_neg <- groupChromPeaks(ms_data_ENDO_neg, param=PeakDensityParam(
  sampleGroups=ms_data_ENDO_neg$sample_group, minFraction=0.7, bw=2.5))

# Get integrated peak intensity per feature/sample
print(head(featureValues(ms_data_ENDO_neg, value="into")))

ppm <- 35  

# missing value imputation, see xcmsSet
#ms_data_ENDO_neg <- fillChromPeaks(ms_data_ENDO_neg, param=FillChromPeaksParam(ppm=ppm, fixedRt=0, expandRt=5))
#ms_data_ENDO_neg

#head(featureValues(ms_data_ENDO_neg))
#head(featureSummary(ms_data_ENDO_neg, group=ms_data_ENDO_neg$sample_group))



# Evaluate grouping
#pdf(file="plots/ENDO_ms1_grouping.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = " endo_neg_1ms2_plots/ENDO_neg_ms_grouping.jpeg", width = 1000, height = 700, quality = 150, bg = "white")
ms_pca_ENDO_neg <- prcomp(t(na.omit(log2(featureValues(ms_data_ENDO_neg, value="into")))), center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms_pca_ENDO_neg$x[, 1], ms_pca_ENDO_neg$x[,2], pch=19, main="PCA: Grouping of samples",
     xlab=paste0("PC1: ", format(summary(ms_pca_ENDO_neg)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms_pca_ENDO_neg)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color, cex=2)
grid()
text(ms_pca_ENDO_neg$x[, 1], ms_pca_ENDO_neg$x[,2], labels=str_sub(ms_data_ENDO_neg$sample_name, - 3, - 1), col=color, pos=3, cex=0.9)
legend("topleft", bty="n", pt.cex=2, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(ms_data_ENDO_neg@phenoData@data[["sample_group"]]))
dev.off()

# broken stick
png(" endo_neg_1ms2_plots/BrokenStick_ENDO_neg_ms_grouping.png", width=10, height=6, units="in", res=100)
evplot = function(ev) {  
  # Broken stick model (MacArthur 1957)  
  n = length(ev)  
  bsm = data.frame(j=seq(1:n), p=0)  
  bsm$p[1] = 1/n  
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))  
  bsm$p = 100*bsm$p/n  
  # Plot eigenvalues and % of variation for each axis  
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))  
  barplot(ev, main="Eigenvalues ENDO MS grouping", col="blue", las=2)  
  abline(h=mean(ev), col="red")  
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")  
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE,   
          main="% variation", col=c("blue",2), las=2)  
  legend("topright", c("% eigenvalue", "Broken stick model"),   
         pch=15, col=c("blue",2), bty="n")  
  par(op)  
} 

ev_pc = ms_pca_ENDO_neg$sdev^2  
evplot(ev_pc)  
dev.off()

# Show peaks
tail(chromPeaks(ms_data_ENDO_neg))
tail(chromPeakData(ms_data_ENDO_neg))

# Show process history
processHistory(ms_data_ENDO_neg)

# save as R object 
MS_ENDO_neg_peak_detection <- ms_data_ENDO_neg




save(MS_ENDO_neg_peak_detection, file = " endo_neg_1ms2_Results/MS_ENDO_neg_peak_detection.RData")

# ---------- Build MS1 feature tables ----------
# Build feature matrix
ms_matrix_ENDO_neg <- featureValues(ms_data_ENDO_neg, method="medret", value="into")

colnames(ms_matrix_ENDO_neg) <- all_files_names
dim(ms_matrix_ENDO_neg)
# transpose feature table
feat_list_ENDO_neg <- t(ms_matrix_ENDO_neg)
feat_list_ENDO_neg

# Build feature summary
ms_summary_ENDO_neg <- featureSummary(ms_data_ENDO_neg)
ms_def_ENDO_neg <- featureDefinitions(ms_data_ENDO_neg)


# Missing value imputation by filling na positions with median of surrounding features
#feat_list_ENDO_neg[is.na(feat_list_ENDO_neg)] <- median(na.omit(as.numeric(unlist(feat_list_ENDO_neg))))

### Transform data
feat_list_ENDO_neg <- log2(feat_list_ENDO_neg)

# change 0 to small value to distinguish between values of 1 and NA
feat_list_ENDO_neg[which(feat_list_ENDO_neg == 0)] <- 0.01

# Missing value imputation
feat_list_ENDO_neg[which(is.na(feat_list_ENDO_neg))] <- 0

# save as csv
write.csv(feat_list_ENDO_neg, file=paste(filename = " endo_neg_1ms2_Results/feature_list_ENDO_neg.csv", sep = ""))

# Plot histogram
#pdf(file="plots/ENDO_feat_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = " endo_neg_1ms2_plots/ENDO_neg_feat_list_hist.jpeg", width = 500, height = 500, quality = 150, bg = "white")
hist(as.numeric(feat_list_ENDO_neg), main="Histogram of feature table")
dev.off()

# PCA of feature table results
#pdf(file="plots/ENDO_ms1_feature_table_pca.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = " endo_neg_1ms2_plots/ENDO_neg_ms_feature_table_pca.jpeg", width = 1000, height = 500, quality = 150, bg = "white")
ms_pca_ENDO_neg <- prcomp(feat_list_ENDO_neg, center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms_pca_ENDO_neg$x[, 1], ms_pca_ENDO_neg$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms_pca_ENDO_neg)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms_pca_ENDO_neg)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color, cex=2)
grid()
text(ms_pca_ENDO_neg$x[, 1], ms_pca_ENDO_neg$x[,2], labels=str_sub(ms_data_ENDO_neg$sample_name, - 3, - 1), col=color, pos=3, cex=0.9)
legend("topleft", bty="n", pt.cex=2, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(ms_data_ENDO_neg@phenoData@data[["sample_group"]]))
dev.off()

# broken stick
jpeg(" endo_neg_1ms2_plots/BrokenStick_ENDO_neg_ms_feature_table_pca.jpeg", width=10, height=6, units="in", res=100)
evplot = function(ev) {  
  # Broken stick model (MacArthur 1957)  
  n = length(ev)  
  bsm = data.frame(j=seq(1:n), p=0)  
  bsm$p[1] = 1/n  
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))  
  bsm$p = 100*bsm$p/n  
  # Plot eigenvalues and % of variation for each axis  
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))  
  barplot(ev, main="Eigenvalues ENDO MS feature table", col="blue", las=2)  
  abline(h=mean(ev), col="red")  
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")  
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE,   
          main="% variation", col=c("blue",2), las=2)  
  legend("topright", c("% eigenvalue", "Broken stick model"),   
         pch=15, col=c("blue",2), bty="n")  
  par(op)  
} 

ev_pc = ms_pca_ENDO_neg$sdev^2  
evplot(ev_pc)  
dev.off()


ms_intensity_cutoff <- 14

# Create single 0/1 matrix
bina_list_ENDO_neg <- t(ms_matrix_ENDO_neg)
bina_list_ENDO_neg[is.na(bina_list_ENDO_neg)] <- 1
bina_list_ENDO_neg <- log2(bina_list_ENDO_neg)
bina_list_ENDO_neg[bina_list_ENDO_neg < ms_intensity_cutoff] <- 0
bina_list_ENDO_neg[bina_list_ENDO_neg != 0] <- 1


# save as csv
write.csv(bina_list_ENDO_neg, file=paste(filename = " endo_neg_1ms2_Results/bina_list_ENDO_neg.csv", sep = ""))

# Only unique compounds in group mzml_pheno$ and not the others
uniq_list_ENDO_neg <- apply(X=bina_list_ENDO_neg, MARGIN=2, FUN=function(x) { if (length(unique(pheno_data_ENDO$sample_group[grepl("1", x)])) == 1) x else rep(0, length(x)) } )
colnames(uniq_list_ENDO_neg) <- colnames(bina_list_ENDO_neg)
rownames(uniq_list_ENDO_neg) <- rownames(bina_list_ENDO_neg)
uniq_list_ENDO_neg

# Create data frame
model_div_ENDO_neg             <- data.frame(features=apply(X=bina_list_ENDO_neg, MARGIN=1, FUN=function(x) { sum(x) } ))
model_div_ENDO_neg$richness    <- apply(X=bina_list_ENDO_neg, MARGIN=1, FUN=function(x) { sum(x) } )
#model_div_ENDO_neg$menhinick   <- apply(X=bina_list_ENDO_neg, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
model_div_ENDO_neg$shannon     <- apply(X=feat_list_ENDO_neg, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
model_div_ENDO_neg$pielou      <- apply(X=scale(feat_list_ENDO_neg, center=F), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#model_div_ENDO_neg$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list_ENDO_neg, species), index="chao")
model_div_ENDO_neg$simpson     <- apply(X=feat_list_ENDO_neg, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
model_div_ENDO_neg$inverse     <- apply(X=feat_list_ENDO_neg, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
model_div_ENDO_neg$fisher      <- apply(X=feat_list_ENDO_neg, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
model_div_ENDO_neg$unique      <- apply(X=uniq_list_ENDO_neg, MARGIN=1, FUN=function(x) { sum(x) })

# Remove NAs if present
model_div_ENDO_neg[is.na(model_div_ENDO_neg)] <- 0

# save as csv
write.csv(model_div_ENDO_neg, file=paste(filename = " endo_neg_1ms2_Results/model_div_ENDO_neg.csv", sep = ""))


# save the objects and tables
save(ms_def_ENDO_neg, file = "endo_neg_1ms2_Results/ms_def_ENDO_neg.RData")
save.image(file = "ENDO_neg_1ms2_Results/ENDO_neg_1MS2_environment.RData")



end.time <- Sys.time()

time.taken <- end.time - start.time
print(time.taken)




