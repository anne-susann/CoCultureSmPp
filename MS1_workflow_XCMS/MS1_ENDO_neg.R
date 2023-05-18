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

###------parallelization----

# library(parallel)
# n.cores <- detectCores()

# # Set up parallel processing using 2 cores
# if (.Platform$OS.type == "unix") {
#   register(bpstart(MulticoreParam(2)))
# } else {
#   register(bpstart(SnowParam(4)))
# }
#R.Version()

########## set directory and list files ########


# set data directory for MS1 data files
input_dir_MS1 <- paste(getwd(), "/MS1/", sep = "")
# set data directory for MS2 data files
input_dir_MS2 <- paste(getwd(), "/MS2/", sep = "")
# file lists
files_MS1 <- list.files(input_dir_MS1)
files_MS2 <- list.files(input_dir_MS2)


###----MS1 ENDO files----
#MS1_ENDO_files <- data.frame(list.files(input_dir_MS1, pattern = "ENDO"))
MS1_ENDO_files <- list.files(input_dir_MS1, pattern = "ENDO")
# basenames of files without path and without extension
#MS1_ENDO_names <- data.frame(str_remove(MS1_ENDO_files[,], ".mzML"))
MS1_ENDO_names <- str_remove(MS1_ENDO_files, ".mzML")


###----MS1 EXO files----
#MS1_EXO_files <- data.frame(list.files(input_dir_MS1, pattern = "EXO"))
#MS1_EXO_files <- list.files(input_dir_MS1, pattern = "EXO")
#MS1_EXO_files_full <- data.frame(list.files(input_dir_MS1, pattern = "EXO", full.names = TRUE))
#MS1_EXO_names <- str_remove(MS1_EXO_files, ".mzML")

###----MS2 files----
#MS2_files <- data.frame(list.files(input_dir_MS2, pattern = ".mzML"))
MS2_files <- list.files(input_dir_MS2, pattern = ".mzML")
# ENDO files
# MS2_ENDO_files <- data.frame(subset(MS2_files, grepl("ENDO", MS2_files[,1]), drop = TRUE))
# MS2_ENDO_neg_files <- data.frame(subset(MS2_ENDO_files, grepl("neg", MS2_ENDO_files[,1]), drop = TRUE))
# MS2_ENDO_neg_files <- data.frame(subset(MS2_ENDO_files, grepl("pos", MS2_ENDO_files[,1]), drop = TRUE))
MS2_ENDO_files <- subset(MS2_files, grepl("ENDO", MS2_files))
MS2_ENDO_neg_files <- subset(MS2_ENDO_files, grepl("neg", MS2_ENDO_files))
#MS2_ENDO_pos_files <- subset(MS2_ENDO_files, grepl("pos", MS2_ENDO_files))

# EXO files
#MS2_EXO_files <- subset(MS2_files, grepl("EXO", MS2_files))
#MS2_EXO_neg_files <- subset(MS2_EXO_files, grepl("neg", MS2_EXO_files))
#MS2_EXO_pos_files <- subset(MS2_EXO_files, grepl("pos", MS2_EXO_files))

###----Separate neg and pos mode in individual files-----
# create result directory
#if (dir.exists(paste(getwd(), "/raw_data/", sep = ""))){
#  raw_data_dir <- paste(getwd(), "/raw_data/", sep = "")
#}  else{
#  dir.create("raw_data")
#  raw_data_dir <- paste(getwd(), "/raw_data/", sep = "")
#}


## separation of ENDO FILES
#for (i in 1:nrow(MS1_ENDO_files)){
#  sps <- Spectra(paste(input_dir_MS1, MS1_ENDO_files[i,1], sep = ""), source = MsBackendMzR())
  
  # set polarity to negative = 0 or positive = 1
#  sps_neg <- sps[sps$polarity==0]
#  sps_pos <- sps[sps$polarity==1] 
  
  # export data in separate files for neg and pos
#  export(sps_neg, backend = MsBackendMzR(), file = paste(raw_data_dir, 
#                                                         gsub("KSS_210324","\\coculture_neg", MS1_ENDO_files[i,1]),
#                                                         sep = ""))
#  export(sps_pos, backend = MsBackendMzR(), file = paste(raw_data_dir, 
#                                                         gsub("KSS_210324","\\coculture_pos", MS1_ENDO_files[i,1]), 
#                                                         sep = ""))
#}

# raw_data_dir <- paste(getwd(), "/MS1_pos_neg/", sep = "")
# setwd(raw_data_dir)
# raw_data_MS1_ENDO <- data.frame(list.files(raw_data_dir, pattern = "ENDO"))
# raw_data_MS1_ENDO_neg <- data.frame(subset(raw_data_MS1_ENDO, grepl("neg", raw_data_MS1_ENDO[,1]), drop = TRUE))
# fls <- raw_data_MS1_ENDO_neg[,1]
# class(fls)


# sps_test <- Spectra(fls, source = MsBackendMzR())
# sps_test
# table(polarity(sps_test))

## separation of EXO FILES
#for (i in 1:nrow(MS1_EXO_files)){
#  sps <- Spectra(paste(input_dir_MS1, MS1_EXO_files[i,1], sep = ""), source = MsBackendMzR())
  
  # set polarity to negative = 0 or positive = 1
#  sps_neg <- sps[sps$polarity==0]
#  sps_pos <- sps[sps$polarity==1] 
  
  # export data in separate files for neg and pos
#  export(sps_neg, backend = MsBackendMzR(), file = paste(raw_data_dir, 
#                                                         gsub("KSS_210324","\\coculture_neg", MS1_EXO_files[i,1]),
#                                                         sep = ""))
#  export(sps_pos, backend = MsBackendMzR(), file = paste(raw_data_dir, 
#                                                         gsub("KSS_210324","\\coculture_pos", MS1_EXO_files[i,1]), 
#                                                         sep = ""))
#}

# raw_data_MS1_EXO <- data.frame(list.files(raw_data_dir, pattern = "EXO"))
# sps_test_x <- Spectra(paste(raw_data_dir, raw_data_MS1_EXO[1,1], sep = ""), source = MsBackendMzR())
# sps_test_x
# table(polarity(sps_test_x))

###----MS1 and MS2 combined-----

# input directories
input_dir_MS1 <- paste(getwd(), "/MS1_pos_neg/", sep = "")
input_dir_MS2 <- paste(getwd(), "/MS2/", sep = "")

# for ENDO_neg files, create list of MS1 files
raw_data_MS1_ENDO_neg <- list.files(input_dir_MS1, pattern = "neg_ENDO")
raw_data_MS1_ENDO_neg_files <- paste(input_dir_MS1, raw_data_MS1_ENDO_neg, sep ="")

# create list of MS2 files
#raw_data_MS2_ENDO_neg <- list.files(input_dir_MS2, pattern = "ENDOneg")
#raw_data_MS2_ENDO_neg_files <- paste(input_dir_MS2, raw_data_MS2_ENDO_neg, sep ="")

# combine all ENDO_neg files of MS1 and MS2
#files <- c(raw_data_MS1_ENDO_neg_files, raw_data_MS2_ENDO_neg_files)
#files_names <- c(raw_data_MS1_ENDO_neg, raw_data_MS2_ENDO_neg)


###----MS1 preparations----
start.time <- Sys.time()

# MS1 variables
# pol <- c(x = polarity, start =0, stop = 0)
ppm <- 35           # needed for missing value imputation
ms1_intensity_cutoff <- 14	          #approx. 0.01%, needed for bina list creation

# mzml_times_ENDO <- NULL

###----Creation of Phenodata/Metadata----
# create vector with sample classes according to culture information sheet
samp_groups <- c("CoCuPp", "CoCuSm", "CoCuSm", "CoCuPp", "CoCuPp", "CoCuSm",
                 rep(x = "Sm", times = 8), 
                 rep(x = "Pp", times = 8),
                 "CoCuSm", "CoCuPp", "MB")

CoCuPp_des <- "Co-culture sample from Prymnesium parvum"
CoCuSm_des <- "Co-culture sample from Skeletonema marinoi"
Sm_des <- "Mono-culture sample from Skeletonema marinoi"
Pp_des <- "Mono-culture sample from Prymnesium parvum"
MB_des <- "Media Blank"

# create vector with sample classes according to culture information sheet
samp_groups_description <- c(CoCuPp_des, CoCuSm_des, CoCuSm_des, CoCuPp_des, CoCuPp_des, CoCuSm_des,
                 rep(x = Sm_des, times = 8), 
                 rep(x = Pp_des, times = 8),
                 CoCuSm_des, CoCuPp_des, MB_des)


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

color <- c(CoCuPp1, CoCuSm1, CoCuPp2, CoCuSm2, Sm, Pp, CoCuSm3, CoCuPp3, MB)


# create phenodata based on culture type
pheno_data_ENDO <- data.frame(sample_name = MS1_ENDO_names, sample_group = samp_groups, samp_groups_description = samp_groups_description)
pheno_color_ENDO <- data.frame(color)

# create variables to store condition name in
#con_1 <- "ENDO"   # ENDO or EXO
#pol_2 <- "neg"    # neg or pos
#ms_3 <- "ms2"     # if ms2 files are present
#all_4 <- "all"    # if ms2 files are present


# Create phenodata base on culture type
#create_pheno_data <- function(sample_name, sample_group, sample_groups_description)
#  if (attr(sample_name, "ENDO") == TRUE) {
#        return(data.frame(sample_name, sample_group, sample_groups_description))
#  }


#stop first
#please check the directory where it should be saved
#dir_MS1_EXO <- "./MAW-Co-culture/Results"
#write.csv(pheno_data_EXO, file = paste(getwd(), "/Results/pheno_data_EXO.csv", sep =""))

# Display MSn levels and check amount of spectra
# mzml_msn_ENDO <- NULL
# for (i in 1:length(MS1_ENDO_files)) {
#   mzml_data_ENDO <- readMSData(files = paste(input_dir_MS1, MS1_ENDO_files[,i], sep = ""), mode="onDisk")
#   mzml_msn_ENDO <- rbind(mzml_msn_ENDO, t(as.matrix(table(msLevel(mzml_data_ENDO)))))
# }
# colnames(mzml_msn_ENDO) <- c("MSn 0", "MSn 1")
# rownames(mzml_msn_ENDO) <- MS1_ENDO_names

# # Plot MSn levels (only if MS1 and MS2 data are in same directory/file, then exchange MSn 0 for MS1 and MS2)
# #pdf(file="plots/ENDOMS1_msn_levels.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
# jpeg(file="plots/ENDOMS1_msn_levels.jpeg", width = 1000, height = 600, quality = 100, bg = "white")
# par(mfrow=c(2,1), mar=c(5,4,4,2), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
# boxplot(mzml_msn_ENDO, main="Number of spectra")

# model_boxplot <- boxplot(t(mzml_msn_ENDO[,2]), main="Number of MS1 spectra per sample", xaxt="n")
# tick <- seq_along(model_boxplot$names)
# axis(1, at=tick, labels=F)
# text(tick, par("usr")[3]-par("usr")[3]/10, model_boxplot$names, adj=0, srt=270, xpd=T)
# dev.off()

###----Import raw DATA ENDO NEGATIVE----
# create plot directory
if (dir.exists(paste(getwd(), "/endo_neg_plots/", sep = ""))){
  print("plots directory already exists")
}  else{
  dir.create("endo_neg_plots")
  print("plots folder has been created")
}

# create Results directory
if (dir.exists(paste(getwd(), "/endo_neg_Results/", sep = ""))){
  print("Results directory already exists")
}  else{
  dir.create("endo_neg_Results")
  print("Results folder has been created")
}


# input directory
input_dir_MS1_polarity <- paste(getwd(), "/MS1_pos_neg/", sep = "")
#ENDO_neg files
MS1_ENDO_neg_files <- list.files(input_dir_MS1_polarity, pattern = "_neg_ENDO")
# Import raw data as MSnbase object OnDiskMsnExp


# only MS1 data
msd <- readMSData(files = paste(input_dir_MS1_polarity, MS1_ENDO_neg_files, sep = ""),
                  pdata = new("NAnnotatedDataFrame",pheno_data_ENDO),
                  mode = "onDisk",
                  centroided = TRUE)

# inspect data 
#table(msLevel(msd))
#head(fData(msd)[, c("scanWindowLowerLimit", "scanWindowUpperLimit",
#                    "originalPeaksCount", "msLevel", 
#                    "polarity", "retentionTime")])

# for MS2 data subset polarity again to pos = 1
#msd <- filterPolarity(msd, polarity = 1)

# Restrict data to 700 seconds
msd <- filterRt(msd, c(0, 700))

# subset data for msLevel = 1 and save raw data
# msd <- filterMsLevel(msd, msLevel = 1)
table(polarity(msd))
write.csv(fData(msd), file=paste(filename = "endo_neg_Results/ENDO_neg_raw_data.csv", sep = ""), row.names=FALSE)

# # Inspect mz values per file
# msd_mz <- mz(msd)
# msd_mz <- split(msd_mz, f = fromFile(msd))
# print(length(msd_mz))

# Get base peak chromatograms
# setwd(input_dir_MS1_polarity)
chromas_ENDO_neg <- chromatogram(msd, 
                             aggregationFun="max", 
                             msLevel = 1, 
                             )

# Plot chromatograms based on phenodata groups
#pdf(file="plots/ENDO_chromas.pdf", encoding="ISOLatin1", pointsize=2, width=6, height=4, family="Helvetica")
jpeg(filename = "endo_neg_plots/ENDO_chromas_pos.jpeg", width = 2000, height = 1200, quality = 150, bg = "white")
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=2, cex.lab=2, cex.main=2)
plot(chromas_ENDO_neg, main="Raw chromatograms", xlab="retention time [s]", ylab="intensity", col = color)
legend("topleft", bty="n", pt.cex=3, cex=1.5, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(msd@phenoData@data[["sample_group"]]))
dev.off()

# Get TICs
#pdf(file="plots/ENDO_tics.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "endo_neg_plots/ENDO_tics_pos.jpeg", width = 2000, height = 1200, quality = 150, bg = "white")
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=2, cex.lab=2, cex.main=2)
tics_ENDO_neg <- split(tic(msd), f=fromFile(msd))
boxplot(tics_ENDO_neg, col=color, ylab="intensity", xlab="sample", main="Total ion current", outline = FALSE)
legend("topleft", bty="n", pt.cex=2, cex=1, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(msd@phenoData@data[["sample_group"]]))
dev.off()

### ---- pre-processing ----

# grouping/binning for peak detection, based on similarity of the base peak chromatogram 
chromas_bin_ENDO_neg <- MSnbase::bin(chromas_ENDO_neg, binSize=2)
chromas_bin_cor_ENDO_neg <- cor(log2(do.call(cbind, lapply(chromas_bin_ENDO_neg, intensity)))) # transformation with log
colnames(chromas_bin_cor_ENDO_neg) <- rownames(chromas_bin_cor_ENDO_neg) <- msd$sample_name
chromas_bin_cor_ENDO_neg[is.na(chromas_bin_cor_ENDO_neg)] <- 0

# representing the data in a heatmap for general overview
#pdf(file="plots/heatmap_chromas_bin_ENDO.pdf", encoding="ISOLatin1", pointsize=10, width=6, 
#    height=4, family="Helvetica")
jpeg(filename = "endo_neg_plots/heatmap_chromas_bin_ENDO_neg.jpeg", width = 500, height = 500, quality = 150, bg = "white")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(7,0,0,7), cex.axis=0.9, cex=0.6)
heatmap(chromas_bin_cor_ENDO_neg)
dev.off()

# Assess retention times and intensities of first file
#head(rtime(chromas_ENDO_neg[1, 1]))
#head(intensity(chromas_ENDO_neg[1, 1]))

# check for polarity
head(fData(msd)[, c("polarity", "filterString", "msLevel", "retentionTime")])
table(polarity(msd))


#--- parameters calculated with IPO ---
# ENDO neg
#ms1_params_ENDO_neg <- CentWaveParam(ppm=39.5, mzCenterFun="wMean", peakwidth=c(14, 59), 
#                                     prefilter=c(3, 140), mzdiff=0.0155, snthresh=7, noise=0, 
#                                     integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
#                                     fitgauss=FALSE, roiList=list(), roiScales=numeric())
# ENDO pos
#ms1_params_ENDO_neg <- CentWaveParam(ppm=9.5, mzCenterFun="wMean", peakwidth=c(12, 51), 
#                                     prefilter=c(4, 60), mzdiff= 0.000099, snthresh=6, noise=0, 
#                                     integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
#                                     fitgauss=FALSE, roiList=list(), roiScales=numeric())
# EXO neg 
#ms1_params_EXO_neg <- CentWaveParam(ppm=24.5, mzCenterFun="wMean", peakwidth=c(12, 53), 
#                                     prefilter=c(4, 60), mzdiff= -0.0032, snthresh=5, noise=0, 
#                                     integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
#                                     fitgauss=FALSE, roiList=list(), roiScales=numeric())
# EXO pos 
#ms1_params_EXO_pos <- CentWaveParam(ppm=15, mzCenterFun="wMean", peakwidth=c(13, 51), 
#                                    prefilter=c(2.6, 41), mzdiff= -0.0098, snthresh=5, noise=0, 
#                                    integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
#                                    fitgauss=FALSE, roiList=list(), roiScales=numeric())

# Peak detection in MS1 data
# define pre-processed data variable
#ms1_data_ENDO_neg <- NULL 
#condition_name <- deparse(substitute(ms1_data_ENDO_neg))
#condition_name

# perform peak detection with corresponding parameters 
#pick_peaks_f <- function(condition_name, MSnobject) {
#if (condition_name == "ms1_data_EXO_neg") {
  # set parameters
  # EXO neg 
#  ms1_params_EXO_neg <- CentWaveParam(ppm=15, mzCenterFun="wMean", peakwidth=c(12, 53), 
#                                      prefilter=c(4, 60), mzdiff= -0.0032, snthresh=5, noise=0, 
#                                      integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
#                                      fitgauss=FALSE, roiList=list(), roiScales=numeric())
  
  # perform peak detection
#  ms1_data_EXO_neg <- findChromPeaks(msd, param=ms1_params_EXO_neg)
#  print("parameters for EXO neg used, peak detection successful")
  
#} else if (condition_name == "ms1_data_EXO_pos") {
  # set parameters
  # EXO pos 
#  ms1_params_EXO_pos <- CentWaveParam(ppm=15, mzCenterFun="wMean", peakwidth=c(13, 51), 
#                                      prefilter=c(2.6, 41), mzdiff= -0.0098, snthresh=5, noise=0, 
#                                      integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
#                                      fitgauss=FALSE, roiList=list(), roiScales=numeric())
  # perform peak detection
#  ms1_data_EXO_pos <- findChromPeaks(msd, param=ms1_params_EXO_pos)
#  print("paramters for EXO pos used, peak detection successful")
  
#} else if (condition_name == "ms1_data_ENDO_neg") {
  # set parameters
  # ENDO neg
  ms1_params_ENDO_neg <- CentWaveParam(ppm=25, mzCenterFun="wMean", peakwidth=c(14, 59), 
                                       prefilter=c(3, 140), mzdiff=0.0155, snthresh=7, noise=0, 
                                       integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
                                       fitgauss=FALSE, roiList=list(), roiScales=numeric())

    # perform peak detection
#  ms1_data_ENDO_neg <- findChromPeaks(msd, param=ms1_params_ENDO_neg)
#  print("parameters for ENDO neg used, peak detection successful")
  
#} else if (condition_name == "ms1_data_ENDO_neg") {
  # set parameters
  # ENDO pos
#  ms1_params_ENDO_neg <- CentWaveParam(ppm=15, mzCenterFun="wMean", peakwidth=c(12, 51), 
#                                       prefilter=c(4, 60), mzdiff= 0.000099, snthresh=6, noise=0, 
#                                       integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
#                                       fitgauss=FALSE, roiList=list(), roiScales=numeric())
  # perform peak detection
#  ms1_data_ENDO_neg <- findChromPeaks(msd, param=ms1_params_ENDO_neg)
#  print("parameters for ENDO pos used, peak detection successful")
  
#} else {
#  print("peak detection failed, check again")
  
#}
#}

# manual 
ms1_data_ENDO_neg <- findChromPeaks(msd, param=ms1_params_ENDO_neg)

# check the detected peaks
#head(chromPeaks(ms1_data_ENDO_neg))
#chromPeakData(ms1_data_ENDO_neg)


# Per file summary
ms1_summary_ENDO_neg <- lapply(split.data.frame(chromPeaks(ms1_data_ENDO_neg), 
                                            f=chromPeaks(ms1_data_ENDO_neg)[, "sample"]), 
                           FUN=function(z) { c(peak_count=nrow(z), rt=quantile(z[, "rtmax"] - z[, "rtmin"])) } )
ms1_summary_ENDO_neg <- do.call(rbind, ms1_summary_ENDO_neg)
rownames(ms1_summary_ENDO_neg) <- basename(fileNames(ms1_data_ENDO_neg))
print(ms1_summary_ENDO_neg)
table(msLevel(ms1_data_ENDO_neg))

write.csv(as.data.frame(table(msLevel(ms1_data_ENDO_neg))), file="endo_neg_Results/ENDO_neg_ms1_data.csv", row.names=FALSE)


# To get a global overview of the peak detection we can plot the frequency of identified peaks per file along the retention time axis. 
# This allows to identify time periods along the MS run in which a higher number of peaks was identified and evaluate whether this is consistent across files.
#pdf(file="plots/ENDO_ms1_data.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "endo_neg_plots/ENDO_neg_ms1_data.jpeg", width = 2000, height = 1200, quality = 150, bg = "white")
par(mfrow=c(1,1), mar=c(5,18,4,1), oma=c(0,0,0,0), cex.axis=1, cex=2, cex.lab=2, cex.main=2)
plotChromPeakImage(ms1_data_ENDO_neg, main="Frequency of identified peaks per RT", binSize = 20)
dev.off()


## Group peaks
ms1_data_ENDO_neg <- groupChromPeaks(ms1_data_ENDO_neg, param=PeakDensityParam(
  sampleGroups=ms1_data_ENDO_neg$sample_group, minFraction=0.7, bw=2.5))

## RT correction
ms1_data_ENDO_neg <- adjustRtime(ms1_data_ENDO_neg, param=PeakGroupsParam(
  minFraction=0.7,smooth="loess",span=0.5,family="gaussian"))


# Plot the difference of raw and adjusted retention times
#pdf(file="plots/ENDO_ms1_raw_adjusted.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
jpeg(filename = "endo_neg_plots/ENDO_neg_ms1_raw_adjusted.jpeg", width = 1000, height = 2000, quality = 100, bg = "white")
par(mfrow=c(2,1), mar=c(5,6,4,1), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(chromas_ENDO_neg, peakType="none", main="Raw chromatograms")
plotAdjustedRtime(ms1_data_ENDO_neg, lwd=2, main="Retention Time correction")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(1,1,0,0), cex.axis=0.9, cex=0.8)
dev.off()

## Group peaks
ms1_data_ENDO_neg <- groupChromPeaks(ms1_data_ENDO_neg, param=PeakDensityParam(
  sampleGroups=ms1_data_ENDO_neg$sample_group, minFraction=0.7, bw=2.5))

# Get integrated peak intensity per feature/sample
print(head(featureValues(ms1_data_ENDO_neg, value="into")))

## Fill peaks
#ms1_data_ENDO_neg <- fillChromPeaks(ms1_data_ENDO_neg, param=FillChromPeaksParam(ppm=ppm, fixedRt=0, expandRt=5))
#head(featureValues(ms1_data_ENDO_neg))
#head(featureSummary(ms1_data_ENDO_neg, group=ms1_data_ENDO_neg$sample_group))

# Evaluate grouping
#pdf(file="plots/ENDO_ms1_grouping.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "endo_neg_plots/ENDO_neg_ms1_grouping.jpeg", width = 1000, height = 700, quality = 150, bg = "white")
ms1_pca_ENDO_neg <- prcomp(t(na.omit(log2(featureValues(ms1_data_ENDO_neg, value="into")))), center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms1_pca_ENDO_neg$x[, 1], ms1_pca_ENDO_neg$x[,2], pch=19, main="PCA: Grouping of samples",
     xlab=paste0("PC1: ", format(summary(ms1_pca_ENDO_neg)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_ENDO_neg)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color, cex=2)
grid()
text(ms1_pca_ENDO_neg$x[, 1], ms1_pca_ENDO_neg$x[,2], labels=str_sub(ms1_data_ENDO_neg$sample_name, - 3, - 1), col=color, pos=3, cex=1.5)
legend("topleft", bty="n", pt.cex=2, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(ms1_data_ENDO_neg@phenoData@data[["sample_group"]]))
dev.off()

# broken stick
png("endo_neg_plots/BrokenStick_ENDO_neg_ms1_grouping.png", width=10, height=6, units="in", res=100)
evplot = function(ev) {  
  # Broken stick model (MacArthur 1957)  
  n = length(ev)  
  bsm = data.frame(j=seq(1:n), p=0)  
  bsm$p[1] = 1/n  
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))  
  bsm$p = 100*bsm$p/n  
  # Plot eigenvalues and % of variation for each axis  
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))  
  barplot(ev, main="Eigenvalues ENDO MS1 grouping", col="blue", las=2)  
  abline(h=mean(ev), col="red")  
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")  
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE,   
          main="% variation", col=c("blue",2), las=2)  
  legend("topright", c("% eigenvalue", "Broken stick model"),   
         pch=15, col=c("blue",2), bty="n")  
  par(op)  
} 

ev_pc = ms1_pca_ENDO_neg$sdev^2  
evplot(ev_pc)  
dev.off()


# Show peaks
#tail(chromPeaks(ms1_data_ENDO_neg))
#tail(chromPeakData(ms1_data_ENDO_neg))

# Show process history
processHistory(ms1_data_ENDO_neg)

# save as R object for later use
MS1_ENDO_neg_peak_detection <- ms1_data_ENDO_neg
save(MS1_ENDO_neg_peak_detection, file = "endo_neg_Results/MS1_ENDO_neg_peak_detection.RData")



# ---------- Build MS1 feature tables ----------
# Build feature matrix
ms1_matrix_ENDO_neg <- featureValues(ms1_data_ENDO_neg, method="medret", value="into")
colnames(ms1_matrix_ENDO_neg) <- MS1_ENDO_names
dim(ms1_matrix_ENDO_neg)
# transpose feature table
feat_list_ENDO_neg <- t(ms1_matrix_ENDO_neg)

# Build feature summary
ms1_summary_ENDO_neg <- featureSummary(ms1_data_ENDO_neg)
ms1_def_ENDO_neg <- featureDefinitions(ms1_data_ENDO_neg)

# Missing value imputation by filling na positions with median of surrounding features
#feat_list_ENDO_neg[is.na(feat_list_ENDO_neg)] <- median(na.omit(as.numeric(unlist(feat_list_ENDO_neg))))

### Transform data
feat_list_ENDO_neg <- log2(feat_list_ENDO_neg)

# change 0 to small value to distinguish between values of 1 and NA
feat_list_ENDO_neg[which(feat_list_ENDO_neg == 0)] <- 0.01

# Missing value imputation
feat_list_ENDO_neg[which(is.na(feat_list_ENDO_neg))] <- 0

# save as csv
write.csv(feat_list_ENDO_neg, file=paste(filename = "endo_neg_Results/feature_list_ENDO_neg.csv", sep = ""))
save(ms1_def_ENDO_neg, file = "endo_neg_Results/ms1_def_ENDO_neg.RData")

# Plot histogram
#pdf(file="endo_neg_plots/ENDO_feat_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "endo_neg_plots/ENDO_neg_feat_list_hist.jpeg", width = 1000, height = 1000, quality = 150, bg = "white")
hist(as.numeric(feat_list_ENDO_neg), main="Histogram of feature table")
dev.off()

# PCA of feature table results
#pdf(file="plots/ENDO_ms1_feature_table_pca.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "endo_neg_plots/ENDO_neg_ms1_feature_table_pca_exc_MB.jpeg", width = 1000, height = 700, quality = 150, bg = "white")
ms1_pca_ENDO_neg <- prcomp(feat_list_ENDO_neg[1:24,], center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms1_pca_ENDO_neg$x[, 1], ms1_pca_ENDO_neg$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms1_pca_ENDO_neg)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_ENDO_neg)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color, cex=2)
grid()
text(ms1_pca_ENDO_neg$x[,1], ms1_pca_ENDO_neg$x[,2], labels=str_sub(ms1_data_ENDO_neg$sample_name[1:24], - 3, - 1), col=color, pos=3, cex=1.5)
legend("topleft", bty="n", pt.cex=2, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(ms1_data_ENDO_neg@phenoData@data[["sample_group"]]))
dev.off()

# PCA of feature table on species level
index_SM <- grep("Sm", pheno_data_ENDO$sample_group)
index_PP <- grep("Pp", pheno_data_ENDO$sample_group)
#test_SM <- data.frame(feat_list_ENDO_pos[index_SM,])
#test_PP <- data.frame(feat_list_ENDO_pos[index_PP,])

# PCA for Sm
jpeg(filename = "endo_neg_plots/ENDO_neg_ms1_feature_table_pca_SM.jpeg", width = 1000, height = 700, quality = 150, bg = "white")
ms1_pca_ENDO_neg_SM <- prcomp(feat_list_ENDO_neg[index_SM,], center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms1_pca_ENDO_neg_SM$x[, 1], ms1_pca_ENDO_neg_SM$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms1_pca_ENDO_neg_SM)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_ENDO_neg_SM)$importance[2, 2] * 100, digits=3), " % variance"),
     col=unique(color[index_SM]), cex=2)
grid()
text(ms1_pca_ENDO_neg_SM$x[,1], ms1_pca_ENDO_neg_SM$x[,2], labels=str_sub(ms1_data_ENDO_neg$sample_name[index_SM], - 3, - 1), col=unique(color[index_SM]), pos=3, cex=1.5)
legend("topleft", bty="n", pt.cex=2, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color[index_SM]), legend= c("CoCuSm", "Sm"))
dev.off()

# PCA for Pp
jpeg(filename = "endo_neg_plots/ENDO_neg_ms1_feature_table_pca_PP.jpeg", width = 1000, height = 700, quality = 150, bg = "white")
ms1_pca_ENDO_neg_PP <- prcomp(feat_list_ENDO_neg[index_PP,], center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms1_pca_ENDO_neg_PP$x[, 1], ms1_pca_ENDO_neg_PP$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms1_pca_ENDO_neg_PP)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_ENDO_neg_PP)$importance[2, 2] * 100, digits=3), " % variance"),
     col=unique(color[index_PP]), cex=2)
grid()
text(ms1_pca_ENDO_neg_PP$x[,1], ms1_pca_ENDO_neg_PP$x[,2], labels=str_sub(ms1_data_ENDO_neg$sample_name[index_PP], - 3, - 1), col=unique(color[index_PP]), pos=3, cex=1.5)
legend("topleft", bty="n", pt.cex=2, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color[index_PP]), legend= c("CoCuPp", "Pp"))
dev.off()


# broken stick
png("endo_neg_plots/BrokenStick_ENDO_neg_ms1_feature_table_pca_exc_MB.png", width=10, height=6, units="in", res=100)
evplot = function(ev) {  
  # Broken stick model (MacArthur 1957)  
  n = length(ev)  
  bsm = data.frame(j=seq(1:n), p=0)  
  bsm$p[1] = 1/n  
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))  
  bsm$p = 100*bsm$p/n  
  # Plot eigenvalues and % of variation for each axis  
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 3, 1, 1))  
  barplot(ev, main="Eigenvalues ENDO MS1 feature table", col="blue", las=2)  
  abline(h=mean(ev), col="red")  
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")  
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE,   
          main="% variation", col=c("blue",2), las=2)  
  legend("topright", c("% eigenvalue", "Broken stick model"),   
         pch=15, col=c("blue",2), bty="n")  
  par(op)  
} 

ev_pc = ms1_pca_ENDO_neg$sdev^2  
evplot(ev_pc)  
dev.off()


ms1_intensity_cutoff <- 14

# Create single 0/1 matrix
bina_list_ENDO_neg <- t(ms1_matrix_ENDO_neg)
bina_list_ENDO_neg[is.na(bina_list_ENDO_neg)] <- 1
bina_list_ENDO_neg <- log2(bina_list_ENDO_neg)
bina_list_ENDO_neg[bina_list_ENDO_neg < ms1_intensity_cutoff] <- 0
bina_list_ENDO_neg[bina_list_ENDO_neg != 0] <- 1

# save as csv
write.csv(bina_list_ENDO_neg, file=paste(filename = "endo_neg_Results/bina_list_ENDO_neg.csv", sep = ""))

# Only unique compounds in group mzml_pheno$ and not the others
uniq_list_ENDO_neg <- apply(X=bina_list_ENDO_neg, MARGIN=2, FUN=function(x) { if (length(unique(pheno_data_ENDO$sample_group[grepl("1", x)])) == 1) x else rep(0, length(x)) } )
colnames(uniq_list_ENDO_neg) <- colnames(bina_list_ENDO_neg)
rownames(uniq_list_ENDO_neg) <- rownames(bina_list_ENDO_neg)


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
# functional hill diversity
#model_div_neg$hillfunc    <- as.numeric(unlist(calcDiv(feat_list_ENDO_neg, compDisMat=scales::rescale(as.matrix(dist(t(feat_list_ENDO_neg)), diag=TRUE, upper=TRUE)), q=1, type="FuncHillDiv")))



# Remove NAs if present
model_div_ENDO_neg[is.na(model_div_ENDO_neg)] <- 0

# save as csv
write.csv(model_div_ENDO_neg, file=paste(filename = "endo_neg_Results/model_div_ENDO_neg.csv", sep = ""))


# save the objects and tables
save.image(file = "endo_neg_Results/ENDO_neg_MS1_environment.RData")
save(ms1_def_ENDO_neg, file = "endo_neg_Results/ms1_def_ENDO_neg.RData")


end.time <- Sys.time()

time.taken <- end.time - start.time
print(time.taken)

# ----- create feature info table -----------
# extract feature names, mz and rt 
feature_info_ENDO_neg <- data.frame(cbind(ms1_data_ENDO_neg@msFeatureData[["featureDefinitions"]]@rownames, 
                                          ms1_data_ENDO_neg@msFeatureData[["featureDefinitions"]]@listData[["mzmed"]], 
                                          ms1_data_ENDO_neg@msFeatureData[["featureDefinitions"]]@listData[["mzmin"]], 
                                          ms1_data_ENDO_neg@msFeatureData[["featureDefinitions"]]@listData[["mzmax"]], 
                                          ms1_data_ENDO_neg@msFeatureData[["featureDefinitions"]]@listData[["rtmed"]], 
                                          ms1_data_ENDO_neg@msFeatureData[["featureDefinitions"]]@listData[["rtmin"]], 
                                          ms1_data_ENDO_neg@msFeatureData[["featureDefinitions"]]@listData[["rtmax"]]))
# add column names
colnames(feature_info_ENDO_neg) <- c("feat_id", "mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax")

# save as csv
write.csv(feature_info_ENDO_neg, file = "endo_neg/endo_neg_Results/feature_info_ENDO_neg.csv", row.names = FALSE)


############# linking MS2 data #################
start.time_linking <- Sys.time()

# --------- preparations -----------
# load object with MS1 and MS2 files preprocessed
load(file = "endo_neg_1ms2_Results/MS_endo_neg_peak_detection.RData")
ms_data_ENDO_neg <- MS_endo_neg_peak_detection
#ms_def_endo_pos
#table(msLevel(ms_data_ENDO_neg))
#table(msLevel(ms1_data_ENDO_neg))

#head(ms1_def_ENDO_neg)


# ---------- MS2 spectra detection ----------
# Estimate precursor intensity
#precursor_intensity_ENDO_neg <- xcms::estimatePrecursorIntensity(ms_data_ENDO_neg)
#print(head(na.omit(precursor_intensity_ENDO_neg)))

# Reconstruct MS2 spectra from MS1 data
ms2_data_ENDO_neg <- chromPeakSpectra(ms_data_ENDO_neg, msLevel=2L, return.type="Spectra")
print(ms2_data_ENDO_neg)
print(length(ms2_data_ENDO_neg$peak_id))
head(ms2_data_ENDO_neg$peak_id)


# Extract all usable MS2 spectra
ms2_spectra_ENDO_neg <- list()
for (i in 1:nrow(ms1_def_ENDO_neg)) {
  #ms2_spectra_ENDO_neg <- foreach(i=1:nrow(ms1_def_ENDO_neg)) %dopar% {
  #print(i)
  # Extract existing MS2 spectra for feature
  feature_of_interest <- ms1_def_ENDO_neg[i, "mzmed"]
  peaks_of_interest <- chromPeaks(ms_data_ENDO_neg, mz=feature_of_interest, ppm=ppm)
  
  # Continue if feature has MS2 peaks
  if (length(which(ms2_data_ENDO_neg$peak_id %in% rownames(peaks_of_interest))) > 0) {
    # Extract spectra
    spectra_of_interest <- ms2_data_ENDO_neg[ms2_data_ENDO_neg$peak_id %in% rownames(peaks_of_interest)]
    combined_spectra_of_interest <- filterIntensity(spectra_of_interest, intensity=c(1,Inf), backend=MsBackendDataFrame)
    combined_spectra_of_interest <- setBackend(combined_spectra_of_interest, backend=MsBackendDataFrame())
    
    # Combine spectra
    combined_spectra_of_interest <- Spectra::combineSpectra(combined_spectra_of_interest, FUN=combinePeaks, ppm=ppm, peaks="union", minProp=0.8, intensityFun=median, mzFun=median, backend=MsBackendDataFrame)#f=rownames(peaks_of_interest))
    
    # Remove noise from spectra
    #combined_spectra_of_interest <- pickPeaks(combined_spectra_of_interest, snr=1.0, method="SuperSmoother") #MAD
    #combined_spectra_of_interest <- Spectra::smooth(combined_spectra_of_interest, method="SavitzkyGolay") #(Weighted)MovingAverage
    
    # Only keep spectral data
    combined_spectra_peaks <- as.data.frame(Spectra::peaksData(combined_spectra_of_interest)[[1]])
    
    # Plot merged spectrum
    #Spectra::plotSpectra(combined_spectra_of_interest)
    #plot(x=combined_spectra_peaks[,1], y=combined_spectra_peaks[,2], type="h", xlab="m/z", ylab="intensity", main=paste("Precursor m/z",combined_spectra_of_interest@backend@spectraData$precursorMz[[1]]))
    #length(spectra_of_interest$peak_id)
    
    ms2_spectra_ENDO_neg[[i]] <- combined_spectra_of_interest
    #return(combined_spectra_of_interest)
  }
}

# Remove empty spectra
names(ms2_spectra_ENDO_neg) <- rownames(ms1_def_ENDO_neg)[1:length(ms2_spectra_ENDO_neg)]
ms2_spectra_ENDO_neg <- ms2_spectra_ENDO_neg[lengths(ms2_spectra_ENDO_neg) != 0]

# Relate all MS2 spectra to MS1 precursors
ms1_def_ENDO_neg$has_ms2 <- as.integer(rownames(ms1_def_ENDO_neg) %in% names(ms2_spectra_ENDO_neg))
print(paste0("Number of MS2 spectra related to precursor: ", length(which(ms1_def_ENDO_neg$has_ms2>0))))

# ADDED
polarity="negative"
pol="neg"

# create a list with file names of feature origin
#ms2_names <- NULL


# extract collision energy
colenergy <- collisionEnergy(ms_data_endo_pos)
head(colenergy)
colenergy <- na.omit(colenergy)


# Save all MS2 spectra in MGF file
mgf_text <- NULL
for (i in names(ms2_spectra_ENDO_neg)) {
  mgf_text <- c(mgf_text, paste0("COM=", i))
  mgf_text <- c(mgf_text, "BEGIN IONS")
  mgf_text <- c(mgf_text, "MSLEVEL=2")
  mgf_text <- c(mgf_text, paste0("TITLE=", i))
  mgf_text <- c(mgf_text, paste0("RTINSECONDS=", ms1_def_ENDO_neg[i, "rtmed"]))
  mgf_text <- c(mgf_text, paste0("PEPMASS=", ms1_def_ENDO_neg[i, "mzmed"]))
  if (polarity == "positive") {
    mgf_text <- c(mgf_text, paste0("CHARGE=", "1"))
  } else {
    mgf_text <- c(mgf_text, paste0("CHARGE=", "0"))
  }
  mgf_text <- c(mgf_text, paste0("COLENERGY=", unique(colenergy)))
  mgf_text <- c(mgf_text, paste(as.data.frame(peaksData(ms2_spectra_ENDO_neg[[i]])[[1]])$mz, as.data.frame(peaksData(ms2_spectra_ENDO_neg[[i]])[[1]])$intensity, sep=" "))
  mgf_text <- c(mgf_text, "END IONS")
  mgf_text <- c(mgf_text, "")
}

# Write MGF file
cat(mgf_text, file="endo_neg_1ms2_Results/ms2_spectra_ENDO_neg.mgf", sep="\n")

end.time_linking <- Sys.time()

time.taken_linking <- end.time_linking - start.time_linking
print(time.taken_linking)


