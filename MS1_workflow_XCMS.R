### Analysis of diatom Mono- and Co-culture MS data
# MS1 data pre-processing and statistical analysis with xcms
# Parameters for peak detection calculated with IPO
# Ms2 data relation to MS1 precursors and preparation for annotation by MAW

###---- library ----


library(Spectra)
library(xcms)
library(mzR)
library(MSnbase)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)
library(SummarizedExperiment)
library(knitr)
library(ggplot2)
library(vegan)
library(stringr)
library(CAMERA)

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

# create plot directory
if (dir.exists(paste(getwd(), "/plots/", sep = ""))){
  print("plots directory already exists")
  start_time <- Sys.time()
}  else{
  dir.create("plots")
  start_time <- Sys.time()
  print("plots folder has been created")
}

###----MS1 ENDO files----
#MS1_ENDO_files <- data.frame(list.files(input_dir_MS1, pattern = "ENDO"))
MS1_ENDO_files <- list.files(input_dir_MS1, pattern = "ENDO")
# basenames of files without path and without extension
#MS1_ENDO_names <- data.frame(str_remove(MS1_ENDO_files[,], ".mzML"))
MS1_ENDO_names <- str_remove(MS1_ENDO_files, ".mzML")


###----MS1 EXO files----
#MS1_EXO_files <- data.frame(list.files(input_dir_MS1, pattern = "EXO"))
MS1_EXO_files <- list.files(input_dir_MS1, pattern = "EXO")
#MS1_EXO_files_full <- data.frame(list.files(input_dir_MS1, pattern = "EXO", full.names = TRUE))
MS1_EXO_names <- str_remove(MS1_EXO_files, ".mzML")

###----MS2 files----
#MS2_files <- data.frame(list.files(input_dir_MS2, pattern = ".mzML"))
MS2_files <- list.files(input_dir_MS2, pattern = ".mzML")
# ENDO files
# MS2_ENDO_files <- data.frame(subset(MS2_files, grepl("ENDO", MS2_files[,1]), drop = TRUE))
# MS2_ENDO_neg_files <- data.frame(subset(MS2_ENDO_files, grepl("neg", MS2_ENDO_files[,1]), drop = TRUE))
# MS2_ENDO_pos_files <- data.frame(subset(MS2_ENDO_files, grepl("pos", MS2_ENDO_files[,1]), drop = TRUE))
MS2_ENDO_files <- subset(MS2_files, grepl("ENDO", MS2_files))
MS2_ENDO_neg_files <- subset(MS2_ENDO_files, grepl("neg", MS2_ENDO_files))
MS2_ENDO_pos_files <- subset(MS2_ENDO_files, grepl("pos", MS2_ENDO_files))

# EXO files
#MS2_EXO_files <- data.frame(subset(MS2_files, grepl("EXO", MS2_files[,1]), drop = TRUE))
#MS2_EXO_neg_files <- data.frame(subset(MS2_EXO_files, grepl("neg", MS2_EXO_files[,1]), drop = TRUE))
#MS2_EXO_pos_files <- data.frame(subset(MS2_EXO_files, grepl("pos", MS2_EXO_files[,1]), drop = TRUE))

###----Separate neg and pos mode in individual files-----
# create result directory
if (dir.exists(paste(getwd(), "/raw_data/", sep = ""))){
  raw_data_dir <- paste(getwd(), "/raw_data/", sep = "")
}  else{
  dir.create("raw_data")
  raw_data_dir <- paste(getwd(), "/raw_data/", sep = "")
}


## separation of ENDO FILES
for (i in 1:nrow(MS1_ENDO_files)){
  sps <- Spectra(paste(input_dir_MS1, MS1_ENDO_files[i,1], sep = ""), source = MsBackendMzR())
  
  # set polarity to negative = 0 or positive = 1
  sps_neg <- sps[sps$polarity==0]
  sps_pos <- sps[sps$polarity==1] 
  
  # export data in separate files for neg and pos
  export(sps_neg, backend = MsBackendMzR(), file = paste(raw_data_dir, 
                                                         gsub("KSS_210324","\\coculture_neg", MS1_ENDO_files[i,1]),
                                                         sep = ""))
  export(sps_pos, backend = MsBackendMzR(), file = paste(raw_data_dir, 
                                                         gsub("KSS_210324","\\coculture_pos", MS1_ENDO_files[i,1]), 
                                                         sep = ""))
}

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
for (i in 1:nrow(MS1_EXO_files)){
  sps <- Spectra(paste(input_dir_MS1, MS1_EXO_files[i,1], sep = ""), source = MsBackendMzR())
  
  # set polarity to negative = 0 or positive = 1
  sps_neg <- sps[sps$polarity==0]
  sps_pos <- sps[sps$polarity==1] 
  
  # export data in separate files for neg and pos
  export(sps_neg, backend = MsBackendMzR(), file = paste(raw_data_dir, 
                                                         gsub("KSS_210324","\\coculture_neg", MS1_EXO_files[i,1]),
                                                         sep = ""))
  export(sps_pos, backend = MsBackendMzR(), file = paste(raw_data_dir, 
                                                         gsub("KSS_210324","\\coculture_pos", MS1_EXO_files[i,1]), 
                                                         sep = ""))
}

# raw_data_MS1_EXO <- data.frame(list.files(raw_data_dir, pattern = "EXO"))
# sps_test_x <- Spectra(paste(raw_data_dir, raw_data_MS1_EXO[1,1], sep = ""), source = MsBackendMzR())
# sps_test_x
# table(polarity(sps_test_x))

###----MS1 and MS2 combined-----

# input directories
input_dir_MS1 <- paste(getwd(), "/MS1_pos_neg/", sep = "")
input_dir_MS2 <- paste(getwd(), "/MS2/", sep = "")

# for ENDO_pos files, create list of MS1 files
raw_data_MS1_ENDO_pos <- list.files(input_dir_MS1, pattern = "pos_ENDO")
raw_data_MS1_ENDO_pos_files <- paste(input_dir_MS1, raw_data_MS1_ENDO_pos, sep ="")

# create list of MS2 files
raw_data_MS2_ENDO_pos <- list.files(input_dir_MS2, pattern = "ENDOpos")
raw_data_MS2_ENDO_pos_files <- paste(input_dir_MS2, raw_data_MS2_ENDO_pos, sep ="")

# combine all ENDO_pos files of MS1 and MS2
all_files <- c(raw_data_MS1_ENDO_pos_files, raw_data_MS2_ENDO_pos_files)
all_files_names <- c(raw_data_MS1_ENDO_pos, raw_data_MS2_ENDO_pos)

# create output directory
outpur_dir <- dir.create("Results_combined")

###----MS1 preparations----
start.time <- Sys.time()
# Analysis of MS1 EXO data only from here on
# MS1 variables
# pol <- c(x = polarity, start =0, stop = 0)
ppm <- 35           # needed for missing value imputation
ms1_intensity_cutoff <- 100	          #approx. 0.01%, needed for bina list creation

# mzml_times_ENDO <- NULL
###----Creation of Phenodata/Metadata----
# create vector with sample classes according to culture information sheet
samp_groups <- c("CoCuPp", "CoCuSm", "CoCuSm", "CoCuPp", "CoCuPp", "CoCuSm",
                 rep(x = "Sm", times = 8), 
                 rep(x = "Pp", times = 8),
                 "CoCuSm", "CoCuPp", "MB")
# vector with MS2 samples
samp_groups_all <- c(samp_groups, rep(x = "MS2", times = 5))

CoCuPp_des <- "Co-culture sample from Prymnesium parvum"
CoCuSm_des <- "Co-culture sample from Skeletonema marinoi"
Sm_des <- "Mono-culture sample from Skeletonema marinoi"
Pp_des <- "Mono-culture sample from Prymnesium parvum"
MB_des <- "Media Blank"
# with MS2 samples
MS2_des <- "MS2 samples from one condition"

# create vector with sample classes according to culture information sheet
samp_groups_description <- c(CoCuPp_des, CoCuSm_des, CoCuSm_des, CoCuPp_des, CoCuPp_des, CoCuSm_des,
                 rep(x = Sm_des, times = 8), 
                 rep(x = Pp_des, times = 8),
                 CoCuSm_des, CoCuPp_des, MB_des)
# with MS2 samples
samp_groups_description_all <- c(samp_groups_description, rep(x = MS2_des, times = 5))

# create vector with colors 
CoCuPp1 <- ("royalblue4")
CoCuSm1 <- rep("darksalmon", 2)
CoCuPp2 <- rep("royalblue4", 2)
CoCuSm2 <- ("darksalmon")
Sm <- rep("violetred", 8)
Pp <- rep("yellow2", 8)
CoCuSm3 <- ("darksalmon")
CoCuPp3 <- ("royalblue4")
MB <- rep("springgreen", 1)
# with MS2 samples
MS2 <- rep("purple", 5)

color <- c(CoCuPp1, CoCuSm1, CoCuPp2, CoCuSm2, Sm, Pp, CoCuSm3, CoCuPp3, MB)
# with MS2 samples
color_all <- c(color, MS2)

# create phenodata based on culture type
pheno_data_EXO <- data.frame(sample_name = MS1_EXO_names, sample_group = samp_groups, samp_groups_description = samp_groups_description)
pheno_color_EXO <- data.frame(color)
# with MS2 samples
pheno_data_all_ENDO <- data.frame(sample_name = all_files_names, sample_group = samp_groups_all, samp_groups_description = samp_groups_description_all)
pheno_color_all_ENDO <- data.frame(color_all)

#stop first
#please check the directory where it should be saved
dir_MS1_EXO <- "./MAW-Co-culture/Results"
write.csv(pheno_data_all_ENDO, file = paste(getwd(), "/Results_combined/pheno_data_all_ENDO.csv", sep =""))

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
# input directory
input_dir_MS1_polarity <- paste(getwd(), "/MS1_pos_neg/", sep = "")
#EXO_pos files
MS1_EXO_pos_files <- list.files(input_dir_MS1_polarity, pattern = "_pos_EXO")
# Import raw data as MSnbase object OnDiskMsnExp

# with MS2 data
msd <- readMSData(files = paste(all_files, sep = ""),
                  pdata = new("NAnnotatedDataFrame",pheno_data_all_ENDO),
                  mode = "onDisk",
                  centroided = TRUE)

# Import as XCMSnExp object for visual analysis, only needed with pos and neg in one file
#msd_XCMS <- readMSData(files = paste(input_dir_MS1_polarity,  MS1_ENDO_neg_files, sep = ""), 
#                       pdata = new("NAnnotatedDataFrame",pheno_data_ENDO),
#                       msLevel = 1,
#                       centroided = TRUE)

# inspect data 
table(msLevel(msd))
head(fData(msd)[, c("scanWindowLowerLimit", "scanWindowUpperLimit",
                    "originalPeaksCount", "msLevel", 
                    "polarity", "retentionTime")])

# for MS2 data subset polarity again to pos = 1
msd <- filterPolarity(msd, polarity = 1)

# Restrict data to 1020 seconds (17 minutes)
msd <- filterRt(msd, c(0, 1020))

# subset data for msLevel = 1 and save raw data
# msd <- filterMsLevel(msd, msLevel = 1)
table(msLevel(msd))
write.csv(fData(msd), file=paste(filename = "Results_combined/all_ENDO_pos_raw_data.csv", sep = ""), row.names=FALSE)

# # Inspect mz values per file
# msd_mz <- mz(msd)
# msd_mz <- split(msd_mz, f = fromFile(msd))
# print(length(msd_mz))

# Get base peak chromatograms
register(bpstart(SnowParam()))
# setwd(input_dir_MS1_polarity)
chromas_all_ENDO_pos <- chromatogram(msd, 
                             aggregationFun="max", 
                             #msLevel = 1, not for all_files
                             BPPARAM = SnowParam(workers = 3))

# Plot chromatograms based on phenodata groups
#pdf(file="plots/ENDO_chromas.pdf", encoding="ISOLatin1", pointsize=2, width=6, height=4, family="Helvetica")
jpeg(filename = "plots_combined/all_ENDO_chromas_pos.jpeg", width = 1000, height = 600, quality = 150, bg = "white")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=1.5)
plot(chromas_all_ENDO_pos, main="Raw chromatograms", xlab="retention time [s]", ylab="intensity", col = color)
legend("topleft", bty="n", pt.cex=2, cex=0.9, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color_all), legend= unique(msd@phenoData@data[["sample_group"]]))
dev.off()

# Get TICs
#pdf(file="plots/ENDO_tics.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots_combined/all_ENDO_tics_pos.jpeg", width = 1000, height = 600, quality = 150, bg = "white")
par(mfrow=c(1,1), mar=c(5,4,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=1.5, cex.lab=2, cex.main=2)
tics_all_ENDO_pos <- split(tic(msd), f=fromFile(msd))
boxplot(tics_all_ENDO_pos, col=color, ylab="intensity", xlab="sample", main="Total ion current", outline = FALSE)
legend("topleft", bty="n", pt.cex=2, cex=0.9, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(msd@phenoData@data[["sample_group"]]))
dev.off()

### ---- pre-processing ----

# grouping/binning for peak detection, based on similarity of the base peak chromatogram 
chromas_bin_all_ENDO_pos <- MSnbase::bin(chromas_all_ENDO_pos, binSize=2)
chromas_bin_cor_all_ENDO_pos <- cor(log2(do.call(cbind, lapply(chromas_bin_all_ENDO_pos, intensity)))) # transformation with log
colnames(chromas_bin_cor_all_ENDO_pos) <- rownames(chromas_bin_cor_all_ENDO_pos) <- msd$sample_name
chromas_bin_cor_all_ENDO_pos[is.na(chromas_bin_cor_all_ENDO_pos)] <- 0

# representing the data in a heatmap for general overview
#pdf(file="plots/heatmap_chromas_bin_ENDO.pdf", encoding="ISOLatin1", pointsize=10, width=6, 
#    height=4, family="Helvetica")
jpeg(filename = "plots_combined/heatmap_chromas_bin_all_ENDO_pos.jpeg", width = 500, height = 500, quality = 150, bg = "white")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
heatmap(chromas_bin_cor_all_ENDO_pos)
dev.off()

# Assess retention times and intensities of first file
head(rtime(chromas_all_ENDO_pos[1, 1]))
head(intensity(chromas_all_ENDO_pos[1, 1]))

# check for polarity
head(fData(msd)[, c("polarity", "filterString", "msLevel", "retentionTime")])
table(polarity(msd))

#--- parameters calculated with IPO ---
# ENDO neg
ms1_params_ENDO_neg <- CentWaveParam(ppm=39.5, mzCenterFun="wMean", peakwidth=c(14, 59), 
                                     prefilter=c(3, 140), mzdiff=0.0155, snthresh=7, noise=0, 
                                     integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
                                     fitgauss=FALSE, roiList=list(), roiScales=numeric())
# ENDO pos
ms1_params_ENDO_pos <- CentWaveParam(ppm=9.5, mzCenterFun="wMean", peakwidth=c(12, 51), 
                                     prefilter=c(4, 60), mzdiff= 0.000099, snthresh=6, noise=0, 
                                     integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
                                     fitgauss=FALSE, roiList=list(), roiScales=numeric())
# EXO neg and pos for now
ms1_params_EXO_neg <- CentWaveParam(ppm=24.5, mzCenterFun="wMean", peakwidth=c(12, 53), 
                                     prefilter=c(4, 60), mzdiff= -0.0032, snthresh=5, noise=0, 
                                     integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
                                     fitgauss=FALSE, roiList=list(), roiScales=numeric())

# Peak detection in MS1 data
# CHOSE THE PARAMETERS 
ms1_data_all_ENDO_pos <- findChromPeaks(msd, param=ms1_params_EXO_neg)

# check the detected peaks
head(chromPeaks(ms1_data_all_ENDO_pos))
chromPeakData(ms1_data_all_ENDO_pos)


# Per file summary
ms1_summary_all_ENDO_pos <- lapply(split.data.frame(chromPeaks(ms1_data_all_ENDO_pos), 
                                            f=chromPeaks(ms1_data_all_ENDO_pos)[, "sample"]), 
                           FUN=function(z) { c(peak_count=nrow(z), rt=quantile(z[, "rtmax"] - z[, "rtmin"])) } )
ms1_summary_all_ENDO_pos <- do.call(rbind, ms1_summary_all_ENDO_pos)
rownames(ms1_summary_all_ENDO_pos) <- basename(fileNames(ms1_data_all_ENDO_pos))
print(ms1_summary_all_ENDO_pos)
table(msLevel(ms1_data_all_ENDO_pos))

write.csv(as.data.frame(table(msLevel(ms1_data_all_ENDO_pos))), file="Results_combined/all_ENDO_pos_ms1_data.csv", row.names=FALSE)


# To get a global overview of the peak detection we can plot the frequency of identified peaks per file along the retention time axis. 
# This allows to identify time periods along the MS run in which a higher number of peaks was identified and evaluate whether this is consistent across files.
#pdf(file="plots/ENDO_ms1_data.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots_combined/all_ENDO_pos_ms1_data.jpeg", width = 1000, height = 600, quality = 150, bg = "white")
par(mfrow=c(1,1), mar=c(4,18,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plotChromPeakImage(ms1_data_all_ENDO_pos, main="Frequency of identified peaks per RT", binSize = 10)
dev.off()


## Group peaks
ms1_data_all_ENDO_pos <- groupChromPeaks(ms1_data_all_ENDO_pos, param=PeakDensityParam(
  sampleGroups=ms1_data_all_ENDO_pos$sample_group, minFraction=0.7, bw=2.5))

## RT correction
ms1_data_all_ENDO_pos <- adjustRtime(ms1_data_all_ENDO_pos, param=PeakGroupsParam(
  minFraction=0.7,smooth="loess",span=0.8,family="gaussian"))

# Plot the difference of raw and adjusted retention times
#pdf(file="plots/ENDO_ms1_raw_adjusted.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
jpeg(filename = "plots_combined/all_ENDO_pos_ms1_raw_adjusted.jpeg", width = 500, height = 1000, quality = 100, bg = "white")
par(mfrow=c(2,1), mar=c(4.5,4.2,4,1), cex=0.8)
plot(chromas_all_ENDO_pos, peakType="none", main="Raw chromatograms")
plotAdjustedRtime(ms1_data_all_ENDO_pos, lwd=2, main="Retention Time correction")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
dev.off()

## Group peaks
ms1_data_all_ENDO_pos <- groupChromPeaks(ms1_data_all_ENDO_pos, param=PeakDensityParam(
  sampleGroups=ms1_data_all_ENDO_pos$sample_group, minFraction=0.7, bw=2.5))

# Get integrated peak intensity per feature/sample
print(head(featureValues(ms1_data_all_ENDO_pos, value="into")))

## Fill peaks
# missing value imputation, see xcmsSet
ms1_data_all_ENDO_pos <- fillChromPeaks(ms1_data_all_ENDO_pos, param=FillChromPeaksParam(ppm=ppm, fixedRt=0, expandRt=5))
head(featureValues(ms1_data_all_ENDO_pos))
head(featureSummary(ms1_data_all_ENDO_pos, group=ms1_data_all_ENDO_pos$sample_group))

# Evaluate grouping
#pdf(file="plots/ENDO_ms1_grouping.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots_combined/all_ENDO_pos_ms1_grouping.jpeg", width = 1000, height = 500, quality = 150, bg = "white")
ms1_pca_all_ENDO_pos <- prcomp(t(na.omit(log2(featureValues(ms1_data_all_ENDO_pos, value="into")))), center=TRUE)
plot(ms1_pca_all_ENDO_pos$x[, 1], ms1_pca_all_ENDO_pos$x[,2], pch=19, main="PCA: Grouping of samples",
     xlab=paste0("PC1: ", format(summary(ms1_pca_all_ENDO_pos)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_all_ENDO_pos)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color, cex=0.8)
grid()
text(ms1_pca_all_ENDO_pos$x[, 1], ms1_pca_all_ENDO_pos$x[,2], labels=ms1_data_all_ENDO_pos$sample_name, col=color, pos=3, cex=0.9)
legend("topleft", bty="n", pt.cex=1, cex=0.8, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(ms1_data_all_ENDO_pos@phenoData@data[["sample_group"]]))
dev.off()

# broken stick
png("plots/BrokenStick_all_ENDO_pos_ms1_grouping.png", width=10, height=6, units="in", res=100)
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

ev_pc = ms1_pca_all_ENDO_pos$sdev^2  
evplot(ev_pc)  
dev.off()


# Show peaks
tail(chromPeaks(ms1_data_all_ENDO_pos))
tail(chromPeakData(ms1_data_all_ENDO_pos))

# Show process history
processHistory(ms1_data_all_ENDO_pos)

# save as R object for later use
MS1_all_ENDO_pos_peak_detection <- ms1_data_all_ENDO_pos
save(MS1_all_ENDO_pos_peak_detection, file = "Results_combined/MS1_all_ENDO_pos_peak_detection.RData")


# ---------- Build MS1 feature tables ----------
# Build feature matrix
ms1_matrix_all_ENDO_pos <- featureValues(ms1_data_all_ENDO_pos, method="medret", value="into")
colnames(ms1_matrix_all_ENDO_pos) <- MS1_EXO_names
dim(ms1_matrix_all_ENDO_pos)
# transpose feature table
feat_list_all_ENDO_pos <- t(ms1_matrix_all_ENDO_pos)

# Build feature summary
ms1_summary_all_ENDO_pos <- featureSummary(ms1_data_all_ENDO_pos)
ms1_def_all_ENDO_pos <- featureDefinitions(ms1_data_all_ENDO_pos)

# Missing value imputation by filling na positions with median of surrounding features
feat_list_all_ENDO_pos[is.na(feat_list_all_ENDO_pos)] <- median(na.omit(as.numeric(unlist(feat_list_all_ENDO_pos))))

### Transform data
feat_list_all_ENDO_pos <- log2(feat_list_all_ENDO_pos)

# Missing value imputation
feat_list_all_ENDO_pos[which(is.na(feat_list_all_ENDO_pos))] <- median(na.omit(as.numeric(unlist(feat_list_all_ENDO_pos))))

# save as csv
write.csv(feat_list_all_ENDO_pos, file=paste(filename = "Results_combined/feature_list_all_ENDO_pos.csv", sep = ""))

# Plot histogram
#pdf(file="plots_combined/ENDO_feat_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots_combined/all_ENDO_pos_feat_list_hist.jpeg", width = 500, height = 500, quality = 150, bg = "white")
hist(as.numeric(feat_list_all_ENDO_pos), main="Histogram of feature table")
dev.off()

# PCA of feature table results
#pdf(file="plots/ENDO_ms1_feature_table_pca.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots_combined/all_ENDO_pos_ms1_feature_table_pca_exc_MB.jpeg", width = 1000, height = 500, quality = 150, bg = "white")
ms1_pca_all_ENDO_pos <- prcomp(feat_list_all_ENDO_pos[1:24,], center=TRUE)
plot(ms1_pca_all_ENDO_pos$x[, 1], ms1_pca_all_ENDO_pos$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms1_pca_all_ENDO_pos)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_all_ENDO_pos)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color, cex=0.8)
grid()
text(ms1_pca_all_ENDO_pos$x[,1], ms1_pca_all_ENDO_pos$x[,2], labels=ms1_data_all_ENDO_pos$sample_name[1:24], col=color, pos=3, cex=0.9)
legend("topleft", bty="n", pt.cex=1, cex=1, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(ms1_data_all_ENDO_pos@phenoData@data[["sample_group"]]))
dev.off()

# broken stick
jpeg("plots_combined/BrokenStick_all_ENDO_pos_ms1_feature_table_pca_exc_MB.jpeg", width=10, height=6, units="in", res=100)
evplot = function(ev) {  
  # Broken stick model (MacArthur 1957)  
  n = length(ev)  
  bsm = data.frame(j=seq(1:n), p=0)  
  bsm$p[1] = 1/n  
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))  
  bsm$p = 100*bsm$p/n  
  # Plot eigenvalues and % of variation for each axis  
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))  
  barplot(ev, main="Eigenvalues ENDO MS1 feature table", col="blue", las=2)  
  abline(h=mean(ev), col="red")  
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")  
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE,   
          main="% variation", col=c("blue",2), las=2)  
  legend("topright", c("% eigenvalue", "Broken stick model"),   
         pch=15, col=c("blue",2), bty="n")  
  par(op)  
} 

ev_pc = ms1_pca_all_ENDO_pos$sdev^2  
evplot(ev_pc)  
dev.off()

# Create single 0/1 matrix
bina_list_all_ENDO_pos <- t(ms1_matrix_all_ENDO_pos)
bina_list_all_ENDO_pos[is.na(bina_list_all_ENDO_pos)] <- 1
bina_list_all_ENDO_pos <- log2(bina_list_all_ENDO_pos)
bina_list_all_ENDO_pos[bina_list_all_ENDO_pos < log2(ms1_intensity_cutoff)] <- 0
bina_list_all_ENDO_pos[bina_list_all_ENDO_pos != 0] <- 1

# save as csv
write.csv(bina_list_all_ENDO_pos, file=paste(filename = "Results_combined/bina_list_all_ENDO_pos.csv", sep = ""))

# Only unique compounds in group mzml_pheno$ and not the others
uniq_list_all_ENDO_pos <- apply(X=bina_list_all_ENDO_pos, MARGIN=2, FUN=function(x) { if (length(unique(pheno_data_EXO$sample_group[grepl("1", x)])) == 1) x else rep(0, length(x)) } )
colnames(uniq_list_all_ENDO_pos) <- colnames(bina_list_all_ENDO_pos)
rownames(uniq_list_all_ENDO_pos) <- rownames(bina_list_all_ENDO_pos)

## unique list is empty, what happend there?

# Create data frame
model_div_all_ENDO_pos             <- data.frame(features=apply(X=bina_list_all_ENDO_pos, MARGIN=1, FUN=function(x) { sum(x) } ))
model_div_all_ENDO_pos$richness    <- apply(X=bina_list_all_ENDO_pos, MARGIN=1, FUN=function(x) { sum(x) } )
#model_div_all_ENDO_pos$menhinick   <- apply(X=bina_list_all_ENDO_pos, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
model_div_all_ENDO_pos$shannon     <- apply(X=feat_list_all_ENDO_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
model_div_all_ENDO_pos$pielou      <- apply(X=scale(feat_list_all_ENDO_pos, center=F), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#model_div_all_ENDO_pos$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list_all_ENDO_pos, species), index="chao")
model_div_all_ENDO_pos$simpson     <- apply(X=feat_list_all_ENDO_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
model_div_all_ENDO_pos$inverse     <- apply(X=feat_list_all_ENDO_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
model_div_all_ENDO_pos$fisher      <- apply(X=feat_list_all_ENDO_pos, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
model_div_all_ENDO_pos$unique      <- apply(X=uniq_list_all_ENDO_pos, MARGIN=1, FUN=function(x) { sum(x) })

# Remove NAs if present
model_div_all_ENDO_pos[is.na(model_div_all_ENDO_pos)] <- 0

# save as csv
write.csv(model_div_all_ENDO_pos, file=paste(filename = "Results_combined/model_div_all_ENDO_pos.csv", sep = ""))


# save the objects and tables



end.time <- Sys.time()

time.taken <- end.time - start.time
print(time.taken)


