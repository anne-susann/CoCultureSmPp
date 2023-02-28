## Co-Culture Analysis: link MS2 data to MS1 data 
## add to inclusion list

## -------library-----

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
library(IPO)

###------parallelization----

library(parallel)
n.cores <- detectCores()

# Set up parallel processing using 2 cores
if (.Platform$OS.type == "unix") {
  register(bpstart(MulticoreParam(2)))
} else {
  register(bpstart(SnowParam(4)))
}




########## set directory and list files #######

# set data directory for MS1 data files
input_dir_MS1 <- paste(getwd(), "/MS1/", sep = "")
input_dir_MS1

# set data directory for MS2 data files
input_dir_MS2 <- paste(getwd(), "/MS2/", sep = "")
input_dir_MS2

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
  print("plots directory is created")
}

###----ENDO files----
MS1_ENDO_files <- data.frame(list.files(input_dir_MS1, pattern = "ENDO"))
MS1_ENDO_files_full <- data.frame(list.files(input_dir_MS1, pattern = "ENDO", full.names = TRUE))

# basenames of files without path and without extension
MS1_ENDO_names <- data.frame(str_remove(MS1_ENDO_files[,], ".mzML"))

###----EXO files----
MS1_EXO_files <- data.frame(list.files(input_dir_MS1, pattern = "EXO"))
MS1_EXO_files_full <- data.frame(list.files(input_dir_MS1, pattern = "EXO", full.names = TRUE))


###----MS2 files----
MS2_files <- data.frame(list.files(input_dir_MS2, pattern = ".mzML"))

# ENDO files
MS2_ENDO_files <- data.frame(subset(MS2_files, grepl("ENDO", MS2_files[,1]), drop = TRUE))
MS2_ENDO_neg_files <- data.frame(subset(MS2_ENDO_files, grepl("neg", MS2_ENDO_files[,1]), drop = TRUE))
MS2_ENDO_pos_files <- data.frame(subset(MS2_ENDO_files, grepl("pos", MS2_ENDO_files[,1]), drop = TRUE))

# EXO files
MS2_EXO_files <- data.frame(subset(MS2_files, grepl("EXO", MS2_files[,1]), drop = TRUE))
MS2_EXO_neg_files <- data.frame(subset(MS2_EXO_files, grepl("neg", MS2_EXO_files[,1]), drop = TRUE))
MS2_EXO_pos_files <- data.frame(subset(MS2_EXO_files, grepl("pos", MS2_EXO_files[,1]), drop = TRUE))


###----MS1 preparations----
# MS1 variables
pol <- c(x = polarity, start =0, stop = 0)
ppm <- 35
ms1_intensity_cutoff <- 100	          #approx. 0.01%

# General variables
# mzml_files_ENDO <- NULL
# mzml_names_neg <- NULL
mzml_times_ENDO <- NULL




###-----data-------
# iESTIMATE code as inspiration (IPB Halle, github)
# https://github.com/ipb-halle/iESTIMATE/blob/main/use-cases/radula-hormones/peak_detection_neg.r

# Analysis of MS1 ENDO data only from here on

# create vector with sample classes according to culture information sheet
samp_groups <- c("CoCuPp", "CoCuSm", "CoCuSm", "CoCuPp", "CoCuPp", "CoCuSm",
                   rep(x = "Sm", times = 8), 
                   rep(x = "Pp", times = 8),
                   "CoCuSm", "CoCuPp", "MB", )

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

col <- c(CoCuPp1, CoCuSm1, CoCuPp2, CoCuSm2, Sm, Pp, CoCuSm3, CoCuPp3, MB)

# create phenodata based on culture type
pheno_data_ENDO <- data.frame(sample_name = MS1_ENDO_names, sample_group = samp_groups)
pheno_col_ENDO <- data.frame(col)
pheno_col_samples_ENDO <- data.frame(cbind(pheno_col_ENDO$col, pheno_data_ENDO$sample_group))

# Save timestamps of samples
# doesn't run
for (i in 1:length(MS1_ENDO_files_full)) {
      fl <- mzR::openMSfile(MS1_ENDO_files_full[i])
      run_info <- mzR::runInfo(fl)
      mzR::close(fl)
      mzml_times_ENDO <- c(mzml_times_ENDO, run_info$startTimeStamp)
}


# Display MSn levels and check amount of spectra
mzml_msn_ENDO <- NULL
for (i in 1:length(MS1_ENDO_files)) {
  mzml_data_ENDO <- readMSData(files = paste(input_dir_MS1, MS1_ENDO_files[,i], sep = ""), mode="onDisk")
  mzml_msn_ENDO <- rbind(mzml_msn_ENDO, t(as.matrix(table(msLevel(mzml_data_ENDO)))))
}
colnames(mzml_msn_ENDO) <- c("MSn 0", "MSn 1")
rownames(mzml_msn_ENDO) <- MS1_ENDO_names

# Plot MSn levels (only if MS1 and MS2 data are in same directory/file, then exchange MSn 0 for MS1 and MS2)
pdf(file="plots/ENDOMS1_msn_levels.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=16, family="Helvetica")
par(mfrow=c(2,1), mar=c(5,4,4,2), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
boxplot(mzml_msn_ENDO, main="Number of spectra")

model_boxplot <- boxplot(t(mzml_msn_ENDO[,2]), main="Number of MS1 spectra per sample", xaxt="n")
tick <- seq_along(model_boxplot$names)
axis(1, at=tick, labels=F)
text(tick, par("usr")[3]-par("usr")[3]/10, model_boxplot$names, adj=0, srt=270, xpd=T)
dev.off()



###----Import raw DATA ENDO NEGATIVE----

# Import raw data as MSnbase object OnDiskMsnExp, segmented for msLevel = 1
msd <- readMSData(files = paste(input_dir_MS1, MS1_ENDO_files[,], sep = ""),
                  pdata = new("NAnnotatedDataFrame",pheno_data_ENDO), 
                  msLevel = 1,
                  mode = "onDisk")

# Import raw data for MS2 in an object, subset for msLevel = 2
#msd_MS2 <- readMSData(files = paste(input_dir_MS2, MS2_ENDO_neg_files[,], sep = ""),
#                      pdata = new("NAnnotatedDataFrame",pheno_data_ENDO), 
#                      msLevel = 2,
#                      mode = "onDisk")
# msd <- msd_MS1


# Import as XCMSnExp object for visual analysis
msd_XCMS <- readMSData(files = paste(input_dir_MS1, MS1_ENDO_files[,], sep = ""), 
                  pdata = new("NAnnotatedDataFrame",pheno_data_ENDO),
                  msLevel = 1)

table(msLevel(msd))
head(fData(msd)[, c("scanWindowLowerLimit", "scanWindowUpperLimit",
                             "originalPeaksCount", "msLevel", 
                    "polarity", "retentionTime")])
write.csv(fData(msd), file="ENDO_raw_data_1.csv", row.names=FALSE)

# Restrict data to 1020 seconds (17 minutes)
msd <- filterRt(msd, c(0, 1020))

# filter for NEGATIVE polarity
msd <- filterPolarity(msd, polarity = 0)

# Inspect mz values per file
msd_mz <- mz(msd)
msd_mz <- split(msd_mz, f = fromFile(msd))
print(length(msd_mz))

# Get base peak chromatograms
register(bpstart(SnowParam()))

chromas_ENDO <- chromatogram(msd_XCMS, 
                             aggregationFun="max", 
                             BPPARAM = SnowParam(workers = 2))


# Plot chromatograms based on phenodata groups
pdf(file="plots/ENDO_chromas.pdf", encoding="ISOLatin1", pointsize=2, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(chromas_ENDO, main="Raw chromatograms", xlab="retention time [s]", ylab="intensity", col = col)
legend("topleft", bty="n", pt.cex=2, cex=1,5, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(col), legend= unique(msd_XCMS@phenoData@data[["sample_group"]]))
dev.off()

# Get TICs
pdf(file="plots/ENDO_tics.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "ENDO_tics.jpeg", width = 1000, height = 600, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(5,4,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=0.4, cex.lab=2, cex.main=2)
tics_ENDO <- split(tic(msd), f=fromFile(msd))
boxplot(tics_ENDO, col=col, ylab="intensity", xlab="sample", main="Total ion current", outline = FALSE)
legend("topleft", bty="n", pt.cex=2, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(col), legend= unique(msd_XCMS@phenoData@data[["sample_group"]]))
dev.off()


### ---- pre-processing ----

# grouping/binning for peak detection, based on similarity of the base peak chromatogram 
chromas_bin_ENDO <- MSnbase::bin(chromas_ENDO, binSize=2)
chromas_bin_cor_ENDO <- cor(log2(do.call(cbind, lapply(chromas_bin_ENDO, intensity))))
colnames(chromas_bin_cor_ENDO) <- rownames(chromas_bin_cor_ENDO) <- msd$sample_name

# representing the data in a heatmap for general overview
pdf(file="plots/heatmap_chromas_bin_ENDO.pdf", encoding="ISOLatin1", pointsize=10, width=6, 
    height=4, family="Helvetica")
jpeg(filename = "plots/heatmap_chromas_bin_ENDO.jpeg", width = 500, height = 500, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
heatmap(chromas_bin_cor_ENDO)
dev.off()

# Assess retention times and intensities of first file
head(rtime(chromas_ENDO[1, 1]))
head(intensity(chromas_ENDO[1, 1]))

# Inspect peak width of standard compound for defining base peakwidth parameter below
# standards stated here are from iESTIMATE project
# check for own standard compounds that are present in dataset 
# to set ppm parameter for subsequent peak detection parameters check maximal mz difference of data points from neighboring scans/spectra
register(bpstart(SnowParam()))
pdf(file="plots/ENDO_chromas_standard_peakwidth_XCMS.pdf", encoding="ISOLatin1", 
    pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(chromatogram(msd_XCMS, mz=c(282, 285), rt=c(545, 555)), col=col, main = "EIC of Biochanin A")
plot(chromatogram(msd, mz=c(244, 247), rt=c(335, 350)), col=col, main = "EIC of N-(3-Indolylacetyl)-L-Ala")
plot(chromatogram(msd, mz=c(294, 297), rt=c(585, 600)), col=col, main = "EIC of Radulanin L")
plot(chromatogram(msd, mz=c(147, 151), rt=c(324, 340)), col=col, main = "EIC of Methionin", BPPARAM = SnowParam(workers = 2))
dev.off()
# first run worked with msd, now neither msd nor msd_XCMS
# with BPPARAM it works, not very well but works with both objects

### needs internal standard data or known compound
rtr <- c(2700, 2900)
mzr <- c(334.9, 335.1)

chr_raw <- chromatogram(msd, mz = mzr, rt = rtr)
plot(chr_raw, col = group_colors[chr_raw$sample_group])


msd %>%
  filterRt(rt = rtr) %>%
  filterMz(mz = mzr) %>%
  plot(type = "XIC")
###

# check for polarity
head(fData(msd)[, c("polarity", "filterString", "msLevel", "retentionTime")])
table(polarity(msd))

# polarity = 0 is negative mode, but distinction here not possible in this way 
# for iESTIMATE, subset data only contains negative data, here neg and pos alternating in data
#if (polarity(msd)==0) {
#  ms1_params_ENDO <- CentWaveParam(ppm=13, mzCenterFun="mean", peakwidth=c(4, 33), prefilter=c(4, 200), 
#                                  mzdiff=0.0023, snthresh=11, noise=0, integrate=1,
#                                  firstBaselineCheck=TRUE, verboseColumns=TRUE, fitgauss=FALSE, 
#                                  roiList=list(), roiScales=numeric())
#} else {
#  ms1_params_ENDO <- CentWaveParam(ppm=32, mzCenterFun="mean", peakwidth=c(4, 32), prefilter=c(2, 100), 
#                                  mzdiff=0.0111, snthresh=8, noise=0, integrate=1,
#                                  firstBaselineCheck=TRUE, verboseColumns=FALSE, fitgauss=FALSE, 
#                                  roiList=list(), roiScales=numeric())
#}


###---- xcms parameter optimization with IPO ----
# test on raw data
datafiles <- list.files(input_dir_MS1, recursive = TRUE, full.names = TRUE)


peakpickingParameters <- getDefaultXcmsSetStartingParams(method = "centWave")
#setting levels for step to 0.2 and 0.3 (hence 0.25 is the center point)
peakpickingParameters$step <- c(0.2, 0.3)
peakpickingParameters$fwhm <- c(40, 50)
peakpickingParameters$quantitative_parameters <- c(2,4)
peakpickingParameters$prefilter <- c(3, 4)
peakpickingParameters$prefilter_value <- c(2, 4)
#setting only one value for steps therefore this parameter is not optimized
peakpickingParameters$steps <- 2
peakpickingParameters$qualitative_parameters <-5
peakpickingParameters$sigma <- 3
peakpickingParameters$max <- 2
peakpickingParameters$index <- FALSE

time.xcmsSet <- system.time({ # measuring time
  resultPeakpicking <- 
    optimizeXcmsSet(files = datafiles[1:2], 
                    params = peakpickingParameters, 
                    BPPARAM = SnowParam(workers = 2),
                    nSlaves = 1, 
                    subdir = NULL,
                    plot = TRUE)
})






###---- resume with iESTIMATE parameters, MS1 data in object ----
# parameters here taken from iESTIMATE, need to be calculated with internal standard
ms1_params_ENDO <- CentWaveParam(ppm=13, mzCenterFun="mean", peakwidth=c(4, 33), prefilter=c(4, 200), 
                                  mzdiff=0.0023, snthresh=11, noise=0, integrate=1,
                                  firstBaselineCheck=TRUE, verboseColumns=TRUE, fitgauss=FALSE, 
                                  roiList=list(), roiScales=numeric())


# Peak detection in MS1 data
ms1_data_ENDO <- findChromPeaks(msd, param=ms1_params_ENDO)
ms1_data_ENDO_default <- findChromPeaks(msd, param = CentWaveParam(snthresh = 2))

# check the detected peaks
head(chromPeaks(ms1_data_ENDO))
chromPeakData(ms1_data_ENDO)


# Per file summary
ms1_summary_ENDO <- lapply(split.data.frame(chromPeaks(ms1_data_ENDO), f=chromPeaks(ms1_data_ENDO)[, "sample"]), FUN=function(z) { c(peak_count=nrow(z), rt=quantile(z[, "rtmax"] - z[, "rtmin"])) } )
ms1_summary_ENDO <- do.call(rbind, ms1_summary_ENDO)
rownames(ms1_summary_ENDO) <- basename(fileNames(ms1_data_ENDO))
print(ms1_summary_ENDO)
table(msLevel(ms1_data_ENDO))

write.csv(as.data.frame(table(msLevel(ms1_data_ENDO))), file="ENDO_ms1_data.csv", row.names=FALSE)



### documented in Master thesis document until here



# To get a global overview of the peak detection we can plot the frequency of identified peaks per file along the retention time axis. 
# This allows to identify time periods along the MS run in which a higher number of peaks was identified and evaluate whether this is consistent across files.
pdf(file="plots/ENDO_ms1_data.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,18,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plotChromPeakImage(ms1_data_ENDO, main="Frequency of identified peaks per RT", binSize = 20)
dev.off()


## Group peaks
# for now use parameters for negative mode

#if (polarity=="negative") {
  ms1_data_ENDO <- groupChromPeaks(ms1_data_ENDO, param=PeakDensityParam(
    sampleGroups=ms1_data_ENDO$sample_group, minFraction=0.7, bw=2.5))
#} else {
#  ms1_data_ENDO <- groupChromPeaks(ms1_data_ENDO, param=PeakDensityParam(sampleGroups=ms1_data_ENDO$sample_group, minFraction=0.7, bw=22))
#}

  
## RT correction
# span parameter in iESTIMATE = 0.2, when applying span = 1 better for our data

#if (polarity=="negative") {
  ms1_data_ENDO <- adjustRtime(ms1_data_ENDO, param=PeakGroupsParam(
    minFraction=0.7,smooth="loess",span=1,family="gaussian"))
#} else {
#  ms1_data_ENDO <- adjustRtime(ms1_data_ENDO, param=PeakGroupsParam(minFraction=0.7,smooth="loess",span=0.2,family="gaussian"))
#}

# Plot the difference of raw and adjusted retention times
pdf(file="plots/ENDO_ms1_raw_adjusted.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
jpeg(filename = "plots/ENDO_ms1_raw_adjusted.jpeg", width = 500, height = 1000, quality = 100, bg = "white")
par(mfrow=c(2,1), mar=c(4.5,4.2,4,1), cex=0.8)
plot(chromas_ENDO, peakType="none", main="Raw chromatograms")
plotAdjustedRtime(ms1_data_ENDO, lwd=2, main="Retention Time correction")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
dev.off()


## Group peaks
# like in xcmsSet, new grouping after retention time correction

#if (polarity=="negative") {
  ms1_data_ENDO <- groupChromPeaks(ms1_data_ENDO, param=PeakDensityParam(
    sampleGroups=ms1_data_ENDO$sample_group, minFraction=0.7, bw=2.5))
#} else {
#  ms1_data_ENDO <- groupChromPeaks(ms1_data_ENDO, param=PeakDensityParam(sampleGroups=ms1_data_ENDO$sample_group, minFraction=0.7, bw=22))
#}

# Get integrated peak intensity per feature/sample
print(head(featureValues(ms1_data_ENDO, value="into")))

## Fill peaks
# missing value imputation, see xcmsSet
ms1_data_ENDO <- fillChromPeaks(ms1_data_ENDO, param=FillChromPeaksParam(ppm=ppm, fixedRt=0, expandRt=5))
head(featureValues(ms1_data_ENDO))
head(featureSummary(ms1_data_ENDO, group=ms1_data_ENDO$sample_group))

# Evaluate grouping
pdf(file="plots/ENDO_ms1_grouping.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
ms1_pca_ENDO <- prcomp(t(na.omit(log2(featureValues(ms1_data_ENDO, value="into")))), center=TRUE)
plot(ms1_pca_ENDO$x[, 1], ms1_pca_ENDO$x[,2], pch=19, main="PCA: Grouping of samples",
     xlab=paste0("PC1: ", format(summary(ms1_pca_ENDO)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_ENDO)$importance[2, 2] * 100, digits=3), " % variance"),
     col=col, cex=0.8)
grid()
text(ms1_pca_ENDO$x[, 1], ms1_pca_ENDO$x[,2], labels=ms1_data_ENDO$sample_name, col=col, pos=3, cex=0.5)
legend("topleft", bty="n", pt.cex=1, cex=0.8, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(col), legend= unique(ms1_data_ENDO@phenoData@data[["sample_group"]]))
dev.off()

# Show peaks
tail(chromPeaks(ms1_data_ENDO))
tail(chromPeakData(ms1_data_ENDO))

# Show process history
processHistory(ms1_data_ENDO)


### makes sense until here

# ---------- Build MS1 feature tables ----------
# Build feature matrix
ms1_matrix_ENDO <- featureValues(ms1_data_ENDO, method="medret", value="into")
colnames(ms1_matrix_ENDO) <- MS1_ENDO_names[,1]
dim(ms1_matrix_ENDO)
feat_list_ENDO <- t(ms1_matrix_ENDO)

# Build feature summary
ms1_summary_ENDO <- featureSummary(ms1_data_ENDO)
ms1_def_ENDO <- featureDefinitions(ms1_data_ENDO)

# Missing value imputation
feat_list_ENDO[is.na(feat_list_ENDO)] <- median(na.omit(as.numeric(unlist(feat_list_ENDO))))

# Transform data
feat_list_ENDO <- log2(feat_list_ENDO)

# Missing value imputation by filling na positions with median of surrounding features
feat_list_ENDO[which(is.na(feat_list_ENDO))] <- median(na.omit(as.numeric(unlist(feat_list_ENDO))))

# Plot histogram
pdf(file="plots/ENDO_feat_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(feat_list_ENDO), main="Histogram of feature table")
dev.off()


# PCA of feature table results
pdf(file="plots/ENDO_ms1_feature_table_pca.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
ms1_pca_ENDO <- prcomp(feat_list_ENDO, center=TRUE)
plot(ms1_pca_ENDO$x[, 1], ms1_pca_ENDO$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms1_pca_ENDO)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_ENDO)$importance[2, 2] * 100, digits=3), " % variance"),
     col=col, cex=0.8)
grid()
text(ms1_pca_ENDO$x[, 1], ms1_pca_ENDO$x[,2], labels=ms1_data_ENDO$sample_name, col=col, pos=3, cex=0.5)
legend("topleft", bty="n", pt.cex=1, cex=0.8, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(col), legend= unique(ms1_data_ENDO@phenoData@data[["sample_group"]]))
dev.off()

# broken stick
png("BrokenStick_ENDO.png", width=10, height=6, units="in", res=100)
evplot = function(ev) {  
  # Broken stick model (MacArthur 1957)  
  n = length(ev)  
  bsm = data.frame(j=seq(1:n), p=0)  
  bsm$p[1] = 1/n  
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))  
  bsm$p = 100*bsm$p/n  
  # Plot eigenvalues and % of variation for each axis  
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))  
  barplot(ev, main="Eigenvalues ENDO", col="blue", las=2)  
  abline(h=mean(ev), col="red")  
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")  
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE,   
          main="% variation", col=c("blue",2), las=2)  
  legend("topright", c("% eigenvalue", "Broken stick model"),   
         pch=15, col=c("blue",2), bty="n")  
  par(op)  
} 

ev_pc = ms1_pca_ENDO$sdev^2  
evplot(ev_pc)  
dev.off()


# Create single 0/1 matrix
bina_list_ENDO <- t(ms1_matrix_ENDO)
bina_list_ENDO[is.na(bina_list_ENDO)] <- 1
bina_list_ENDO <- log2(bina_list_ENDO)
bina_list_ENDO[bina_list_ENDO < log2(ms1_intensity_cutoff)] <- 0
bina_list_ENDO[bina_list_ENDO != 0] <- 1

# Only unique compounds in group mzml_pheno$ and not the others
uniq_list_ENDO <- apply(X=bina_list_ENDO, MARGIN=2, FUN=function(x) { if (length(unique(pheno_data_ENDO$sample_group[grepl("1", x)])) == 1) x else rep(0, length(x)) } )
colnames(uniq_list_ENDO) <- colnames(bina_list_ENDO)
rownames(uniq_list_ENDO) <- rownames(bina_list_ENDO)

## unique list is empty, what happend there?

# Create data frame
model_div_ENDO             <- data.frame(features=apply(X=bina_list_ENDO, MARGIN=1, FUN=function(x) { sum(x) } ))
model_div_ENDO$richness    <- apply(X=bina_list_ENDO, MARGIN=1, FUN=function(x) { sum(x) } )
#model_div_ENDO$menhinick   <- apply(X=bina_list_ENDO, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
model_div_ENDO$shannon     <- apply(X=feat_list_ENDO, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
model_div_ENDO$pielou      <- apply(X=scale(feat_list_ENDO, center=F), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#model_div_ENDO$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list_ENDO, species), index="chao")
model_div_ENDO$simpson     <- apply(X=feat_list_ENDO, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
model_div_ENDO$inverse     <- apply(X=feat_list_ENDO, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
model_div_ENDO$fisher      <- apply(X=feat_list_ENDO, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
model_div_ENDO$unique      <- apply(X=uniq_list_ENDO, MARGIN=1, FUN=function(x) { sum(x) })

# Remove NAs if present
model_div_ENDO[is.na(model_div_ENDO)] <- 0




# ############################## MS2 ##############################

### ---- ENDO -----

inclusion_list <- read.table(paste(input_dir_MS2, str_remove("/KSS_210324_ENDO.txt", "."), sep = ""), 
                             header = TRUE, sep = "\t", dec = ".", fill = TRUE)
inc_list_ENDO <- data.frame(inclusion_list$Mass..m.z., inclusion_list$CS..z., 
                       inclusion_list$Polarity, inclusion_list$Start..min., 
                       inclusion_list$End..min.)
colnames(inc_list_ENDO) <- c("mz", "CS", "polarity", "start_sec", "end_sec")

# convert retention time from minutes to seconds
inc_list_ENDO$start_sec <- inc_list_ENDO$start_sec*60
inc_list_ENDO$end_sec <- inc_list_ENDO$end_sec*60

# ---------- MS2 spectra detection ----------
# ms1_data_ENDO is feature table


# first OnDiskMsnExp object, like with MS1 data subset for msLevel = 2
mzml_data_ENDO_MS2 <- readMSData(files = paste(input_dir_MS2, MS2_ENDO_files[,], sep = ""), 
                  pdata = NULL, 
                  mode = "onDisk")

msd_MS2 <- readMSData(files = paste(input_dir_MS2, MS2_ENDO_files[,], sep = ""), 
                      pdata = NULL, 
                      msLevel = 2,
                      mode = "onDisk")

# check for msLevel
table(msLevel(mzml_data_ENDO_MS2))
mzml_data_ENDO_MS2 <- filterMsLevel(mzml_data_ENDO_MS2, msLevel = "1")
table(msLevel(msd_MS2))
head(fData(msd_MS2)[, c("precursorMZ", "precursorCharge", "precursorIntensity",
                         "precursorScanNum", "originalPeaksCount",
                        "scanWindowLowerLimit", "scanWindowUpperLimit", "msLevel", "retentionTime")])

# save fData of MS2 spectra
write.csv(fData(msd_MS2), file="ENDO_raw_data_MS2.csv", row.names=FALSE)


### extract precursorMz with forward pipe through msLevel
# from https://www.bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms-lcms-ms.html
library(magrittr)

# precursorMz
mzml_data_ENDO_MS2 %>%
  filterMsLevel(2L) %>%
  precursorMz() %>%
  head()

# precursorIntensity
mzml_data_ENDO_MS2 %>%
  filterMsLevel(2L) %>%
  precursorIntensity() %>%
  head()

cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                       peakwidth = c(3, 30))
mzml_data_ENDO_MS2 <- findChromPeaks(mzml_data_ENDO_MS2, param = cwp)
nrow(chromPeaks(mzml_data_ENDO_MS2))

library(Spectra)
mzml_data_ENDO_MS2_spectra <- chromPeakSpectra(
  mzml_data_ENDO_MS2, msLevel = 2L, return.type = "Spectra")
mzml_data_ENDO_MS2_spectra$peak_id


ex_mz <- 644
chromPeaks(mzml_data_ENDO_MS2, mz = ex_mz, ppm = 20)

###

library(foreach)

ms2_data_ENDO <- chromPeakSpectra(mzml_data_ENDO_MS2, msLevel=2L, return.type="Spectra")

# Extract all usable MS2 spectra
ms2_spectra_ENDO <- list()
#for (i in 1:nrow(ms1_def_ENDO)) {
ms2_spectra_ENDO <- foreach(i=1:nrow(ms1_def_ENDO)) %dopar% {
  # Extract existing MS2 spectra for feature
  feature_of_interest <- ms1_def_ENDO[i, "mzmed"]
  peaks_of_interest <- chromPeaks(ms1_data_ENDO, mz=feature_of_interest, ppm=ppm)
  if (length(which(ms2_data_ENDO$peak_id %in% rownames(peaks_of_interest))) > 0) {
    spectra_of_interest <- ms2_data_ENDO[ms2_data_ENDO$peak_id %in% rownames(peaks_of_interest)]
    combined_spectra_of_interest <- filterIntensity(spectra_of_interest, intensity=c(1,Inf), backend=MsBackendDataFrame)
    combined_spectra_of_interest <- setBackend(combined_spectra_of_interest, backend=MsBackendDataFrame())
    combined_spectra_of_interest <- Spectra::combineSpectra(combined_spectra_of_interest, 
                                                            FUN=combinePeaks, ppm=ppm, peaks="union", minProp=0.8, intensityFun=median, 
                                                            mzFun=median, backend=MsBackendDataFrame)#f=rownames(peaks_of_interest))
    combined_spectra_peaks <- as.data.frame(peaksData(combined_spectra_of_interest)[[1]])
    #Spectra::plotSpectra(combined_spectra_of_interest)
    #plot(x=combined_spectra_peaks[,1], y=combined_spectra_peaks[,2], type="h", xlab="m/z", ylab="intensity", main=paste("Precursor m/z",combined_spectra_of_interest@backend@spectraData$precursorMz[[1]]))
    #length(spectra_of_interest$peak_id)
    
    #ms2_spectra_ENDO[[i]] <- combined_spectra_of_interest
    return(combined_spectra_of_interest)
  }
}
###

# Remove empty spectra
names(ms2_spectra_ENDO) <- rownames(ms1_def_ENDO)
ms2_spectra_ENDO <- ms2_spectra_ENDO[lengths(ms2_spectra_ENDO) != 0]


###----Link MS2 and MS1-----

# Relate all MS2 spectra to MS1 precursors
ms1_def_ENDO$has_ms2 <- as.integer(rownames(ms1_def_ENDO) %in% names(ms2_spectra_ENDO))


# Take data from peaks table for MS1, or potentially data from ENDO_raw_data.csv/table is enough
# for MS2 data, ENDO_raw_data_MS2_2.csv/table generated from msd_MS2 object fData is the necessary information
# contain the precursorMz info and sample/scan/file names --> link should be possible







###----with MS-Dial----
# C:\Users\abela\Downloads\MSDIAL ver.4.9.221218 Windowsx64\MSDIAL ver.4.9.221218 Windowsx64
# MS-Dial