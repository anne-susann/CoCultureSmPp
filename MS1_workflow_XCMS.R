### Analysis of diatom Mono- and Co-culture MS data
# MS1 data pre-processing and statistical analysis with xcms
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

library(parallel)
n.cores <- detectCores()

# Set up parallel processing using 2 cores
if (.Platform$OS.type == "unix") {
  register(bpstart(MulticoreParam(2)))
} else {
  register(bpstart(SnowParam(4)))
}


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
MS1_ENDO_files <- data.frame(list.files(input_dir_MS1, pattern = "ENDO"))
# basenames of files without path and without extension
MS1_ENDO_names <- data.frame(str_remove(MS1_ENDO_files[,], ".mzML"))


###----MS1 EXO files----
#MS1_EXO_files <- data.frame(list.files(input_dir_MS1, pattern = "EXO"))
#MS1_EXO_files_full <- data.frame(list.files(input_dir_MS1, pattern = "EXO", full.names = TRUE))


###----MS2 files----
MS2_files <- data.frame(list.files(input_dir_MS2, pattern = ".mzML"))

# ENDO files
MS2_ENDO_files <- data.frame(subset(MS2_files, grepl("ENDO", MS2_files[,1]), drop = TRUE))
MS2_ENDO_neg_files <- data.frame(subset(MS2_ENDO_files, grepl("neg", MS2_ENDO_files[,1]), drop = TRUE))
MS2_ENDO_pos_files <- data.frame(subset(MS2_ENDO_files, grepl("pos", MS2_ENDO_files[,1]), drop = TRUE))

# EXO files
#MS2_EXO_files <- data.frame(subset(MS2_files, grepl("EXO", MS2_files[,1]), drop = TRUE))
#MS2_EXO_neg_files <- data.frame(subset(MS2_EXO_files, grepl("neg", MS2_EXO_files[,1]), drop = TRUE))
#MS2_EXO_pos_files <- data.frame(subset(MS2_EXO_files, grepl("pos", MS2_EXO_files[,1]), drop = TRUE))

###----MS1 preparations----
# Analysis of MS1 ENDO data only from here on
# MS1 variables
pol <- c(x = polarity, start =0, stop = 0)
ppm <- 35
ms1_intensity_cutoff <- 100	          #approx. 0.01%

mzml_times_ENDO <- NULL

# create vector with sample classes according to culture information sheet
samp_groups <- c("CoCuPp", "CoCuSm", "CoCuSm", "CoCuPp", "CoCuPp", "CoCuSm",
                 rep(x = "Sm", times = 8), 
                 rep(x = "Pp", times = 8),
                 "CoCuSm", "CoCuPp", "MB")

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

# Display MSn levels and check amount of spectra
mzml_msn_ENDO <- NULL
for (i in 1:length(MS1_ENDO_files)) {
  mzml_data_ENDO <- readMSData(files = paste(input_dir_MS1, MS1_ENDO_files[,i], sep = ""), mode="onDisk")
  mzml_msn_ENDO <- rbind(mzml_msn_ENDO, t(as.matrix(table(msLevel(mzml_data_ENDO)))))
}
colnames(mzml_msn_ENDO) <- c("MSn 0", "MSn 1")
rownames(mzml_msn_ENDO) <- MS1_ENDO_names

# Plot MSn levels (only if MS1 and MS2 data are in same directory/file, then exchange MSn 0 for MS1 and MS2)
#pdf(file="plots/ENDOMS1_msn_levels.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(file="plots/ENDOMS1_msn_levels.jpeg", width = 1000, height = 600, quality = 100, bg = "white")
par(mfrow=c(2,1), mar=c(5,4,4,2), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
boxplot(mzml_msn_ENDO, main="Number of spectra")

model_boxplot <- boxplot(t(mzml_msn_ENDO[,2]), main="Number of MS1 spectra per sample", xaxt="n")
tick <- seq_along(model_boxplot$names)
axis(1, at=tick, labels=F)
text(tick, par("usr")[3]-par("usr")[3]/10, model_boxplot$names, adj=0, srt=270, xpd=T)
dev.off()

###----Import raw DATA ENDO NEGATIVE----

# Import raw data as MSnbase object OnDiskMsnExp
msd <- readMSData(files = paste(input_dir_MS1, MS1_ENDO_files[,], sep = ""),
                  pdata = new("NAnnotatedDataFrame",pheno_data_ENDO),
                  mode = "onDisk",
                  centroided = TRUE)

# Import as XCMSnExp object for visual analysis
msd_XCMS <- readMSData(files = paste(input_dir_MS1, MS1_ENDO_files[,], sep = ""), 
                       pdata = new("NAnnotatedDataFrame",pheno_data_ENDO),
                       msLevel = 1,
                       centroided = TRUE)

# inspect data 
table(msLevel(msd))
head(fData(msd)[, c("scanWindowLowerLimit", "scanWindowUpperLimit",
                    "originalPeaksCount", "msLevel", 
                    "polarity", "retentionTime")])

# Restrict data to 1020 seconds (17 minutes)
msd <- filterRt(msd, c(0, 1020))

# subset data for msLevel = 1 and save raw data
msd <- filterMsLevel(msd, msLevel = 1)
table(msLevel(msd))
write.csv(fData(msd), file="ENDO_raw_data.csv", row.names=FALSE)

# Inspect mz values per file
msd_mz <- mz(msd)
msd_mz <- split(msd_mz, f = fromFile(msd))
print(length(msd_mz))

# Get base peak chromatograms
register(bpstart(SnowParam()))

chromas_ENDO <- chromatogram(msd_XCMS, 
                             aggregationFun="max", 
                             msLevel = 1,
                             BPPARAM = SnowParam(workers = 3))

# Plot chromatograms based on phenodata groups
#pdf(file="plots/ENDO_chromas.pdf", encoding="ISOLatin1", pointsize=2, width=6, height=4, family="Helvetica")
jpeg(filename = "plots/ENDO_chromas.jpeg", width = 1000, height = 600, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(chromas_ENDO, main="Raw chromatograms", xlab="retention time [s]", ylab="intensity", col = col)
legend("topleft", bty="n", pt.cex=2, cex=1,5, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(col), legend= unique(msd_XCMS@phenoData@data[["sample_group"]]))
dev.off()

# Get TICs
#pdf(file="plots/ENDO_tics.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots/ENDO_tics.jpeg", width = 1000, height = 600, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(5,4,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=0.4, cex.lab=2, cex.main=2)
tics_ENDO <- split(tic(msd), f=fromFile(msd))
boxplot(tics_ENDO, col=col, ylab="intensity", xlab="sample", main="Total ion current", outline = FALSE)
legend("topleft", bty="n", pt.cex=2, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(col), legend= unique(msd_XCMS@phenoData@data[["sample_group"]]))
dev.off()

### ---- pre-processing ----

# grouping/binning for peak detection, based on similarity of the base peak chromatogram 
chromas_bin_ENDO <- MSnbase::bin(chromas_ENDO, binSize=2)
chromas_bin_cor_ENDO <- cor(log2(do.call(cbind, lapply(chromas_bin_ENDO, intensity)))) # transformation with log
colnames(chromas_bin_cor_ENDO) <- rownames(chromas_bin_cor_ENDO) <- msd$sample_name

# representing the data in a heatmap for general overview
#pdf(file="plots/heatmap_chromas_bin_ENDO.pdf", encoding="ISOLatin1", pointsize=10, width=6, 
#    height=4, family="Helvetica")
jpeg(filename = "plots/heatmap_chromas_bin_ENDO.jpeg", width = 500, height = 500, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
heatmap(chromas_bin_cor_ENDO)
dev.off()

# Assess retention times and intensities of first file
head(rtime(chromas_ENDO[1, 1]))
head(intensity(chromas_ENDO[1, 1]))

# check for polarity
head(fData(msd)[, c("polarity", "filterString", "msLevel", "retentionTime")])
table(polarity(msd))

# parameters here taken from iESTIMATE, need to be calculated with internal standard or IPO
ms1_params_ENDO <- CentWaveParam(ppm=13, mzCenterFun="mean", peakwidth=c(4, 33), prefilter=c(4, 200), 
                                 mzdiff=0.0023, snthresh=11, noise=0, integrate=1,
                                 firstBaselineCheck=TRUE, verboseColumns=TRUE, fitgauss=FALSE, 
                                 roiList=list(), roiScales=numeric())

# Peak detection in MS1 data
ms1_data_ENDO <- findChromPeaks(msd, param=ms1_params_ENDO)

# check the detected peaks
head(chromPeaks(ms1_data_ENDO))
chromPeakData(ms1_data_ENDO)


# Per file summary
ms1_summary_ENDO <- lapply(split.data.frame(chromPeaks(ms1_data_ENDO), 
                                            f=chromPeaks(ms1_data_ENDO)[, "sample"]), 
                           FUN=function(z) { c(peak_count=nrow(z), rt=quantile(z[, "rtmax"] - z[, "rtmin"])) } )
ms1_summary_ENDO <- do.call(rbind, ms1_summary_ENDO)
rownames(ms1_summary_ENDO) <- basename(fileNames(ms1_data_ENDO))
print(ms1_summary_ENDO)
table(msLevel(ms1_data_ENDO))

write.csv(as.data.frame(table(msLevel(ms1_data_ENDO))), file="ENDO_ms1_data.csv", row.names=FALSE)


# To get a global overview of the peak detection we can plot the frequency of identified peaks per file along the retention time axis. 
# This allows to identify time periods along the MS run in which a higher number of peaks was identified and evaluate whether this is consistent across files.
#pdf(file="plots/ENDO_ms1_data.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots/ENDO_ms1_data.jpeg", width = 1000, height = 600, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(4,18,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plotChromPeakImage(ms1_data_ENDO, main="Frequency of identified peaks per RT", binSize = 20)
dev.off()


## Group peaks
ms1_data_ENDO <- groupChromPeaks(ms1_data_ENDO, param=PeakDensityParam(
  sampleGroups=ms1_data_ENDO$sample_group, minFraction=0.7, bw=2.5))

## RT correction
ms1_data_ENDO <- adjustRtime(ms1_data_ENDO, param=PeakGroupsParam(
  minFraction=0.7,smooth="loess",span=0.8,family="gaussian"))

# Plot the difference of raw and adjusted retention times
#pdf(file="plots/ENDO_ms1_raw_adjusted.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
jpeg(filename = "plots/ENDO_ms1_raw_adjusted.jpeg", width = 500, height = 1000, quality = 100, bg = "white")
par(mfrow=c(2,1), mar=c(4.5,4.2,4,1), cex=0.8)
plot(chromas_ENDO, peakType="none", main="Raw chromatograms")
plotAdjustedRtime(ms1_data_ENDO, lwd=2, main="Retention Time correction")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
dev.off()

## Group peaks
ms1_data_ENDO <- groupChromPeaks(ms1_data_ENDO, param=PeakDensityParam(
  sampleGroups=ms1_data_ENDO$sample_group, minFraction=0.7, bw=2.5))

# Get integrated peak intensity per feature/sample
print(head(featureValues(ms1_data_ENDO, value="into")))

## Fill peaks
# missing value imputation, see xcmsSet
ms1_data_ENDO <- fillChromPeaks(ms1_data_ENDO, param=FillChromPeaksParam(ppm=ppm, fixedRt=0, expandRt=5))
head(featureValues(ms1_data_ENDO))
head(featureSummary(ms1_data_ENDO, group=ms1_data_ENDO$sample_group))

# Evaluate grouping
#pdf(file="plots/ENDO_ms1_grouping.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots/ENDO_ms1_grouping.jpeg", width = 1000, height = 500, quality = 100, bg = "white")
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

# broken stick
png("BrokenStick_ENDO_ms1_grouping.png", width=10, height=6, units="in", res=100)
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

ev_pc = ms1_pca_ENDO$sdev^2  
evplot(ev_pc)  
dev.off()


# Show peaks
tail(chromPeaks(ms1_data_ENDO))
tail(chromPeakData(ms1_data_ENDO))

# Show process history
processHistory(ms1_data_ENDO)


# ---------- Build MS1 feature tables ----------
# Build feature matrix
ms1_matrix_ENDO <- featureValues(ms1_data_ENDO, method="medret", value="into")
colnames(ms1_matrix_ENDO) <- MS1_ENDO_names[,1]
dim(ms1_matrix_ENDO)
# transpose feature table
feat_list_ENDO <- t(ms1_matrix_ENDO)

# Build feature summary
ms1_summary_ENDO <- featureSummary(ms1_data_ENDO)
ms1_def_ENDO <- featureDefinitions(ms1_data_ENDO)

# Missing value imputation by filling na positions with median of surrounding features
feat_list_ENDO[is.na(feat_list_ENDO)] <- median(na.omit(as.numeric(unlist(feat_list_ENDO))))

### Transform data
feat_list_ENDO <- log2(feat_list_ENDO)

# Missing value imputation
feat_list_ENDO[which(is.na(feat_list_ENDO))] <- median(na.omit(as.numeric(unlist(feat_list_ENDO))))

# Plot histogram
#pdf(file="plots/ENDO_feat_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots/ENDO_feat_list_hist.jpeg", width = 500, height = 500, quality = 100, bg = "white")
hist(as.numeric(feat_list_ENDO), main="Histogram of feature table")
dev.off()

# PCA of feature table results
#pdf(file="plots/ENDO_ms1_feature_table_pca.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots/ENDO_ms1_feature_table_pca.jpeg", width = 1000, height = 500, quality = 100, bg = "white")
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
jpeg("plots/BrokenStick_ENDO_ms1_feature_table_pca.jpeg", width=10, height=6, units="in", res=100)
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







