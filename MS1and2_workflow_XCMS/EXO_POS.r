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

start.time <- Sys.time()

########## set directory and list files ########

# set data directory for MS1 data files
input_dir_MS1 <- "/Users/mahnoorzulfiqar/OneDriveUNI/CoCulture_mzML2/MS1_pos_neg"

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
if (dir.exists(paste(getwd(), "/plots/", sep = ""))){
  print("plots directory already exists")
  start_time <- Sys.time()
}  else{
  dir.create("plots")
  start_time <- Sys.time()
  print("plots folder has been created")
}

raw_data_MS1_EXO_pos <- list.files(input_dir_MS1, pattern = "pos_EXO")
raw_data_MS1_EXO_pos

# type(paste(input_dir_MS1, "/", raw_data_MS1_EXO, sep = ""))



# sps_all <- Spectra(paste(input_dir_MS1, "/", raw_data_MS1_EXO_pos, sep = ""), source = MsBackendMzR())
# sps_all

raw_data_MS1_EXO_pos <- list.files(input_dir_MS1, pattern = "pos_EXO")
raw_data_MS1_EXO_pos_files <- paste(input_dir_MS1, "/", raw_data_MS1_EXO_pos, sep ="")
raw_data_MS1_EXO_pos_files

raw_data_MS2_EXO_pos <- list.files(input_dir_MS2, pattern = "EXOpos")
raw_data_MS2_EXO_pos_files <- paste(input_dir_MS2, raw_data_MS2_EXO_pos, sep ="")
raw_data_MS2_EXO_pos_files

all_files <- c(raw_data_MS1_EXO_pos_files, raw_data_MS2_EXO_pos_files)
length(all_files)

36-11

all_files_names <- str_remove_all(all_files, ".mzML")
all_files_names <- str_remove_all(all_files_names,"/Users/mahnoorzulfiqar/OneDriveUNI/CoCulture2/MS1_pos_neg/")

all_files_names[26:36]<- str_remove_all(all_files_names[26:36], '/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/CoCulture2/MS2/')

all_files_names

# create vector with sample classes according to culture information sheet
samp_groups <- c("CoCuPp", "CoCuSm", "CoCuSm", "CoCuPp", "CoCuPp", "CoCuSm",
                 rep(x = "Sm", times = 8), 
                 rep(x = "Pp", times = 8),
                 "CoCuSm", "CoCuPp", "MB", rep(x = "MS2", times = 11))

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
                 CoCuSm1, CoCuPp1, MB1, rep(x = ms2, times = 11))

samp_groups_description

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
ms2<- rep("aquamarine", 11)
color <- c(CoCuPp1, CoCuSm1, CoCuPp2, CoCuSm2, Sm, Pp, CoCuSm3, CoCuPp3, MB, ms2)



all_files



# create phenodata based on culture type
pheno_data_EXO <- data.frame(sample_name = all_files_names, sample_group = samp_groups, samp_groups_description=samp_groups_description)
pheno_col_EXO <- data.frame(color)


pheno_data_EXO

msd <- readMSData(files = all_files,
                  pdata = new("NAnnotatedDataFrame",pheno_data_EXO),
                  mode = "onDisk",
                  centroided = TRUE)
msd 

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
write.csv(fData(msd), file=paste(filename = "Results/EXO_pos_raw_data.csv", sep = ""), row.names=FALSE)



# Get base peak chromatograms
register(bpstart(SnowParam()))
# setwd(input_dir_MS1_polarity)

chromas_EXO_pos <- chromatogram(msd, 
                             aggregationFun="max", 
                             msLevel = 1,
                             BPPARAM = SnowParam(workers = 3))
chromas_EXO_pos

# Plot chromatograms based on phenodata groups
#pdf(file="plots/EXO_chromas.pdf", encoding="ISOLatin1", pointsize=2, width=6, height=4, family="Helvetica")
jpeg(filename = "plots/EXO_chromas_pos.jpeg", width = 1000, height = 600, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(chromas_EXO_pos, main="Raw chromatograms", xlab="retention time [s]", ylab="intensity", col = color)
legend("topleft", bty="n", pt.cex=2, cex=1,5, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(msd@phenoData@data[["sample_group"]]))
dev.off()


# Get TICs
#pdf(file="plots/EXO_tics.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots/EXO_tics_pos.jpeg", width = 1000, height = 600, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(5,4,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=0.4, cex.lab=2, cex.main=2)
tics_EXO_pos <- split(tic(msd), f=fromFile(msd))
boxplot(tics_EXO_pos, col=color, ylab="intensity", xlab="sample", main="Total ion current", outline = FALSE)
legend("topleft", bty="n", pt.cex=2, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(msd@phenoData@data[["sample_group"]]))
dev.off()

### ---- pre-processing ----

# grouping/binning for peak detection, based on similarity of the base peak chromatogram 
chromas_bin_EXO_pos <- MSnbase::bin(chromas_EXO_pos, binSize=2)
chromas_bin_EXO_pos

chromas_bin_cor_EXO_pos <- cor(log2(do.call(cbind, lapply(chromas_bin_EXO_pos, intensity)))) # transformation with log
chromas_bin_cor_EXO_pos



colnames(chromas_bin_cor_EXO_pos) <- rownames(chromas_bin_cor_EXO_pos) <- msd$sample_name


chromas_bin_cor_EXO_pos[is.na(chromas_bin_cor_EXO_pos)] <- 0


jpeg(filename = "plots/heatmap_chromas_bin_EXO_pos.jpeg", width = 500, height = 500, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
heatmap(chromas_bin_cor_EXO_pos)
dev.off()


# Assess retention times and intensities of first file
head(rtime(chromas_EXO_pos[1, 1]))
head(intensity(chromas_EXO_pos[1, 1]))

# check for polarity
head(fData(msd)[, c("polarity", "filterString", "msLevel", "retentionTime")])
table(polarity(msd))


ms_params_EXO_pos <- CentWaveParam(ppm=24.5, mzCenterFun="wMean", peakwidth=c(12, 53), 
                                     prefilter=c(4, 60), mzdiff= -0.0032, snthresh=5, noise=0, 
                                     integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
                                     fitgauss=FALSE, roiList=list(), roiScales=numeric())
ms_params_EXO_pos

# CHOSE THE PARAMETERS 
ms_data_EXO_pos <- findChromPeaks(msd, param=ms_params_EXO_pos)
ms_data_EXO_pos

# check the detected peaks
head(chromPeaks(ms_data_EXO_pos))
chromPeakData(ms_data_EXO_pos)

ms_summary_EXO_pos <- lapply(split.data.frame(chromPeaks(ms_data_EXO_pos), 
                                            f=chromPeaks(ms_data_EXO_pos)[, "sample"]), 
                           FUN=function(z) { c(peak_count=nrow(z), rt=quantile(z[, "rtmax"] - z[, "rtmin"])) } )

ms_summary_EXO_pos <- do.call(rbind, ms_summary_EXO_pos)

rownames(ms_summary_EXO_pos) <- basename(fileNames(ms_data_EXO_pos))
rownames(ms_summary_EXO_pos)

print(ms_summary_EXO_pos)

table(msLevel(ms_data_EXO_pos))



write.csv(as.data.frame(table(msLevel(ms_data_EXO_pos))), file="Results/EXO_pos_ms_data.csv", row.names=FALSE)

jpeg(filename = "plots/EXO_ms_data.jpeg", width = 1000, height = 600, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(4,18,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plotChromPeakImage(ms_data_EXO_pos, main="Frequency of identified peaks per RT", binSize = 20)
#dev.off()

plotChromPeakImage(ms_data_EXO_pos, main="Frequency of identified peaks per RT", binSize = 20)

## Group peaks
ms_data_EXO_pos <- groupChromPeaks(ms_data_EXO_pos, param=PeakDensityParam(
  sampleGroups=ms_data_EXO_pos$sample_group, minFraction=0.7, bw=2.5))
ms_data_EXO_pos

## RT correction
ms_data_EXO_pos <- adjustRtime(ms_data_EXO_pos, param=PeakGroupsParam(
  minFraction=0.7,smooth="loess",span=0.8,family="gaussian"))

# Plot the difference of raw and adjusted retention times
#pdf(file="plots/EXO_ms1_raw_adjusted.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
jpeg(filename = "plots/EXO_ms_raw_adjusted.jpeg", width = 500, height = 1000, quality = 100, bg = "white")
par(mfrow=c(2,1), mar=c(4.5,4.2,4,1), cex=0.8)
plot(chromas_EXO_pos, peakType="none", main="Raw chromatograms")
plotAdjustedRtime(ms_data_EXO_pos, lwd=2, main="Retention Time correction")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
dev.off()

## Group peaks
ms_data_EXO_pos <- groupChromPeaks(ms_data_EXO_pos, param=PeakDensityParam(
  sampleGroups=ms_data_EXO_pos$sample_group, minFraction=0.7, bw=2.5))

# Get integrated peak intensity per feature/sample
print(head(featureValues(ms_data_EXO_pos, value="into")))

ppm <- 35  

# missing value imputation, see xcmsSet
ms_data_EXO_pos <- fillChromPeaks(ms_data_EXO_pos, param=FillChromPeaksParam(ppm=ppm, fixedRt=0, expandRt=5))
ms_data_EXO_pos

head(featureValues(ms_data_EXO_pos))
head(featureSummary(ms_data_EXO_pos, group=ms_data_EXO_pos$sample_group))



# Evaluate grouping
#pdf(file="plots/EXO_ms1_grouping.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots/EXO_pos_ms_grouping.jpeg", width = 1000, height = 500, quality = 150, bg = "white")
ms_pca_EXO_pos <- prcomp(t(na.omit(log2(featureValues(ms_data_EXO_pos, value="into")))), center=TRUE)
plot(ms_pca_EXO_pos$x[, 1], ms_pca_EXO_pos$x[,2], pch=19, main="PCA: Grouping of samples",
     xlab=paste0("PC1: ", format(summary(ms_pca_EXO_pos)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms_pca_EXO_pos)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color, cex=0.8)
grid()
text(ms_pca_EXO_pos$x[, 1], ms_pca_EXO_pos$x[,2], labels=ms_data_EXO_pos$sample_name, col=color, pos=3, cex=0.9)
legend("topleft", bty="n", pt.cex=1, cex=0.8, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(ms_data_EXO_pos@phenoData@data[["sample_group"]]))
dev.off()

# broken stick
png("plots/BrokenStick_EXO_pos_ms_grouping.png", width=10, height=6, units="in", res=100)
evplot = function(ev) {  
  # Broken stick model (MacArthur 1957)  
  n = length(ev)  
  bsm = data.frame(j=seq(1:n), p=0)  
  bsm$p[1] = 1/n  
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))  
  bsm$p = 100*bsm$p/n  
  # Plot eigenvalues and % of variation for each axis  
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))  
  barplot(ev, main="Eigenvalues EXO MS grouping", col="blue", las=2)  
  abline(h=mean(ev), col="red")  
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")  
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE,   
          main="% variation", col=c("blue",2), las=2)  
  legend("topright", c("% eigenvalue", "Broken stick model"),   
         pch=15, col=c("blue",2), bty="n")  
  par(op)  
} 

ev_pc = ms_pca_EXO_pos$sdev^2  
evplot(ev_pc)  
dev.off()

# Show peaks
tail(chromPeaks(ms_data_EXO_pos))
tail(chromPeakData(ms_data_EXO_pos))

# Show process history
processHistory(ms_data_EXO_pos)

# save as R object for use on Monday
MS_EXO_pos_peak_detection <- ms_data_EXO_pos




save(MS_EXO_pos_peak_detection, file = "Results/MS_EXO_pos_peak_detection.RData")

# ---------- Build MS1 feature tables ----------
# Build feature matrix
ms_matrix_EXO_pos <- featureValues(ms_data_EXO_pos, method="medret", value="into")

colnames(ms_matrix_EXO_pos) <- all_files_names
dim(ms_matrix_EXO_pos)
# transpose feature table
feat_list_EXO_pos <- t(ms_matrix_EXO_pos)
feat_list_EXO_pos

# Build feature summary
ms_summary_EXO_pos <- featureSummary(ms_data_EXO_pos)
ms_def_EXO_pos <- featureDefinitions(ms_data_EXO_pos)


# Missing value imputation by filling na positions with median of surrounding features
feat_list_EXO_pos[is.na(feat_list_EXO_pos)] <- median(na.omit(as.numeric(unlist(feat_list_EXO_pos))))
### Transform data
feat_list_EXO_pos <- log2(feat_list_EXO_pos)

# Missing value imputation
feat_list_EXO_pos[which(is.na(feat_list_EXO_pos))] <- median(na.omit(as.numeric(unlist(feat_list_EXO_pos))))

# save as csv
write.csv(feat_list_EXO_pos, file=paste(filename = "Results/feature_list_EXO_pos.csv", sep = ""))

# Plot histogram
#pdf(file="plots/EXO_feat_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots/EXO_pos_feat_list_hist.jpeg", width = 500, height = 500, quality = 150, bg = "white")
hist(as.numeric(feat_list_EXO_pos), main="Histogram of feature table")
dev.off()

# PCA of feature table results
#pdf(file="plots/EXO_ms1_feature_table_pca.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "plots/EXO_pos_ms_feature_table_pca.jpeg", width = 1000, height = 500, quality = 150, bg = "white")
ms_pca_EXO_pos <- prcomp(feat_list_EXO_pos, center=TRUE)
plot(ms_pca_EXO_pos$x[, 1], ms_pca_EXO_pos$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms_pca_EXO_pos)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms_pca_EXO_pos)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color, cex=0.8)
grid()
text(ms_pca_EXO_pos$x[, 1], ms_pca_EXO_pos$x[,2], labels=ms_data_EXO_pos$sample_name, col=color, pos=3, cex=0.9)
legend("topleft", bty="n", pt.cex=1, cex=1, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(ms_data_EXO_pos@phenoData@data[["sample_group"]]))
dev.off()

# broken stick
jpeg("plots/BrokenStick_EXO_pos_ms_feature_table_pca.jpeg", width=10, height=6, units="in", res=100)
evplot = function(ev) {  
  # Broken stick model (MacArthur 1957)  
  n = length(ev)  
  bsm = data.frame(j=seq(1:n), p=0)  
  bsm$p[1] = 1/n  
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))  
  bsm$p = 100*bsm$p/n  
  # Plot eigenvalues and % of variation for each axis  
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))  
  barplot(ev, main="Eigenvalues EXO MS feature table", col="blue", las=2)  
  abline(h=mean(ev), col="red")  
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")  
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE,   
          main="% variation", col=c("blue",2), las=2)  
  legend("topright", c("% eigenvalue", "Broken stick model"),   
         pch=15, col=c("blue",2), bty="n")  
  par(op)  
} 

ev_pc = ms_pca_EXO_pos$sdev^2  
evplot(ev_pc)  
dev.off()


ms_intensity_cutoff <- 100

# Create single 0/1 matrix
bina_list_EXO_pos <- t(ms_matrix_EXO_pos)
bina_list_EXO_pos[is.na(bina_list_EXO_pos)] <- 1
bina_list_EXO_pos <- log2(bina_list_EXO_pos)
bina_list_EXO_pos[bina_list_EXO_pos < log2(ms_intensity_cutoff)] <- 0
bina_list_EXO_pos[bina_list_EXO_pos != 0] <- 1


# save as csv
write.csv(bina_list_EXO_pos, file=paste(filename = "Results/bina_list_EXO_pos.csv", sep = ""))

# Only unique compounds in group mzml_pheno$ and not the others
uniq_list_EXO_pos <- apply(X=bina_list_EXO_pos, MARGIN=2, FUN=function(x) { if (length(unique(pheno_data_EXO$sample_group[grepl("1", x)])) == 1) x else rep(0, length(x)) } )
colnames(uniq_list_EXO_pos) <- colnames(bina_list_EXO_pos)
rownames(uniq_list_EXO_pos) <- rownames(bina_list_EXO_pos)
uniq_list_EXO_pos

# Create data frame
model_div_EXO_pos             <- data.frame(features=apply(X=bina_list_EXO_pos, MARGIN=1, FUN=function(x) { sum(x) } ))
model_div_EXO_pos$richness    <- apply(X=bina_list_EXO_pos, MARGIN=1, FUN=function(x) { sum(x) } )
#model_div_EXO_pos$menhinick   <- apply(X=bina_list_EXO_pos, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
model_div_EXO_pos$shannon     <- apply(X=feat_list_EXO_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
model_div_EXO_pos$pielou      <- apply(X=scale(feat_list_EXO_pos, center=F), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#model_div_EXO_pos$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list_EXO_pos, species), index="chao")
model_div_EXO_pos$simpson     <- apply(X=feat_list_EXO_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
model_div_EXO_pos$inverse     <- apply(X=feat_list_EXO_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
model_div_EXO_pos$fisher      <- apply(X=feat_list_EXO_pos, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
model_div_EXO_pos$unique      <- apply(X=uniq_list_EXO_pos, MARGIN=1, FUN=function(x) { sum(x) })

# Remove NAs if present
model_div_EXO_pos[is.na(model_div_EXO_pos)] <- 0

# save as csv
write.csv(model_div_EXO_pos, file=paste(filename = "Results/model_div_EXO_pos.csv", sep = ""))


# save the objects and tables



end.time <- Sys.time()

time.taken <- end.time - start.time
print(time.taken)






