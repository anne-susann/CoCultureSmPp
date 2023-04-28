# ---------- Preparations ----------
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
if (dir.exists(paste(getwd(), "/endo_pos_1ms2_plots/", sep = ""))){
	print("plots directory already exists")
	start_time <- Sys.time()
}  else{
	dir.create("endo_pos_1ms2_plots")
	start_time <- Sys.time()
	print("plots folder has been created")
}

raw_data_MS1_endo_pos <- list.files(input_dir_MS1, pattern = "pos_ENDO")
raw_data_MS1_endo_pos

raw_data_MS1_endo_pos <- list.files(input_dir_MS1, pattern = "pos_ENDO")
raw_data_MS1_endo_pos_files <- paste(input_dir_MS1, raw_data_MS1_endo_pos, sep ="")
raw_data_MS1_endo_pos_files

raw_data_MS2_endo_pos <- list.files(input_dir_MS2, pattern = "ENDOpos")
raw_data_MS2_endo_pos_files <- paste(input_dir_MS2, raw_data_MS2_endo_pos, sep ="")
raw_data_MS2_endo_pos_files

all_files <- c(raw_data_MS1_endo_pos_files, raw_data_MS2_endo_pos_files)
all_files

all_files_names <- str_remove(all_files, ".mzML")
all_files_names <- str_remove_all(all_files_names, input_dir_MS1)
all_files_names[26:(25 + length(raw_data_MS2_endo_pos))]<- str_remove_all(all_files_names[26:(25 + length(raw_data_MS2_endo_pos))], input_dir_MS2)
all_files_names



# create vector with sample classes according to culture information sheet
samp_groups <- c("CoCuPp", "CoCuSm", "CoCuSm", "CoCuPp", "CoCuPp", "CoCuSm",
				 rep(x = "Sm", times = 8), 
				 rep(x = "Pp", times = 8),
				 "CoCuSm", "CoCuPp", "MB", rep("MS2", length(raw_data_MS2_endo_pos)))

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
							 CoCuSm1, CoCuPp1, MB1, rep(ms2, length(raw_data_MS2_endo_pos)))

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
ms2<- rep("aquamarine", length(raw_data_MS2_endo_pos))
color <- c(CoCuPp1, CoCuSm1, CoCuPp2, CoCuSm2, Sm, Pp, CoCuSm3, CoCuPp3, MB, ms2)



all_files



# create phenodata based on culture type
pheno_data_endo_pos <- data.frame(sample_name = all_files_names, sample_group = samp_groups, samp_groups_description=samp_groups_description)
pheno_col_endo_pos <- data.frame(color)


pheno_data_endo_pos

msd <- readMSData(files = all_files,
				  pdata = new("NAnnotatedDataFrame",pheno_data_endo_pos),
				  mode = "onDisk",
				  centroided = TRUE)
msd 

# inspect data 
table(msLevel(msd))
table(polarity(msd))
head(fData(msd)[, c("scanWindowLowerLimit", "scanWindowUpperLimit",
					"originalPeaksCount", "msLevel", 
					"polarity", "retentionTime")])

# Restrict data to 1020 seconds (17 minutes)
msd <- filterRt(msd, c(0, 1020))

# subset data for msLevel = 1 and save raw data
#msd <- filterMsLevel(msd, msLevel = 1)
#table(msLevel(msd))
# filter polarity 
#msd <- filterPolarity(msd, polarity = 1)

# create result directory
if (dir.exists(paste(getwd(), "/endo_pos_1ms2_Results/", sep = ""))){
	print("plots directory already exists")
	start_time <- Sys.time()
}  else{
	dir.create("endo_pos_1ms2_Results")
	start_time <- Sys.time()
	print("results folder has been created")
}
write.csv(fData(msd), file=paste(filename = "endo_pos_1ms2_Results/endo_pos_1ms2_raw_data.csv", sep = ""), row.names=FALSE)



# Get base peak chromatograms
#register(bpstart(SnowParam()))
# setwd(input_dir_MS1_polarity)

chromas_endo_pos <- chromatogram(msd, 
								 aggregationFun="max", 
								 msLevel = 1 
								 )
chromas_endo_pos

# Plot chromatograms based on phenodata groups
#pdf(file="plots/ENDO_chromas.pdf", encoding="ISOLatin1", pointsize=2, width=6, height=4, family="Helvetica")
jpeg(filename = "endo_pos_1ms2_plots/endo_pos_chromas.jpeg", width = 2000, height = 1200, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=2, cex.lab=2, cex.main=2)
plot(chromas_endo_pos, main="Raw chromatograms", xlab="retention time [s]", ylab="intensity", col = color)
legend("topleft", bty="n", pt.cex=3, cex=1.5, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(msd@phenoData@data[["sample_group"]]))
dev.off()


# Get TICs
#pdf(file="plots/ENDO_tics.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "endo_pos_1ms2_plots/endo_pos_tics.jpeg", width = 2000, height = 1200, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=2, cex.lab=2, cex.main=2)
tics_endo_pos <- split(tic(msd), f=fromFile(msd))
boxplot(tics_endo_pos, col=color, ylab="intensity", xlab="sample", main="Total ion current", outline = FALSE)
legend("topleft", bty="n", pt.cex=2, cex=1, y.intersp=0.7, text.width=0.5, pch=20, 
	   col= unique(color), legend= unique(msd@phenoData@data[["sample_group"]]))
dev.off()

### ---- pre-processing ----

# grouping/binning for peak detection, based on similarity of the base peak chromatogram 
chromas_bin_endo_pos <- MSnbase::bin(chromas_endo_pos, binSize=2)
head(chromas_bin_endo_pos)

chromas_bin_cor_endo_pos <- cor(log2(do.call(cbind, lapply(chromas_bin_endo_pos, intensity)))) # transformation with log
head(chromas_bin_cor_endo_pos)



colnames(chromas_bin_cor_endo_pos) <- rownames(chromas_bin_cor_endo_pos) <- msd$sample_name


chromas_bin_cor_endo_pos[is.na(chromas_bin_cor_endo_pos)] <- 0


jpeg(filename = "endo_pos_1ms2_plots/heatmap_chromas_bin_endo_pos.jpeg", width = 1000, height = 1000, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(7,0,0,7), cex.axis=0.9, cex=0.6)
heatmap(chromas_bin_cor_endo_pos)
dev.off()


# Assess retention times and intensities of first file
head(rtime(chromas_endo_pos[1, 1]))
head(intensity(chromas_endo_pos[1, 1]))

# check for polarity
head(fData(msd)[, c("polarity", "filterString", "msLevel", "retentionTime")])
table(polarity(msd))

ms_params_endo_pos <- CentWaveParam(ppm=25, mzCenterFun="wMean", peakwidth=c(12, 51), 
                                    prefilter=c(4, 60), mzdiff= 0.000099, snthresh=6, noise=0, 
                                    integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
                                    fitgauss=FALSE, roiList=list(), roiScales=numeric())
ms_params_endo_pos

# CHOSE THE PARAMETERS 
ms_data_endo_pos <- findChromPeaks(msd, param=ms_params_endo_pos)
backup <- ms_data_endo_pos

# check the detected peaks
head(chromPeaks(ms_data_endo_pos))
chromPeakData(ms_data_endo_pos)

ms_summary_endo_pos <- lapply(split.data.frame(chromPeaks(ms_data_endo_pos), 
											   f=chromPeaks(ms_data_endo_pos)[, "sample"]), 
							  FUN=function(z) { c(peak_count=nrow(z), rt=quantile(z[, "rtmax"] - z[, "rtmin"])) } )

ms_summary_endo_pos <- do.call(rbind, ms_summary_endo_pos)

rownames(ms_summary_endo_pos) <- basename(fileNames(ms_data_endo_pos))
rownames(ms_summary_endo_pos)

print(ms_summary_endo_pos)

table(msLevel(ms_data_endo_pos))


write.csv(as.data.frame(table(msLevel(ms_data_endo_pos))), file="endo_pos_1ms2_Results/endo_pos_1ms2_data.csv", row.names=FALSE)

jpeg(filename = "endo_pos_1ms2_plots/endo_pos_ms_data.jpeg", width = 2000, height = 1200, quality = 100, bg = "white")
par(mfrow=c(1,1), mar=c(4,18,4,1), oma=c(0,0,0,0), cex.axis=1, cex=2, cex.lab=2, cex.main=2)
plotChromPeakImage(ms_data_endo_pos, main="Frequency of identified peaks per RT", binSize = 20)
dev.off()

plotChromPeakImage(ms_data_endo_pos, main="Frequency of identified peaks per RT", binSize = 20)


## Group peaks
ms_data_endo_pos <- groupChromPeaks(ms_data_endo_pos, param=PeakDensityParam(
	sampleGroups=ms_data_endo_pos$sample_group, minFraction=0.7, bw=2.5))
ms_data_endo_pos

## RT correction
ms_data_endo_pos <- adjustRtime(ms_data_endo_pos, param=PeakGroupsParam(
	minFraction=0.7,smooth="loess",span=0.5,family="gaussian"))

# Plot the difference of raw and adjusted retention times
#pdf(file="plots/ENDO_ms1_raw_adjusted.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
jpeg(filename = "endo_pos_1ms2_plots/endo_pos_ms_raw_adjusted.jpeg", width = 500, height = 1000, quality = 100, bg = "white")
par(mfrow=c(2,1), mar=c(4.5,4.2,4,1), cex=0.8)
plot(chromas_endo_pos, peakType="none", main="Raw chromatograms")
plotAdjustedRtime(filterRt(ms_data_endo_pos, rt = c(0, 670)), lwd=2, main="Retention Time correction")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
text(labels=ms_data_endo_pos$sample_name, pos=3, cex=1)
dev.off()

## Group peaks
ms_data_endo_pos <- groupChromPeaks(ms_data_endo_pos, param=PeakDensityParam(
	sampleGroups=ms_data_endo_pos$sample_group, minFraction=0.7, bw=2.5))

# Get integrated peak intensity per feature/sample
print(head(featureValues(ms_data_endo_pos, value="into")))

ppm <- 25  

# missing value imputation, see xcmsSet
#ms_data_endo_pos <- fillChromPeaks(ms_data_endo_pos, param=FillChromPeaksParam(ppm=ppm, fixedRt=0, expandRt=5))
ms_data_endo_pos

head(featureValues(ms_data_endo_pos))
head(featureSummary(ms_data_endo_pos, group=ms_data_endo_pos$sample_group))



# Evaluate grouping
#pdf(file="plots/ENDO_ms1_grouping.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "endo_pos_1ms2_plots/endo_pos_ms_grouping.jpeg", width = 1000, height = 700, quality = 150, bg = "white")
ms_pca_endo_pos <- prcomp(t(na.omit(log2(featureValues(ms_data_endo_pos, value="into")))), center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms_pca_endo_pos$x[, 1], ms_pca_endo_pos$x[,2], pch=19, main="PCA: Grouping of samples",
	 xlab=paste0("PC1: ", format(summary(ms_pca_endo_pos)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(ms_pca_endo_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=color, cex=2)
grid()
text(ms_pca_endo_pos$x[, 1], ms_pca_endo_pos$x[,2], labels=str_sub(ms_data_endo_pos$sample_name, - 3, - 1), col=color, pos=3, cex=0.9)
legend("topleft", bty="n", pt.cex=2, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
	   col= unique(color), legend= unique(ms_data_endo_pos@phenoData@data[["sample_group"]]))
dev.off()

# broken stick
png("endo_pos_1ms2_plots/BrokenStick_endo_pos_ms_grouping.png", width=10, height=6, units="in", res=100)
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

ev_pc = ms_pca_endo_pos$sdev^2  
evplot(ev_pc)  
dev.off()

# Show peaks
tail(chromPeaks(ms_data_endo_pos))
tail(chromPeakData(ms_data_endo_pos))

# Show process history
processHistory(ms_data_endo_pos)

# save as R object for use on Monday
MS_endo_pos_peak_detection <- ms_data_endo_pos


save(MS_endo_pos_peak_detection, file = "endo_pos_1ms2_Results/MS_endo_pos_peak_detection.RData")

# ---------- Build MS1 feature tables ----------
# Build feature matrix
ms_matrix_endo_pos <- featureValues(ms_data_endo_pos, method="medret", value="into")

colnames(ms_matrix_endo_pos) <- all_files_names
dim(ms_matrix_endo_pos)
# transpose feature table
feat_list_endo_pos <- t(ms_matrix_endo_pos)
feat_list_endo_pos

# Build feature summary
ms_summary_endo_pos <- featureSummary(ms_data_endo_pos)
ms_def_endo_pos <- featureDefinitions(ms_data_endo_pos)


# Missing value imputation by filling na positions with median of surrounding features
#feat_list_endo_pos[is.na(feat_list_endo_pos)] <- median(na.omit(as.numeric(unlist(feat_list_endo_pos))))

### Transform data
feat_list_endo_pos <- log2(feat_list_endo_pos)

# change 0 to small value to distinguish between values of 1 and NA
feat_list_endo_pos[which(feat_list_endo_pos == 0)] <- 0.01

# Missing value imputation
feat_list_endo_pos[which(is.na(feat_list_endo_pos))] <- 0

# save as csv
write.csv(feat_list_endo_pos, file=paste(filename = "endo_pos_1ms2_Results/feature_list_endo_pos.csv", sep = ""))

# Plot histogram
#pdf(file="plots/ENDO_feat_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "endo_pos_1ms2_plots/endo_pos_feat_list_hist.jpeg", width = 500, height = 500, quality = 150, bg = "white")
hist(as.numeric(feat_list_endo_pos), main="Histogram of feature table")
dev.off()

# PCA of feature table results
#pdf(file="plots/ENDO_ms1_feature_table_pca.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
jpeg(filename = "endo_pos_1ms2_plots/endo_pos_ms_feature_table_pca.jpeg", width = 1000, height = 700, quality = 150, bg = "white")
ms_pca_endo_pos <- prcomp(feat_list_endo_pos, center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms_pca_endo_pos$x[, 1], ms_pca_endo_pos$x[,2], pch=19, main="PCA of feature table",
	 xlab=paste0("PC1: ", format(summary(ms_pca_endo_pos)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(ms_pca_endo_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=color, cex=2)
grid()
text(ms_pca_endo_pos$x[, 1], ms_pca_endo_pos$x[,2], labels=str_sub(ms_data_endo_pos$sample_name, - 3, - 1), col=color, pos=3, cex=0.9)
legend("topleft", bty="n", pt.cex=2, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
	   col= unique(color), legend= unique(ms_data_endo_pos@phenoData@data[["sample_group"]]))
dev.off()

# broken stick
jpeg("endo_pos_1ms2_plots/BrokenStick_endo_pos_ms_feature_table_pca.jpeg", width=10, height=6, units="in", res=100)
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

ev_pc = ms_pca_endo_pos$sdev^2  
evplot(ev_pc)  
dev.off()


ms_intensity_cutoff <- 14

# Create single 0/1 matrix
bina_list_endo_pos <- t(ms_matrix_endo_pos)
bina_list_endo_pos[is.na(bina_list_endo_pos)] <- 1
bina_list_endo_pos <- log2(bina_list_endo_pos)
bina_list_endo_pos[bina_list_endo_pos < ms_intensity_cutoff] <- 0
bina_list_endo_pos[bina_list_endo_pos != 0] <- 1


# save as csv
write.csv(bina_list_endo_pos, file=paste(filename = "endo_pos_1ms2_Results/bina_list_endo_pos.csv", sep = ""))

# Only unique compounds in group mzml_pheno$ and not the others
uniq_list_endo_pos <- apply(X=bina_list_endo_pos, MARGIN=2, FUN=function(x) { if (length(unique(pheno_data_ENDO$sample_group[grepl("1", x)])) == 1) x else rep(0, length(x)) } )
colnames(uniq_list_endo_pos) <- colnames(bina_list_endo_pos)
rownames(uniq_list_endo_pos) <- rownames(bina_list_endo_pos)
uniq_list_endo_pos

# Create data frame
model_div_endo_pos             <- data.frame(features=apply(X=bina_list_endo_pos, MARGIN=1, FUN=function(x) { sum(x) } ))
model_div_endo_pos$richness    <- apply(X=bina_list_endo_pos, MARGIN=1, FUN=function(x) { sum(x) } )
#model_div_endo_pos$menhinick   <- apply(X=bina_list_endo_pos, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
model_div_endo_pos$shannon     <- apply(X=feat_list_endo_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
model_div_endo_pos$pielou      <- apply(X=scale(feat_list_endo_pos, center=F), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#model_div_endo_pos$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list_endo_pos, species), index="chao")
model_div_endo_pos$simpson     <- apply(X=feat_list_endo_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
model_div_endo_pos$inverse     <- apply(X=feat_list_endo_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
model_div_endo_pos$fisher      <- apply(X=feat_list_endo_pos, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
model_div_endo_pos$unique      <- apply(X=uniq_list_endo_pos, MARGIN=1, FUN=function(x) { sum(x) })
# funtional hill
#model_div_pos$hillfunc    <- as.numeric(unlist(calcDiv(feat_list_endo_pos, compDisMat=scales::rescale(as.matrix(dist(t(feat_list_endo_pos)), diag=TRUE, upper=TRUE)), q=1, type="FuncHillDiv")))


# Remove NAs if present
model_div_endo_pos[is.na(model_div_endo_pos)] <- 0

# save as csv
write.csv(model_div_endo_pos, file=paste(filename = "endo_pos_1ms2_Results/model_div_endo_pos.csv", sep = ""))


# save the objects and tables
save(ms1_def_ENDO_pos, file = "endo_pos_Results_MS1_Anne/ms1_def_ENDO_pos.RData")
save.image(file = "endo_pos_1ms2_Results/endo_pos_ms_environment.RData")


# --------- preparations linking MS2 data -----------
# object with MS1 and MS2 files preprocessed
#load()

# MS1 and MS2 files 
ms_data_endo_pos # needed for linking, created in this script
ms_def_endo_pos
table(msLevel(ms_data_endo_pos))


#table(msLevel(ms12_data_endo_pos))

# MS1 files
load("C:/Users/abela/Documents/Uni_Jena/Masterarbeit/MAW-Co-culture/endo_pos_Results_MS1_Anne/ms1_def_ENDO_pos.RData")
#ms1_data_endo_pos
#ms1_def_ENDO_pos # most important file for linking, created in MS1 script
ms1_def_ENDO_pos <- ms1_def_ENDO_pos
table(msLevel(ms1_data_endo_pos))


# ---------- MS2 spectra detection ----------
# Estimate precursor intensity
#precursor_intensity_endo_pos <- xcms::estimatePrecursorIntensity(ms_data_endo_pos)
#print(head(na.omit(precursor_intensity_endo_pos)))

# Reconstruct MS2 spectra from MS1 data
ms2_data_endo_pos <- chromPeakSpectra(ms_data_endo_pos, msLevel=2L, return.type="Spectra")
print(ms2_data_endo_pos)
print(length(ms2_data_endo_pos$peak_id))
head(ms2_data_endo_pos$peak_id)

# Extract all usable MS2 spectra
ms2_spectra_endo_pos <- list()
for (i in 1:nrow(ms1_def_endo_pos)) {
#ms2_spectra_endo_pos <- foreach(i=1:nrow(ms1_def_endo_pos)) %dopar% {
	#print(i)
	# Extract existing MS2 spectra for feature
	feature_of_interest <- ms1_def_endo_pos[i, "mzmed"]
	peaks_of_interest <- chromPeaks(ms_data_endo_pos, mz=feature_of_interest, ppm=ppm)
	
	# Continue if feature has MS2 peaks
	if (length(which(ms2_data_endo_pos$peak_id %in% rownames(peaks_of_interest))) > 0) {
		# Extract spectra
		spectra_of_interest <- ms2_data_endo_pos[ms2_data_endo_pos$peak_id %in% rownames(peaks_of_interest)]
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
		
		ms2_spectra_endo_pos[[i]] <- combined_spectra_of_interest
		#return(combined_spectra_of_interest)
	}
}

# Remove empty spectra
names(ms2_spectra_endo_pos) <- rownames(ms1_def_endo_pos)[1:length(ms2_spectra_endo_pos)]
ms2_spectra_endo_pos <- ms2_spectra_endo_pos[lengths(ms2_spectra_endo_pos) != 0]

# Relate all MS2 spectra to MS1 precursors
ms1_def_endo_pos$has_ms2 <- as.integer(rownames(ms1_def_endo_pos) %in% names(ms2_spectra_endo_pos))
print(paste0("Number of MS2 spectra related to precursor: ", length(which(ms1_def_endo_pos$has_ms2>0))))

# ADDED
polarity="positive"
pol="pos"

# Save all MS2 spectra in MGF file
mgf_text <- NULL
for (i in names(ms2_spectra_endo_pos)) {
	mgf_text <- c(mgf_text, paste0("COM=", i))
	mgf_text <- c(mgf_text, "BEGIN IONS")
	mgf_text <- c(mgf_text, "MSLEVEL=2")
	mgf_text <- c(mgf_text, paste0("TITLE=", i))
	mgf_text <- c(mgf_text, paste0("RTINSECONDS=", ms1_def_endo_pos[i, "rtmed"]))
	mgf_text <- c(mgf_text, paste0("PEPMASS=", ms1_def_endo_pos[i, "mzmed"]))
	if (polarity == "positive") {
		mgf_text <- c(mgf_text, paste0("CHARGE=", "1+"))
	} else {
		mgf_text <- c(mgf_text, paste0("CHARGE=", "1-"))
	}
	mgf_text <- c(mgf_text, paste(as.data.frame(peaksData(ms2_spectra_endo_pos[[i]])[[1]])$mz, as.data.frame(peaksData(ms2_spectra_endo_pos[[i]])[[1]])$intensity, sep=" "))
	mgf_text <- c(mgf_text, "END IONS")
	mgf_text <- c(mgf_text, "")
}

# Write MGF file
cat(mgf_text, file="ms2_spectra_endo_pos.mgf", sep="\n")



################# MAW #######################

# Fix some naming
mzml_names_endo_pos <- all_files_names




# ---------- Annotate MS2 spectra with SIRIUS ----------
# Apply annotated compounds onto feature table
ms1_def_endo_pos$has_id <- ""
ms1_def_endo_pos$smiles <- ""
ms1_def_endo_pos$name <- ""

# Read SIRIUS annotation
annotator_sirius_endo_pos <- read.table(file=paste0("compound_identifications.tsv"), header=TRUE, sep="\t", quote="\"", comment.char="", fill=TRUE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
annotator_sirius_endo_pos$Metabolite.name <- gsub(x=annotator_sirius_endo_pos$id, pattern=".*_", replacement="")

for (j in unique(annotator_sirius_endo_pos$Metabolite.name)) {
	obj <- annotator_sirius_endo_pos[annotator_sirius_endo_pos$Metabolite.name %in% j, "smiles"]
	ms1_def_endo_pos[j, "smiles"] <- obj[1]
	
	obj <- annotator_sirius_endo_pos[annotator_sirius_endo_pos$Metabolite.name %in% j, "name"]
	if (obj[1] != "null" ) {
		ms1_def_endo_pos[j, "name"] <- obj[1]
	}
}

# Annotation
sirius_compounds_endo_pos <- annotator_sirius_endo_pos$smiles
sirius_compound_ids_endo_pos <- gsub(x=annotator_sirius_endo_pos$id, pattern=".*_", replacement="")

# Calculate molecular descriptors with CDK
cdk_descriptors_endo_pos <- NULL
for (i in sirius_compounds_endo_pos) {
	# Get Structure from SMILES
	cdk_mol = parse.smiles(i)[[1]]
	
	# Get simple measures
	cdk_atoms = get.atoms(cdk_mol)
	cdk_bonds = get.bonds(cdk_mol)
	
	# Calculate simple measures
	MolWeight = get.exact.mass(cdk_mol)
	nAtoms = get.atom.count(cdk_mol)
	cdk_num_atoms = as.factor(unlist(lapply(cdk_atoms, get.symbol)))
	cdk_num_atoms = tapply(cdk_num_atoms, cdk_num_atoms, length)
	numC = as.numeric(cdk_num_atoms["C"])
	numN = as.numeric(cdk_num_atoms["N"])
	numP = as.numeric(cdk_num_atoms["P"])
	numO = as.numeric(cdk_num_atoms["O"])
	numHydrogen = get.total.hydrogen.count(cdk_mol)
	CNRatio = as.numeric(numC / numN)
	XLogP = get.xlogp(cdk_mol)
	
	# Calculate descriptors and restrict to only "constitutional" and "topological"
	cdk_mol_des_cats = get.desc.categories()
	cdk_mol_des_names = c(get.desc.names(cdk_mol_des_cats[3]), get.desc.names(cdk_mol_des_cats[4]))
	cdk_mol_des = as.data.frame(eval.desc(cdk_mol, cdk_mol_des_names))
	cdk_descriptors_endo_pos <- plyr::rbind.fill(cdk_descriptors_endo_pos, cbind(data.frame(MolWeight=MolWeight, nAtoms=nAtoms, numC=numC, numN=numN, numP=numP, numO=numO, numHydrogen=numHydrogen, CNRatio=CNRatio, XLogP=XLogP), cdk_mol_des))
}
rownames(cdk_descriptors_endo_pos) <- paste0(sirius_compound_ids_endo_pos,"_",pol)

# Properly format NAs and convert to numeric
for (i in 1:nrow(cdk_descriptors_endo_pos)) {
	cdk_descriptors_endo_pos[i, as.character(cdk_descriptors_endo_pos[i,]) %in% 'list(NULL)'] <- 0
	cdk_descriptors_endo_pos[i, as.character(cdk_descriptors_endo_pos[i,]) %in% 'NA'] <- 0
	cdk_descriptors_endo_pos[i, as.character(cdk_descriptors_endo_pos[i,]) %in% 'NULL'] <- 0
	cdk_descriptors_endo_pos[i, as.character(cdk_descriptors_endo_pos[i,]) %in% 'NaN'] <- 0
	cdk_descriptors_endo_pos[i, ] <- as.numeric(unlist(cdk_descriptors_endo_pos[i,]))
}

# Bug: Remove left-over variables
rm(CNRatio, MolWeight, XLogP, nAtoms, numC, numHydrogen, numN, numO, numP)



# ---------- Classify MS2 spectra with CANOPUS ----------
# Apply classified classes onto feature table
ms1_def_endo_pos$primary_class <- ""
ms_def_endo_pos$alternative_classes <- ""
ms_def_endo_pos$pathway <- ""
SIRIUS_VERSION <- 5

# Read SIRIUS/CANOPUS classifier
if (SIRIUS_VERSION == 4) {
	classifier_canopus_endo_pos <- read.table(file=paste0("canopus_compound_summary.tsv"), header=TRUE, sep="\t", quote="\"", fill=FALSE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
	classifier_canopus_endo_pos$Metabolite.name <- gsub(x=classifier_canopus_endo_pos$name, pattern=".*_", replacement="")
	classifier_canopus_endo_pos$primary_class <- paste("Organic compounds", classifier_canopus_endo_pos$superclass, classifier_canopus_endo_pos$class, classifier_canopus_endo_pos$subclass, classifier_canopus_endo_pos$level.5, sep="; ")
	classifier_canopus_endo_pos$primary_class <- gsub(x=classifier_canopus_endo_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_endo_pos$primary_class <- gsub(x=classifier_canopus_endo_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_endo_pos$primary_class <- gsub(x=classifier_canopus_endo_pos$primary_class, pattern="; $", replacement="", perl=TRUE)
	classifier_canopus_endo_pos$primary_class <- gsub(x=classifier_canopus_endo_pos$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	classifier_canopus_endo_pos$Annotation..putative. <- classifier_canopus_endo_pos$primary_class
	classifier_canopus_endo_pos$alternative_classes <- classifier_canopus_endo_pos$all.classifications
} else {
	classifier_canopus_endo_pos <- read.table(file=paste0("canopus_compound_summary.tsv"), header=TRUE, sep="\t", quote="\"", comment.char="", fill=TRUE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
	classifier_canopus_endo_pos$Metabolite.name <- ""
	for (i in 1:length(classifier_canopus_endo_pos$id)) {
		x = which(gsub(x=classifier_canopus_endo_pos$id, pattern=".*?_", replacement="", perl=TRUE) %in% gsub(x=classifier_canopus_endo_pos$id[i], pattern=".*?_", replacement="", perl=TRUE))
		if (length(x) > 0) classifier_canopus_endo_pos$Metabolite.name[i] <- gsub(x=classifier_canopus_endo_pos$id[x], pattern=".*?_", replacement="", perl=TRUE)
	}
	classifier_canopus_endo_pos$Metabolite.name[classifier_canopus_endo_pos$Metabolite.name == "null"] <- ""
	classifier_canopus_endo_pos$primary_class <- paste("Organic compounds", classifier_canopus_endo_pos$ClassyFire.superclass, classifier_canopus_endo_pos$ClassyFire.class, classifier_canopus_endo_pos$ClassyFire.subclass, classifier_canopus_endo_pos$ClassyFire.level.5, sep="; ")
	classifier_canopus_endo_pos$primary_class <- gsub(x=classifier_canopus_endo_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_endo_pos$primary_class <- gsub(x=classifier_canopus_endo_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_endo_pos$primary_class <- gsub(x=classifier_canopus_endo_pos$primary_class, pattern="; $", replacement="", perl=TRUE)
	classifier_canopus_endo_pos$primary_class <- gsub(x=classifier_canopus_endo_pos$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	classifier_canopus_endo_pos$Annotation..putative. <- classifier_canopus_endo_pos$primary_class
	classifier_canopus_endo_pos$alternative_classes <- classifier_canopus_endo_pos$all.classifications
}

for (j in unique(classifier_canopus_endo_pos$Metabolite.name)) {
	# CHEMONT
	obj <- classifier_canopus_endo_pos[classifier_canopus_endo_pos$Metabolite.name %in% j, "Annotation..putative."]
	ms1_def_endo_pos[j, "primary_class"] <- obj[1]
	
	# NPClassifier
	if (SIRIUS_VERSION > 4) {
		obj <- classifier_canopus_endo_pos[classifier_canopus_endo_pos$Metabolite.name %in% j, "NPC.pathway"]
		ms1_def_endo_pos[j, "pathway"] <- obj[1]
	}
}



# ---------- Diversity of MS2 classes ----------
# Create CANOPUS classifier object for each sample
classifiers_endo_pos <- list()
for (i in mzml_names_endo_pos) {
	obj <- names(which(bina_list_endo_pos[rownames(bina_list_endo_pos)==i, colnames(bina_list_endo_pos) %in% classifier_canopus_endo_pos$Metabolite.name] > 0))
	classifiers_endo_pos[[i]] <- classifier_canopus_endo_pos[classifier_canopus_endo_pos$Metabolite.name %in% obj, ]
}

# Diversity of classes per sample
div_classes_samples_endo_pos <- NULL
for (i in mzml_names_endo_pos) {
	obj <- table(classifiers_endo_pos[[i]][,"Annotation..putative."])
	obj <- data.frame(classes=names(obj), frequency=as.numeric(obj))
	if (is.null(div_classes_samples_endo_pos)) {
		div_classes_samples_endo_pos <- obj
	} else {
		div_classes_samples_endo_pos <- merge(div_classes_samples_endo_pos, obj, by="classes", all.x=TRUE, all.y=TRUE)
	}
}
rownames(div_classes_samples_endo_pos) <- div_classes_samples_endo_pos$classes
div_classes_samples_endo_pos <- div_classes_samples_endo_pos[, -which(colnames(div_classes_samples_endo_pos)=="classes")]
colnames(div_classes_samples_endo_pos) <- mzml_names_endo_pos

# Diversity of classes
div_classes_endo_pos <- div_classes_samples_endo_pos
div_classes_endo_pos[is.na(div_classes_endo_pos)] <- 0
div_classes_endo_pos <- apply(X=div_classes_endo_pos, MARGIN=1, FUN=function(x) { sum(x) })
div_classes_endo_pos <- data.frame(row.names=names(div_classes_endo_pos), frequency=as.numeric(div_classes_endo_pos))

# Plot diversity of classes in all samples
pdf(file="plots/pos_ms2_classes_diversity.pdf", encoding="ISOLatin1", pointsize=8, width=6, height=14, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,15,4,1), oma=c(0,0,0,0), cex.axis=0.8, cex=0.9)
barplot(div_classes_endo_pos$frequency, names.arg=gsub('.*; ','',rownames(div_classes_endo_pos)), las=1, horiz=TRUE, xlab="frequency", main="Diversity of compound classes", col=rainbow(n=nrow(div_classes_endo_pos), alpha=0.6))
dev.off()

# Sunburst plot of classes of all samples
pdf(file="plots/pos_ms2_classes_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
sunBurstPlotFromSubstanceClasses(rownames(div_classes_endo_pos), div_classes_endo_pos$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_classes_endo_pos), "Number of spectra"=div_classes_endo_pos$frequency), file="plots/pos_ms2_classes_sunburst.csv", row.names=FALSE)

# Classes
classes_endo_pos <- rownames(div_classes_endo_pos)
classes_endo_pos <- classes_endo_pos[which(grepl(pattern="^Organic compounds", x=classes_endo_pos))]
classes_endo_pos <- gsub(x=classes_endo_pos, pattern=":", replacement="")
classes_endo_pos <- gsub(x=classes_endo_pos, pattern="/", replacement="; ")



# ---------- Build diversity objects ----------
# Imputation of NA with zeros
div_classes_endo_pos[is.na(div_classes_endo_pos)] <- 0
div_classes_samples_endo_pos[is.na(div_classes_samples_endo_pos)] <- 0

# Classification list for statistics
class_list_endo_pos <- as.data.frame(t(div_classes_samples_endo_pos))
class_list_endo_pos[is.na(class_list_endo_pos)] <- 0

# Log Transform
#class_list_endo_pos <- log2(class_list_endo_pos + 1)

# Only keep class names
colnames(class_list_endo_pos) <- gsub(x=colnames(class_list_endo_pos), pattern='.*; ', replacement='')

# Generate class_int_list_endo_pos with abundances instead of counts
class_int_list_endo_pos <- class_list_endo_pos

for (i in 1:nrow(class_list_endo_pos)) {
	samp <- rownames(class_list_endo_pos)[i]
	for (j in 1:ncol(class_list_endo_pos)) {
		cl <- colnames(class_list_endo_pos)[j]
		ft <- classifier_canopus_endo_pos$Metabolite.name[which(gsub(x=classifier_canopus_endo_pos$primary_class, pattern='.*; ', replacement='') == cl)]
		ints <- as.numeric(feat_list_endo_pos[i, which(colnames(feat_list_endo_pos) %in% ft)])
		class_int_list_endo_pos[i, j] <- sum(ints) * as.numeric(class_list_endo_pos[i, j])
	}
}



# ---------- Classification at CHEMONT level of classes ----------
# Make superclasses at CHEMONT level 2
superclass_level_endo_pos <- 2
div_superclasses_samples_names_endo_pos <- NULL
for (i in c(1:superclass_level_endo_pos)) {
	div_superclasses_samples_names_endo_pos <- c(div_superclasses_samples_names_endo_pos, lapply(X=strsplit(rownames(div_classes_samples_endo_pos), '; '), FUN=function(x) { gsub(x=paste(x[1:i],sep='',collapse='; '),pattern='; NA',replacement='') }))
}
div_superclasses_samples_endo_pos <- data.frame()
for (i in c(1:ncol(div_classes_samples_endo_pos))) div_superclasses_samples_endo_pos <- rbind(div_superclasses_samples_endo_pos, rep(0, length(unique(div_superclasses_samples_names_endo_pos))))
div_superclasses_samples_endo_pos <- t(div_superclasses_samples_endo_pos)
colnames(div_superclasses_samples_endo_pos) <- colnames(div_classes_samples_endo_pos)
rownames(div_superclasses_samples_endo_pos) <- unique(div_superclasses_samples_names_endo_pos)
for (i in rownames(div_superclasses_samples_endo_pos)) {
	for (j in c(1:ncol(div_classes_samples_endo_pos))) {
		div_superclasses_samples_endo_pos[rownames(div_superclasses_samples_endo_pos)==i, j] <- sum(div_classes_samples_endo_pos[grep(x=rownames(div_classes_samples_endo_pos), pattern=i), j])
	}
}

# Diversity of superclasses
div_superclasses_endo_pos <- div_superclasses_samples_endo_pos
div_superclasses_endo_pos[is.na(div_superclasses_endo_pos)] <- 0
div_superclasses_endo_pos <- apply(X=div_superclasses_endo_pos, MARGIN=1, FUN=function(x) { sum(x) })
div_superclasses_endo_pos <- data.frame(row.names=names(div_superclasses_endo_pos), frequency=as.numeric(div_superclasses_endo_pos))

# Sunburst plot
pdf(file="plots/pos_ms2_superclasses_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(div_superclasses_endo_pos), div_superclasses_endo_pos$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_superclasses_endo_pos), "Number of spectra"=div_superclasses_endo_pos$frequency), file="plots/pos_ms2_superclasses_sunburst.csv", row.names=FALSE)

# Imputation of NA with zeros
div_superclasses_endo_pos[is.na(div_superclasses_endo_pos)] <- 0
div_superclasses_samples_endo_pos[is.na(div_superclasses_samples_endo_pos)] <- 0

# Classification list for statistics
superclass_list_endo_pos <- as.data.frame(t(div_superclasses_samples_endo_pos))
superclass_list_endo_pos[is.na(superclass_list_endo_pos)] <- 0

# Log Transform
#superclass_list_endo_pos <- log2(superclass_list_endo_pos + 1)

# Only keep superclass names
colnames(superclass_list_endo_pos) <- gsub(x=colnames(superclass_list_endo_pos), pattern='.*; ', replacement='')

# Generate superclass_int_list_endo_pos with abundances instead of counts
superclass_int_list_endo_pos <- superclass_list_endo_pos

for (i in 1:nrow(superclass_list_endo_pos)) {
	samp <- rownames(superclass_list_endo_pos)[i]
	for (j in 1:ncol(superclass_list_endo_pos)) {
		cl <- colnames(superclass_list_endo_pos)[j]
		ft <- classifier_canopus_endo_pos$Metabolite.name[which(classifier_canopus_endo_pos$primary_class %in% classifier_canopus_endo_pos$primary_class[grep(x=classifier_canopus_endo_pos$primary_class, pattern=cl)]) ]
		ints <- as.numeric(feat_list_endo_pos[i, which(colnames(feat_list_endo_pos) %in% ft)])
		superclass_int_list_endo_pos[i, j] <- sum(ints) * as.numeric(superclass_list_endo_pos[i, j])
	}
}



# ---------- Classification at CHEMONT level of subclasses ----------
# Make subclasses at CHEMONT level 3
subclass_level_endo_pos <- 3
div_subclasses_samples_names_endo_pos <- NULL
for (i in c(1:subclass_level_endo_pos)) {
	div_subclasses_samples_names_endo_pos <- c(div_subclasses_samples_names_endo_pos, lapply(X=strsplit(rownames(div_classes_samples_endo_pos), '; '), FUN=function(x) { gsub(x=paste(x[1:i],sep='',collapse='; '),pattern='; NA',replacement='') }))
}
div_subclasses_samples_endo_pos <- data.frame()
for (i in c(1:ncol(div_classes_samples_endo_pos))) div_subclasses_samples_endo_pos <- rbind(div_subclasses_samples_endo_pos, rep(0, length(unique(div_subclasses_samples_names_endo_pos))))
div_subclasses_samples_endo_pos <- t(div_subclasses_samples_endo_pos)
colnames(div_subclasses_samples_endo_pos) <- colnames(div_classes_samples_endo_pos)
rownames(div_subclasses_samples_endo_pos) <- unique(div_subclasses_samples_names_endo_pos)
for (i in rownames(div_subclasses_samples_endo_pos)) {
	for (j in c(1:ncol(div_classes_samples_endo_pos))) {
		div_subclasses_samples_endo_pos[rownames(div_subclasses_samples_endo_pos)==i, j] <- sum(div_classes_samples_endo_pos[grep(x=rownames(div_classes_samples_endo_pos), pattern=i), j])
	}
}

# Diversity of subclasses
div_subclasses_endo_pos <- div_subclasses_samples_endo_pos
div_subclasses_endo_pos[is.na(div_subclasses_endo_pos)] <- 0
div_subclasses_endo_pos <- apply(X=div_subclasses_endo_pos, MARGIN=1, FUN=function(x) { sum(x) })
div_subclasses_endo_pos <- data.frame(row.names=names(div_subclasses_endo_pos), frequency=as.numeric(div_subclasses_endo_pos))

# Sunburst plot
pdf(file="plots/pos_ms2_subclasses_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(div_subclasses_endo_pos), div_subclasses_endo_pos$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_subclasses_endo_pos), "Number of spectra"=div_subclasses_endo_pos$frequency), file="plots/pos_ms2_subclasses_sunburst.csv", row.names=FALSE)

# Imputation of NA with zeros
div_subclasses_endo_pos[is.na(div_subclasses_endo_pos)] <- 0
div_subclasses_samples_endo_pos[is.na(div_subclasses_samples_endo_pos)] <- 0

# Classification list for statistics
subclass_list_endo_pos <- as.data.frame(t(div_subclasses_samples_endo_pos))
subclass_list_endo_pos[is.na(subclass_list_endo_pos)] <- 0

# Log Transform
#subclass_list_endo_pos <- log2(subclass_list_endo_pos + 1)

# Only keep subclass names
colnames(subclass_list_endo_pos) <- gsub(x=colnames(subclass_list_endo_pos), pattern='.*; ', replacement='')

# Generate subclass_int_list_endo_pos with abundances instead of counts
subclass_int_list_endo_pos <- subclass_list_endo_pos

for (i in 1:nrow(subclass_list_endo_pos)) {
	samp <- rownames(subclass_list_endo_pos)[i]
	for (j in 1:ncol(subclass_list_endo_pos)) {
		cl <- colnames(subclass_list_endo_pos)[j]
		ft <- classifier_canopus_endo_pos$Metabolite.name[which(classifier_canopus_endo_pos$primary_class %in% classifier_canopus_endo_pos$primary_class[grep(x=classifier_canopus_endo_pos$primary_class, pattern=cl)]) ]
		ints <- as.numeric(feat_list_endo_pos[i, which(colnames(feat_list_endo_pos) %in% ft)])
		subclass_int_list_endo_pos[i, j] <- sum(ints) * as.numeric(subclass_list_endo_pos[i, j])
	}
}




# ############################## MS1 statistics examples ##############################

mzml_pheno_samples <- samp_groups_description
mzml_pheno_colors <- color
principal_components <- 5

# comp_list
comp_list <- feat_list_endo_pos[, c(rownames(ms1_def_endo_pos)[which(ms1_def_endo_pos$has_ms2==1)])]
colnames(comp_list) <- paste0(rownames(ms1_def_endo_pos)[which(ms1_def_endo_pos$has_ms2==1)], "_pos")

# bina_list
bina_list <- feat_list_endo_pos[, c(rownames(ms1_def_endo_pos)[which(ms1_def_endo_pos$has_ms2==1)])]
colnames(bina_list) <- paste0(rownames(ms1_def_endo_pos)[which(ms1_def_endo_pos$has_ms2==1)], "_pos")

# uniq_list
uniq_list <- uniq_list_endo_pos[, c(rownames(ms1_def_endo_pos)[which(ms1_def_endo_pos$has_ms2==1)])]
rownames(uniq_list_endo_pos) <- gsub(x=rownames(uniq_list_endo_pos), pattern="\\.auto.*", replacement="")
colnames(uniq_list_endo_pos) <- paste0(rownames(ms1_def_endo_pos)[which(ms1_def_endo_pos$has_ms2==1)], "_pos")

# Merge pos and pos div_classes
div_classes <- div_classes_endo_pos
div_classes_samples <- div_classes_samples_endo_pos

# Merge pos and pos div_superclasses
div_superclasses <- div_superclasses_endo_pos
div_superclasses_samples <- div_superclasses_samples_endo_pos

# Merge pos and pos div_subclasses
div_subclasses <- div_subclasses_endo_pos
div_subclasses_samples <- div_subclasses_samples_endo_pos

# class_list
class_list <- class_list_endo_pos
rownames(class_list) <- gsub(x=rownames(class_list_endo_pos), pattern="\\.auto.*", replacement="")

# class_int_list
class_int_list <- class_int_list_endo_pos
rownames(class_int_list) <- gsub(x=rownames(class_int_list_endo_pos), pattern="\\.auto.*", replacement="")

# superclass_list
superclass_list <- superclass_list_endo_pos
rownames(superclass_list) <- gsub(x=rownames(superclass_list_endo_pos), pattern="\\.auto.*", replacement="")

# superclass_int_list
superclass_int_list <- superclass_int_list_endo_pos
rownames(superclass_int_list) <- gsub(x=rownames(superclass_int_list_endo_pos), pattern="\\.auto.*", replacement="")

# subclass_list
subclass_list <- subclass_list_endo_pos
rownames(subclass_list) <- gsub(x=rownames(subclass_list_endo_pos), pattern="\\.auto.*", replacement="")

# subclass_int_list
subclass_int_list <- subclass_int_list_endo_pos
rownames(subclass_int_list) <- gsub(x=rownames(subclass_int_list_endo_pos), pattern="\\.auto.*", replacement="")

cdk_descriptors <- cdk_descriptors_endo_pos

# Create data frame
model_div             <- data.frame(features=apply(X=bina_list, MARGIN=1, FUN=function(x) { sum(x) } ))
model_div$richness    <- apply(X=bina_list, MARGIN=1, FUN=function(x) { sum(x) } )
model_div$menhinick   <- apply(X=bina_list, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
model_div$shannon     <- apply(X=comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
model_div$pielou      <- apply(X=scale(comp_list, center=FALSE, scale=FALSE), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#model_div$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list, species), index="chao")
model_div$simpson     <- apply(X=comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
model_div$inverse     <- apply(X=comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
model_div$fisher      <- apply(X=comp_list, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
model_div$unique      <- apply(X=uniq_list, MARGIN=1, FUN=function(x) { sum(x) })
model_div$hillfunc    <- as.numeric(unlist(calcDiv(comp_list, compDisMat=scales::rescale(as.matrix(dist(t(comp_list)), diag=TRUE, upper=TRUE)), q=1, type="FuncHillDiv")))

# Plot Shannon index
pdf(paste("plots/ms1_comp_list_diversity_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$shannon ~ mzml_pheno_samples, col=mzml_pheno_colors, names=NA, main="Shannon diversity (H\')", xlab="treatment", ylab="Shannon diversity index (H\')")
text(1:length(levels(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=levels(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$shannon, term=as.factor(mzml_pheno_samples))
text(1:length(levels(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# PLS
sel_pls_comp_list <- f.select_features_pls(feat_matrix=comp_list, sel_factor=mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=principal_components, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/ms1_comp_list_select_pls_roc.pdf")
print(paste("Number of selected variables:", f.count.selected_features(sel_feat=sel_pls_comp_list$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=comp_list, sel_feat=sel_pls_comp_list$`_selected_variables_`, sel_names=paste0("         ",sel_pls_comp_list$`_selected_variables_`), sample_colors=mzml_pheno_colors, plot_width=7, plot_height=7, cex_col=0.5, cex_row=0.4, filename="plots/ms1_comp_list_select_pls.pdf", main="PLS")
heatmaply(scale(comp_list[, which(colnames(comp_list) %in% sel_pls_comp_list$`_selected_variables_`)]), k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), file="plots/ms1_comp_list_select_pls.html", selfcontained=TRUE)
sel_pls_comp_list$`_selected_variables_`
sel_pls_comp_list$`_model_r2_`
sel_pls_comp_list$`_multiclass_metrics_`



# ############################## MS2 statistics examples ##############################

# PLS
sel_pls_class_list <- f.select_features_pls(feat_matrix=class_list, sel_factor=mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=principal_components, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/ms2_class_list_select_pls_roc.pdf")
print(paste("Number of selected variables:", f.count.selected_features(sel_feat=sel_pls_class_list$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=class_list, sel_feat=sel_pls_class_list$`_selected_variables_`, sample_colors=mzml_pheno_colors, plot_width=7, plot_height=7, cex_col=0.5, cex_row=0.4, filename="plots/ms2_class_list_select_pls.pdf", main="PLS")
heatmaply(scale(class_list[, which(colnames(class_list) %in% sel_pls_class_list$`_selected_variables_`)]), k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), file="plots/ms2_class_list_select_pls.html", selfcontained=TRUE)
sel_pls_class_list$`_multiclass_metrics_`
sel_pls_class_list$`_model_r2_`



# ############################## Molecular descriptors examples ##############################


# Table L: samples x metabolites
mdes_tab_l <- bina_list
mdes_tab_l <- bina_list[, which(colnames(bina_list) %in% paste0(rownames(ms1_def_endo_pos)[which(ms1_def_endo_pos$smiles != "")], "_pos"))]
mdes_tab_l <- as.data.frame(mdes_tab_l)

# Table R: samples x species
mdes_tab_r <- as.data.frame.matrix(table(rownames(mdes_tab_l), mzml_pheno_samples))
rownames(mdes_tab_r) <- rownames(mdes_tab_l)

# Table Q: metabolites x traits
mdes_tab_q <- cdk_descriptors
mdes_tab_q[is.na(mdes_tab_q)] <- 0

# Perform matrix operation
mdes_list <- as.data.frame(as.matrix(mdes_tab_l) %*% as.matrix(mdes_tab_q))

# PLS
sel_pls_mdes_list <- f.select_features_pls(feat_matrix=mdes_list, sel_factor=mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=principal_components, tune_length=10, quantile_threshold=0.995, plot_roc_filename="plots/descriptors_bina_list_select_pls_roc.pdf")
print(paste("Number of selected descriptors:", f.count.selected_features(sel_feat=sel_pls_mdes_list$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=mdes_list, sel_feat=sel_pls_mdes_list$`_selected_variables_`, sel_names=paste0("        ",sel_pls_mdes_list$`_selected_variables_`), sample_colors=mzml_pheno_colors, plot_width=7, plot_height=7, cex_col=0.5, cex_row=0.4, filename="plots/descriptors_bina_list_select_pls.pdf", main="PLS")
heatmaply(scale(mdes_list[, which(colnames(mdes_list) %in% sel_pls_mdes_list$`_selected_variables_`)]), k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), file="plots/descriptors_bina_list_select_pls.html", selfcontained=TRUE)
sel_pls_mdes_list$`_multiclass_metrics_`
sel_pls_mdes_list$`_model_r2_`
sel_pls_mdes_list$`_selected_variables_`



