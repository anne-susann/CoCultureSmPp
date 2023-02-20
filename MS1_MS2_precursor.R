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
MS2_ENDO_files <- data.frame(subset(MS2_files, grepl("ENDO", MS2_files[,1]), drop = TRUE))
MS2_EXO_files <- data.frame(subset(MS2_files, grepl("EXO", MS2_files[,1]), drop = TRUE))




###----MS1 preparations----
# MS1 variables
# pol <- substr(x = polarity, start = 1, stop = 3)
ppm <- 35
ms1_intensity_cutoff <- 100	          #approx. 0.01%

# General variables
# mzml_files_ENDO <- NULL
# mzml_names_neg <- NULL
mzml_times_ENDO <- NULL




###-----data-------
# iESTIMATE code as inspiration (IPB Halle, github)
# https://github.com/ipb-halle/iESTIMATE/blob/main/use-cases/radula-hormones/peak_detection_neg.r

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
for (i in 1:length(MS1_ENDO_files_full)) {
  mzml_data_ENDO <- readMSData(files = paste(input_dir_MS1, MS1_ENDO_files[,i], sep = ""), mode="onDisk")
  mzml_msn_ENDO <- rbind(mzml_msn_ENDO, t(as.matrix(table(msLevel(mzml_data_ENDO)))))
}
colnames(mzml_msn_ENDO) <- c("MSn 0", "MSn 1")
rownames(mzml_msn_ENDO) <- mzml_names_ENDO

# Plot MSn levels (only if MS1 and MS2 data are in same directory/file, then exchange MSn 0 for MS1 and MS2)
pdf(file="plots/ENDOMS1_msn_levels.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=16, family="Helvetica")
par(mfrow=c(2,1), mar=c(5,4,4,2), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
boxplot(mzml_msn_ENDO, main="Number of spectra")

model_boxplot <- boxplot(t(mzml_msn_ENDO[,2]), main="Number of MS1 spectra per sample", xaxt="n")
tick <- seq_along(model_boxplot$names)
axis(1, at=tick, labels=F)
text(tick, par("usr")[3]-par("usr")[3]/10, model_boxplot$names, adj=0, srt=270, xpd=T)
dev.off()



###----Import raw data----

# Import raw data as MSnbase object OnDiskMsnExp, segement for msLevel = 1
msd <- readMSData(files = paste(input_dir_MS1, MS1_ENDO_files[,], sep = ""), 
                  pdata = new("NAnnotatedDataFrame",pheno_data_ENDO), 
                  msLevel = 1,
                  mode = "onDisk")

# Import as XCMSnExp object for visual analysis
msd_XCMS <- readMSData(files = paste(input_dir_MS1, MS1_ENDO_files[,], sep = ""), 
                  pdata = new("NAnnotatedDataFrame",pheno_data_ENDO),
                  msLevel = 1)

table(msLevel(msd))
head(fData(msd)[, c("scanWindowLowerLimit", "scanWindowUpperLimit",
                             "originalPeaksCount", "msLevel", "retentionTime")])
write.csv(fData(msd), file="ENDO_raw_data_1.csv", row.names=FALSE)

# Restrict data to 1020 seconds (17 minutes)
msd <- filterRt(msd, c(0, 1020))

# Inspect mz values per file
msd_mz <- mz(msd)
msd_mz <- split(msd_mz, f = fromFile(msd))
print(length(msd_mz))

# Get base peak chromatograms
register(bpstart(SnowParam()))
# register(SerialParam())
chromas_ENDO <- chromatogram(msd_XCMS, 
                             aggregationFun="max", 
                             BPPARAM = SnowParam(workers = 2))


# Plot chromatograms based on phenodata groups
pdf(file="ENDO_chromas.pdf", encoding="ISOLatin1", pointsize=2, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(chromas_ENDO, main="Raw chromatograms", xlab="retention time [s]", ylab="intensity", col = col)
legend("topleft", bty="n", pt.cex=2, cex=1,5, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(col), legend= unique(msd_XCMS@phenoData@data[["sample_group"]]))
dev.off()



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
# Estimate precursor intensity
# ms1_data_neg is feature table, try with feature table from xcmsSet
# when take MS2 data??? not here for estimating precursor insensity?
# maybe create new XCMSnExp object for MS2 data later

#msd_XCMS_MS2 <- readMSData(files = paste(input_dir_MS2, MS2_ENDO_files[,], sep = ""), 
#                       pdata = NULL,
#                       msLevel = 2)


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
table(msLevel(msd_MS2))
head(fData(msd_MS2)[, c("precursorMZ", "precursorCharge", "precursorIntensity",
                         "precursorScanNum", "originalPeaksCount",
                        "scanWindowLowerLimit", "scanWindowUpperLimit", "msLevel", "retentionTime")])
write.csv(fData(msd_MS2), file="ENDO_raw_data_MS2.csv", row.names=FALSE)



# normal script from here on
precursor_intensity_ENDO <- estimatePrecursorIntensity(msd_MS2_2, method = "interpolation", 
                                                       BPPARAM = SnowParam(workers = 2))
print(head(na.omit(precursor_intensity_ENDO)))
# has only NAs??


# Reconstruct MS2 spectra from MS1 data
ms2_data_ENDO <- chromPeakSpectra(msd_XCMS_MS2, msLevel=2L, return.type="Spectra")
print(ms2_data_ENDO)
print(length(ms2_data_ENDO$peak_id))

# Extract all usable MS2 spectra
ms2_spectra_ENDO <- list()
#for (i in 1:nrow(ms1_def_neg)) {
ms2_spectra_neg <- foreach(i=1:nrow(ms1_def_neg)) %dopar% {
  # Extract existing MS2 spectra for feature
  feature_of_interest <- ms1_def_neg[i, "mzmed"]
  peaks_of_interest <- chromPeaks(ms1_data_neg, mz=feature_of_interest, ppm=ppm)
  if (length(which(ms2_data_neg$peak_id %in% rownames(peaks_of_interest))) > 0) {
    spectra_of_interest <- ms2_data_neg[ms2_data_neg$peak_id %in% rownames(peaks_of_interest)]
    combined_spectra_of_interest <- filterIntensity(spectra_of_interest, intensity=c(1,Inf), backend=MsBackendDataFrame)
    combined_spectra_of_interest <- setBackend(combined_spectra_of_interest, backend=MsBackendDataFrame())
    combined_spectra_of_interest <- Spectra::combineSpectra(combined_spectra_of_interest, FUN=combinePeaks, ppm=ppm, peaks="union", minProp=0.8, intensityFun=median, mzFun=median, backend=MsBackendDataFrame)#f=rownames(peaks_of_interest))
    combined_spectra_peaks <- as.data.frame(peaksData(combined_spectra_of_interest)[[1]])
    #Spectra::plotSpectra(combined_spectra_of_interest)
    #plot(x=combined_spectra_peaks[,1], y=combined_spectra_peaks[,2], type="h", xlab="m/z", ylab="intensity", main=paste("Precursor m/z",combined_spectra_of_interest@backend@spectraData$precursorMz[[1]]))
    #length(spectra_of_interest$peak_id)
    
    #ms2_spectra_neg[[i]] <- combined_spectra_of_interest
    return(combined_spectra_of_interest)
  }
}

# Remove empty spectra
names(ms2_spectra_neg) <- rownames(ms1_def_neg)
ms2_spectra_neg <- ms2_spectra_neg[lengths(ms2_spectra_neg) != 0]


###----Link MS2 and MS1-----

# Relate all MS2 spectra to MS1 precursors
ms1_def_neg$has_ms2 <- as.integer(rownames(ms1_def_neg) %in% names(ms2_spectra_neg))







###----with MS-Dial----
# C:\Users\abela\Downloads\MSDIAL ver.4.9.221218 Windowsx64\MSDIAL ver.4.9.221218 Windowsx64
# MS-Dial