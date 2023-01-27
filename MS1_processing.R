## Co-Culture Analysis: MS1 data pre-processing
## XCMS package for analysis
## xcmsSet object

###-------library-----

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

###------parallelization----

library(parallel)
n.cores <- detectCores()

# Set up parallel processing using 2 cores
if (.Platform$OS.type == "unix") {
  register(bpstart(MulticoreParam(2)))
} else {
  register(bpstart(SnowParam(2)))
}



###------set directory----

# set data directory for MS1 data files
input_dir <- paste(getwd(), "/MS1/", sep = "")
input_dir

# create a list of all MS1 ENDO metabolome files in data folder
MS1_ENDO_files <- data.frame(list.files(input_dir, pattern = "ENDO"))
View(MS1_ENDO_files)

###------pre-processing----

## create xcmsSet object 

# for phenodata would make sense to add column with class name
pd <- data.frame(file = basename(paste(input_dir, MS1_ENDO_files[,], sep = "")))
pd$class <- c("NA")

# creating an xcmsSet object, including peak identification with given parameters
# xcmsSet class does peak detection, peak grouping, nonlinear retention time correction
# missing value imputation, can generate extracted ion chromatograms
xs <- xcmsSet(files = paste(input_dir, MS1_ENDO_files[,], sep = ""), 
              phenoData = pd, 
              method = "matchedFilter",
              fwhm = 20, 
              max = 50, 
              snthresh = 4, 
              step = 0.05,
              BPPARAM = SnowParam(workers = 2))

# Need to define sampclass to identify the different classes of data, eg Pp, Sm, CoCu
# loop runs but doesn't add vector/name to correct file, names all of the the same

for (i in 1:length(phenoData(xs))){
  if (grepl("1A|1b|2a|2b|3a|3b|4a|4b", phenoData(xs)) == TRUE) {
    sampclass(xs) <- c("Sm")
  } else if (grepl("5a|5b|6a|6b|7a|7b|8a|8b", phenoData(xs)) == TRUE){
    sampclass(xs) <- c("Pp")
  } else if (grepl("MB", phenoData(xs)) == TRUE){
    sampclass(xs) <- c("MB")
  } else {
    sampclass(xs) <- c("CoCu")
  }
}

# until loop runs for class naming in either xcmsSet or pd1, do it manuallly
sampclass(xs) <- c(rep(x = "CoCu", times = 6), 
                   rep(x = "Sm", times = 8), 
                   rep(x = "Pp", times = 8),
                   rep(x = "CoCu", times = 2),
                   "MB")

# summary of xcmsSet
xs

# xs peak table
xs@peaks
# better view
View(xs@peaks)

# overview: plot a spectra
# relative intensity and mz or rt
plot(x = xs@peaks[1:1928,"rt"], y = xs@peaks[1:1928,"intf"], col=xs@peaks[,"sample"])

## peaks are now detected by matched filter function and listed in xs@peaks list
## sample says which origin file/sample the peak comes from

###----post-processing------
## Assignment in groups, now peaks that represent the same metabolite need to be grouped together

# peaks are grouped together to features according to m/z ranges and retention time
# peak detection method was density
xsg <- group(object = xs)

# summary
xsg
xsg@groups[1:5, ]

# rt correction using Loess filter
xsg <- retcor(object = xsg)

# re-grouping
xsg <- group(object = xsg)
xsg@groups[1:5, ]

## Missing value imputation
# First test by outputting first values of dataset
groupval(object = xsg, value = "into")[1:5, 1:10]

# now missing value imputation to remove NA values, check with groupval() for NA
xsgf <- fillPeaks(object = xsg, BPPARAM = MulticoreParam(workers = 2))
# groupval of xsgf feature list of the peaks grouped across the samples n = nrow(input_files)
# the samples are the colums 
groupval(object = xsgf, value = "into")[1:5, 1:25]

# Quality control of data to see major outliers or gaussian shape of distribution
# save as picture
plotQC(xsgf, what = "mzdevhist")
plotQC(xsgf, what = "rtdevhist")

# Get peak intensity matrix, groupvalue shows for each sample each group as rt/mz correlation
results <- groupval(xsgf, "medret", "into")
resultsx <- groupval(xsgfx, "medret", "into")


###----save and re-import calculated data-----


# Save results as .csv file 
write.csv(results, file = "MS1_ENDO_results.csv", 
          quote = TRUE, row.names = TRUE, col.names = TRUE)

write.csv(resultsx, file = "MS1_EXO_results.csv")


# save peak tables as csv
write.csv(peakTable(xsgf), file = "ENDO_peakTable.csv")
write.csv(peakTable(xsgfx), file = "EXO_peakTable.csv")



###----PCA----


# get intensity values of peak groups 
# same data as results table
data.pca <- groupval(object = xsgf, value = "into")
# defines the symbols to be drawn in the plot
symb <- c(rep(0, 6), rep(1, 8), rep(2, 8), rep(0, 2), rep(3, 1))
# execute PCA with scaled data
pc <- prcomp(x = t(na.omit(data.pca)), scale. = T)

# define colors for classes
CoCu <- rep("royalblue4", 6)
Sm <- rep("violetred", 8)
Pp <- rep("yellow2", 8)
CoCu1 <- rep("royalblue4", 2)
MB <- rep("springgreen", 1)

col <- c(CoCu, Sm, Pp, CoCu1, MB)

# check loadings and score data
scree.data <- as.data.frame(pc$importance)
score.data <- as.data.frame(pc$x)
loadings.data <- as.data.frame(pc$rotation)

# collect results for Plot
pcSummary <- summary(pc)

pdf("PCA_ENDO.pdf")
plot(x = -1*pc$x[, 1], y = pc$x[,2], pch = symb, main = "PCA: S. marinoi, P. parvum, Co-Culture", 
     xlab = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100,
                                   digits = 3), " % variance"),
     ylab = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100,
                                   digits = 3), " % variance"),
     col = col, bg = col, cex = 1.5, xlim = c(-40, 40), ylim = c(-40, 40))
abline(h = 0, v = 0, col = "black")
legend("topleft", col = unique(col), legend = levels(sampclass(xsgf)), 
       pch = unique(symb))
dev.off()



