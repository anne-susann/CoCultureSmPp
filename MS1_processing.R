## Co-Culture Analysis: MS1 data pre-processing
## XCMS package for analysis
## xcmsSet object

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

###----ENDO files----
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

#for (i in 1:length(phenoData(xs))){
#  if (grepl("1A|1b|2a|2b|3a|3b|4a|4b", phenoData(xs)) == TRUE) {
#     sampclass(xs) <- c("Sm")
#  } else if (grepl("5a|5b|6a|6b|7a|7b|8a|8b", phenoData(xs)) == TRUE){
#    sampclass(xs) <- c("Pp")
#  } else if (grepl("MB", phenoData(xs)) == TRUE){
#    sampclass(xs) <- c("MB")
#  } else {
#    sampclass(xs) <- c("CoCu")
#  }
#}

# until loop runs for class naming in either xcmsSet or pd1, do it manuallly
sampclass(xs) <- c("CoCuPp", "CoCuSm", "CoCuSm", "CoCuPp", "CoCuPp", "CoCuSm",
                   rep(x = "Sm", times = 8), 
                   rep(x = "Pp", times = 8),
                   "CoCuSm", "CoCuPp", "MB")

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

###----EXO files-----
# create a list of all MS1 EXO metabolome files in data folder
MS1_EXO_files <- data.frame(list.files(input_dir, pattern = "EXO"))
View(MS1_EXO_files)

###------pre-processing----

## create xcmsSet object 

# for phenodata would make sense to add column with class name
pdx <- data.frame(file = basename(paste(input_dir, MS1_EXO_files[,], sep = "")))
pdx$class <- c("NA")

# creating an xcmsSet object, including peak identification with given parameters
# xcmsSet class does peak detection, peak grouping, nonlinear retention time correction
# missing value imputation, can generate extracted ion chromatograms
xsx <- xcmsSet(files = paste(input_dir, MS1_EXO_files[,], sep = ""), 
              phenoData = pdx, 
              method = "matchedFilter",
              fwhm = 20, 
              max = 50, 
              snthresh = 4, 
              step = 0.05,
              BPPARAM = SnowParam(workers = 2))

# Need to define sampclass to identify the different classes of data, eg Pp, Sm, CoCu
# loop runs but doesn't add vector/name to correct file, names all of the the same

#for (i in 1:length(phenoData(xsx))){
#  if (grepl("1A|1b|2a|2b|3a|3b|4a|4b", phenoData(xsx)) == TRUE) {
#     sampclass(xsx) <- c("Sm")
#  } else if (grepl("5a|5b|6a|6b|7a|7b|8a|8b", phenoData(xsx)) == TRUE){
#    sampclass(xsx) <- c("Pp")
#  } else if (grepl("MB", phenoData(xsx)) == TRUE){
#    sampclass(xsx) <- c("MB")
#  } else {
#    sampclass(xsx) <- c("CoCu")
#  }
#}

# until loop runs for class naming in either xcmsSet or pd1, do it manuallly
sampclass(xsx) <- c("CoCuPp", "CoCuSm", "CoCuSm", "CoCuPp", "CoCuPp", "CoCuSm",
                   rep(x = "Sm", times = 8), 
                   rep(x = "Pp", times = 8),
                   "CoCuSm", "CoCuPp", "MB")

# summary of xcmsSet
xsx

# xs peak table
xsx@peaks
# better view
View(xsx@peaks)

# overview: plot a spectra
# relative intensity and mz or rt
plot(x = xsx@peaks[1:1928,"rt"], y = xsx@peaks[1:1928,"intf"], col=xsx@peaks[,"sample"])

## peaks are now detected by matched filter function and listed in xs@peaks list
## sample says which origin file/sample the peak comes from

###----post-processing------
## Assignment in groups, now peaks that represent the same metabolite need to be grouped together

# peaks are grouped together to features according to m/z ranges and retention time
# peak detection method was density
xsgx <- group(object = xsx)

# summary
xsgx
xsgx@groups[1:5, ]

# rt correction using Loess filter
xsgx <- retcor(object = xsgx)

# re-grouping
xsgx <- group(object = xsgx)
xsgx@groups[1:5, ]


## Missing value imputation
# First test by outputting first values of dataset
groupval(object = xsgx, value = "into")[1:5, 1:10]

# now missing value imputation to remove NA values, check with groupval() for NA
xsgfx <- fillPeaks(object = xsgx, BPPARAM = MulticoreParam(workers = 2))
# groupval of xsgf feature list of the peaks grouped across the samples n = nrow(input_files)
# the samples are the colums 
groupval(object = xsgfx, value = "into")[1:5, 1:25]

# Quality control of data to see major outliers or gaussian shape of distribution
# save as picture
plotQC(xsgfx, what = "mzdevhist")
plotQC(xsgfx, what = "rtdevhist")


resultsx <- groupval(xsgfx, "medret", "into")


###----save and re-import calculated data-----


# Save results as .csv file 
write.csv(results, file = "MS1_ENDO_results.csv") 
write.csv(resultsx, file = "MS1_EXO_results.csv")


# save peak tables as csv
write.csv(peakTable(xsgf), file = "ENDO_peakTable.csv")
write.csv(peakTable(xsgfx), file = "EXO_peakTable.csv")


# ENDO: re-import data from results.csv file for statistical analysis
data.pca <- read.csv("MS1_ENDO_results.csv")

# identify duplicated values and rename 
dup <- anyDuplicated(data.pca[,1])
data.pca[dup,1] <- paste(data.pca[dup,1], "x", sep = "")

# set first row as rownames
# subset without first row
# set column names to numbers 1-25
rownames(data.pca) <- data.pca[,1]
data.pca <- data.pca[, -1]
colnames(data.pca) <- seq.int(1, 25)


# EXO: re-import data from results.csv
data.pcax <- read.csv("MS1_EXO_results.csv")

# identify duplicated values and rename 
dupx <- anyDuplicated(data.pcax[,1])
data.pcax[dupx,1] <- paste(data.pcax[dupx,1], "x", sep = "")

# set first row as rownames
# subset without first row
# set column names to numbers 1-25
rownames(data.pcax) <- data.pcax[,1]
data.pcax <- data.pcax[, -1]
colnames(data.pcax) <- seq.int(1, 25)



###----PCA----


# get intensity values of peak groups 
# same data as results table
# either import from csv or recalculate
data.pca <- groupval(object = xsgfx, value = "into")
# defines the symbols to be drawn in the plot
symb <- c(rep(0,1), rep(1,2), rep(0,2), rep(1,1), 
          rep(2,8), rep(3,8), 
          rep(1,1), rep(0,1))
# execute PCA with scaled data
pc <- prcomp(x = t(na.omit(data.pca)), scale. = T)

# define colors for classes
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

# check loadings and score data
scree.data <- as.data.frame(pc$importance)
score.data <- as.data.frame(pc$x)
loadings.data <- as.data.frame(pc$rotation)

# collect results for Plot
pcSummary <- summary(pc)

png("PCA_ENDO.png", width=10, height=6, units="in", res=100)
plot(x = -1*pc$x[, 1], y = pc$x[,2], pch = symb, main = "PCA: S. marinoi, P. parvum, Co-Culture", 
     xlab = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100,
                                   digits = 3), " % variance"),
     ylab = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100,
                                   digits = 3), " % variance"),
     col = col, bg = col, cex = 1.5, xlim = c(-40, 40), ylim = c(-40, 40))
abline(h = 0, v = 0, col = "black")
legend("topright", col = unique(col), legend = levels(sampclass(xsgf)), 
       pch = unique(symb))
dev.off()

###----Significance of features----

# perform broken stick test to evaluate which components are statistically important
# modified after this code https://github.com/mohanwugupta/Machine-Learning-Neuroimaging-Tutuorial/blob/master/Broken_Stick.R


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

ev_pc = pc$sdev^2  
evplot(ev_pc)  
dev.off()



###----other analysis----
# model_div_ENDO <- data.frame(rownames(data.pca))
model_div_ENDO <- data.pca
# shannon diversity index
model_div_ENDO$shannon <- apply(X=data.pca, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })























