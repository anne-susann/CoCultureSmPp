## Co-Culture Analysis: MS1 data pre-processing
## XCMS package for analysis
## testing different objects and methods

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

## ------parallelization----

library(parallel)
n.cores <- detectCores()

# Set up parallel processing using 2 cores
if (.Platform$OS.type == "unix") {
  register(bpstart(MulticoreParam(2)))
} else {
  register(bpstart(SnowParam(2)))
}



## ------set directory----

# set data directory for MS1 data files
input_dir <- paste(getwd(), "/MS1/", sep = "")
input_dir

# create a list of all MS1 ENDO metabolome files in data folder
MS1_ENDO_files <- data.frame(list.files(input_dir, pattern = "ENDO"))
View(MS1_ENDO_files)


## ------pre-processing----

# store all raw data for MS1 ENDO files in one OnDiskMsnExp object (smaller memory usage) 
raw_data <- readMSData(files = paste(input_dir, MS1_ENDO_files[,], sep = ""),
                       pdata = NULL,
                       mode = "onDisk")

# store mz values for all spectra in mzs, independent of source file/sample all together
mzs <- mz(raw_data)

# Split the list by file
mzs_by_file <- split(mzs, f = fromFile(raw_data))
length(mzs_by_file)

###----Testing with Tutorials-----

### works until here, then chromatogram function fails, follows this tutorial
### https://www.bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html 

# ## Get the base peak chromatograms. This reads data from the files.
bpis <- chromatogram(raw_data_1, aggregationFun = "max")
## Define colors for the two groups
group_colors <- paste0(brewer.pal(3, "Set1")[1:2], "60")
names(group_colors) <- c("KO", "WT")

## Plot all chromatograms.
plot(bpis, col = group_colors[raw_data_1$sample_group])

bpi_1 <- bpis[1, 1]
head(rtime(bpi_1))

head(intensity(bpi_1))

## Get the total ion current by file
tc <- split(tic(raw_data_1), f = fromFile(raw_data_1))
boxplot(tc, col = group_colors[raw_data_1$sample_group],
        ylab = "intensity", main = "Total ion current")


## Bin the BPC
bpis_bin <- MSnbase::bin(bpis, binSize = 1)

## Calculate correlation on the log2 transformed base peak intensities
cormat <- cor(log2(do.call(cbind, lapply(bpis_bin, intensity))))
colnames(cormat) <- rownames(cormat) <- raw_data_1$sample_name

## Define which phenodata columns should be highlighted in the plot
ann <- data.frame(group = raw_data_1$sample_group) # does not work somehow
rownames(ann) <- raw_data_1$sample_name

## Perform the cluster analysis
pheatmap(cormat, annotation = raw_data_1$sample_group,
         annotation_color = list(group = group_colors))




### trying another tutorial 
### https://jorainer.github.io/metabolomics2018/xcms-preprocessing.html 


## first read data in (also following Blockpraktikum Script)
# pb is phennodata, additional data supplied in the files
pd <- data.frame(file = basename(paste(input_dir, MS1_ENDO_files[,], sep = "")))
                 
data <- readMSData(files = paste(input_dir, MS1_ENDO_files[,], sep = ""),
                   pdata = new("NAnnotatedDataFrame", pd),
                   mode = "onDisk") 

# retention time splitted by files
rts <- split(rtime(data), fromFile(data))
length(rts)

# subset the data with filter functions, eg filterFile, filterRtime
#' Use %>% to avoid nested function calls
#' spectra creates a list of Spectrum objects, not sorted by file
sps <- data %>%
  spectra 

length(sps)

plot(sps[[1000]]) # shows spectra of specific rt, in increasing order 

pData(data)

chr <- chromatogram(data)
chr

group_colors <- paste0(brewer.pal(3, "Set1")[1:2], "60")
names(group_colors) <- c("KO", "WT")
# Plot all chromatograms in one plot
plot(chr)

# extract total ion intensities from the TIC
ints <-  intensity(chr[1,1])
head(ints)

# also still contains all phenotype information like original OnDiskMSnExp object
pData(chr)

###----Pre-processing Blockpraktikum Protocol-----

### try following Blockpraktikum ### so far so good, works
## create xcmsSet object 

# for phenodata would make sense to add column with class name
pd1 <- data.frame(file = basename(paste(input_dir, MS1_ENDO_files[,], sep = "")))
pd1$class <- c("NA")

# add species name to corresponding row of file name, but how????
for (i in 1:nrow(pd1)){
  
  if (grepl("1A|1b|2a|2b|3a|3b|4a|4b", pd1[i,]) == TRUE) {
    pd1[i,2] <- "Sm"
  } else if (grepl("5a|5b|6a|6b|7a|7b|8a|8b", pd1[i,]) == TRUE){
    pd1[i,"class"] <- "Pp"
  } else if (grepl("MB", pd1[i,]) == TRUE){
    pd1[i, 2] <- "MB"
  } else {
    pd1[i, 2] <- "CoCu"
  }
}



# store pheno data of the assay and the sample names
pd <- read.AnnotatedDataFrame(MS1_ENDO_files[1,], 
                              input_dir, 
                              row.names = NULL, blank.lines.skip = TRUE, 
                              fill = TRUE, varMetadata.char = "$", quote = "\"")
#sampleNames(pd) <- pd$col_names

# creating an xcmsSet object, including peak identification with given parameters
# xcmsSet class does peak detection, peak grouping, nonlinear retention time correction
# missing value imputation, can generate extracted ion chromatograms
xs <- xcmsSet(files = paste(input_dir, MS1_ENDO_files[,], sep = ""), 
              phenoData = pd1, 
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
plot(x = xs@peaks[1:1928,"rt"], y = xs@peaks[1:1928,"mz"], col=xs@peaks[,"sample"])

xs@rt[["raw"]][[1]]

## peaks are now detected by matched filter function and listed in xs@peaks list
## sample says which origin file/sample the peak comes from

###----Post-processing------
## Assignment in groups, now peaks that represent the same metabolite need to be grouped together

# peaks are grouped together to features according to m/z ranges and retention time
# peak detection method was density
xsg <- group(object = xs)

# try peak detection method nearest neighbor
xsg_1 <- group(object = xs, method = "nearest")

# summary
xsg
xsg_1

# xsg_1 has triple the amount of peaks gouped together, better or worse?

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
groupval(object = xsgf, value = "into")[1:5, 1:10]

# Step 6. Get peak intensity matrix, groupvalue shows for each sample each group as rt/mz correlation
results <- groupval(xsgf, "medret", "into")
resultsx <- groupval(xsgfx, "medret", "into")


# Step 7. Save results as .csv file 
write.csv(results, file="MS1_ENDO_results.csv")

### important to reformat into table to export and import data


# Quality control of data to see major outliers or gaussian shape of distribution
# save as picture
plotQC(xsgf, what = "mzdevhist")
plotQC(xsgf, what = "rtdevhist")

# mit farbe funktioniert irgendwie nicht so 
MB <- rep("royalblue4", 8)
Sm <- rep("violetred", 8)
Pp <- rep("yellow2", 8)
CoCu <- rep("springgreen", 8)

col <- c(CoCu, Sm, Pp, CoCu)
plotQC(xsgf, sampColors = col, what = "mzdevsample")
plotQC(xsgf, sampColors = col, what = "rtdevsample")

###----saving pre-processed data-----
# save plots in working directory
pdf(file = "rtQC_MS1.pdf")
plotQC(xsgf, what = "rtdevhist")
dev.off()

# read groupval results in, no need for new calculations
# for statistical analysis just import results file to not repeat calculations
t <- read.csv(file = "MS1_ENDO_results.csv")
rnam <- c(ENDO$X)
resultsx <- ENDO[, -1]
rownames(ENDO, do.NULL = TRUE, prefix = rnam)

resultsx <- read.table(file = "MS1_EXO_resultsx.txt", 
                       header = TRUE, sep = "\t", dec = ".", fill = FALSE,
                       row.names = rnam)

# how to set first row as row names and column names as simple numbers
rownames(resultsx) <- resultsx[, 1]  ## set rownames
resultsx <- resultsx[, -1]           ## remove the first variable

# set column names
colnames(resultsx) <- seq.int(1, 25)






# save xcmsSet xsgf as mzTab format, gives an overview of the results
mzt <- data.frame(character(0))
mzt <- xcms:::mzTabHeader(mzt,
                          version="1.1.0", mode="Complete", type="Quantification",
                          description="MS1_ENDO",
                          xset=xsgf)
mzt <- xcms:::mzTabAddSME(mzt, xsgf)
xcms:::writeMzTab(mzt, "MS1_ENDO_processed.mzTab")

# doesn't save group information and/or mz information
# also present in MetaboAnalystR package


# save peakTable for ENDO and EXO with groups


write.csv(peakTable(xsgf), file = "ENDO_peakTable.csv")
write.csv(peakTable(xsgfx), file = "EXO_peakTable.csv")


chr_endo <- chromatogram(xsgf) # not possible

xset <- as(xsgf, "XCMSnExp") # not possible

# problem with importing data before was a duplicated rt/mz value
# could not be set as row because double name

ENDO_1 <- read.csv("MS1_ENDO_results_test.csv")

# identify duplicated value and rename 
dup <- anyDuplicated(ENDO_1[,1])
ENDO_1[dup,1] <- paste(ENDO_1[dup,1], "x", sep = "")

# set first row as rownames
# subset without first row
# set column names to numbers 1-25
rownames(ENDO_1) <- ENDO_1[,1]
ENDO_1 <- ENDO_1[, -1]
colnames(ENDO_1) <- seq.int(1, 25)





###----Try XCMSnExp object instead of xcmsSet----
### change workflow from xcmsSet object to XCMSnExp object (follow up version)
# https://bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcms.html
# all these methods only work for XCMSnExp objects and not for xcmsSet

getdata2(
  path = input_dir,
  index = F,
  snames = NULL,
  sclass = pd$class,
  phenoData = pd,
  BPPARAM = BiocParallel::SnowParam(workers = 2),
  mode = "onDisk",
  ppp = xcms::CentWaveParam(ppm = 5, peakwidth = c(5, 25), prefilter = c(3, 5000)),
  rtp = xcms::ObiwarpParam(binSize = 1),
  gpp = xcms::PeakDensityParam(sampleGroups = 1, minFraction = 0.67, bw = 2, binSize = 0.025),
  fpp = xcms::FillChromPeaksParam()
)

# ppp use MatchedFilterParam
# rtp use Loess function
# gpp is same, density function

# run failed again, looks like error from OnDiskMsnExp??






###----PCA----
prcomp()


###----1. approach to PCA----

# Processing and Visualization of Metabolomics Data Using R http://dx.doi.org/10.5772/65405

## descriptive statistics for the variables of the data
sumstats <- function(z) { 
  Mean <- apply(z, 1, mean) 
  Median <- apply(z, 1, median) 
  SD <- apply(z, 1, sd) 
  SE <- apply(z, 1, function(x) sd(x)/sqrt(length(x))) 
  CV <- apply(z, 1, function(x) sd(x)/mean(x)) 
  result <- data.frame(Mean, Median, SD, SE, CV) 
  return(result) 
  }
sumstats(results)

replacezero <- function(x) "[<-"(x, !x | is.na(x), min(x[x > 0], 
  na.rm = TRUE) / 2) 
newdata <- apply(results, 1, replacezero)
newdata_xsgf <- apply(peakTable(xsgf), 1, replacezero)
peakTable(xsgf)


paretoscale <- function(z) {
  rowmean <- apply(z, 1, mean) # row means
  rowsd <- apply(z, 1, sd) # row standard deviation
  rowsqrtsd <- sqrt(rowsd) # sqrt of sd
  rv <- sweep(z, 1, rowmean,"-") # mean center
  rv <- sweep(rv, 1, rowsqrtsd, "/") # divide by sqrtsd
  return(rv)
}

# why do we use this function? Normalization and scaling?
logdata <- log(newdata, 2)
pareto.logdata <- paretoscale(logdata)
# here rows equal samples, but where is group data?

# normalization and scaling has already been done, so set to false here
pca <- prcomp(t(pareto.logdata), center = FALSE, scale = FALSE)
pcaresults <- summary(pca)
pcaresults


scree.data <- as.data.frame(pcaresults$importance)
score.data <- as.data.frame(pcaresults$x)
loadings.data <- as.data.frame(pcaresults$rotation)
write.csv(scree.data, "pca_scree.csv") 
write.csv(score.data, "pca_scores.csv") 
write.csv(loadings.data, "pca_loadings.csv")

datapca <- read.csv("pca_scores.csv", header=T) 
datapca3<- datapca[, c(1:3)] # subset columns 1-3

# not the right association, just 
groupvec <- c(rep(x = c("Sm", "Pp", "CoCu", "MB"), times = nrow(datapca3)))

datapca3 <- cbind(datapca3, groupvec)
# how add names of samples???? information lost somewhere between xcmsSet and here
# and how to do this with other object types???


ggplot(datapca3, aes(PC1, PC2)) + 
  geom_point(aes(shape=groupvec)) + 
  geom_text(aes(label=datapca3$X)) + 
  stat_ellipse(aes(fill=groupvec))








###----2. approach to PCA-----

## follow practical course PCA and statistical analysis
# n times p dimensional space, n = 25 samples, p = 2777 peak groups (features)
# also using prcomp and xcmsSet object, so similar to 1. approach



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



###----Significance of features----

# using a diffreport to show the biggest differences between two groups
# Creation of the diffreport regarding control vs. group 1 
# and additional saving in diffreport_XCMS.tsv
dr <- diffreport(object = xsgf, class1 = "Sm", class2 = "CoCu", 
                 filebase = "diffreport_XCMS", sortpval = FALSE)

BiocManager::install("multtest")

# first 10 significant values 
pval<- dr$pvalue(which(dr$pvalue < 0.05))
pval [1:10]

# 10 most significant features
drpv<- dr[order(dr$pvalue,decreasing = FALSE) , ]
drpv[1:10,]

## diffreport
dr <- diffreport(object = xsgfx, class1 = "Pp", class2 = "CoCuPp", 
                 filebase = "diffreport_XCMS", sortpval = FALSE)

# check for signifantly different features
sigxPp <- which(dr$pvalue<=0.05)
length(sigPp) # Pp against Co-culture
length(sigSm) # Sm against Co-culture
length(sigSP) # Pp against Sm

length(sigxPp)


## perform broken stick test to evaluate which components are statistically important

png("BrokenStick_EXO.png", width=10, height=6, units="in", res=100)
evplot = function(ev) {  
  # Broken stick model (MacArthur 1957)  
  n = length(ev)  
  bsm = data.frame(j=seq(1:n), p=0)  
  bsm$p[1] = 1/n  
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))  
  bsm$p = 100*bsm$p/n  
  # Plot eigenvalues and % of variation for each axis  
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))  
  barplot(ev, main="Eigenvalues EXO", col="blue", las=2)  
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

# --> works, 3 PCs seem to be statistically relevant
# pc$sdev is the standard 


# broken stick determine statistically relevant pcs
brokenStick()

library(PCDimension)
AuerGervini(pc$sdev, dd=NULL, epsilon = 2e-16)

brokenStick(pc$sdev[1:25], length(pc$sdev))
# bad values error??


## varpart
var_pc <- vegan::varpart(pc)

data(mite)
data(mite.env)
data(mite.pcnm)

## See detailed documentation:
vegandocs("partition")

# Two explanatory matrices -- Hellinger-transform Y
# Formula shortcut "~ ." means: use all variables in 'data'.
mod <- varpart(mite, ~ ., mite.pcnm, data=mite.env, transfo="hel")
mod
showvarparts(2)
plot(mod)
# Alterna














## random forest
library(randomForest)
str(data.pcax)

data <- data.pcax


set.seed(222)
ind <- sample(2, nrow(data.pcax), replace = TRUE, prob = c(0.7, 0.3))
train <- data.pcax[ind==1,]
test <- data.pcax[ind==2,]

rf <- randomForest(data.pcax$`3`~., data=train, proximity=TRUE) 
print(rf)
#Call:
#  randomForest(formula = Species ~ ., data = train)
#Type of random forest: classification
#Number of trees: 500
#No. of variables tried at each split: 2
#OOB estimate of  error rate: 2.83%





###----diversity shannon-----

# feature list
View(data.pca)
# model_div_ENDO <- data.frame(rownames(data.pca))
model_div_ENDO <- data.pca
# shannon diversity index
model_div_ENDO$shannon <- apply(X=data.pca, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })













###----INCLUSION LIST WITH MS2-----
## MS DIAL



spss <- Spectra(paste(input_dir, MS1_ENDO_files[1,], sep = ""), backend = MsBackendMzR())
spss











