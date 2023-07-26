###---- library ----
library(RColorBrewer)           # For colors
library(MSnbase)                # MS features
library(xcms)                   # Swiss army knife for metabolomics
library(CAMERA)                 # Metabolite Profile Annotation
library(Spectra)                # Spectra package needed for XCMS3
library(vegan)                  # For shannon diversity
library(multcomp)               # For Tukey test
library(Hmisc)                  # For correlation test
library(gplots)                 # For fancy heatmaps
library(circlize)               # For sunburst plot
library(plotrix)                # For sunburst plot
library(caret)                  # Swiss-army knife for statistics
library(pROC)                   # Evaluation metrics
library(PRROC)                  # Evaluation metrics
library(multiROC)               # Evaluation metrics
library(plotly)                 # For creating html plots
library(htmlwidgets)            # For creating html plots
library(stringr)              
source("https://raw.githubusercontent.com/ipb-halle/iESTIMATE/main/R/_functions.r")


# ---- data preparation ----
data_preparation <- function(files, phenodata, result_dir_name, plots_dir_name, chrom_run_sec, msLevel) {
  # files = list of mzML files
  # phenodata = csv file of phenodata, including (1) sample_name, (2) sample_group, (3) sample_description
  # result_dir_name = name for directory to store all result files and objects
  # plots_dir_name = name for directory to store all plots 
  # chrom_run_sec = length of the chromatography run in seconds
  # msLevel = level of MS for files, default is ms level 1
  
  # create results directory
  if (dir.exists(paste(getwd(), result_dir_name, sep = ""))){
    print("results directory already exists")
  }  else{
    dir.create(result_dir_name)
    print("results folder has been created")
  }
  
  # create plot directory
  if (dir.exists(paste(getwd(), plots_dir_name, sep = ""))){
    print("plots directory already exists")
  }  else{
    dir.create(plots_dir_name)
    print("plots folder has been created")
  }
  
  # read MS data from mzML file list
  msd <- readMSData(files = files,
                    pdata = new("NAnnotatedDataFrame", phenodata),
                    mode = "onDisk",
                    centroided = TRUE)
  
  # restrict data to length of chromatography run
  msd <- filterRt(msd, c(0, chrom_run_sec))
  
  # restrict data to MS1
  msd <- filterMsLevel(msd, msLevel = msLevel)
  
  # create csv of raw data
  write.csv(fData(msd), file=paste(result_dir_name, filename = "_raw_data.csv", sep = ""), row.names=FALSE)
  
  ### create quality control
  # create color palette
  color_pal <- brewer.pal(n = 9, name = "Set1")  
  
  # create color vector from phenodata
  uniq_sample_group <- data.frame(unique(phenodata$sample_group))
  uniq_sample_group$color <- color_pal[1:nrow(uniq_sample_group)]
  
  for (i in 1:nrow(uniq_sample_group)) {
    for (j in 1:nrow(phenodata)) {
      if (uniq_sample_group[i,1] %in% phenodata$sample_group[j]){
        phenodata$color[j] <- uniq_sample_group$color[i]
        
      }
    }
  }
  
  color <- phenodata$color
  
  # # create chromatogram object for quality control plots
  # chromas <- chromatogram(msd, 
  # aggregationFun="max", 
  # msLevel = 1)
  # 
  # # Plot base peak chromatograms based on phenodata groups
  # jpeg(filename = paste(plots_dir_name, "/chromas_bpc.jpeg", sep = ""), width = 2000, height = 1200, quality = 100, bg = "white")
  # par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=2, cex.lab=2, cex.main=2)
  # plot(chromas, main="Raw chromatograms", xlab="retention time [s]", ylab="intensity", col = color)
  # legend("topleft", bty="n", pt.cex=3, cex=1.5, y.intersp=0.7, text.width=0.5, pch=20, 
  #        col= unique(color), legend= unique(msd@phenoData@data[["sample_group"]]))
  # dev.off()
  # 
  # # Plot total ion current chromatogram
  # jpeg(filename = paste(plots_dir_name, "/chromas_tic.jpeg", sep = ""), width = 2000, height = 1200, quality = 100, bg = "white")
  # par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=2, cex.lab=2, cex.main=2)
  # tics <- split(tic(msd), f=fromFile(msd))
  # boxplot(tics, col=color, ylab="intensity", xlab="sample", main="Total ion current", outline = FALSE)
  # legend("topleft", bty="n", pt.cex=2, cex=1, y.intersp=0.7, text.width=0.5, pch=20, 
  #        col= unique(color), legend= unique(msd@phenoData@data[["sample_group"]]))
  # dev.off()
  # 
  # # Qugrouping/binning for peak detection, based on similarity of the base peak chromatogram 
  # chromas_bin <- MSnbase::bin(chromas, binSize=2)
  # chromas_bin_cor <- cor(log2(do.call(cbind, lapply(chromas, intensity)))) # transformation with log
  # colnames(chromas_bin_cor) <- rownames(chromas_bin_cor) <- msd$sample_name
  # chromas_bin_cor[is.na(chromas_bin_cor)] <- 0
  # 
  # # representing the data in a heatmap for general overview
  # jpeg(filename = paste(plots_dir_name, "/heatmap_chromas_bin.jpeg", sep = ""), width = 500, height = 500, quality = 100, bg = "white")
  # par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(7,0,0,7), cex.axis=0.9, cex=0.6)
  # heatmap(chromas_bin_cor)
  # dev.off()
  
  save(msd, file = paste(result_dir_name, "/msd.RData", sep = ""))
  return(msd)
  
}

### trial run data_preparation
# files = raw_data_MS1_ENDO_neg
raw_data_MS1_ENDO_neg <- list.files("C:/Users/abela/Documents/Uni_Jena/Masterarbeit/MAW-Co-culture/MS1_pos_neg", pattern = "neg_ENDO", full.names = TRUE)
# phenodata = phenodata
phenodata <- read.csv("phenodata_ENDO_neg.csv")
# result_dir_name = ENDO_neg_testing_results
# plots_dir_name = ENDO_neg_testing_plots
# chrom_run_sec = 700
# msLevel = 1


data_preparation(files = raw_data_MS1_ENDO_neg, phenodata = phenodata, result_dir_name = "ENDO_neg_testing_results", 
                 plots_dir_name = "ENDO_neg_testing_plots", chrom_run_sec = 700, msLevel = 1)

# ---- peak_detection -----
peak_detection <- function(msd, CentWaveParam) {
  # msd = prepared ms object from data_preparation function
  # CentWaveParam = CentWaVeParam class for parameters for peak detection algorithm
  # 
  
  # define peak detection parameters, if default
  
  # perform peak detection
  
  # create summary file .csv
  
  # create quality control plots
  
  # save detected object
  
  return(ms1_data)
}


# ----grouping_1 ----
grouping_1 <- function(ms1_data) {
  
}

# ----rt_correction ----
rt_correction <- function(ms1_data) {
  
  
}

# ----grouping_2----
grouping_2 <- function(ms1_data) {
  

}

# ----feature_extraction ----
feature_extraction <- function(ms1_data) {
  
}

# ----feature_transformation----
feature_transformation <- function(ms1_data) {
  
}

# ----bina_list_creation-----
bina_list_creation <- function(feature_table) {
  
}