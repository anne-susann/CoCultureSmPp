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

### MS1_workflow_functions
source("C:/Users/abela/Documents/GitHub/CoCultureSmPp/MS1_workflow_XCMS/MS1_workflow_functions.R")

# ---- analysis for ENDO neg ----
# defining variable for functions to use
raw_data_MS1_ENDO_neg <- list.files("C:/Users/abela/Documents/Uni_Jena/Masterarbeit/MAW-Co-culture/MS1_pos_neg", pattern = "neg_ENDO", full.names = TRUE)
phenodata <- read.csv("phenodata_ENDO_neg.csv")
result_dir_name <- "ENDO_neg_testing_results"
plots_dir_name <- "ENDO_neg_testing_plots"
ms1_params_ENDO_neg <- CentWaveParam(ppm=25, mzCenterFun="wMean", peakwidth=c(14, 59), 
                                     prefilter=c(3, 140), mzdiff=0.0155, snthresh=7, noise=0, 
                                     integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
                                     fitgauss=FALSE, roiList=list(), roiScales=numeric())

### trial run data_preparation
msd <- data_preparation(files = raw_data_MS1_ENDO_neg, phenodata = phenodata, result_dir_name = "ENDO_neg_testing_results", 
                        plots_dir_name = "ENDO_neg_testing_plots", chrom_run_sec = 700, msLevel = 1)

### trial run chromatogram
chrom_msd <- chromatogram_qc(msd, phenodata, result_dir_name, plots_dir_name)

### trial run peak_detection
ms1_data <- peak_detection(msd, CentWaveParam = ms1_params_ENDO_neg, result_dir_name, plots_dir_name)

### trial run grouping_1
ms1_data <- grouping_1(ms1_data, PeakDensityParam_gr = NULL, result_dir_name)

### trial run rt_correction
ms1_data <- rt_correction(ms1_data, plots_dir_name, result_dir_name, PeakGroupsParam_rt = NULL, chromas_msd = chrom_msd)

### trial run grouping_2
ms1_data <- grouping_2(ms1_data, PeakDensityParam_gr = NULL, result_dir_name)

### trial run feature_extraction
feature_list <- feature_extraction(ms1_data, result_dir_name)

### trial run feature_transformation
feature_list <- feature_transformation(feature_list, result_dir_name, plots_dir_name)

### trial run bina_list_creation
bina_list <- bina_list_creation(ms1_matrix = NULL, intensity_cutoff = 14, result_dir_name)

# ---- analysis for ENDO pos ----
raw_data_MS1_ENDO_pos <- list.files("C:/Users/abela/Documents/Uni_Jena/Masterarbeit/MAW-Co-culture/MS1_pos_neg", pattern = "pos_ENDO", full.names = TRUE)
phenodata <- read.csv("phenodata_ENDO_pos.csv")
result_dir_name <- "ENDO_pos_testing_results"
plots_dir_name <- "ENDO_pos_testing_plots"
ms1_params_ENDO_pos <- CentWaveParam(ppm=25, mzCenterFun="wMean", peakwidth=c(12, 51), 
                                     prefilter=c(4, 60), mzdiff= 0.000099, snthresh=6, noise=0, 
                                     integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
                                     fitgauss=FALSE, roiList=list(), roiScales=numeric())


msd <- data_preparation(files = raw_data_MS1_ENDO_pos, phenodata = phenodata, result_dir_name, plots_dir_name, chrom_run_sec = 700, msLevel = 1)

chrom_msd <- chromatogram_qc(msd, phenodata, result_dir_name, plots_dir_name)

ms1_data <- peak_detection(msd, CentWaveParam = ms1_params_ENDO_pos, result_dir_name, plots_dir_name)

ms1_data <- grouping_1(ms1_data, PeakDensityParam_gr = NULL, result_dir_name)

ms1_data <- rt_correction(ms1_data, plots_dir_name, result_dir_name, PeakGroupsParam_rt = NULL, chromas_msd = chrom_msd)

ms1_data <- grouping_2(ms1_data, PeakDensityParam_gr = NULL, result_dir_name)

feature_list <- feature_extraction(ms1_data, result_dir_name)

feature_list <- feature_transformation(feature_list, result_dir_name, plots_dir_name)

bina_list <- bina_list_creation(ms1_matrix = NULL, intensity_cutoff = 14, result_dir_name)

