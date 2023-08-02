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

# ------ data_preparation_stat --------
data_preparation_stat <- function()