### Statistics analysis combining the feature tables of EXO and ENDO
# Data from feat_list_ENDO and feat_list_EXO from the individual statistics files

###---- library ----
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
library(chemodiv)               # Chemodiversity (Petren 2022)
#library(rcdk)                   # CDK
#library(rinchi)                 # Converting SMILES to InchiKey
library(plotly)                 # For creating html plots
library(htmlwidgets)            # For creating html plots
#library(shiny)                  # HTML in R
#library(sunburstR)              # HTML-sunburst plots
library(heatmaply)              # HTML heatmaps
library(stringr)  
library(M3C)                    # t-SNE
library(randomForest)           # Random forest
library(MetaboAnalystR)         # Random forest
library(gplots)
#library(iESTIMATE)
source("https://raw.githubusercontent.com/ipb-halle/iESTIMATE/main/R/_functions.r")

########## Import objects ########

# import feature tables
feat_list_ENDO <- read.csv("ENDO_stats_Results/feat_list_ENDO.csv", row.names = 1)
feat_list_EXO <- read.csv("exo_stats_Results/feat_list_EXO.csv", row.names = 1)

# import binary tables
bina_list_ENDO <- read.csv("ENDO_stats_Results/bina_list_ENDO.csv", row.names = 1)
bina_list_EXO <- read.csv("exo_stats_Results/bina_list_EXO.csv", row.names = 1)


# name the features per polarity
colnames(feat_list_ENDO) <- paste(colnames(feat_list_ENDO), "endo", sep = "_")
colnames(feat_list_EXO) <- paste(colnames(feat_list_EXO), "exo", sep = "_")

# name the features per polarity
colnames(bina_list_ENDO) <- paste(colnames(bina_list_ENDO), "endo", sep = "_")
colnames(bina_list_EXO) <- paste(colnames(bina_list_EXO), "exo", sep = "_")

# combine the feature tables  into one
feat_list_COMBINED <- cbind(feat_list_ENDO, feat_list_EXO)

# combine the feature tables  into one
bina_list_COMBINED <- cbind(bina_list_ENDO, bina_list_EXO)


# remove MB
#feat_list_COMBINED <- feat_list_COMBINED[,1:24]

# rename the rows
rownames(feat_list_COMBINED) <- gsub("ENDO", "COMBINED", rownames(feat_list_COMBINED))

# rename the rows
rownames(bina_list_COMBINED) <- gsub("ENDO", "COMBINED", rownames(bina_list_COMBINED))


###----Creation of Phenodata/Metadata----
# create vector with sample classes according to culture information sheet
samp_groups <- c("CoCuPp", "CoCuSm", "CoCuSm", "CoCuPp", "CoCuPp", "CoCuSm",
                 rep(x = "Sm", times = 8), 
                 rep(x = "Pp", times = 8),
                 "CoCuSm", "CoCuPp", "MB")

CoCuPp_des <- "Co-culture sample from Prymnesium parvum"
CoCuSm_des <- "Co-culture sample from Skeletonema marinoi"
Sm_des <- "Mono-culture sample from Skeletonema marinoi"
Pp_des <- "Mono-culture sample from Prymnesium parvum"
MB_des <- "Media Blank"

# create vector with sample classes according to culture information sheet
samp_groups_description <- c(CoCuPp_des, CoCuSm_des, CoCuSm_des, CoCuPp_des, CoCuPp_des, CoCuSm_des,
                             rep(x = Sm_des, times = 8), 
                             rep(x = Pp_des, times = 8),
                             CoCuSm_des, CoCuPp_des, MB_des)
# create vector for colors
CoCuPp1 <- ("#3E134F")
CoCuSm1 <- rep("#F36E35", 2)
CoCuPp2 <- rep("#3E134F", 2)
CoCuSm2 <- ("#F36E35")
Sm <- rep("#F8B83C", 8)
Pp <- rep("#C53270", 8)
CoCuSm3 <- ("#F36E35")
CoCuPp3 <- ("#3E134F")
MB <- rep("#040404", 1)


color <- c(CoCuPp1, CoCuSm1, CoCuPp2, CoCuSm2, Sm, Pp, CoCuSm3, CoCuPp3, MB)

MS1_COMBINED_names <- rownames(feat_list_COMBINED)


# create phenodata based on culture type
pheno_data <- data.frame(sample_name = MS1_COMBINED_names, sample_group = samp_groups[1:24], samp_groups_description = samp_groups_description[1:24])
pheno_color <- data.frame(color)


###---- divide into species level ----
# index for the species types
index_SM <- grep("Sm", pheno_data$sample_group)
index_PP <- grep("Pp", pheno_data$sample_group)

# divide the feature table by species type
feat_list_COMBINED_Sm <- feat_list_COMBINED[index_SM,]
feat_list_COMBINED_Pp <- feat_list_COMBINED[index_PP,]

# divide the feature table by species type
bina_list_COMBINED_Sm <- bina_list_COMBINED[index_SM,]
bina_list_COMBINED_Pp <- bina_list_COMBINED[index_PP,]




###----Preparations----
# create plot directory
if (dir.exists(paste(getwd(), "/COMBINED_stats_plots/", sep = ""))){
  print("plots directory already exists")
}  else{
  dir.create("COMBINED_stats_plots")
  print("plots folder has been created")
}

# create Results directory
if (dir.exists(paste(getwd(), "/COMBINED_stats_Results/", sep = ""))){
  print("Results directory already exists")
}  else{
  dir.create("COMBINED_stats_Results")
  print("Results folder has been created")
}


# define variables by which to analyse
mzml_pheno_samples <- samp_groups_description[1:24]
mzml_pheno_colors <- color[1:24]
# By Culture type and Species
mzml_pheno_samples_type <- samp_groups[1:24]
# By Culture type
mzml_pheno_origin_samples <- as.factor(origin_samp_groups <- c("CoCu", "CoCu", "CoCu", "CoCu", "CoCu", "CoCu",
                                                               rep(x = "MonoCu", times = 8), rep(x = "MonoCu", times = 8),
                                                               "CoCu", "CoCu"))
# By Species
mzml_pheno_origin_species <- as.factor(species_samp_groups <- c("P. parvum", "S. marinoi", "S. marinoi", "P. parvum", "P. parvum", "S. marinoi",
                                                                rep(x = "S. marinoi", times = 8), rep(x = "P. parvum", times = 8),
                                                                "S. marinoi", "P. parvum"))
# For plot legends
mzml_pheno_legend <-  c("Co-culture P.p", "Co-culture S.m", "Co-culture S.m", "Co-culture P.p", "Co-culture P.p", "Co-culture S.m",
                        rep(x = "Mono-culture S.m", times = 8), rep(x = "Mono-culture P.p", times = 8),
                        "Co-culture S.m", "Co-culture P.p")

# overview dataframe
mzml_pheno <- data.frame(mzml_pheno_samples_type, mzml_pheno_origin_samples, mzml_pheno_origin_species, mzml_pheno_legend)


# comp list
comp_list <- feat_list_COMBINED
comp_list_Sm <- feat_list_COMBINED_Sm
comp_list_Pp <- feat_list_COMBINED_Pp

# bina list
bina_list <- bina_list_COMBINED
bina_list_Sm <- bina_list_COMBINED_Sm
bina_list_Pp <- bina_list_COMBINED_Pp

# ############################## MS1 statistics (ALL FEATURES) ##############################


# ---------- Diversity  ----------

# Create data frame
model_div             <- data.frame(features=apply(X=bina_list , MARGIN=1, FUN=function(x) { sum(x) } ))
model_div$richness    <- apply(X=bina_list , MARGIN=1, FUN=function(x) { sum(x) } )
model_div$menhinick   <- apply(X=bina_list , MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
model_div$shannon     <- apply(X=comp_list , MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
model_div$pielou      <- apply(X=scale(comp_list , center=FALSE, scale=FALSE), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#model_div$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list, species), index="chao")
model_div$simpson     <- apply(X=comp_list , MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
model_div$inverse     <- apply(X=comp_list , MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
model_div$fisher      <- apply(X=comp_list , MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
#model_div$unique      <- apply(X=uniq_list, MARGIN=1, FUN=function(x) { sum(x) })
#model_div$hillfunc    <- as.numeric(unlist(calcDiv(comp_list, compDisMat=scales::rescale(as.matrix(dist(t(comp_list)), diag=TRUE, upper=TRUE)), q=1, type="FuncHillDiv")))
#rownames(model_div) <- paste(pheno_data_COMBINED[1:24,1], pheno_data_COMBINED[1:24,2], sep = "_" )

# test for functional hill
# subset comp_list for the first 20 features
#comp_list_hillfunc <- comp_list[,1:1000]
#model_div$hillfunc    <- as.numeric(unlist(calcDiv(comp_list_hillfunc, compDisMat=scales::rescale(as.matrix(dist(t(comp_list_hillfunc)), diag=TRUE, upper=TRUE)), q=1, type="FuncHillDiv")))

# save as model_div as csv
write.csv(model_div, file=paste(filename = "COMBINED_stats_Results/model_div_COMBINED.csv", sep = ""))

# color for shannon plot
col_shan <- c(unique(CoCuPp1), unique(CoCuSm1), unique(CoCuPp1), unique(CoCuSm1))

# Plot Shannon index
pdf(paste("COMBINED_stats_plots/ms1_comp_list_diversity_shannon_Culture_1.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=10, height=10, family="Helvetica")
boxplot(model_div$shannon ~ mzml_pheno_samples_type, col=col_shan, main="Shannon diversity (H\')", xlab="culture types", ylab="Shannon diversity index (H\')")
text(1:length(levels(mzml_pheno_samples_type)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=levels(mzml_pheno_samples_type), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$shannon, term=as.factor(mzml_pheno_samples_type))
#text(1:length(levels(mzml_pheno_samples_type)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
legend("bottomleft", bty="n", pt.cex=5, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(col_shan), legend= c("P.parvum", "S.marinoi"))
dev.off()

# ---------- PLS ----------
principal_components <- 5
# # Cross validation for number of components
# library(pls)
# set.seed(123)
# 
# # Define the number of PLS components to test
# ncomp <- 1:10
# 
# # Define the number of folds for cross-validation
# nfolds <- 10
# 
# # Create a matrix with the response variable (Class) and the selected features from the PLS analysis
# pls_data <- cbind(Class = mzml_pheno_origin_samples, comp_list[, sel_pls_comp_list$`_selected_variables_`])
# 
# # Perform cross-validation for each number of PLS components
# cv_results <- list()
# for (i in ncomp) {
#   # Create a PLS model with i components
#   pls_model <- pls(Class ~ ., data = pls_data, ncomp = i)
#   
#   # Perform cross-validation
#   cv <- cv.pls(pls_model, ncomp = i, method = "crossval", nfolds = nfolds, verbose = FALSE)
#   
#   # Store the cross-validation results
#   cv_results[[i]] <- cv$mserror
# }
# 
# # Plot the cross-validation results
# plot(ncomp, unlist(cv_results), type = "b", xlab = "Number of PLS components", ylab = "Cross-validation error")


# PLS according to Culture type on Species level Skeletonema marinoi
sel_pls_comp_list <- f.select_features_pls(feat_matrix=comp_list_Sm, 
                                           sel_factor=mzml_pheno_origin_samples[index_SM], 
                                           sel_colors=mzml_pheno_colors[index_SM], 
                                           components=principal_components, tune_length=10, 
                                           quantile_threshold=0.95, plot_roc_filename="COMBINED_stats_plots/ms1_comp_list_select_pls_roc_Sm.pdf")
print(paste("Number of selected variables:", f.count.selected_features(sel_feat=sel_pls_comp_list$`_selected_variables_`)))
jpeg(filename = "COMBINED_stats_plots/COMBINED_comp_list_pls_Sm.jpeg", width = 1000, height = 1800, quality = 100, bg = "white")
par(mar=c(8,6,4,3), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=2, cex.main=2)
f.heatmap.selected_features(feat_list=comp_list_Sm, 
                            sel_feat=sel_pls_comp_list$`_selected_variables_`, 
                            sel_names=paste0("",sel_pls_comp_list$`_selected_variables_`), 
                            sample_colors=mzml_pheno_colors[index_SM], plot_width=7, plot_height=7, 
                            #cex_col=0.5, cex_row=0.4, 
                            filename=NULL, main="PLS")
#text(ms1_pca_COMBINED$x[,1], ms1_pca_COMBINED$x[,2], labels=str_sub(ms1_data_COMBINED_pos$sample_name[1:24], - 3, - 1), col=color, pos=3, cex=1.5)
dev.off()
heatmaply(scale(comp_list_Sm[, which(colnames(comp_list_Sm) %in% sel_pls_comp_list$`_selected_variables_`)]), 
          k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256),
          #label = comp_list[25,],
          
          file="COMBINED_stats_plots/ms1_comp_list_select_pls_Sm.html", selfcontained=TRUE)
sel_pls_comp_list$`_selected_variables_`
sel_pls_comp_list$`_model_r2_`
sel_pls_comp_list$`_multiclass_metrics_`

save(sel_pls_comp_list, file = "COMBINED_stats_Results/sel_pls_comp_list_culture_Sm.RData")

# PLS according to Culture type on Species level Prymnesium parvum
sel_pls_comp_list <- f.select_features_pls(feat_matrix=comp_list_Pp, 
                                           sel_factor=mzml_pheno_origin_samples[index_PP], 
                                           sel_colors=mzml_pheno_colors[index_PP], 
                                           components=principal_components, tune_length=10, 
                                           quantile_threshold=0.95, plot_roc_filename="COMBINED_stats_plots/ms1_comp_list_select_pls_roc_Pp.pdf")
print(paste("Number of selected variables:", f.count.selected_features(sel_feat=sel_pls_comp_list$`_selected_variables_`)))
jpeg(filename = "COMBINED_stats_plots/COMBINED_comp_list_pls_Pp.jpeg", width = 1000, height = 1800, quality = 100, bg = "white")
par(mar=c(8,6,4,3), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=2, cex.main=2)
f.heatmap.selected_features(feat_list=comp_list_Pp, 
                            sel_feat=sel_pls_comp_list$`_selected_variables_`, 
                            sel_names=paste0("",sel_pls_comp_list$`_selected_variables_`), 
                            sample_colors=mzml_pheno_colors[index_PP], plot_width=7, plot_height=7, 
                            #cex_col=0.5, cex_row=0.4, 
                            filename=NULL, main="PLS")
#text(ms1_pca_COMBINED$x[,1], ms1_pca_COMBINED$x[,2], labels=str_sub(ms1_data_COMBINED_pos$sample_name[1:24], - 3, - 1), col=color, pos=3, cex=1.5)
dev.off()
heatmaply(scale(comp_list_Pp[, which(colnames(comp_list_Pp) %in% sel_pls_comp_list$`_selected_variables_`)]), 
          k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256),
          #label = comp_list[25,],
          
          file="COMBINED_stats_plots/ms1_comp_list_select_pls_Pp_complete_comp_list.html", selfcontained=TRUE)
sel_pls_comp_list$`_selected_variables_`
sel_pls_comp_list$`_model_r2_`
sel_pls_comp_list$`_multiclass_metrics_`

save(sel_pls_comp_list, file = "COMBINED_stats_Results/sel_pls_comp_list_culture_Pp.RData")

#check if selected features are annotated
# pls_feature_neg <- c("FT18839", "FT18078", "FT17890", "FT13643", "FT04103", "FT09428", "FT14727", "FT14411", "FT07576", "FT02408")
# matching_features <- match(pls_feature_neg, COMBINED_neg_annotated$ft_id_x)


# PLS including both factors Species and Culture type
# sel_pls_comp_list <- f.select_features_pls(feat_matrix=comp_list, 
#                                            sel_factor=as.factor(c(mzml_pheno_origin_samples, mzml_pheno_origin_species)), # include species parameter
#                                            sel_colors=mzml_pheno_colors,
#                                            components=principal_components, tune_length=10, 
#                                            quantile_threshold=0.95, plot_roc_filename="COMBINED_stats_plots/ms1_comp_list_select_pls_roc_SPCU.pdf")
# print(paste("Number of selected variables:", f.count.selected_features(sel_feat=sel_pls_comp_list$`_selected_variables_`)))
# f.heatmap.selected_features(feat_list=comp_list, 
#                             sel_feat=sel_pls_comp_list$`_selected_variables_`, 
#                             sel_names=paste0("",sel_pls_comp_list$`_selected_variables_`), 
#                             sample_colors=mzml_pheno_colors, plot_width=7, plot_height=7, cex_col=0.5, cex_row=0.4, filename=NULL, main="PLS")
# heatmaply(scale(comp_list[, which(colnames(comp_list) %in% sel_pls_comp_list$`_selected_variables_`)]), 
#           k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), 
#           file="COMBINED_stats_plots/ms1_comp_list_select_pls_SPCU.html", selfcontained=TRUE)
# sel_pls_comp_list$`_selected_variables_`
# sel_pls_comp_list$`_model_r2_`
# sel_pls_comp_list$`_multiclass_metrics_`
# 
# save(sel_pls_comp_list, file = "COMBINED_stats_Results/sel_pls_comp_list_SPCU.RData")

# ---------- PCA  ----------
# combined feature table pos and neg
jpeg(filename = "COMBINED_stats_plots/COMBINED_ms1_feature_table_pca.jpeg", width = 1000, height = 700, quality = 100, bg = "white")
ms1_pca_COMBINED <- prcomp(comp_list[1:24,], center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms1_pca_COMBINED$x[, 1], ms1_pca_COMBINED$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms1_pca_COMBINED)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_COMBINED)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color, cex=2)
grid()
text(ms1_pca_COMBINED$x[,1], ms1_pca_COMBINED$x[,2], col=color, pos=3, cex=1.5)
legend("topleft", bty="n", pt.cex=5, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(mzml_pheno$mzml_pheno_samples_type))
dev.off()

# Skeletonema marinoi
jpeg(filename = "COMBINED_stats_plots/COMBINED_ms1_feature_table_pca_Sm.jpeg", width = 1000, height = 700, quality = 100, bg = "white")
ms1_pca_COMBINED <- prcomp(comp_list_Sm, center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms1_pca_COMBINED$x[, 1], ms1_pca_COMBINED$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms1_pca_COMBINED)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_COMBINED)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color[index_SM], cex=2)
grid()
text(ms1_pca_COMBINED$x[,1], ms1_pca_COMBINED$x[,2], col=color[index_SM], labels=str_sub(MS1_COMBINED_names[1:24], - 3, - 1), pos=3, cex=1.5)
legend("topleft", bty="n", pt.cex=5, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color[index_SM]), legend= unique(mzml_pheno$mzml_pheno_samples_type[index_SM]))
dev.off()

# Prymnesium parvum
jpeg(filename = "COMBINED_stats_plots/COMBINED_ms1_feature_table_pca_Pp.jpeg", width = 1000, height = 700, quality = 100, bg = "white")
ms1_pca_COMBINED <- prcomp(comp_list_Pp, center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms1_pca_COMBINED$x[, 1], ms1_pca_COMBINED$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms1_pca_COMBINED)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_COMBINED)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color[index_PP], cex=2)
grid()
text(ms1_pca_COMBINED$x[,1], ms1_pca_COMBINED$x[,2], col=color[index_PP], labels=str_sub(ms1_data_ENDO_neg$sample_name[1:24], - 3, - 1), pos=3, cex=1.5)
legend("topleft", bty="n", pt.cex=5, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color[index_PP]), legend= unique(mzml_pheno$mzml_pheno_samples_type[index_PP]))
dev.off()



