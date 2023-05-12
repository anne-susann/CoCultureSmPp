### Statistics analysis combining the feature tables of EXO pos and EXO neg
# Data from MS1_EXO_neg and MS1_EXO_pos scripts

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
library(randomForest)           # Random forest
library(MetaboAnalystR)         # Random forest
library(gplots)
#library(iESTIMATE)
source("https://raw.githubusercontent.com/ipb-halle/iESTIMATE/main/R/_functions.r")

########## Import objects ########

# import feature tables
feature_table_pos <- read.csv("exo_pos_Results/feature_list_EXO_pos.csv", row.names = 1)
# load(file = "exo_stats_Results/MS_exo_pos_peak_detection.RData")
feature_table_neg <- read.csv("exo_neg/exo_neg_Results/feature_list_EXO_neg.csv", row.names = 1)
# load(file = "exo_neg_Results/MS_exo_pos_peak_detection.RData")

# name the features per polarity
colnames(feature_table_pos) <- paste(colnames(feature_table_pos), "pos", sep = "_")
colnames(feature_table_neg) <- paste(colnames(feature_table_neg), "neg", sep = "_")

# combine the feature tables of neg and pos into one
feat_list_EXO <- cbind(feature_table_pos, feature_table_neg)

# import binary tables
binary_table_pos <- read.csv("exo_pos_Results/bina_list_EXO_pos.csv", row.names = 1)
binary_table_neg <- read.csv("exo_neg/exo_neg_Results/bina_list_EXO_neg.csv", row.names = 1)

# name the features per polarity
colnames(binary_table_pos) <- paste(colnames(binary_table_pos), "pos", sep = "_")
colnames(binary_table_neg) <- paste(colnames(binary_table_neg), "neg", sep = "_")

# combine the feature tables of neg and pos into one
bina_list_EXO <- cbind(binary_table_pos, binary_table_neg)

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


# create vector with colors 
colors <- rainbow(4)
colors

CoCuPp1 <- ("#FF0000")
CoCuSm1 <- rep("#009900", 2)
CoCuPp2 <- rep("#FF0000", 2)
CoCuSm2 <- ("#009900")
Sm <- rep("#00CCFF", 8)
Pp <- rep("#8000FF", 8)
CoCuSm3 <- ("#009900")
CoCuPp3 <- ("#FF0000")
MB <- rep("#040404", 1)

color <- c(CoCuPp1, CoCuSm1, CoCuPp2, CoCuSm2, Sm, Pp, CoCuSm3, CoCuPp3, MB)

# create vector with symbols for plotting
CoCuPp1_s <- (15)
CoCuSm1_s <- rep(17, 2)
CoCuPp2_s <- rep(15, 2)
CoCuSm2_s <- (17)
Sm_s <- rep(18, 8)
Pp_s <- rep(19, 8)
CoCuSm3_s <- (17)
CoCuPp3_s <- (15)
MB_s <- rep(4, 1)

pch_s <- c(CoCuPp1_s, CoCuSm1_s, CoCuPp2_s, CoCuSm2_s, Sm_s, Pp_s, CoCuSm3_s, CoCuPp3_s, MB_s)


# file names
MS1_EXO_names <- rownames(feat_list_EXO)

# create phenodata based on culture type
pheno_data_EXO <- data.frame(sample_name = MS1_EXO_names, sample_group = samp_groups, samp_groups_description = samp_groups_description)
pheno_color_EXO <- data.frame(color)

###----Preparations----
# create plot directory
if (dir.exists(paste(getwd(), "/exo_stats_plots/", sep = ""))){
  print("plots directory already exists")
}  else{
  dir.create("exo_stats_plots")
  print("plots folder has been created")
}

# create Results directory
if (dir.exists(paste(getwd(), "/exo_stats_Results/", sep = ""))){
  print("Results directory already exists")
}  else{
  dir.create("exo_stats_Results")
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


# remove media blank sample from analysis
feat_list_EXO <- feat_list_EXO[-25,]
bina_list_EXO <- bina_list_EXO[-25,]

# create objets for Statistical analysis
# comp_list 
comp_list <- feat_list_EXO
comp_list <- comp_list[] # removes missing values
rownames(comp_list) <- paste(pheno_data_EXO[1:24,1], pheno_data_EXO[1:24,2], sep = "_" )
row_names <- rownames(comp_list)
new_row_names <- gsub("KSS_210324_", "", row_names)
rownames(comp_list) <- new_row_names
#colnames(comp_list) <- paste0(rownames(ms1_def_EXO_pos)[which(ms1_def_EXO_pos$has_ms2==1)], "_pos")

# bina_list
bina_list <- bina_list_EXO
rownames(bina_list) <- paste(pheno_data_EXO[1:24,1], pheno_data_EXO[1:24,2], sep = "_" )
#colnames(bina_list) <- paste0(rownames(ms1_def_EXO_pos)[which(ms1_def_EXO_pos$has_ms2==1)], "_pos")

# uniq_list
#uniq_list <- uniq_list
#rownames(uniq_list) <- paste(pheno_data_EXO[1:24,1], pheno_data_EXO[1:24,2], sep = "_" )
#length(table(which(uniq_list == 1)))

# ############################## Integrate annotation in feature table ##############################
# load annotated MS2 features
EXO_neg_annotated <- data.frame(read.csv(paste(getwd(), "/exo_neg/mergedResults-with-one-Candidates.csv", sep = ""), row.names = 1))

# remove NA, meaning remove any unannotated features
EXO_neg_annotated <- EXO_neg_annotated[EXO_neg_annotated$AnnotationCount >= 1,]
EXO_neg_annotated <- EXO_neg_annotated[rowSums(is.na(EXO_neg_annotated)) != ncol(EXO_neg_annotated),]

# add by feature name to featuer table
feature_table_neg <- feature_table_neg[-25,]

# add row for annotation
annotation <- c(rep.int(0, length(feature_table_neg)))
feature_table_neg <- rbind(feature_table_neg, annotation)

# search for feature id from annotation file in colnames of feature table
# if they match, add annotation info to additonal row
matching_col <- match(EXO_neg_annotated$ft_id_x, colnames(feature_table_neg))
matching_vals <- EXO_neg_annotated$IUPAC
for (i in seq_along(matching_col)) {
  if (!is.na(matching_col[i])) {
    feature_table_neg[25, matching_col[i]] <- matching_vals[i]
  }
}
print(feature_table_neg[25, matching_col[1:5]])
# look 
#feature_table_neg$FT00621
#feature_table_neg$FT00628

#comp_list <- feature_table_neg

# ############################## MS1 statistics (ALL FEATURES) ##############################


# ---------- Diversity  ----------

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
#model_div$unique      <- apply(X=uniq_list, MARGIN=1, FUN=function(x) { sum(x) })
#model_div$hillfunc    <- as.numeric(unlist(calcDiv(comp_list, compDisMat=scales::rescale(as.matrix(dist(t(comp_list)), diag=TRUE, upper=TRUE)), q=1, type="FuncHillDiv")))
#rownames(model_div) <- paste(pheno_data_EXO[1:24,1], pheno_data_EXO[1:24,2], sep = "_" )

# test for functional hill
# subset comp_list for the first 20 features
comp_list_hillfunc <- comp_list[,1:1000]
model_div$hillfunc    <- as.numeric(unlist(calcDiv(comp_list_hillfunc, compDisMat=scales::rescale(as.matrix(dist(t(comp_list_hillfunc)), diag=TRUE, upper=TRUE)), q=1, type="FuncHillDiv")))


# Plot Shannon index
pdf(paste("exo_stats_plots/ms1_comp_list_diversity_shannon_Culture.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$shannon ~ mzml_pheno_origin_samples, col=mzml_pheno_colors, main="Shannon diversity (H\')", xlab="culture types", ylab="Shannon diversity index (H\')")
text(1:length(levels(mzml_pheno_origin_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=levels(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$shannon, term=as.factor(mzml_pheno_origin_samples))
text(1:length(levels(mzml_pheno_origin_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
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


# PLS according to Culture type
sel_pls_comp_list <- f.select_features_pls(feat_matrix=comp_list, 
                                           sel_factor=mzml_pheno_origin_samples, 
                                           sel_colors=mzml_pheno_colors, 
                                           components=principal_components, tune_length=10, 
                                           quantile_threshold=0.95, plot_roc_filename="exo_stats_plots/ms1_comp_list_select_pls_roc_culture_2.pdf")
print(paste("Number of selected variables:", f.count.selected_features(sel_feat=sel_pls_comp_list$`_selected_variables_`)))
jpeg(filename = "exo_stats_plots/EXO_comp_list_pls.jpeg", width = 1000, height = 1800, quality = 100, bg = "white")
par(mar=c(8,6,4,3), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=2, cex.main=2)
f.heatmap.selected_features(feat_list=comp_list, 
                            sel_feat=sel_pls_comp_list$`_selected_variables_`, 
                            sel_names=paste0("",sel_pls_comp_list$`_selected_variables_`), 
                            sample_colors=mzml_pheno_colors, plot_width=7, plot_height=7, 
                            #cex_col=0.5, cex_row=0.4, 
                            filename=NULL, main="PLS")
#text(ms1_pca_EXO$x[,1], ms1_pca_EXO$x[,2], labels=str_sub(ms1_data_EXO_pos$sample_name[1:24], - 3, - 1), col=color, pos=3, cex=1.5)
dev.off()
heatmaply(scale(comp_list[, which(colnames(comp_list) %in% sel_pls_comp_list$`_selected_variables_`)]), 
          k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256),
          #label = comp_list[25,],
          
          file="exo_stats_plots/ms1_comp_list_select_pls_cuture_2.html", selfcontained=TRUE)


heatmaply(scale(comp_list[, which(colnames(comp_list) %in% sel_pls_comp_list$`_selected_variables_`)]), 
          k_row=1, k_col=1, 
          colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256),
          file="exo_stats_plots/ms1_comp_list_select_pls_cuture_3.html", 
          selfcontained=TRUE,
          layout_opts = list(xaxis = list(tickfont = list(size = 20)), 
                             yaxis = list(tickfont = list(size = 20, color = ms1_data_EXO_pos$sample_color),
                                          titlefont = list(size = 25, color = ms1_data_EXO_pos$sample_color))))


sel_pls_comp_list$`_selected_variables_`
sel_pls_comp_list$`_model_r2_`
sel_pls_comp_list$`_multiclass_metrics_`



save(sel_pls_comp_list, file = "exo_stats_Results/sel_pls_comp_list_culture.RData")

#check if selected features are annotated
pls_feature_neg <- c("FT18839", "FT18078", "FT17890", "FT13643", "FT04103", "FT09428", "FT14727", "FT14411", "FT07576", "FT02408")
matching_features <- match(pls_feature_neg, EXO_neg_annotated$ft_id_x)


# PLS including both factors Species and Culture type
# sel_pls_comp_list <- f.select_features_pls(feat_matrix=comp_list, 
#                                            sel_factor=as.factor(c(mzml_pheno_origin_samples, mzml_pheno_origin_species)), # include species parameter
#                                            sel_colors=mzml_pheno_colors,
#                                            components=principal_components, tune_length=10, 
#                                            quantile_threshold=0.95, plot_roc_filename="exo_stats_plots/ms1_comp_list_select_pls_roc_SPCU.pdf")
# print(paste("Number of selected variables:", f.count.selected_features(sel_feat=sel_pls_comp_list$`_selected_variables_`)))
# f.heatmap.selected_features(feat_list=comp_list, 
#                             sel_feat=sel_pls_comp_list$`_selected_variables_`, 
#                             sel_names=paste0("",sel_pls_comp_list$`_selected_variables_`), 
#                             sample_colors=mzml_pheno_colors, plot_width=7, plot_height=7, cex_col=0.5, cex_row=0.4, filename=NULL, main="PLS")
# heatmaply(scale(comp_list[, which(colnames(comp_list) %in% sel_pls_comp_list$`_selected_variables_`)]), 
#           k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), 
#           file="exo_stats_plots/ms1_comp_list_select_pls_SPCU.html", selfcontained=TRUE)
# sel_pls_comp_list$`_selected_variables_`
# sel_pls_comp_list$`_model_r2_`
# sel_pls_comp_list$`_multiclass_metrics_`
# 
# save(sel_pls_comp_list, file = "exo_stats_Results/sel_pls_comp_list_SPCU.RData")

# ---------- PCA for combined feature tables ----------
jpeg(filename = "exo_stats_plots/EXO_ms1_feature_table_pca_exc_MB.jpeg", width = 1000, height = 700, quality = 100, bg = "white")
ms1_pca_EXO <- prcomp(comp_list[1:24,], center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms1_pca_EXO$x[, 1], ms1_pca_EXO$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms1_pca_EXO)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_EXO)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color, cex=4)
grid()
#text(ms1_pca_EXO$x[,1], ms1_pca_EXO$x[,2], col=color, pos=3, cex=1.5)
legend("bottomleft", bty="n", pt.cex=5, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(mzml_pheno$mzml_pheno_legend))
dev.off()


# ---------- Variation partitioning ----------
# 
library(vegan)
model_varpart  <- varpart(scale(comp_list), ~ mzml_pheno_samples , ~ mzml_pheno_samples )
length(comp_list)
nrow(comp_list)
summary(mzml_pheno_samples)
summary(mzml_pheno_origin_samples)
summary(mzml_pheno_origin_species)



# Plot results
pdf(file="exo_stats_plots/pos_ms1_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart , Xnames=c("samples","culture type"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()
 
# --------- ANOVA ----------
model_anova <- aov(comp_list, ~ mzml_pheno_origin_samples , ~ mzml_pheno_origin_species, data = comp_list )

# -------- Random forest -----
# Load the data
comp_list

# Split the data into training and testing sets
set.seed(123) # for reproducibility
train_idx <- sample(1:nrow(comp_list), nrow(comp_list) * 0.8)
train <- comp_list[train_idx, ]
test <- comp_list[-train_idx, ]

# Fit a random forest model
rf_model <- randomForest(Class ~ ., data = train, ntree = 500, importance = TRUE)

# Evaluate the model on the test set
rf_pred <- predict(rf_model, test, type = "class")
conf_mat <- table(test$Class, rf_pred)
accuracy <- sum(diag(conf_mat)) / sum(conf_mat)

# Assess feature importance
var_imp <- rf_model$importance

# Plot variable importance
var_imp_rank <- sort(var_imp, decreasing = TRUE)
jpeg(filename = "exo_stats_plots/EXO _ms1_random_forest.jpeg", width = 1000, height = 700, quality = 150, bg = "white")
barplot(var_imp_rank, las = 2, cex.names = 0.7, main = "Random Forest Variable Importance")
dev.off()



############################## MS1 statistics (containing MS2 spectra) ##############################
# 
# mzml_pheno_samples <- samp_groups_description
# mzml_pheno_colors <- color
# principal_components <- 5
# 
# # remove MB from lists for analysis
# feat_list_nMB <- feat_list_EXO_pos[-(25),]
# bina_list_nMB <- bina_list_EXO_pos[-(25),]
# 
# # creates lists of MS1 features that have ms2 spectra 
# # comp_list
# comp_list <- feat_list_EXO_pos_nMB[, c(rownames(ms1_def_EXO_pos)[which(ms1_def_EXO_pos$has_ms2==1)])]
# comp_list <- comp_list[]
# colnames(comp_list) <- paste0(rownames(ms1_def_EXO_pos)[which(ms1_def_EXO_pos$has_ms2==1)], "_pos")
# 
# # bina_list
# bina_list <- feat_list_EXO_pos_nMB[, c(rownames(ms1_def_EXO_pos)[which(ms1_def_EXO_pos$has_ms2==1)])]
# colnames(bina_list) <- paste0(rownames(ms1_def_EXO_pos)[which(ms1_def_EXO_pos$has_ms2==1)], "_pos")
# 
# # uniq_list
# uniq_list <- uniq_list_EXO_pos[, c(rownames(ms1_def_EXO_pos)[which(ms1_def_EXO_pos$has_ms2==1)])]
# rownames(uniq_list_EXO_pos) <- gsub(x=rownames(uniq_list_EXO_pos), pattern="\\.auto.*", replacement="")
# colnames(uniq_list_EXO_pos) <- paste0(rownames(ms1_def_EXO_pos)[which(ms1_def_EXO_pos$has_ms2==1)], "_pos")
# 
# # Merge pos and neg div_classes
# #div_classes <- div_classes_exo_pos
# #div_classes_samples <- div_classes_samples_exo_pos
# 
# # Merge pos and neg div_superclasses
# #div_superclasses <- div_superclasses_exo_pos
# #div_superclasses_samples <- div_superclasses_samples_exo_pos
# 
# # Merge pos and neg div_subclasses
# #div_subclasses <- div_subclasses_exo_pos
# #div_subclasses_samples <- div_subclasses_samples_exo_pos
# 
# # class_list
# #class_list <- class_list_exo_pos
# #rownames(class_list) <- gsub(x=rownames(class_list_exo_pos), pattern="\\.auto.*", replacement="")
# 
# # class_int_list
# #class_int_list <- class_int_list_exo_pos
# #rownames(class_int_list) <- gsub(x=rownames(class_int_list_exo_pos), pattern="\\.auto.*", replacement="")
# 
# # superclass_list
# #superclass_list <- superclass_list_exo_pos
# #rownames(superclass_list) <- gsub(x=rownames(superclass_list_exo_pos), pattern="\\.auto.*", replacement="")
# 
# # superclass_int_list
# #superclass_int_list <- superclass_int_list_exo_pos
# #rownames(superclass_int_list) <- gsub(x=rownames(superclass_int_list_exo_pos), pattern="\\.auto.*", replacement="")
# 
# # subclass_list
# #subclass_list <- subclass_list_exo_pos
# #rownames(subclass_list) <- gsub(x=rownames(subclass_list_exo_pos), pattern="\\.auto.*", replacement="")
# 
# # subclass_int_list
# #subclass_int_list <- subclass_int_list_exo_pos
# #rownames(subclass_int_list) <- gsub(x=rownames(subclass_int_list_exo_pos), pattern="\\.auto.*", replacement="")
# 
# #cdk_descriptors <- cdk_descriptors_exo_pos
# 
# 
# ## only on features with MS2 spectra found
# # Create data frame
# model_div             <- data.frame(features=apply(X=bina_list, MARGIN=1, FUN=function(x) { sum(x) } ))
# model_div$richness    <- apply(X=bina_list, MARGIN=1, FUN=function(x) { sum(x) } )
# model_div$menhinick   <- apply(X=bina_list, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
# model_div$shannon     <- apply(X=comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
# model_div$pielou      <- apply(X=scale(comp_list, center=FALSE, scale=FALSE), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
# #model_div$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list, species), index="chao")
# model_div$simpson     <- apply(X=comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
# model_div$inverse     <- apply(X=comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
# model_div$fisher      <- apply(X=comp_list, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
# model_div$unique      <- apply(X=uniq_list, MARGIN=1, FUN=function(x) { sum(x) })
# model_div$hillfunc    <- as.numeric(unlist(calcDiv(comp_list, compDisMat=scales::rescale(as.matrix(dist(t(comp_list)), diag=TRUE, upper=TRUE)), q=1, type="FuncHillDiv")))
# 
# # Plot Shannon index
# pdf(paste("exo_stats_plots/ms1_comp_list_diversity_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
# boxplot(model_div$shannon ~ mzml_pheno_samples, col=mzml_pheno_colors, names=NA, main="Shannon diversity (H\')", xlab="treatment", ylab="Shannon diversity index (H\')")
# text(1:length(levels(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=levels(mzml_pheno_samples), xpd=TRUE, cex=0.9)
# div_tukey <- tukey.test(response=model_div$shannon, term=as.factor(mzml_pheno_samples))
# text(1:length(levels(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
# dev.off()
# 
# # PLS
# sel_pls_comp_list <- f.select_features_pls(feat_matrix=comp_list, sel_factor=mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=principal_components, tune_length=10, quantile_threshold=0.95, plot_roc_filename="exo_stats_plots/ms1_comp_list_select_pls_roc.pdf")
# print(paste("Number of selected variables:", f.count.selected_features(sel_feat=sel_pls_comp_list$`_selected_variables_`)))
# f.heatmap.selected_features(feat_list=comp_list, sel_feat=sel_pls_comp_list$`_selected_variables_`, sel_names=paste0("         ",sel_pls_comp_list$`_selected_variables_`), sample_colors=mzml_pheno_colors, plot_width=7, plot_height=7, cex_col=0.5, cex_row=0.4, filename=NULL, main="PLS")
# heatmaply(scale(comp_list[, which(colnames(comp_list) %in% sel_pls_comp_list$`_selected_variables_`)]), k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), file="exo_stats_plots/ms1_comp_list_select_pls.html", selfcontained=TRUE)
# sel_pls_comp_list$`_selected_variables_`
# sel_pls_comp_list$`_model_r2_`
# sel_pls_comp_list$`_multiclass_metrics_`
# 
# save(sel_pls_comp_list, file = "exo_stats_Results/sel_pls_comp_list.RData")

# # ---------- Variation partitioning ----------
# # sample dependent
# mzml_pheno_samples_type <- samp_groups
# mzml_pheno_samples
# # mono-culture or co-culture
# mzml_pheno_origin_samples_pos <- as.factor(origin_samp_groups <- c("CoCu", "CoCu", "CoCu", "CoCu", "CoCu", "CoCu",
#                                                             rep(x = "MonoCu", times = 8), 
#                                                             rep(x = "MonoCu", times = 8),
#                                                             "CoCu", "CoCu", "MB")
# )
#   
# model_varpart_pos <- varpart(scale(comp_list), ~ mzml_pheno_samples, ~ mzml_pheno_origin_samples_pos )
# 
# # Plot results
# pdf(file="exo_stats_plots/pos_ms1_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
# plot(model_varpart_pos, Xnames=c("samples","culture type"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
# dev.off()


# ############################## MS2 statistics examples ##############################
# 
# # PLS
# sel_pls_class_list <- f.select_features_pls(feat_matrix=class_list, sel_factor=mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=principal_components, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/ms2_class_list_select_pls_roc.pdf")
# print(paste("Number of selected variables:", f.count.selected_features(sel_feat=sel_pls_class_list$`_selected_variables_`)))
# f.heatmap.selected_features(feat_list=class_list, sel_feat=sel_pls_class_list$`_selected_variables_`, sample_colors=mzml_pheno_colors, plot_width=7, plot_height=7, cex_col=0.5, cex_row=0.4, filename="plots/ms2_class_list_select_pls.pdf", main="PLS")
# heatmaply(scale(class_list[, which(colnames(class_list) %in% sel_pls_class_list$`_selected_variables_`)]), k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), file="plots/ms2_class_list_select_pls.html", selfcontained=TRUE)
# sel_pls_class_list$`_multiclass_metrics_`
# sel_pls_class_list$`_model_r2_`
# 


# ############################## Molecular descriptors examples ##############################
# 
# 
# # Table L: samples x metabolites
# mdes_tab_l <- bina_list
# mdes_tab_l <- bina_list[, which(colnames(bina_list) %in% paste0(rownames(ms1_def_EXO_pos)[which(ms1_def_EXO_pos$smiles != "")], "_neg"))]
# mdes_tab_l <- as.data.frame(mdes_tab_l)
# 
# # Table R: samples x species
# mdes_tab_r <- as.data.frame.matrix(table(rownames(mdes_tab_l), mzml_pheno_samples))
# rownames(mdes_tab_r) <- rownames(mdes_tab_l)
# 
# # Table Q: metabolites x traits
# mdes_tab_q <- cdk_descriptors
# mdes_tab_q[is.na(mdes_tab_q)] <- 0
# 
# # Perform matrix operation
# mdes_list <- as.data.frame(as.matrix(mdes_tab_l) %*% as.matrix(mdes_tab_q))
# 
# # PLS
# sel_pls_mdes_list <- f.select_features_pls(feat_matrix=mdes_list, sel_factor=mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=principal_components, tune_length=10, quantile_threshold=0.995, plot_roc_filename="plots/descriptors_bina_list_select_pls_roc.pdf")
# print(paste("Number of selected descriptors:", f.count.selected_features(sel_feat=sel_pls_mdes_list$`_selected_variables_`)))
# f.heatmap.selected_features(feat_list=mdes_list, sel_feat=sel_pls_mdes_list$`_selected_variables_`, sel_names=paste0("        ",sel_pls_mdes_list$`_selected_variables_`), sample_colors=mzml_pheno_colors, plot_width=7, plot_height=7, cex_col=0.5, cex_row=0.4, filename="plots/descriptors_bina_list_select_pls.pdf", main="PLS")
# heatmaply(scale(mdes_list[, which(colnames(mdes_list) %in% sel_pls_mdes_list$`_selected_variables_`)]), k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), file="plots/descriptors_bina_list_select_pls.html", selfcontained=TRUE)
# sel_pls_mdes_list$`_multiclass_metrics_`
# sel_pls_mdes_list$`_model_r2_`
# sel_pls_mdes_list$`_selected_variables_`
# 

