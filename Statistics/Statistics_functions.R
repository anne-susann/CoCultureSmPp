###---- library ----
library(RColorBrewer)           # For colors
library(MSnbase)                # MS features
library(xcms)                   # Swiss army knife for metabolomics
library(CAMERA)                 # Metabolite Profile Annotation
library(Spectra)                # Spectra package needed for XCMS3
library(vegan)                  # For Statistics/Varpart
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
library(plotly)                 # For creating html plots
library(htmlwidgets)            # For creating html plots
library(heatmaply)              # HTML heatmaps
library(stringr)              
library(randomForest)           # Random forest
library(MetaboAnalystR)         # Random forest
#library(iESTIMATE)
source("https://raw.githubusercontent.com/ipb-halle/iESTIMATE/main/R/_functions.r")


# statistics is performed on species level
# if the data is separated by polarity for the pre-processing part, the first step in the analysis is to create a feature table 
# containing both polarities. This combination can only be done when the data distribution in the histogram look similar for the polarities
# ------ pre_prep_stat --------
pre_prep_stat <- function(feat_table_pos, feat_table_neg, result_dir_name){
  # feat_table1 = feat_table of polarity a, can also be binary list
  # feat_table2 = feat_table of polarity b, can also be binary list
  
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
  
  
  # name the features per polarity
  colnames(feat_table_pos) <- paste(colnames(feat_table_pos), "pos", sep = "_")
  colnames(feat_table_neg) <- paste(colnames(feat_table_neg), "neg", sep = "_")
  
  # combine the feature tables of neg and pos into one
  feat_list <- cbind(feat_table_pos, feat_table_neg)
  
  # save
  write.csv(feat_list, file=paste(result_dir_name, "/feat_list_combined.csv", sep = ""), row.names=FALSE)
  return(feat_list)
  
}


### trial run pre_prep_stat
# import feature tables
feature_table_pos <- read.csv("endo_pos/endo_pos_Results/feature_list_ENDO_pos.csv", row.names = 1)
feature_table_neg <- read.csv("endo_neg/endo_neg_Results/feature_list_ENDO_neg.csv", row.names = 1)


# for statistical anaylses, any blank needs to be excluded from the feature list
# ------ phenodata_preparation_stat --------
phenodata_preparation_stat <- function( phenodata, result_dir_name, plots_dir_name){
  # feat_list = combined or simple feature list from pre-processing workflow
  # phenodata = csv table of phenodata including (1) sample_name, (2) sample_group, (3) sample_description 
  
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
  
  return(phenodata)
 
}


data_preparation_stats <- function(feat_list){
 # set rownames of feature table
  for (i in nrow(feat_list)) {
    rownames(feat_list)[i] <- phenodata$sample_name[i]
  }
  
  return(feat_list)
}


# ---- diversity_index ----
diversity_index <- function(){
  
}

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
#rownames(model_div) <- paste(pheno_data_ENDO[1:24,1], pheno_data_ENDO[1:24,2], sep = "_" )

# test for functional hill
# subset comp_list for the first 20 features
#comp_list_hillfunc <- comp_list[,1:1000]
#model_div$hillfunc    <- as.numeric(unlist(calcDiv(comp_list_hillfunc, compDisMat=scales::rescale(as.matrix(dist(t(comp_list_hillfunc)), diag=TRUE, upper=TRUE)), q=1, type="FuncHillDiv")))

# color for shannon plot
col_shan <- c(unique(CoCuPp1), unique(CoCuSm1), unique(CoCuPp1), unique(CoCuSm1))

# Plot Shannon index
pdf(paste("ENDO_stats_plots/ms1_comp_list_diversity_shannon_Culture.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
#jpeg(filename = "ENDO_stats_plots/ms1_comp_list_diversity_shannon_Culture.jpeg", width = 1000, height = 700, quality = 100, bg = "white")
boxplot(model_div$shannon ~ mzml_pheno_samples_type, col= col_shan, main="Shannon diversity (H\')", xlab="culture types", ylab="Shannon diversity index (H\')")
text(1:length(levels(mzml_pheno_samples_type)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=levels(mzml_pheno_samples_type), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$shannon, term=as.factor(mzml_pheno_samples_type))
text(1:length(levels(mzml_pheno_samples_type)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
legend("bottomleft", bty="n", pt.cex=5, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(col_shan), legend= c("P.parvum", "S.marinoi"))
dev.off()


# ---- PCA_analysis----
PCA_analysis <- function(feat_list, phenodata, plots_dir_name){
  jpeg(filename = paste(plots_dir_name, "/ENDO_ms1_feature_table_pca.jpeg", sep = ""), width = 1000, height = 700, quality = 100, bg = "white")
  ms1_pca <- prcomp(feat_list, center=TRUE)
  par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
  plot(ms1_pca$x[, 1], ms1_pca$x[,2], pch=19, main="PCA of feature table",
       xlab=paste0("PC1: ", format(summary(ms1_pca_ENDO)$importance[2, 1] * 100, digits=3), " % variance"),
       ylab=paste0("PC2: ", format(summary(ms1_pca_ENDO)$importance[2, 2] * 100, digits=3), " % variance"),
       col=color, cex=2)
  grid()
  text(ms1_pca$x[,1], ms1_pca$x[,2], labels=str_sub(rownames(feat_list), - 3, - 1), col=color, pos=3, cex=1.5)
  legend("topleft", bty="n", pt.cex=5, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
         col= unique(phenodata$color), legend= unique(phenodata$sample_description))
  dev.off()
  
}


# ---- PLS_analysis ----
PLS_analysis <- function(feat_list, principal_components, sel_factor, sel_color){
  principal_components <- 5
# PLS according to culture type by species Pp
sel_pls_comp_list <- f.select_features_pls(feat_matrix=bina_list_Pp, 
                                           sel_factor=mzml_pheno_origin_samples[index_PP], 
                                           sel_colors=mzml_pheno_colors[index_PP], 
                                           components=principal_components, tune_length=10, 
                                           quantile_threshold=0.95, plot_roc_filename="/ms1_comp_list_select_pls_roc_Pp_bina.pdf")
print(paste("Number of selected variables:", f.count.selected_features(sel_feat=sel_pls_comp_list$`_selected_variables_`)))
heatmaply(scale(comp_list_Pp[, which(colnames(comp_list_Pp) %in% sel_pls_comp_list$`_selected_variables_`)]), 
          k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), 
          file="/ms1_comp_list_select_pls_Pp_bina.html", 
          selfcontained=TRUE,
          column_text_angle = 90, 
          showticklabels = c(TRUE, TRUE),
          #titleX = FALSE,
          branches_lwd = 0.3,
          #width = 1000, heigth = 800,
          #colorbar_thickness = 50,
          #colorbar_len = 1,
          cexRow = 2, cexCol = 1.5, cexColorbar = 2)
}




sel_pls_comp_list$`_selected_variables_`
sel_pls_comp_list$`_model_r2_`
sel_pls_comp_list$`_multiclass_metrics_`

# selected features from PLS 
features_Pp <- sel_pls_comp_list$`_selected_variables_`
features_Pp <- sort(features_Pp, decreasing = FALSE)

#save(sel_pls_comp_list, file = "ENDO_stats_Results/sel_pls_comp_list_culture_Pp_bina.RData")
load("ENDO_stats_Results/sel_pls_comp_list_culture_Pp_bina.RData")
# which features in which condition P.PARVUM
bina_list_Pp_features <- bina_list_Pp[,features_Pp]
bina_list_Pp_features_Co <- bina_list_Pp_features[c(1,2,3,12),]
bina_list_Pp_features_Co <- bina_list_Pp_features_Co[,colSums(bina_list_Pp_features_Co)>0]
bina_features_Pp_Co <- data.frame(colnames(bina_list_Pp_features_Co))
bina_features_Pp_Co$Description <- "present in Pp Co-culture"
bina_features_Pp_Co$Condition <- "Co"
colnames(bina_features_Pp_Co) <- c("Feature", "Description","Condition")

bina_list_Pp_features_Mono <- bina_list_Pp_features[c(4:11),]
bina_list_Pp_features_Mono <- bina_list_Pp_features_Mono[,colSums(bina_list_Pp_features_Mono)>0]
bina_features_Pp_Mono <- data.frame(colnames(bina_list_Pp_features_Mono))
bina_features_Pp_Mono$Description <- "present in Pp Mono-culture"
bina_features_Pp_Mono$Condition <- "Mono"
colnames(bina_features_Pp_Mono) <- c("Feature", "Description","Condition")


sel_pls_features_Pp_bina <- data.frame(rbind(bina_features_Pp_Co, bina_features_Pp_Mono))
sel_pls_features_Pp_bina$Feature_id <- gsub("_neg","", sel_pls_features_Pp_bina$Feature)
sel_pls_features_Pp_bina$Feature_id <- gsub("_pos","", sel_pls_features_Pp_bina$Feature_id)
write.csv(sel_pls_features_Pp_bina, file = "ENDO_stats_Results/manuscript/sel_pls_features_Pp_bina.csv")


