#### First link the feature table with inclusion list

##### separate inclusion lists ENDO
endo_inc <- read.csv("/KSS_210324_ENDO.csv", sep =";")
colnames(endo_inc) <- c("Mass",	"Formula",	"Formula_type",	"Species",	"CS",	"Polarity",	"Start",	"End",	"CE",	"Comment")
endo_inc["Start"] <- endo_inc["Start"]*60
endo_inc["End"] <- endo_inc["End"]*60
neg <- endo_inc[endo_inc["Polarity"]=="Negative" , ]
pos <- endo_inc[endo_inc["Polarity"]=="Positive" , ]
write.csv(neg, "KSS_210324_ENDO_neg.csv")
write.csv(pos, "KSS_210324_ENDO_pos.csv")

##### separate inclusion lists EXO
exo_inc <- read.csv("/KSS_210324_EXO.csv", sep =";")
colnames(exo_inc) <- c("Mass",	"Formula",	"Formula_type",	"Species",	"CS",	"Polarity",	"Start",	"End",	"CE",	"Comment")
exo_inc["Start"] <- exo_inc["Start"]*60
exo_inc["End"] <- exo_inc["End"]*60
neg <- exo_inc[exo_inc["Polarity"]=="Negative" , ]
pos <- exo_inc[exo_inc["Polarity"]=="Positive" , ]
write.csv(neg, "KSS_210324_EXO_neg.csv")
write.csv(pos, "KSS_210324_EXO_pos.csv")

#---------------------------------------------------------------------------------------------------------------------------

##### Linking #####


linking_with_inc_list <- function(ft_csv, inc_csv, cond, mass_col, rtmin_col, rtmax_col){
  # make feature table smaller by excluding features not in MS2 inclusion list
  ft_tbl <- read.csv(ft_csv)
  inc_tbl <- read.csv(inc_csv)

  inc_tbl['feat_inc'] <- NA
  inc_tbl["Comment"]<- NA
  
  for (i in 1:nrow(inc_tbl)){

    for (j in 1:nrow(ft_tbl)){
      if (inc_tbl[i, mass_col]>ft_tbl[j, "mzmin"] && ft_tbl[j, "mzmax"]>inc_tbl[i, mass_col]){
        if (ft_tbl[j, "rtmed"]>inc_tbl[i, rtmin_col] && inc_tbl[i, rtmax_col]>ft_tbl[j, "rtmed"]){
          inc_tbl[i, 'feat_inc'] <- ft_tbl[j, "feat_id"]
        }
      }
      else if(round(inc_tbl[i, mass_col],4)>=round(ft_tbl[j, "mzmin"],4) && round(ft_tbl[j, "mzmax"],4)>=round(inc_tbl[i, mass_col],4)){
        if (ft_tbl[j, "rtmed"]>inc_tbl[i, rtmin_col] && inc_tbl[i, rtmax_col]>ft_tbl[j, "rtmed"]){
          inc_tbl[i, 'feat_inc'] <- ft_tbl[j, "feat_id"]
        }
      }
      else if(round(inc_tbl[i, mass_col],3)>=round(ft_tbl[j, "mzmin"],3) && round(ft_tbl[j, "mzmax"],3)>=round(inc_tbl[i, mass_col],3)){
        if (ft_tbl[j, "rtmed"]>inc_tbl[i, rtmin_col] && inc_tbl[i, rtmax_col]>ft_tbl[j, "rtmed"]){
          inc_tbl[i, 'feat_inc'] <- ft_tbl[j, "feat_id"]
        }
      }
    }
  }
  
  
  for (i in 1:nrow(inc_tbl)){
    if (is.na(inc_tbl[i, "feat_inc"])){
      ft_tbl_ids <- c()
      ft_rt_med <- c()
      for (j in 1:nrow(ft_tbl)){
        if(round(inc_tbl[i, mass_col],2)>=round(ft_tbl[j, "mzmin"],2) && round(ft_tbl[j, "mzmax"],2)>=round(inc_tbl[i, mass_col],2)){
          ft_tbl_ids <- c(ft_tbl_ids, ft_tbl[j, "feat_id"])
          ft_rt_med <- c(ft_rt_med, ft_tbl[j, "rtmed"])
        }
      }
      ft_tbl_rt_sh <- cbind(ft_tbl_ids, ft_rt_med)
      print(nrow(ft_tbl_rt_sh))
      if (!(is.null(ft_tbl_rt_sh))){
        if (nrow(ft_tbl_rt_sh) == 1){
          inc_tbl[i, 'feat_inc'] <- ft_tbl_ids
          inc_tbl[i, "Comment"] <- paste("The ONLY rtmed for this premz is:", ft_tbl_rt_sh[1,"ft_rt_med"])
        }
        else if (nrow(ft_tbl_rt_sh)>1){
          # Calculate the absolute differences between each number and the range boundaries
          differences <- abs(as.numeric(ft_tbl_rt_sh[,"ft_rt_med"]) - inc_tbl[i, rtmin_col])
          differences <- pmin(differences, abs(as.numeric(ft_tbl_rt_sh[,"ft_rt_med"]) - inc_tbl[i, rtmax_col]))
          
          # Find the index of the number with the minimum difference
          closest_index <- which.min(differences)
          inc_tbl[i, 'feat_inc'] <- ft_tbl_rt_sh[closest_index, "ft_tbl_ids"]
          # Get the closest number
          closest_number <- as.numeric(ft_tbl_rt_sh[,"ft_rt_med"])[closest_index]
          
          # Print the result
          inc_tbl[i, "Comment"] <- paste("The rtmed closest to the range is:", closest_number)
        }
      }
    }
  }
  
  write.csv(inc_tbl,paste("ms_tbl_", cond, ".csv", sep = ""))
  return(inc_tbl)
}

link_ann_inc <- function(inc_tbl, ann_tbl, cond){
  ann_tbl <- read.csv(ann_tbl) # results from MAW
  inc_tbl <- read.csv(inc_tbl) # result from previous function
  ann_tbl['feat_inc'] <- NA
  ann_tbl["Comment"]<- NA
  for (i in 1:nrow(ann_tbl)) {
    for (j in 1:nrow(inc_tbl)) {
      if (round(ann_tbl[i, "premz"], 4) == round(inc_tbl[j, "Mass"], 4)){
        if (ann_tbl[i, "rtmed"]>inc_tbl[j, "Start"] && inc_tbl[j,"End"]>ann_tbl[i, "rtmed"]){
          x <- x + 1
          ann_tbl[i, "Comment"] <- inc_tbl[j, "Comment"]
          ann_tbl[i, "feat_inc"] <- inc_tbl[j, "feat_inc"]
        }
      }
      else if (round(ann_tbl[i, "premz"], 3) == round(inc_tbl[j, "Mass"], 3)){
        if (ann_tbl[i, "rtmed"]>inc_tbl[j, "Start"] && inc_tbl[j,"End"]>ann_tbl[i, "rtmed"]){
          x <- x + 1
          ann_tbl[i, "Comment"] <- inc_tbl[j, "Comment"]
          ann_tbl[i, "feat_inc"] <- inc_tbl[j, "feat_inc"]
        }
      }
    }
  }

  for (i in 1:nrow(ann_tbl)){
    if (!(grepl("FT", ann_tbl[i, "feat_inc"]))){
      ft_tbl_ids <- c()
      ft_rt_min <- c()
      ft_rt_max <- c()
      for (j in 1:nrow(inc_tbl)){
        if(round(ann_tbl[i, "premz"], 3) == round(inc_tbl[j, "Mass"], 3)){
          ft_tbl_ids <- c(ft_tbl_ids, inc_tbl[j, "feat_inc"])
          ft_rt_min <- c(ft_rt_min, inc_tbl[j, "Start"])
          ft_rt_max <- c(ft_rt_max, inc_tbl[j, "End"])
        }
      }
      ft_tbl_rt_sh <- cbind(ft_tbl_ids, ft_rt_min, ft_rt_max)
      print(nrow(ft_tbl_rt_sh))
      if (nrow(ft_tbl_rt_sh) == 1){
        ann_tbl[i, 'feat_inc'] <- ft_tbl_ids
        ann_tbl[i, "Comment"] <- paste("The ONLY range for this ann_premz is:", ft_tbl_rt_sh[1,"ft_rt_min"], "-", ft_tbl_rt_sh[1,"ft_rt_max"], sep = "")
      }
      else if (nrow(ft_tbl_rt_sh)>1){
        # Calculate the absolute differences between each number and the range boundaries
        differences <- abs(as.numeric(ann_tbl[i,"rtmed"]) - as.numeric(ft_tbl_rt_sh[, "ft_rt_min"]))
        differences <- pmin(differences, abs(as.numeric(ann_tbl[i,"rtmed"]) - as.numeric(ft_tbl_rt_sh[, "ft_rt_max"])))
        # Find the index of the number with the minimum difference
        closest_index <- which.min(differences)
        ann_tbl[i, 'feat_inc'] <- ft_tbl_rt_sh[closest_index, "ft_tbl_ids"]
        # Get the closest number
        closest_min <- as.numeric(ft_tbl_rt_sh[,"ft_rt_min"])[closest_index]
        closest_max <- as.numeric(ft_tbl_rt_sh[,"ft_rt_max"])[closest_index]
        # Print the result
        ann_tbl[i, "Comment"] <- paste("The rtmed closest to range is:", closest_min, "-", closest_min)
      }
    }
  }


  write.csv(ann_tbl, paste(cond, "_ann.csv", sep = ""))
  return(ann_tbl)
}

assign_sig_feat <- function(sig_feat, ann_tbl_modified, cond, pol, feat_col, origin_col){
  
  sig_feat <- read.csv(sig_feat)
  ann_tbl <- read.csv(ann_tbl_modified)
  sig_feat_pos <- sig_feat[grepl("_pos", sig_feat[,feat_col]), ]
  
  sig_feat_pos["sig_feat_ids"] <-gsub('_pos', '', sig_feat_pos[,feat_col])
  sig_feat_pos["sig_feat_ids"]
  
  ann_tbl["origin"] <- NA
  ann_tbl
  colnames(sig_feat_pos)
  for (i in 1:nrow(ann_tbl)){
    if (!(is.na(as.character(ann_tbl[i, "feat_inc"])))){
      
      if (as.character(ann_tbl[i, "feat_inc"]) == as.character(sig_feat_pos[j, "sig_feat_ids"])){
        ann_tbl[i, "origin"] <- sig_feat_pos[j, origin_col]
      }
    }
  }
  
  write.csv(ann_tbl, paste(cond, "_mergedResults-with-one-Candidates_sig_feat.csv"), sep = "")
  return(ann_tbl)
}



#---------------------------------------------------------------------------------------------------------------------------

##### script #####



setwd("/any/dir")

ft_csv = "feature_info_EXO_pos.csv"
inc_csv = "KSS_210324_EXO_pos.csv"
cond = "exo_pos"
mass_col = "Mass"
rtmin_col = "Start"
rtmax_col = "End"
ann_tbl <- "exo_pos_mergedResults-with-one-Candidates.csv"
ann_tbl_modified <- paste(cond, "_ann.csv", sep ="")
sig_feat <- "sel_pls_feat_ENDO.csv"
feat_col <- "sel_pls_features_EXO"
origin <- "originEXO"
pol = "_pos" # or "_neg"

inc_tbl1 <- linking_with_inc_list(ft_csv, inc_csv, cond, mass_col, rtmin_col, rtmax_col)
inc_tbl1 

ann_tbl1 <- link_ann_inc(inc_tbl = paste("ms_tbl_", cond, ".csv", sep = ""), ann_tbl, cond)
ann_tbl1

significant_annotated_feats <- assign_sig_feat(sig_feat, ann_tbl_modified, cond, pol, feat_col, origin_col)
significant_annotated_feats
