# CoCultureSmPp
Analysis of the MS1 and MS2 data of Co-Culture experiment with _S. marinoi_ and _P. parvum_

## HOW TO:
1. Pre-process **ONLY MS1 data** with scripts from _/MS1_workflow_XCMS/_ directory. The MS1 data is pre-processed and analysed with specific parameters for each condition (see ENDO_pos, EXO_neg etc.). Run code until ### linking MS2 data ###. Until here you get an overview over the data and first insights into its nature with PCA plots. Created results will be:
	- ms1_data_xxx_pol: contains all detected peaks and additional information like polarity, msLevel, filterString etc. and contains the further pre-processing 		(saved as XXX_pol_raw_data.csv and MS1_XXX_pol_peak_detection.RData)
	- feat_list_XXX_pol: contains the features (saved as feature_list_XXX_pol.csv)
	- bina_list_XXX_pol: a binary list of absence/presence of peaks across the samples (saved as bina_list_XXX_pol.csv)
	- model_div_XXX_pol: a list with different diversity measures for the samples
	- ms1_def_XXX_pol: data frame containing the feature definitions


2. Run the pre-processing on **MS1 AND MS2 data COMBINED** with the same parameters using the scripts in the _/MS1and2_workflow_XCMS/_ directory. The necessary objects are saved as .RData.

3. For the linking of MS1 and MS2 data, return to the script in _/MS1_workflow_XCMS/_ and follow from ### linking MS2 data ### on, loading the combined object created in the other directory *ms_data_xxx_pol*. The script establishes a connection between the precursor masses and the origin files using the feature definitions stored in *ms1_def_XXX_pol*, giving the MS2 spectra as a .mgf file to be used with MAW for feature annotation.

4. For in-depth statistical analysis, the _/Statistics/_ directory contains RScripts for several statistical methods. The objects used are the previously produced .csv tables of feat_list and bina_list. The negative and positive feature tables are combined and analysis is performed within the conditions (EXO and ENDO) combined and on species level. Statistical methods used include:
	- PCA
	- Diversity measures (Shannon diversity index, ...)
	- PLS 
	- t-SNE
	- Variation partitioning
	-~~ANOVA~~
	- Random forest 

5. ~~The information from the annotation can be integrated with the statistical analysis, giving annotation to the siginifcant features determined in the statistics~~



### Additional files and info
xxx = stand-in for conditions (ENDO, EXO)
pol = stand-in for polarity (pos, neg)
related code: [iESTIMATE](https://github.com/ipb-halle/iESTIMATE/blob/main/use-cases/radula-hormones/peak_detection_neg.r)


coculture_modified_endo_pos.R:
Is the MS1 AND MS2 COMBINED script with the linking included


MS1_processing.R: 
Pre-processing with old XCMS functionalities (xcmsSet) and some statistical analysis of the data (PCA, broken stick)


testing_MS1_processing.R: 
My personal notes for code ideas to be implemented in MS1_processing



MS1_MS2_precursor.R:
Find the MS2 precursor mass in the MS1 data by following the iESTIMATE code of ipb-halle


MS1_workflow_XCMS.R:
Clearer version of MS1_MS2_precursor.R, is the pre-processing with newest XCMS functionalities 
