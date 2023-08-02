# Co-culture workflow
This workflow can be used to pre-process and analyse MS1 data stemming from a two species, untargeted, comparative metabolomics liquid-chromatography tandem mass spectrometry approach. It takes .mzML files as input and gives processed feature lists and statistical analyses of the MS level 1 data.

The workflow was developed using a co-culture experiment between the two microalgal species _S. marinoi_ and _P. parvum_

## HOW TO:
1. 

2. 

3. 

4. For in-depth statistical analysis, the _/Statistics/_ directory contains RScripts for several statistical methods. The objects used are the previously produced .csv tables of feat_list and bina_list. The negative and positive feature tables are combined and analysis is performed within the conditions (EXO and ENDO) combined and on species level. Statistical methods used include:
	- PCA
	- Diversity measures (Shannon diversity index, ...)
	- PLS 
	- t-SNE
	- Variation partitioning 

5. By following the _linking.R_ script, MS1 precursor can be linked to MS2 fragementation data. The MS2 fragmentation can be annotated using [MAW](https://github.com/zmahnoor14/MAW)



### Additional files and info
some statistical function are from the [iESTIMATE](https://github.com/ipb-halle/iESTIMATE/blob/main/use-cases/radula-hormones/peak_detection_neg.r) repository


testing_MS1_processing.R: 
My personal notes for code ideas to be implemented in MS1_processing

MS1_MS2_precursor.R:
Find the MS2 precursor mass in the MS1 data by following the iESTIMATE code of ipb-halle

MS1_workflow_XCMS.R:
Clearer version of MS1_MS2_precursor.R, is the pre-processing with newest XCMS functionalities 
