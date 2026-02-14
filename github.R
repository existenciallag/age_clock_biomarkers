##This serves as an example of training biological aging measures using the NHANES III (1988 - 1994) and projecting into NHANES IV (1999 - 2018) dataset. It also provides documentation for fit parameters contained in the BioAge package. The cleaned NHANES dataset is loaded as the dataset NHANES3 and NHANES4. The original KDM bioage and phenoage values are saved as kdm0 and phenoage0 as part of NHANES dataset.

library(BioAge) #topic of example
library(dplyr)


##I train in the NHANES III and project biological aging measures into the NHANES IV by using the hd_nhanes, kdm_nhanes, and phenoage_nhanes function of the BioAge package.




#HD using NHANES (separate training for men and women)
hd = hd_nhanes(biomarkers=c("albumin","alp","lncrp","totchol","lncreat","hba1c","sbp","bun","uap","lymph","mcv","wbc"))

#KDM bioage using NHANES (separate training for men and women)
kdm = kdm_nhanes(biomarkers=c("albumin","alp","lncrp","totchol","lncreat","hba1c","sbp","bun","uap","lymph","mcv","wbc"))

#phenoage using NHANES
phenoage = phenoage_nhanes(biomarkers=c("albumin_gL","alp","lncrp","totchol","lncreat_umol","hba1c","sbp","bun","uap","lymph","mcv","wbc"))


##The projected data and estimated models are saved as part of the list structure. The dataset can be drawn by typing data. The model can be drawn by typing fit.


#assemble NHANES IV dataset with projected biological aging measures for analysis
data = merge(hd$data, kdm$data) %>% merge(., phenoage$data)

##In the figure below, the graphs titled “KDM Biological Age” and “Levine Phenotypic Age” show measures based on the original biomarker sets published in Levine 2013 J Geron A and Levine et al. 2018 AGING. The remaining graphs shows the new measures computed with the biomarker set specified within this code.



#select biological age variables
agevar = c("kdm0","phenoage0","kdm","phenoage","hd","hd_log")

#prepare labels
label = c("KDM\nBiological Age",
          "Levine\nPhenotypic Age",
          "Modified-KDM\nBiological Age",
          "Modified-Levine\nPhenotypic Age",
          "Homeostatic\nDysregulation",
          "Log\nHomeostatic\nDysregulation")

#plot age vs bioage
plot_ba(data, agevar, label)


##The figure plots associations among the different biological aging measures. Cells below the diagonal show scatter plots of the measures listed above the cell (x-axis) and to the right (y-axis). Cells above the diagonal show the Pearson correlations for the measures listed below the cell and to the left. For this analysis, KDM Biological Age and Levine Phenotypic Age measures are differenced from chronological age (i.e. plotted values = BA-CA).


#select biological age variables
agevar = c("kdm_advance0","phenoage_advance0","kdm_advance","phenoage_advance","hd","hd_log")

#prepare labels
#values should be formatted for displaying along diagonal of the plot
#names should be used to match variables and order is preserved
label = c(
  "kdm_advance0"="KDM\nBiological Age\nAdvancement",
  "phenoage_advance0"="Levine\nPhenotypic Age\nAdvancement",
  "kdm_advance"="Modified-KDM\nBiological Age\nAdvancement",
  "phenoage_advance"="Modified-Levine\nPhenotypic Age\nAdvancement",
  "hd" = "Homeostatic\nDysregulation",
  "hd_log" = "Log\nHomeostatic\nDysregulation")

#use variable name to define the axis type ("int" or "float")
axis_type = c(
  "kdm_advance0"="float",
  "phenoage_advance0"="float",
  "kdm_advance"="float",
  "phenoage_advance"="flot",
  "hd"="flot",
  "hd_log"="float")

#plot BAA corplot
plot_baa(data,agevar,label,axis_type)

table_surv(data, agevar, label)

table2 = table_health(data,agevar,outcome = c("health","adl","lnwalk","grip_scaled"), label)

#pull table
table2$table




