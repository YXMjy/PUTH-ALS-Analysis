# PUTH-ALS-Analysis 

*** This is a series of codes on an epigenome-wide association study of ALS cases in the Chinese population. 

*** Firstly, we identified a number of differentially methylated positions (DMPs), 5 of which were of high confidence hypermethylated biomarkers 
in the PUTH-ALS cohort (Peking University Third Hospital ALS cohort).
   
    reference file "1-age_corplot.R" (Figure S1),
                   "1-PC_plot.R" (Figure 1),
                   "2-DMPs identification.R" (Table 2/Table S1),
                   "2-GC_correct.R" (Figure 2/Table S1),
                   "2-go_plot.R" (Figure S3/Table S3/Table S4),
                   "2-pattern_plot.R" (Figure S3),
                   "2-volcano_plot.R" (Figure 3).


*** Second, we conducted an external validation of these DMPs and 5 biomarkers based on ALS postmortem expression profiles from the NYGC-ALS cohort 
(coordinated by the New York Genomic Center; Access link:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153960; Access ID：GSE153960).
    
    reference file "3-gene exp.R","3-pvalue-corplot.R" (Figure 3/Table S2).


*** Thirdly, we also identified 3 DMPs associated with survival time of ALS cases in the PUTH-ALS cohort.
    
    reference file "4-survival.R" (Figure 4)

*** Fourth, we identified a 27-loci signature using machine learning methods, and calculated the methylation risk score for each individual in the 
PUTH-ALS cohort. Considering the possible interaction of genetic and epigenetic factors with ALS, we further generated the genetic-profile scores 
(the weighted sum of the effect allele count) based on summary statistics from an EU-ALS cohort (derived from a genome-wide association study (GWAS) 
based on large European population cohorts consisting of 20,806 ALS cases and 59,804 controls; Access link:http://als.umassmed.edu; reference research: https://doi.org/10.1016/j.neuron.2018.02.027).
    
    reference file "5-MS model.R","5-MS_rocplot.R","5-PRS_MS.R","5-PRS_plink.txt","5-PRS_snp-match.R" (Figure 5).
