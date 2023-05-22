# Differential gene expression analysis

## Description of analysis

Differential gene expression analysis was performed using the QuasiSeq v1.0-11-0 R package. Read counts are normalized for library size using upper-quartile normalization, and the negative binomial QLSpline method, with fixed effects for replicates and strain x treatment combinations, is used to estimate the mean Ln(fold change) in expression and calculate the associated p-value for each gene for each comparison. The model uses the read count data for each gene across all 66 samples to estimate gene-specific error variances, while differential expression is subsequently assessed only for specified sample pairs of interest. To account for multiple testing, all *p*-values are converted to *q*-values in QuasiSeq using histogram-based false discovery estimation. 

QuasiSeq can be installed from source at: <br/> https://cran.r-project.org/src/contrib/Archive/QuasiSeq/.  
<br/>

## Input file

* B278aLight_read_counts.csv  
<br/>

## Analysis script

* Quasi_analysis_B278aLight.R

*__Note:__ A* ``` figures/ ```  *subdirectory should be created within the working directory before running the analysis. If this subdirectory does not already exist when the script is run, the p-value histograms and volcano plots will not be saved.*  
<br/>

## Notes regarding file formatting

1.  The input .csv file is the table of read counts for each feature (output from HTSeq-count) but with:
    * headers removed
    * columns sorted so that the three replicates for each strain x treatment combination are grouped in adjacent columns  
 <br/>   

2. The R script uses two-character treatment codes (e.g., 1A, 2A, 4B, etc.) to define columns (i.e., samples) and specify which samples are to be compared for assessing differential expression. The treatment codes correspond to the following strain x treatment combinations:
<p align="center">
<img src="https://user-images.githubusercontent.com/128737867/236014974-670872d6-ac6c-47b4-a8c1-5b64d742f2c0.png" width=350>
</p>
