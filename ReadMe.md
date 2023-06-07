# Custom R scripts for the manuscript '*Light cues induce protective anticipation of environmental water loss in terrestrial bacteria*'

**Authors:** Bridget M. Hatfield, Breah LaSarre, Meiling Liu, Haili Dong, Dan Nettleton, Gwyn A. Beattie

**Manuscript Status:** submitted

**Contact:** gbeattie@iastate.edu  
<br /> 


## Differential gene expression analysis

[ReadMe](DifferentialExpression/)

### *Files*

* __[Quasi_analysis_B278aLight.R](DifferentialExpression/Quasi_analysis_B278aLight.R)__ : R script for differential expression analysis of RNAseq data from *Pseudomonas syringae* exposed to distinct light wavelengths using QuasiSeq 

* __[B278aLight_read_counts.csv](DifferentialExpression/B278aLight_read_counts.csv)__ : input data for QuasiSeq analysis; table contains the number of reads aligning to each gene in the genome in each sample  
<br /> 

## Functional over-representation analysis

[ReadMe](OverRepresentation/)

### *Files*

* __[FishersExact_Enrichment_jabes.R](OverRepresentation/FishersExact_Enrichment_jabes.R)__ : R script for analyzing functional overrepresentation among BphP1-dependent light-responsive genes using a Fisher's exact test

* __[jabes.txt](OverRepresentation/jabes.txt)__ : text file containing the jabes.q function for histogram-based multiple-comparison correction

* __[BphP1DepDEGs_function_counts.csv](OverRepresentation/BphP1DepDEGs_function_counts.csv)__ : input data for overrepresentation analysis; table contains the number of total genes and number of BphP1-dependent light-responsive genes assigned to different functional categories  
<br />

## Source data and statistics

[ReadMe](SourceData/)

### *Files*

* __[Source Data Fig 3.xlsx](SourceData/Source%20Data%20Fig%203.xlsx)__ : Excel file containing the source data and full statistical analysis results for plots in Fig. 3

* __[Source Data Fig 4.xlsx](SourceData/Source%20Data%20Fig%204.xlsx)__ : Excel file containing the source data and full statistical analysis results for plots in Fig. 4

* __[Source Data Fig S4.xlsx](SourceData/Source%20Data%20Fig%20S4.xlsx)__ : Excel file containing the source data and full statistical analysis results for plots in *SI Appendix*, Fig. S4

* __[Source Data Fig S5.xlsx](SourceData/Source%20Data%20Fig%20S5.xlsx)__ : Excel file containing the source data and full statistical analysis results for plots in *SI Appendix*, Fig. S5

* __[Source Data Fig S6.xlsx](SourceData/Source%20Data%20Fig%20S6.xlsx)__ : Excel file containing the source data and full statistical analysis results for plots in *SI Appendix*, Fig. S6
