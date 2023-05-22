# Functional over-representation analysis

## Description of analysis

Functional over-representation analysis among BphP1-dependent photoresponsive genes in *Pss* B728a was conducted using a Fisher’s exact test with histogram-based correction for multiple testing.  
<br/>

## Input file

* BphP1DepDEGs_function_counts.csv  
<br/>

## Analysis script

* FishersExact_Enrichment_jabes.R

__Note:__ The file ``` jabes.txt ``` is required as input during analysis and should be saved in the working directory.  
<br/>

## Notes regarding input file formatting

1.  The input .csv file is a table with headers containing the following:
    * __Column A:__ total number of genes assigned to each functional category
    * __Column B:__ titles of functional categories (*Note: "All genes" indicates gene counts regardless of functional categorization*)
    * __Columns C-E:__ number of genes that were upregulated in a BphP1-dependent manner in response to blue, red, or far-red (FR) light
    * __Columns F-H:__ number of genes that were downregulated in a BphP1-dependent manner in response to blue, red, or far-red (FR) light
