# PRANA

PRANA: A Pseudo-Value Regression Approach for Differential Network Analysis of Co-Expression Data <br />
Submitted to *Bioinformatics* for consideration.
</br>


### Usage
This repository contains clinical and expression data, saved as **combinedCOPD_RelatedGenesOnly.RDS**. Of note, this looks at 28 COPD-related genes that were highlighted in a recent genome-wide association study (Sakornsakolpat *et al.*, 2019). The full data is available from the Gene Expression Database with accession number GSE158699 (Wang *et al.*, 2021).

* **TotalConnectivity.R** is to calculate the total connectivity of each gene. 
* **EmpiricalBayes_Datta_2005.R** is to compute the adjusted p-values via empirical Bayes approach (Datta and Datta, 2005)
* **PRANA_main.R** is the primary code to analyze the data. This will load the first two codes, so please download all three R codes to apply our method. The code is heavily annotated for the ease of implementation.
</br>

### References
Sakornsakolpat, P., Prokopenko, D., Lamontagne, M., and et al. (2019). Genetic landscape of chronic obstructive pulmonary disease identifies heterogeneous cell- type and phenotype associations. *Nature Genetics*, **51**(3), 494–505. </br></br>
Wang, Z., Masoomi, A., Xu, Z., Boueiz, A., Lee, S., Zhao, T., Bowler, R., Cho, M., Silverman, E., Hersh, C., Dy, J., and Castaldi, P. (2021). Improved prediction of smoking status via isoform-aware RNA-seq deep learning models. *PLoS Computational Biology*, **17**(10). </br></br>
Datta, S. and Datta, S. (2005). Empirical Bayes screening of many p-values with applications to microarray studies. *Bioinformatics*, **21**(9), 1987–1994.
</br>


### Contact
Please feel free to email me (Seungjun Ahn) at sahn1@ufl.edu if you need any further information or suggestions to improve our work.

