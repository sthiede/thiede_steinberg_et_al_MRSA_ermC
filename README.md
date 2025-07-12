# Genomic epidemiology analysis to detect the role of antibiotic use on the spread of MRSA in an urban jail and surrounding communities 

Select code for data analysis. 

# Preprocessing steps:
Data not provided on Github for these scripts because input fastas are too large. Provided to show how calculations were completed, but not runnable out of the box. 

- get_pw_distances_CCJ: Calculate pairwise distances between sequences and output summary stats
  
- Bowers_et_al_ST8_classifier: Codified logic from Bowers et al paper (mSphere 2018 May 2;3(3):e00464-17. doi: 10.1128/mSphere.00464-17) to characterize CC8 MRSA 
  
# GWAS
Raw data provided on Github

- GWAS: CCJ and CCH datasets; A Fisher’s exact test was conducted to assess the association between transmission and gene using the R package exact2x2 v1.6.5. We made the intentional decision not to explicitly control for population structure, as the clustering of strains is the signal we were trying to detect a genetic basis for. The relationship between odds ratio (OR) and p value was plotted and colored by gene. Significance was assessed with a Bonferonni-adjusted p value. We made the intentional decision not to explicitly control for population structure in our GWAS, as the clustering of strains is the signal we were trying to identify a genetic basis for. For example, using regression-based approaches that account for population structure by including a distance-matrix as a random effect, would control for the genetic clustering we are trying to detect.

# Phylogenetic analysis
Raw data not provided on Github

- phylogenetic_analysis: Evaluate phylogenetics of ermC plasmid in MRSA dataset
