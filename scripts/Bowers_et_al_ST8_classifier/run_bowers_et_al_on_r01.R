# canSNPs SCRIPT
source('Bowers_et_al_ST8_classification.R')
library(ape)

# RAN ON A CLUSTER BECAUSE ALIGNMENT FILE IS LARGE

# READ IN ALIGNMENT 
path_to_dna = '../processed_data/2023_12_11_13_35_06_MRSA_USA_300_MRSA_CCH_QC_Passed_ST8_genomes_aln_w_alt_allele_unmapped.fa'
dna = read.dna(path_to_dna, format = 'fasta')
row.names(dna) = gsub('_$','',row.names(dna))

# Bowers classification scheme for ST8 
class = sapply(labels(dna), function(i){
  where_are_the_probes(dna[i,])
})

print('saving')
save(class, file = '../processed_data/r01_probe_matches.RData')

# Classify based on probes found 
final_class = data.frame(bowers_class(class))
