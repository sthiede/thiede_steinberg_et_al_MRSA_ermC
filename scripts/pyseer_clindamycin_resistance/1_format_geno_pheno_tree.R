# PHENO ---- 
# Goal: format phenotype to work with pyseer 
# Phenotype will be a .tsv with 2 columns & column names (no row names): 
# samples penicillin
# 6925_1#49       0
# 6925_1#50       0
# 6925_1#51       0

library(readxl)
library(ape)
# METADATA
cdc_new_metadata = read.csv('../2019_Project_BAA_MRSA_CO_HA/2021-05-06_make_master_metadata_updated/data/2021-05-06_master_metadata.csv')
cdc_new_metadata = cdc_new_metadata[cdc_new_metadata$molecular_type == "USA300-NAE" & cdc_new_metadata$seq_id != "NULL",]
cdc_new_metadata = cdc_new_metadata[!duplicated(cdc_new_metadata$study_id),]
cdc_new_metadata = cdc_new_metadata[cdc_new_metadata$Cd != 'NULL',]
cdc_new_metadata$Cd[cdc_new_metadata$Cd == 'I'] = 'R'

# TREE 
tree = read.tree('../2019_Project_BAA_MRSA_CO_HA/data/2020_02_04_14_32_53_MRSA_USA_300_genome_aln_w_alt_allele_unmapped_exclude_NEG_controls.treefile')
tree$tip.label = gsub('_$','',tree$tip.label)
sum(!cdc_new_metadata$seq_id %in% tree$tip.label)
tree = keep.tip(tree, cdc_new_metadata$seq_id)

cdc_new_metadata$Cd[cdc_new_metadata$Cd == 'R'] = 1
cdc_new_metadata$Cd[cdc_new_metadata$Cd == 'S'] = 0

pheno = data.frame(samples = cdc_new_metadata$seq_id, 
                   clinda = as.numeric(cdc_new_metadata$Cd))



# GENO ---- 


# SUBSET GENO TO SAMPLES 
# REMOVE INVARIANT SITES 

# Read in panaroo
geno = read.table('2021-02-05_run_panaroo_cdc_genomes/data/moderate/gene_presence_absence.Rtab', 
                  header = TRUE, 
                  row.names = 1)

# intersect of geno and pheno
samp_in_both = intersect(pheno$samples, colnames(geno))

geno = geno[,samp_in_both]

# remove invariant rows
geno =  geno[rowSums(geno) != ncol(geno),] # all 1s
geno = geno[rowSums(geno == 0) != ncol(geno),] # all 0s 
write.table(geno, file = '2021-05-10_pyseer_cdc_pangenome_vs_clindaR/data/geno.tsv', row.names = TRUE, sep = '\t')



# PHENO AGAIN ---- 
# subset pheno further
pheno = pheno[pheno$samples %in% samp_in_both,]
write.table(pheno, file = '2021-05-10_pyseer_cdc_pangenome_vs_clindaR/data/pheno.tsv', row.names = FALSE, sep = '\t')



# TREE ----
# ROOT TO USA500 isolate 
is.rooted(tree) # FALSE

# DROP SAMPLES NOT IN TREE
tree = keep.tip(tree, pheno$samples)


# SAVE TREE
write.tree(tree, file = '2021-05-10_pyseer_cdc_pangenome_vs_clindaR/data/cdc_tree_rooted.tree')

