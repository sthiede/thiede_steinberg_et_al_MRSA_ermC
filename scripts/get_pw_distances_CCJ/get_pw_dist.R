# LIBRARIES ---- 
library(openxlsx)
library(ape)

# FUNCTION ----
# dna_obj = DNAbin object 
# prefix = string you want to incoporate in the output file name

# Output: 
# QC figures 
# pairwise distance matrices RData 
# summary about core genome size, number of Ns and dashes
# RData of DNAbin with only variant positions 

get_pw_dist_and_core_size <- function(dna_obj, prefix){
  
  # Ns and dashes per sample 
  Ns_and_dash_per_sample = apply(dna_obj, 1, function(sample){
    char = ape::as.character.DNAbin(sample)
    c(N = sum(char == 'n'), 
      dash = sum(char == '-'))
  })
  write.table(Ns_and_dash_per_sample, file = paste0('2020-05-13_manuscript_figures/qc_figures/', prefix, '_Ns_and_dash_per_sample.txt'))

  # Ns and dashes per position 
  
  Ns_and_dash_per_position = apply(dna_obj, 2, function(pos){
    char = ape::as.character.DNAbin(pos)
    c(N = sum(char == 'n'), 
      dash = sum(char == '-'))
  })
  colnames(Ns_and_dash_per_position) = paste0('pos', 1:ncol(Ns_and_dash_per_position))
  write.table(Ns_and_dash_per_position, file = paste0('2020-05-13_manuscript_figures/qc_figures/', prefix, '_Ns_and_dash_per_position.txt'))
 
  # Get variant positions 
  variant = apply(dna_obj, 2, function(pos){
    char = ape::as.character.DNAbin(pos)
    as.numeric(sum(unique(char) %in% c('a', 'c', 't', 'g')) > 1)
  })
  
  # Start writing summary report 
  sink(paste0('2020-05-13_manuscript_figures/qc_figures/', prefix, '_summary.txt'))
  print(paste('The genomes are', ncol(dna_obj), 'bp long (including Ns and dashes)'))
  print(paste('There are', nrow(dna_obj), 'genomes'))
  print(paste('There are', sum(variant), 'variant sites (defined as a position with more than one A, C, T or G'))
  print(paste('There are', sum(colSums(Ns_and_dash_per_position) > 0), 'positions with at least one N or dash in at least one sample'))
  print(paste('Thus the core genome (positions with no Ns or dashes) size is', sum(colSums(Ns_and_dash_per_position) == 0)))
  print(paste('In the core genome, there are', sum(variant[colSums(Ns_and_dash_per_position) == 0]), 'variant positions'))

  print('Number of Ns per sample')
  print(summary(Ns_and_dash_per_sample['N',]))

  print('Number of dashes per sample')
  print(summary(Ns_and_dash_per_sample['dash',]))
  sink()
  
  # Plot number of Ns and dashes per sample  
  pdf(paste0('2020-05-13_manuscript_figures/qc_figures/', prefix, '_distribution_of_N_and_dash_per_sample.pdf'))
  hist(Ns_and_dash_per_sample['N',],100, main = 'Number of Ns per sample', xlab = 'Number of Ns')
  hist(Ns_and_dash_per_sample['dash',],100, main = 'Number of dash per sample', xlab = 'Number of dashes')
  dev.off()

  # Subset alignment on variant positions & save 
  dna_var = dna_obj[,variant == 1]
  save(dna_var, file = paste0('2020-05-13_manuscript_figures/processed_data/dna_var_', prefix, '.RData'))

  # Core SNVs
  pw_dist_N_core = dist.dna(x = dna_var, model = 'N', as.matrix = TRUE, pairwise.deletion = FALSE)
  write.table(pw_dist_N_core, file = paste0('2020-05-13_manuscript_figures/processed_data/pw_dist_N_core_', prefix, '.csv'))
  
  # Noncore SNVs - pariwise comparisons 
  pw_dist_N_noncore = dist.dna(x = dna_var, model = 'N', as.matrix = TRUE, pairwise.deletion = TRUE)
  write.table(pw_dist_N_noncore, file = paste0('2020-05-13_manuscript_figures/processed_data/pw_dist_N_noncore_', prefix, '.csv'))

  # Plot core vs. noncore pairwise SNVs 
  pdf(paste0('2020-05-13_manuscript_figures/qc_figures/', prefix, '_core_vs_noncore_pw_dist.pdf'))
  hist(pw_dist_N_noncore[lower.tri(pw_dist_N_noncore, diag = FALSE)],
       breaks = seq(0,max(pw_dist_N_noncore)+10, 10),
       col = rgb(1,0,0,0.5), 
       main = 'Core vs. Noncore', 
       xlab = 'SNV dist')
  hist(pw_dist_N_core[lower.tri(pw_dist_N_core, diag = FALSE)], 
       col = rgb(0,0,1,0.5),
       breaks = seq(0,max(pw_dist_N_core)+10, 10),
       add = TRUE)
  legend('topright', fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), 
         legend = c('pwdel = T:noncore', 'pwdel = F:core'), bty = 'n')
  dev.off()
  
  
}

# ANALYSIS ----
# WHAT SAMPLES ARE USA300 
molecular_type = read.table('2020-05-13_manuscript_figures/processed_data/Bowers_class_probes_and_mlst.csv', 
           stringsAsFactors = FALSE,
           header = TRUE)
USA300_seq_num = molecular_type$SeqNum[molecular_type$molecular_type == 'USA300']

# WHOLE GENOME ALIGNMENT 
path_to_aln = 'data/2019_11_21_12_33_02_MRSA_USA_300_genome_aln_w_alt_allele_unmapped.fa'

# READ IN DNA 
dna = read.dna(path_to_aln,
               format = 'fasta')
row.names(dna) = gsub('_$','',row.names(dna))

# REMOVE DUE TO QC: 
# 192 has lots of dashes, 216 & 566 are flagged from earlier QC script
# 522 has MLST of ND and is 0 from the reference genome 
to_remove = c('MRSA_jail_216', 'MRSA_jail_566', 'MRSA_jail_192', 'MRSA_jail_522')
filepath = paste0('2020-05-13_manuscript_figures/qc_figures/', 'README_removed_genomes')
write(x = paste('Removed genomes include:', to_remove), file = filepath)
write(x = 'Flagged by Evans script, too many Ns or dashes, or MLST = ND. See code for more comments.', file = filepath, append = TRUE)
dna = dna[!row.names(dna) %in% to_remove,]

# SAVE 
save(dna, file = '2020-05-13_manuscript_figures/processed_data/dna_whole.RData')

# SUBSET USA300 and SAVE 
dna_usa300 = dna[row.names(dna) %in% USA300_seq_num, ]
save(dna_usa300, file = '2020-05-13_manuscript_figures/processed_data/dna_usa300_whole.RData')

# GET PAIRWISE DISTANCE, VARIANT POSITIONS, CORE GENOME SIZE 
get_pw_dist_and_core_size(dna, 'all')
get_pw_dist_and_core_size(dna_usa300, 'USA300')






