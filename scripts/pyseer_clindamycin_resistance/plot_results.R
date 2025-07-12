results = read.delim('2021-05-10_pyseer_cdc_pangenome_vs_clindaR/data/model_results/LMM_panaroo_vs_clindaR.txt')
results$OR = exp(results$beta)
results$log10p = -1 * log10(results$lrt.pvalue)
results$sig_gene = results$variant
results$sig_gene[results$log10p  < 50]= ''

library(ggplot2)

pdf('2021-05-10_pyseer_cdc_pangenome_vs_clindaR/figures/cdc_pyseer_clindaR_volcano_plot.pdf', width = 8, height = 5)
ggplot(results, aes(x =OR, y = log10p)) + 
  geom_point(alpha = 0.5) +  
  theme_bw() + 
  theme(text = element_text(size = 16)) + 
  geom_text(aes(label = sig_gene), nudge_x = -0.15) + 
  ylab('-log10 p value') + 
  xlab('OR') 
dev.off()

