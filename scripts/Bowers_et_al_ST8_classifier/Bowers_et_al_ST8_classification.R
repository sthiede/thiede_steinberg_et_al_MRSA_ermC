
#### DEFINE FUNCTIONS: ####

isEmpty <- function(x) {
  return(identical(unlist(x), numeric(0)))
}

# Input string with no spaces
# Returns complement of nucleotide sequence 
comp <- function(string){
  complement = paste(sapply(unlist(strsplit(string,'')),  switch,  "A"="T", "T"="A","G"="C","C"="G"), collapse = '')
  return(complement)
  }

# DEFINE PROBES 
# As described in: 
# Bowers JR, Driebe EM, Albrecht V, McDougal LK, Granade M, Roe CC, Lemmer D, 
# Rasheed JK, Engelthaler DM, Keim P, Limbago BM. Improved Subtyping of 
# Staphylococcus aureus Clonal Complex 8 Strains Based on Whole-Genome Phylogenetic Analysis. 
# mSphere. 2018 May 2;3(3):e00464-17. 
# doi: 10.1128/mSphere.00464-17. PMID: 29720527; PMCID: PMC5932376.

where_are_the_probes <- function(dna){
  
  #### DEFINE PROBES ####
  CladeCC8 = 'ACCTATACCTGAACGTCAA'
  nonCladeCC8 = 'CTATACCTGAGCGTCAAA'
  
  innerCC8 = 'ATCGGACCCGGTAACC'
  noninnerCC8 = 'TAATCGGACCTGGTAACC'
  
  CC8a = 'ATTACTGTAGCAGGGCTG'
  nonCC8a = 'CTGTAGCAGGGTTGC'
  
  CC8b = 'AAGCTAACAAAATCACCTACTG'
  nonCC8b = 'CAAAGCTAACAAAATTACCTAC'
  
  Newlber = 'TGCACTTACATATCATCCAT'
  nonNewlber = 'CACTTACATACCATCCATC'
  
  CC8e = 'TATTAGATGAAGGCCTCAATA'
  nonCC8e = 'TTTATTAGATGAAGGCTTCAATA'
  
  CC8f = 'TAAACGTCGTAAAGTAGAACAA'
  nonCC8f = 'ACGTAAACGTCGTAAAGAAGAAC'
  
  ST239 = 'TACGACTGACCTGATGC'
  nonST239 = 'CGACTGACTTGATGCC'
  
  #### CLASSIFY ISOLATE ####
  
  
  store = c()
  # CC8 
  if (isEmpty(where(dna, CladeCC8)) == FALSE & isEmpty(where(dna, nonCladeCC8)) == TRUE){
    #print('CladeCC8')
    store = c(store, 'CladeCC8')
  }
  # CC8 complemented 
  if (isEmpty(where(dna, comp(CladeCC8))) == FALSE & isEmpty(where(dna, comp(nonCladeCC8))) == TRUE){
    #print('CladeCC8')
    store = c(store, 'CladeCC8')
  }
  # innerCC8
  if (isEmpty(where(dna, innerCC8)) == FALSE & isEmpty(where(dna, noninnerCC8)) == TRUE){
    #print('innerCC8')
    store = c(store, 'innerCC8')
  }
  # inner CC8 complemented 
  if (isEmpty(where(dna, comp(innerCC8))) == FALSE & isEmpty(where(dna, comp(noninnerCC8))) == TRUE){
    #print('innerCC8')
    store = c(store, 'innerCC8')
  }
  
  # CC8a
  if (isEmpty(where(dna, CC8a)) == FALSE & isEmpty(where(dna, nonCC8a)) == TRUE){
    #print('CC8a')
    store = c(store, 'CC8a')
  }
  # CC8a complemented 
  if (isEmpty(where(dna, comp(CC8a))) == FALSE & isEmpty(where(dna, comp(nonCC8a))) == TRUE){
    #print('CC8a')
    store = c(store, 'CC8a')
  }
  
  #CC8b
  if (isEmpty(where(dna, CC8b)) == FALSE & isEmpty(where(dna, nonCC8b)) == TRUE){
    #print('CC8b')
    store = c(store, 'CC8b')
  }
  #CC8b complemented 
  if (isEmpty(where(dna, comp(CC8b))) == FALSE & isEmpty(where(dna, comp(nonCC8b))) == TRUE){
    #print('CC8b')
    store = c(store, 'CC8b')
  }
  
  # Newlber
  if (isEmpty(where(dna, Newlber)) == FALSE & isEmpty(where(dna, nonNewlber)) == TRUE){
    #print('Newlber')
    store = c(store, 'Newlber')
  }
  # Newlber complemented
  if (isEmpty(where(dna, comp(Newlber))) == FALSE & isEmpty(where(dna, comp(nonNewlber))) == TRUE){
    #print('Newlber')
    store = c(store, 'Newlber')
  }
  
  # ST39 
  if (isEmpty(where(dna, ST239)) == FALSE & isEmpty(where(dna, nonST239)) == TRUE){
    #print('ST239')
    store = c(store, 'ST239')
  }
  # ST239 complemented
  if (isEmpty(where(dna, comp(ST239))) == FALSE & isEmpty(where(dna, comp(nonST239))) == TRUE){
    #print('ST239')
    store = c(store, 'ST239')
  }
  
  #CC8f
  if (isEmpty(where(dna, CC8f)) == FALSE & isEmpty(where(dna, nonCC8f)) == TRUE){
    #print('CC8f')
    store = c(store, 'CC8f')
  }
  #CC8f complemented 
  if (isEmpty(where(dna, comp(CC8f))) == FALSE & isEmpty(where(dna, comp(nonCC8f))) == TRUE){
    #print('CC8f')
    store = c(store, 'CC8f')
  }
  
  # CC8e
  if (isEmpty(where(dna, CC8e)) == FALSE & isEmpty(where(dna, nonCC8e)) == TRUE){
    #print('CC8e')
    store = c(store, 'CC8e')
  }
  # CC8e complemented
  if (isEmpty(where(dna, comp(CC8e))) == FALSE & isEmpty(where(dna, comp(nonCC8e))) == TRUE){
    #print('CC8e')
    store = c(store, 'CC8e')
  }
  
  return(store)
}


# BOWERS CLASSIFICATION 
# input is the output of where_are_the_probes
bowers_class <- function(seq_classification){
  
seq_classification_final = sapply(seq_classification, function(s){
  
  store = c()
  if(sum(c('innerCC8', 'CC8f', 'CC8e') %in% s) == 3 & length(s)==3){
    store = c(store, 'USA300_NAE')
  }else if(sum(c('innerCC8', 'CC8e') %in% s) == 2 & length(s)==2){
    store = c(store, 'USA500-CC8e')
  }else if(sum(c('innerCC8', 'Newlber') %in% s) == 2 & length(s)==2){
    store = c(store, 'USA500-CC8c-Newlber')
  }else if(sum(c('CC8b') %in% s) == 1 & length(s)==1){
    store = c(store, 'non-CC8')
  }else if(sum(c('CladeCC8', 'innerCC8', 'CC8b') %in% s) == 3 & length(s)==3){
    store = c(store, 'CC8b')
  }else{
    store = c(store, 'unclassified')
  }
  
  return(store)
  
  
})

return(seq_classification_final)
} 
