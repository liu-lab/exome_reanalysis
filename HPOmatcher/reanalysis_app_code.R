library(tidyverse)
library(parallel)
library(biomaRt)
library(magrittr)

select = dplyr::select


####get gene-disease-HPO from the HPO website
genes_to_diseases = read_tsv('data/genes_to_diseases.txt',
                             skip = 1,
                             col_names = c('entrez_gene_id', 'entrez_gene_symbol', 'disease_id'))
disease_to_hpo = read_tsv('data/phenotype_annotation_hpoteam.tab', 
                          col_names = c('db',
                                        'db_object_id',
                                        'db_name',
                                        'qualifier',
                                        'hpo_id',
                                        'db_reference',
                                        'evidence_code',
                                        'onset_modifier',
                                        'frequency_modifier',
                                        'with',
                                        'aspect',
                                        'synonym',
                                        'date',
                                        'assigned_by',
                                        'frequency_description')) 

####get the link for gene symbol-entrez gene id 
entrez_to_hgnc = read_tsv('data/entrez_to_hgnc.tsv') %>% select(`Approved Symbol`, `Entrez Gene ID`)
names(entrez_to_hgnc) = names(entrez_to_hgnc) %>% map_chr(~gsub(' ', '_', .x))

####functions###
load('data/MICA_mat.RData')

compare_term_sets = function(annotated_set, disease_set){
  
  if (length(disease_set) == 0) {return(0)}
  
  sub_mica_1 = MICA_mat[annotated_set, disease_set, drop = FALSE]
  resnik_1 = mean(apply(sub_mica_1, 1, max))
  
  sub_mica_2 = MICA_mat[disease_set, annotated_set, drop = FALSE]
  resnik_2 = mean(apply(sub_mica_2, 1, max))
  
  round(mean(c(resnik_1, resnik_2)),digits = 3)
}

rank_patients_diseases_by_gene = function(pt_data){
  pt_hpo_annotations = strsplit(as.character(pt_data$annotated_terms), ";") %>% unlist %>% unique
  
  if (length(pt_hpo_annotations) == 0) {return(data_frame(db_reference = vector(mode = 'character',
                                                                                length = 0),
                                                          annotation_disease_similarity = vector('numeric',
                                                                                                 length = 0)))}
  pt_data %>% 
    select(Gene_name, disease_terms) %>% 
    unnest() %>% 
    group_by(Gene_name, db_reference) %>% 
    summarise(annotation_disease_similarity = compare_term_sets(pt_hpo_annotations, hpo_id)) %>% 
    arrange(desc(annotation_disease_similarity)) %>% 
    ungroup
  
  
}


hpo_terms_from_diseases = function(gene_diseases_tbl){
  disease_to_hpo %>% 
    filter(db_reference %in% gene_diseases_tbl$disease_id) %>% 
    unique
}




