
#############################################################################################################
##                                                                                                        ##
## The code below is developed by graduate student Andrew Ghazi Andrew.Ghazi@bcm.edu and his thesis       ##
## supervisor Dr. Chad Shaw cashaw@bcm.edu;  Mr. Ghazi is a student in Dr. Shaw’s laboratory, and the     ##
## code for this work is a re-implementation and extension of the algorithm from a previous publication   ##
## developed in Dr. Shaw’s laboratory at BCM [James et al. A visual and curatorial approach to clinical   ##
## variant prioritization and disease gene discovery in genome-wide diagnostics.(2016) Genome Med. 8:13.] ##
##                                                                                                        ##
############################################################################################################

## This module of the codes require much longer compute time compared to the other module of the codes.  
## One can generate a MICA matrix file and use this file for a period of time for the downstream analyses. 
## It is only necessary to re-run this code when one desires to incorportae updates from disesase databases. 
##
## The two files generated from this analysis below, MICA_mat.RData and pairs_MICA.RData are too large 
##   to be uploaded to GitHub. They can be downloaded in the following link: https://bcm.box.com/v/exome-reanalysis

setwd("~/exome_reanalysis/")


###############################################################
#################   Prepare MICA matrix    ####################
###############################################################

library(tidyverse)
library(ontologyIndex)
library(parallel)
library(magrittr)

############1. get hpo terms and the ancestors for each hpo term
hpo = get_ontology('data/input/hp.obo',extract_tags = 'everything') ##data-version: releases/2017-12-12

id_to_term = data_frame(id = hpo$id,
                        term = hpo$name)
save(id_to_term, file = 'data/output/hpo_id_to_term.RData')

hpo_ancestry = hpo$ancestors %>%
  data_frame(ID = names(.),
             ancestors = .)
save(hpo_ancestry, file = 'data/output/hpo_ancestry.RData')

############2. get the association between disease and hpo terms
############http://human-phenotype-ontology.github.io/documentation.html

disease_to_hpo = read_tsv('data/input/phenotype_annotation_hpoteam.tab',
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

############3. calculate the information content (IC) for each hpo ancestry

####build a n*m matrix, where each of n is a unique disease and each of m is a hpo ancestry
disease_term_mat = disease_to_hpo %>%
  dplyr::select(db_name, hpo_id) %>%
  unique %>%
  left_join(hpo_ancestry, by = c('hpo_id' = 'ID')) %>%
  dplyr::select(-hpo_id) %>%
  unnest %>%
  mutate(value = 1) %>%
  unique %>%
  spread(ancestors, value, fill = 0)%>%
  dplyr::select(-db_name) %>%
  as.matrix

####calculate IC using IC(m)= -log(f(m)), where m is a hpo ancestry,
####f(m) is the mean value of the frequency for m being obeserved across all diseaases
IC_vec = -log(colMeans(disease_term_mat))
save(IC_vec, file = 'data/output/IC_vec.RData')


############4. Calculate the most informative common ancestor (MICA) for each pair of HPO terms

###function to get all the common ancestors for a pair of hpo terms

get_ancestor_overlap = function(id_1, id_2){
  # get_max_ic now only calls this function if id_1 != id_2
  
  anc_list = hpo_ancestry %>% 
    filter(ID %in% c(id_1, id_2)) %>% 
    .$ancestors
  
  intersect(anc_list[[1]], anc_list[[2]])
  
}

###function to calculate MICA using the maximum value of IC for all common ancestors


get_max_unannotated_ic = function(unannotated_term){
  # for diagonal elements where the term isn't annotated, just return the
  # maximum IC of any ancestor
  
  term_ancestors = hpo_ancestry %>% 
    filter(ID == unannotated_term) %>% 
    pull(ancestors) %>% 
    .[[1]]
  
  max(IC_vec[term_ancestors], na.rm = TRUE)
}

get_max_ic = function(id_1, id_2) {
  
  # each term is an ancestor of itself, so if the terms are the same (the
  # diagonal) just return the IC of the term
  if (id_1 == id_2) {
    if (id_1 %in% names(IC_vec)) {
      return(IC_vec[id_1]) 
    } else {
      # if the term isn't in the IC vector (because it didn't show up in the
      # disease database, return the maximum IC of any of its ancestors)
      return(get_max_unannotated_ic(id_1))
    }
    
  }
  
  overlapping_terms = get_ancestor_overlap(id_1, id_2)
  
  if (length(overlapping_terms) == 0) return(0)
  max(IC_vec[names(IC_vec) %in% overlapping_terms])
}

try_get_max_ic = function(id_1, id_2) {
  try(get_max_ic(id_1, id_2))
}


####generate a matrix with all pairs of hpo ancestors, each row is a pair
ancestor_pairs = combn(hpo_ancestry$ID,2) %>%
  t %>%
  as_tibble

####calculate the MICA for all pairs of hpo terms and ancestors
pairs_MICA = ancestor_pairs %>%
  mutate(max_shared_ancestor_IC = mcmapply(try_get_max_ic, 
                                           V1, V2, 
                                           mc.cores = 20, 
                                           SIMPLIFY = TRUE))

diag_pairs = data_frame(V1 = hpo$id,
                        V2 = hpo$id) %>% 
  mutate(max_shared_ancestor_IC = mcmapply(try_get_max_ic, 
                                           V1, V2, 
                                           mc.cores = 20, 
                                           SIMPLIFY = TRUE))

pairs_MICA = pairs_MICA %>%
  bind_rows(diag_pairs)

save(pairs_MICA, file = 'data/output/pairs_MICA.RData')

MICA_mat = pairs_MICA %>%
  spread(V2, value = max_shared_ancestor_IC, fill = NA) %>%
  as.data.frame()

rownames(MICA_mat) = MICA_mat$V1
MICA_mat$V1 = NULL
MICA_mat %<>% as.matrix
MICA_mat[lower.tri(MICA_mat)] = t(MICA_mat)[lower.tri(MICA_mat)]

save(MICA_mat,
     file = 'data/output/MICA_mat.RData')

####get disease-HPO from the HPO website
genes_to_diseases = read_tsv('data/input/genes_to_diseases.txt',
                            skip = 1,
                             col_names = c('entrez_gene_id', 'entrez_gene_symbol', 'disease_id'))

####get the link for gene symbol-entrez gene id 
entrez_to_hgnc = read_tsv('data/input/entrez_to_hgnc.tsv') %>% 
  dplyr::select(`Approved Symbol`, `Entrez Gene ID`)   # https://www.genenames.org/cgi-bin/download
names(entrez_to_hgnc) = names(entrez_to_hgnc) %>% map_chr(~gsub(' ', '_', .x))

dz_gene_db_ref_hpo <- merge(genes_to_diseases, disease_to_hpo, by.x = "disease_id", by.y = "db_reference")
dz_gene_db_ref_hpo2 <- unique(dz_gene_db_ref_hpo[,c(1:3,6,8)]) %>%
group_by(disease_id, entrez_gene_id, entrez_gene_symbol,db_name) %>% 
summarise(hpo_id = list(unique(hpo_id)))

dz_gene_db_ref_hpo<- dz_gene_db_ref_hpo2
save(dz_gene_db_ref_hpo, file = "data/output/dz_gene_db_ref_hpo.RData")

all_hpo_term <- unique(unlist(disease_to_hpo$hpo_id))
save(all_hpo_term, file="data/output/all_hpo_term.RData")
