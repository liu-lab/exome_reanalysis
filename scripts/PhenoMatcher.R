setwd("~/exome_reanalysis/")

library(tidyverse)
library(parallel)
library(biomaRt)
library(magrittr)

select = dplyr::select

# loading MICA matrix
load("data/output/all_hpo_term.RData")
load("data/output/dz_gene_db_ref_hpo.RData")
load('data/output/MICA_mat.RData')

#### Function to load patient phenotype data. 
#   The data input format should be patient ID and HPO numbers separated by a tab. The HPO numbers are deliminated by semicolons.
#   For example:
#   100001  HP:0001290;HP:0001250;HP:0000253;HP:0001298;HP:0002126
#   100001 is the patient ID. The patient HPOs are HP:0001290;HP:0001250;HP:0000253;HP:0001298;HP:0002126. 
#   Alternatively, the HPO data can be entered into the PhenoMatcher website to generate the same results. 
LoadPhenotypeData = function(Patient_HPO_path)
{
  ptdata <- read_tsv(Patient_HPO_path,col_names = F)
  ptHPOinput <- strsplit(ptdata[[2]],";")[[1]]
  ptID <- ptdata[[1]]
  HPO_valid <- intersect(ptHPOinput, all_hpo_term)
  list(HPO_valid, ptID)
}

# Given a set of terms annotated to a patient and a set associated with a disease, compute the Resnik similarity between the two sets.
compare_term_sets = function(annotated_set, disease_set)
{
  if (length(disease_set) == 0) {return(0)}
  
  sub_mica_1 = MICA_mat[annotated_set, disease_set, drop = FALSE]
  resnik_1 = mean(apply(sub_mica_1, 1, max))
  
  sub_mica_2 = MICA_mat[disease_set, annotated_set, drop = FALSE]
  resnik_2 = mean(apply(sub_mica_2, 1, max))
  
  round(mean(c(resnik_1, resnik_2)),digits = 3)
}


# Generate semantic matching score for all disease genes in relation to the patient's phenotypes. 
GeneratePenotypeScores = function(Pt_processed_HPO_data)
{
  Pt_HPO_processed <- Pt_processed_HPO_data[1]
  Pt_ID <- Pt_processed_HPO_data[2]
  PtHPO_DzHPO_compare <- dz_gene_db_ref_hpo %>%
                        rowwise() %>%
                          mutate(
                                  Pt_HPO_valid= Pt_HPO_processed,
                                  PhenoMatch_score = compare_term_sets(hpo_id,Pt_HPO_valid )) %>%
                          select(-hpo_id, -Pt_HPO_valid) %>%
                          group_by(entrez_gene_symbol) %>%
                          arrange(desc(PhenoMatch_score)) %>%
                          mutate(disease_id_max = dplyr::first(disease_id),
                                 PhenoMatch_score_max = dplyr::first(PhenoMatch_score),
                                disease_name_max = dplyr::first(db_name)) %>%
                           group_by( entrez_gene_symbol, disease_id_max, PhenoMatch_score_max ) %>%
                           dplyr::summarise(
                                  dz_ID_all = paste(disease_id, collapse  = ";"),
                                  scores = paste(PhenoMatch_score, collapse  = ";"))  %>%
                          arrange(desc(PhenoMatch_score_max)) %>% 
                           mutate(ID = unlist(Pt_ID),
                                  Patient_HPO = paste(unlist(Pt_HPO_processed), collapse = ";") )
  
      # the output phenotype matching file is in the same format of the file generated from the PhenoMatcher website. 
      write.csv(PtHPO_DzHPO_compare, 
                file=paste0('data/test_data/output/',Pt_ID, '_PhenoMatcher_output.csv'), row.names = FALSE)

}


LoadPhenotypeData("data/test_data/input/Phenotype_100001.txt") %>% GeneratePenotypeScores(.)
LoadPhenotypeData("data/test_data/input/Phenotype_100002.txt") %>% GeneratePenotypeScores(.)
LoadPhenotypeData("data/test_data/input/Phenotype_100003.txt") %>% GeneratePenotypeScores(.)

