setwd("/mnt/bigData2/collaborators/udn/pengfeil/exome_reanalysis/")

library(tidyverse)
library(parallel)
library(biomaRt)
library(magrittr)
require(VariantAnnotation)

select = dplyr::select

########## FUNCTIONS ################

#### Function to load and transform vcf file into a table
MakeTableFromVCF = function(vcf_path)
{
  # Load annotated filtered VCF file from patient, with the readVcf function from the VariantAnnotation R package
  vcf_file <- readVcf(vcf_path)
  vcf_table <- as.tibble(data.frame(rowRanges(expand(vcf_file)), geno(vcf_file)$GT, 
                                  geno(vcf_file)$VR, geno(vcf_file)$DP,info(expand(vcf_file)), stringsAsFactors = F))

  # Cleaning up and formatting the loaded variant file
  vcf_table$GI  <- sapply(vcf_table$GI, function(x) {paste(x,collapse = ",")})
  vcf_table <- vcf_table %>% mutate_at(.vars = vars(GN,DT,NS,DS,IHM,IHB,PLP,HGMD,HGMD_same_codon,LMDDD), 
                                     .funs = unlist)
  vcf_table <- vcf_table %>% mutate_at(.vars = vars(EXAC_AC_HET,EXAC_AC_HMZ), .funs = as.integer)
  vcf_table <- vcf_table %>% mutate_at(.vars = vars(PLI,MSZ), .funs = as.numeric)

  PtVrnt_anno <- vcf_table %>% dplyr::select(-c(end,width,strand,paramRangeID))
  colnames(PtVrnt_anno) <- c("CHROM","POS","REF","ALT","QUAL","FILTER","Zygosity","Variant_read","Total_read",
                           "Gene_name","Gene_information","Deleterious","Novel_seen","Disease",
                           "Monoallelic_inheritance", "Biallelic_inheritance",
                           "EXAC_AC_HET","EXAC_AC_HMZ","Previous_pathogenic","HGMD","HGMD_same_codon",
                           "pLI","Mis_Z","Monoallelic_LoF_DDD")
  PtVrnt_anno <-  PtVrnt_anno %>% 
                    mutate(PtID= samples(header(vcf_file))) # write patient ID into the file
  PtVrnt_anno
}

###### Function to load phenotype matching data. 
# This can either be outputs generated from the PhenoMatcher script module or outputs generated from the website.
MergePhenotype = function(Variant_table, Phenotype_file_path)
{
  PtHPO_DzHPO_compare <- read_csv(Phenotype_file_path)
  # joining the phenotype file and the variant table
  Variant_table_phenotype_merged <- Variant_table %>% 
            left_join(PtHPO_DzHPO_compare, by=c('Gene_name' = 'entrez_gene_symbol'))  # Merging the variant table with phenotype matching table. 
  Variant_table_phenotype_merged
}

##### Function to compute monoallelic and biallelic hypothesis and prioritize variants
# The function input "Parental_sample_availability" asks for a TURE or FALSE logical vector. 
# This logical vector informs prioritization rules, as is explained below. 


Compute_Hypothesis = function(Variant_table_phenotype_merged, Parental_sample_availability)
{
  # The Monoallelic hypothesis filters for variants that are
  # 1. located in genes with known monoallelic disease inheritance, and 
  # 2. rare in the ExAC database (<=5 heterozygous counts), or not rare in ExAC (>5 heterozygous counts) but previously classified as pathogenic or likely pathogenic in our internal database, and
  # 3. fulfilling one of the four following rules
  #     a. previously reported as pathogenic or likely pathogenic internally
  #     b. truncating variants with likely haploinsufficiency disease mechanism (defined as a stringent pLI cutoff of >0.99 or classified as loss-of-function disease mechanism by the DDD group)
  #     c. variants that are not seen in the control database but seen in patients from the literature (or allelic variants seen in patients)
  #     d. variants that are not seen in either control databases or in patients; however, the gene is associated with likely haploinsufficiency disease mechanism (same as mentioned above) 
  #         or predicted to be intolerant to missense variations (missense intolerant Z score from the ExAC database > 3). To limit the number of variants included from this subgroup, the phenotype 
  #         matching score has to be higher than a certain cutoff (set as 0.3 in this study); genes without HPO annotation and therefore no phenotype matching score are not excluded. This subgroup of 
  #         variants, based on the ACMG/AMP variant classification guideline, likely needs extra evidence to be classified as pathogenic or likely pathogenic. The extra piece of evidence usually is that 
  #         the variant is found to be of de novo origin after parental studies, which requires availability of parental samples for proband only exome analysis. Therefore, an additional logic is added 
  #         here: this subset of variants will only be retained if both parental samples are available for subsequent targeted variant study. 
  
  vrnt_df_mono_filt <- Variant_table_phenotype_merged %>% 
    filter( Monoallelic_inheritance==TRUE & 
              (EXAC_AC_HET<=5 | Previous_pathogenic==TRUE) &
              (Previous_pathogenic==TRUE |
                 (Deleterious=="deleterious" & (pLI > 0.99 | Monoallelic_LoF_DDD==TRUE) ) |
                 (Novel_seen=="novel" & (HGMD==TRUE | HGMD_same_codon==TRUE)) | 
                 (Novel_seen=="novel" & 
                    ((pLI > 0.99 | Monoallelic_LoF_DDD==TRUE) | Mis_Z > 3) & 
                    (PhenoMatch_score_max>0.3 | is.na(PhenoMatch_score_max) ) & Parental_sample_availability)) )  
  
  # The Biallelic hypothesis filters for variants that are 
  #       located in genes with known biallelic disease inheritance, and 
  #       exclude variants that are located in genes with low phenotype matching scores. 
  #   It then filters variants that are
  #       a. previously classified as pathogenic or likely pathogenic internally
  #       b. truncating variants 
  #       c. rare in the ExAC database 
  #       d. potentially biallelic. The current pipeline is designed for proband only exome sequencing. So, a "potential biallelic" status is defined as more than one variant appearing in the same gene 
  #         or variants being hemizygous or homozygous.
  
  vrnt_df_bi_filt <- Variant_table_phenotype_merged %>% 
    filter(Biallelic_inheritance==TRUE & 
             (PhenoMatch_score_max>0.3 | is.na(PhenoMatch_score_max) ) & 
             ( Previous_pathogenic==TRUE |  Deleterious=="deleterious" | EXAC_AC_HET<=5) )  %>% 
    group_by(Gene_name) %>% filter( n()>1 | Zygosity=="1/1") %>% ungroup
  
  
  Variants_prioritized <- bind_rows(vrnt_df_mono_filt, vrnt_df_bi_filt) %>%
    distinct(PtID, CHROM, POS, REF, ALT, .keep_all = TRUE)
  
  Variants_prioritized
}


########### Reanalysis pipeline ###############

MakeTableFromVCF("data/test_data/input/100001.vcf") %>%
  MergePhenotype("data/test_data/output/100001_PhenoMatcher_output.csv") %>%
  Compute_Hypothesis(Parental_sample_availability=FALSE)


MakeTableFromVCF("data/test_data/input/100002.vcf") %>%
  MergePhenotype("data/test_data/output/100002_PhenoMatcher_output.csv") %>%
  Compute_Hypothesis(Parental_sample_availability=TRUE)


MakeTableFromVCF("data/test_data/input/100003.vcf") %>%
  MergePhenotype("data/test_data/output/100003_PhenoMatcher_output.csv") %>%
  Compute_Hypothesis(Parental_sample_availability=FALSE)
