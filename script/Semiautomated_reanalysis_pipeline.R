
##############################################################################################################
##                                                                                                          ##
## Section II of the code is developed by graduate student Andrew Ghazi Andrew.Ghazi@bcm.edu and his        ##
## thesis supervisor Dr. Chad Shaw cashaw@bcm.edu;  Mr. Ghazi is a student in Dr. Shaw’s laboratory, and    ##
## the code for this work is a re-implementation and extension of the algorithm from a previous publication ##
## developed in Dr. Shaw’s laboratory at BCM [James et al. A visual and curatorial approach to clinical     ##
## variant prioritization and disease gene discovery in genome-wide diagnostics.(2016) Genome Med. 8:13.]   ##
##                                                                                                          ##
##############################################################################################################



####################################################
##############I. Load input files ##################
####################################################

setwd("~/exome_reanalysis/")

####load phenotypes and variants from mock patient data
####three patient mock data are available
####The scripts below processes data from Pt100001
load("data/test_data/MockPtData.RData")

PtVrnt_anno <- Pt100001
#PtVrnt_anno <- Pt100002
#PtVrnt_anno <- Pt100003



#######resource data
Constraint <- read.table("data/input/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt",header=T,stringsAsFactors=F)
Constraint <-Constraint[,c("gene","pLI","mis_z")]
HI <- read.csv("data/input/HI_Predictions_Version3.csv",header=T,stringsAsFactors=F)
colnames(HI) <- c("Gene","HI_score")
DDD <- read.csv("data/input/DDG2P_1_9_2017.csv",header=T,stringsAsFactors=F)


###############################################################
##############II. Annotate a HPO match score ##################
###############################################################

library(tidyverse)
library(ontologyIndex)
library(parallel)
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

phenotypes = read_tsv('data/input/phenotype_annotation_hpoteam.tab',
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
disease_term_mat = phenotypes %>%
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
    anc_list = hpo_ancestry[c(id_1, id_2), ]$ancestors
    intersect(anc_list[[1]], anc_list[[2]])
}

###function to calculate MICA using the maximum value of IC for all common ancestors
get_max_ic = function(id_1, id_2) {
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
mutate(max_shared_ancestor_IC = mcmapply(try_get_max_ic, V1, V2, mc.cores = 20))

diag_pairs = data_frame(V1 = hpo$id,
V2 = hpo$id) %>%
mutate(max_shared_ancestor_IC = mcmapply(try_get_max_ic, V1, V2, mc.cores = 20))

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

############5. Semantic matching

library(tidyverse)
library(parallel)
library(biomaRt)
library(magrittr)

select = dplyr::select
setwd("~/exome_reanalysis/")

####get gene-disease-HPO from the HPO website
genes_to_diseases = read_tsv('data/input/genes_to_diseases.txt',
skip = 1,
col_names = c('entrez_gene_id', 'entrez_gene_symbol', 'disease_id'))
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

####get the link for gene symbol-entrez gene id
entrez_to_hgnc = read_tsv('data/input/entrez_to_hgnc.tsv') %>% select(`Approved Symbol`, `Entrez Gene ID`)
names(entrez_to_hgnc) = names(entrez_to_hgnc) %>% map_chr(~gsub(' ', '_', .x))

####functions###
load('data/output/MICA_mat.RData')

compare_term_sets = function(annotated_set, disease_set){
    
    if (length(disease_set) == 0) {return(0)}
    
    sub_mica_1 = MICA_mat[annotated_set, disease_set, drop = FALSE]
    resnik_1 = mean(apply(sub_mica_1, 1, max))
    
    sub_mica_2 = MICA_mat[disease_set, annotated_set, drop = FALSE]
    resnik_2 = mean(apply(sub_mica_2, 1, max))
    
    round(mean(c(resnik_1, resnik_2)),digits = 3)
}

rank_patients_diseases_by_gene = function(pt_data){
  pt_hpo_annotations = pt_data$annotated_terms %>% unlist %>% unique
  
  pt_data %>% 
    select(Gene_name, disease_terms) %>% 
    unnest() %>% 
    group_by(Gene_name, db_reference) %>% 
    summarise(annotation_disease_similarity = ifelse(length(pt_hpo_annotations)>0,compare_term_sets(pt_hpo_annotations, hpo_id),NA)) %>% 
    arrange(desc(annotation_disease_similarity)) %>% 
    ungroup
  
  
}


####annotate disease id, diease hpo terms, patient hpo terms for each variant gene

hpo_terms_from_diseases = function(gene_diseases_tbl){
    disease_to_hpo %>%
    filter(db_reference %in% gene_diseases_tbl$disease_id) %>%
    unique
}


variants_with_diseases = PtVrnt_anno %>%
left_join(entrez_to_hgnc, by = c('Gene_name' = 'Approved_Symbol')) %>%
mutate(gene_diseases = map(Entrez_Gene_ID, ~filter(genes_to_diseases, entrez_gene_id == .x))) %>%
mutate(disease_terms = mclapply(gene_diseases, hpo_terms_from_diseases, mc.cores = 20))%>%
select(PtID,Gene_name, Entrez_Gene_ID, gene_diseases, disease_terms) %>%
left_join(ptHPOdata, by = c("PtID"="dnaID")) %>%
rename(annotated_terms = all_terms)


####rank genes based on the HPO match score

disease_rankings_by_patient_by_gene = variants_with_diseases  %>%
filter(map_lgl(disease_terms, ~nrow(.x) > 0)) %>%
select(-gene_diseases) %>%
group_by(PtID) %>%
nest(.key = 'patient_data') %>%
mutate(ranked_diseases = mclapply(patient_data, rank_patients_diseases_by_gene, mc.cores = 20)) %>%
unnest(ranked_diseases) %>%
group_by(PtID, Gene_name) %>%
arrange(desc(annotation_disease_similarity)) %>%
mutate(dz_ID_max = dplyr::first(db_reference),
HPOmatch_score_max = dplyr::first(annotation_disease_similarity)) %>%
group_by(PtID, Gene_name, dz_ID_max,  HPOmatch_score_max ) %>%
dplyr::summarise(
dz_ID_all = paste(db_reference, collapse  = ";"),
HPOmatch_scores = paste(annotation_disease_similarity, collapse  = ";")) %>%
right_join(PtVrnt_anno, by = c("PtID", "Gene_name"))

save(disease_rankings_by_patient_by_gene,
file = 'data/test_data/pt_disease_rankings_by_patient.RData')

######################################################
##############III. VARIANT FILTERING #################
######################################################
setwd("~/exome_reanalysis/")
load('data/test_data/pt_disease_rankings_by_patient.RData')
vrnt_df <- unique(disease_rankings_by_patient_by_gene)
rm(disease_rankings_by_patient_by_gene)

############1. Exclude variants

####1.0 exclude variants that are likely benign based on an in-house scoring. 
#### filter for variants with previous reports, deleterious effect, being potential compound heterozygous, or novel. 
vrnt_df <- vrnt_df %>%
  filter(CommonVrnt.FP!="TRUE" & Disease=="TRUE" & 
           (PMID=="TRUE" | PMID_same_codon=="TRUE" | Potential.cmphet=="TRUE" | Novel_seen=="novel" | Deleterious=="deleterious")
  ) 

####1.1 exclude all variants with any homo ExAC count
Exclude.Idx.with.hmzExAC <- which(vrnt_df$Counts.hmzExAC>0)

####1.2 exclude all variants with >5 het ExAC count and not being compound het or hom or hem
Exclude.Idx.with.hetExAC <- which(vrnt_df$Counts.hetExAC>5 &
vrnt_df$Potential.cmphet==FALSE &
!vrnt_df$Zygosity %in% c("Hom", "Hem"))

####1.3 exclude all variants with AR inheritance but with only one hit

Exclude.Idx.with.AR <- which(vrnt_df$Inheritance %in% c("AR") &
vrnt_df$Potential.cmphet==FALSE &
!vrnt_df$Zygosity %in% c("Hom", "Hem") )

####1.4 exclude all variants previously scored as likely benign or benign for more than one time, or known VUS for three or more times or novel for five or more times
Exclude.Idx.previous.scores <- unique(which(vrnt_df$Counts.previous.B.LB>=1|vrnt_df$Counts.previous.VUS_seen>=3|vrnt_df$Counts.previous.VUS_novel>=5))

####1.5 retrive the previously reported as pathogenic or likely pathogenic
Exclude.Idx <- unique(c(Exclude.Idx.with.hmzExAC,Exclude.Idx.with.hetExAC, Exclude.Idx.with.AR, Exclude.Idx.previous.scores))

Filt.Idx <- c(1:nrow(vrnt_df))[-setdiff(Exclude.Idx, which(vrnt_df$Previous.P.LP==TRUE))]


vrnt_df.filt <- vrnt_df[Filt.Idx,]

############2. Asign an inclusion criteria for each vriant

Idx.Previous.P.LP <- which(vrnt_df.filt$Previous.P.LP==TRUE)
#### cluster variants based on variant type-- truncating variants or non-truncating variants
Idx.trunc <- which(vrnt_df.filt$Deleterious=="deleterious")

Idx.nontrunc <- which(is.na(vrnt_df.filt$Deleterious) | vrnt_df.filt$Deleterious== "VUS")

#### cluster variants based on disease inheritance model
Idx.dom <- which(vrnt_df.filt$Inheritance %in% c("AD","\\?AD", "AD\\?", "XL","XLD", "X-linked", "X-linked dominant", "XLR", "YL", "Epigenetic", "Epigenetic/imprinted", "IP", "AD/IM") )
Idx.domrec <- which( vrnt_df.filt$Inheritance %in% c("AD","\\?AD", "AD\\?", "XL","XLD", "X-linked", "X-linked dominant", "XLR", "YL", "Epigenetic", "Epigenetic/imprinted", "IP", "AD/AR",
"AD/AR/Complex", "AD/AR/Digenic", "AR/AD","AD/digenic", "AD/Complex", "AD/IM", "?AD/AR", "AD/AR/DD" ) )

Idx.dom.trunc <- intersect(Idx.dom,Idx.trunc)
Idx.dom.nontrunc <- intersect(Idx.dom,Idx.nontrunc)
Idx.domrec.nontrunc <- intersect(Idx.domrec, Idx.nontrunc)

#### cluster variants based on novelty, using control databases, inhouse-notes, associated PubMed id and HGMD
Idx.novel <- which(vrnt_df.filt$Novel_seen== "novel")

Idx.PMID <- which(vrnt_df.filt$PMID==TRUE)

Idx.PMIDsameCodon <- which(vrnt_df.filt$PMID_same_codon == TRUE & vrnt_df.filt$PMID==FALSE) ## A different variant affecting the identical amino acid with literature report

Idx.noPMID <- which(vrnt_df.filt$PMID_same_codon == FALSE & vrnt_df.filt$PMID==FALSE)


#### cluster variants based on the availibility of a parent-child trio
Idx.withparents <- which(vrnt_df.filt$withparents==TRUE)

#### cluster variants based on the HPO similarity between the disease gene associated phenotypes and the patient pheynotypes
# 0.3 is the lowest score for a diagnositic vrnt from the cohort1
Idx.SemanticSim.pass <- which(vrnt_df.filt$HPOmatch_score_max>0.3 | is.na(vrnt_df.filt$HPOmatch_score_max))
Idx.SemanticSim.high <- which(vrnt_df.filt$HPOmatch_score_max>1.8) # 1.8 is the 99% percentile of all variants' score

#### clustering variants based on pLI, mis_z, HI and DDD data

CheckDDD <- function(x){
    if (x %in% DDD$gene.symbol){
        if("loss of function" %in% DDD[which(DDD$gene.symbol==x & DDD$allelic.requirement %in% c("monoallelic", "x-linked dominant", "mosaic", "x-linked over-dominance", "hemizygous")), "mutation.consequence"] )
        {TRUE} else {FALSE}
    }	else {FALSE}
}


Idx.pLI <- which( (vrnt_df.filt$pLI>0.95 & vrnt_df.filt$HI_score<0.1) |
vrnt_df.filt$pLI>0.99 |
do.call(c,lapply(vrnt_df.filt %>% .$Gene_name, CheckDDD)))

Idx.highZ <- which(vrnt_df.filt$mis_z >3)

#### assign an initial value of the inclusion criteria
vrnt_df.filt$Inclusion.criteria<-0

#### assign 99 for variants with a high score of HPO match
if(length(Idx.SemanticSim.high)>0){
    vrnt_df.filt$Inclusion.criteria[Idx.SemanticSim.high] <- 99}


#### asign 11 for variants that have been reported as pathogenic or likely-pathogenic
if(length(Idx.Previous.P.LP)>0){
    vrnt_df.filt$Inclusion.criteria[Idx.Previous.P.LP] <- 11}

#### assign 12 for truncating variants that have a promising score
Idx.trunc.pLI<- intersect(Idx.trunc, Idx.pLI)
if(length(Idx.trunc.pLI )>0){
    vrnt_df.filt$Inclusion.criteria[Idx.trunc.pLI] <- 12}

#### assign 13 for variants that are dominant, non-truncating, novel, with PMID
Idx.domrec.nontrunc.PMID.novel <- Reduce(intersect, list(Idx.domrec.nontrunc , Idx.PMID, Idx.novel))
if(length(Idx.domrec.nontrunc.PMID.novel)>0){
    vrnt_df.filt$Inclusion.criteria[Idx.domrec.nontrunc.PMID.novel] <- 13}

#### assign 14 for variants that are dominant, non-truncating, novel, no PMID but with previously reported same codon variant
Idx.domrec.nontrunc.noPMID.novel.PMIDsameCodon <- Reduce(intersect, list(Idx.domrec.nontrunc , Idx.PMIDsameCodon, Idx.novel))
if(length(Idx.domrec.nontrunc.noPMID.novel.PMIDsameCodon)>0){
    vrnt_df.filt$Inclusion.criteria[Idx.domrec.nontrunc.noPMID.novel.PMIDsameCodon] <- 14}

####assign 15 for variants that are missense pLI high
Idx.nontrunc.IDwithparents.novel.highpLI <- Reduce(intersect, list(Idx.nontrunc, Idx.withparents,
Idx.novel,Idx.SemanticSim.pass, Idx.pLI))
if(length(Idx.nontrunc.IDwithparents.novel.highpLI)>0){
    vrnt_df.filt$Inclusion.criteria[Idx.nontrunc.IDwithparents.novel.highpLI] <- 15}

#### assign 16 for varaitns that are missense pLI low, high missense Z score
Idx.nontrunc.IDwithparents.lowpLI.highZ <- Reduce(intersect, list(Idx.nontrunc,Idx.withparents,Idx.novel,
Idx.SemanticSim.pass, Idx.highZ))
if(length(Idx.nontrunc.IDwithparents.lowpLI.highZ)>0){
    vrnt_df.filt$Inclusion.criteria[Idx.nontrunc.IDwithparents.lowpLI.highZ] <- 16}

############# recessive model

Inheritance.AR <- DDD$gene.symbol[DDD$allelic.requirement %in% c("biallelic", "digenic")]

vrnt_df.filt.unique.rec.gene <- vrnt_df.filt %>% 
filter(Gene_name %in% Inheritance.AR) %>%
.$Gene_name %>% unique

Idx.rec <- which(vrnt_df.filt$Gene_name %in% Inheritance.AR)

Idx.rec.mix <- if( length(Idx.rec) >0 ){
    do.call(c, lapply(vrnt_df.filt.unique.rec.gene, function(x) {
        
        Idx.AR <- which(vrnt_df.filt$Gene_name == x)
        
        Idx.AR.deleterious <- intersect(Idx.AR, Idx.trunc)
        
        Idx.AR.PMID <- intersect(Idx.AR, Idx.PMID)
        
        Idx.AR.PMIDsameCodon <- intersect(Idx.AR, Idx.PMIDsameCodon)
        
        Idx.AR.novelPMID <- Reduce(intersect,list(Idx.AR, Idx.novel, Idx.PMID))
        
        Idx.AR.novelPMIDsameCodon <- Reduce(intersect,list(Idx.AR, Idx.novel, Idx.PMIDsameCodon))
        
        Idx.AR.novelnoPMID <-  Reduce(intersect,list(Idx.AR, Idx.novel, Idx.noPMID))
        
        Idx.AR.exaclow <- setdiff(Idx.AR, unique(c(Idx.AR.deleterious,Idx.AR.PMIDsameCodon,Idx.AR.novelPMID,Idx.AR.novelPMIDsameCodon,Idx.AR.novelnoPMID)))
        
        Idx.rec.twotrunc<- NULL
        Idx.rec.twonovelPMID<- NULL
        Idx.rec.twoPMID<- NULL
        Idx.rec.novelPMIDsameCodon<- NULL
        Idx.rec.onetrunconenovelPMID<- NULL
        Idx.rec.onetrunconePMID<- NULL
        Idx.rec.onetrunconenovelPMIDsameCodon<- NULL
        Idx.rec.onetrunconenovel<- NULL
        Idx.rec.onenovelPMIDonePMID<- NULL
        Idx.rec.onenovelPMIDonenovelPMIDsameCodon<- NULL
        Idx.rec.onePMIDonenovelPMIDsameCodon<- NULL
        Idx.rec.onenovelPMIDonenovel<- NULL
        Idx.rec.onePMIDonenovel<- NULL
        Idx.rec.onePMIDsameCodononenovel<- NULL
        Idx.rec.twonovel<- NULL
        Idx.rec.onenoveloneexaclow <- NULL
        Idx.rec.twoexaclow <- NULL
        
        #### assign 21 for variants that are recessive, two truncating
        if(length(Idx.AR.deleterious)>0){
            if( length(Idx.AR.deleterious)>1 | length(intersect(vrnt_df.filt$Zygosity[Idx.AR.deleterious], c("Hom", "Hem")))>0 )
            {
                Idx.rec.twotrunc <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.twotrunc] <<- 21
            }
        }
        
        
        ##### assign 23 for variants that are recessive, with PMID
        if(length(Idx.AR.PMID)>0)
        {
            if( length(Idx.AR.PMID)>1 |  length(intersect(vrnt_df.filt$Zygosity[Idx.AR.PMID], c("Hom", "Hem")))>0 )
            {
                Idx.rec.twoPMID <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.twoPMID] <<- 23
            }
        }
        
        #### assign 22 for variants that are recessive, novel, with PMID
        if(length(Idx.AR.novelPMID)>0)
        {
            if( length(Idx.AR.novelPMID)>1 | length(intersect(vrnt_df.filt$Zygosity[Idx.AR.novelPMID], c("Hom", "Hem")))>0 )
            {
                Idx.rec.twonovelPMID <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.twonovelPMID] <<- 22
            }
        }
        
        
        #### assign 24 for variants that are recessive, two PMID novel related
        if(length(Idx.AR.novelPMIDsameCodon>0))
        {
            if( length(Idx.AR.novelPMIDsameCodon)>1 |  length(intersect(vrnt_df.filt$Zygosity[Idx.AR.novelPMIDsameCodon], c("Hom", "Hem")))>0 )
            {
                Idx.rec.novelPMIDsameCodon <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.novelPMIDsameCodon] <<- 24
            }
        }
        
        #### assign 31 for variants that are recessive, one truncating, one PMID novel
        if ( !identical(Idx.AR.deleterious, integer(0)) & !identical(Idx.AR.novelPMID,integer(0)) )
        {
            if( length(Idx.AR.deleterious)>0 & length(Idx.AR.novelPMID)>0 & !identical(Idx.AR.deleterious,Idx.AR.novelPMID))
            {
                Idx.rec.onetrunconenovelPMID <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.onetrunconenovelPMID] <<- 31
            }
        }
        
        #### assign 32 for variants that are compound het, one truncating, one PMID not novel
        if( !identical(Idx.AR.deleterious, integer(0)) & !identical(Idx.AR.PMID,integer(0)) )
        {
            if( length(Idx.AR.deleterious)>0 & length(Idx.AR.PMID)>0 & !identical(Idx.AR.deleterious, Idx.AR.PMID))
            {
                Idx.rec.onetrunconePMID <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.onetrunconePMID] <<- 32
            }
        }
        
        #### assign 33 for variants that are compound het, one truncating, one novel PMID related
        if( !identical(Idx.AR.deleterious, integer(0)) & !identical(Idx.AR.novelPMIDsameCodon,integer(0)) )
        {
            if( length(Idx.AR.deleterious)>0 & length(Idx.AR.novelPMIDsameCodon)>0 & !identical(Idx.AR.deleterious, Idx.AR.novelPMIDsameCodon))
            {
                Idx.rec.onetrunconenovelPMIDsameCodon <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.onetrunconenovelPMIDsameCodon] <<- 33
            }
        }
        
        #### assign 34 for variants that are compound het, one truncating, one novel
        if( !identical(Idx.AR.deleterious, integer(0)) & !identical(Idx.AR.novelnoPMID,integer(0)) )
        {
            if(length(Idx.AR.deleterious)>0 & length(Idx.AR.novelnoPMID)>0 & !identical(Idx.AR.deleterious, Idx.AR.novelnoPMID))
            {
                Idx.rec.onetrunconenovel <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.onetrunconenovel] <<- 34
            }
        }
        
        #### assign 35 for variants that are recessive, one novel PMID, one PMID not novel
        if( !identical(Idx.AR.novelPMID, integer(0)) & !identical(Idx.AR.PMID,integer(0)) )
        {
            if( length(Idx.AR.novelPMID)>0 & length(Idx.AR.PMID)>0 & !identical(Idx.AR.novelPMID, Idx.AR.PMID))
            {
                Idx.rec.onenovelPMIDonePMID <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.onenovelPMIDonePMID] <<- 35
            }
        }
        
        #### assign 36 for variants that are recessive, one novel PMID, one novel PMID related
        if( !identical(Idx.AR.novelPMID, integer(0)) & !identical(Idx.AR.novelPMIDsameCodon,integer(0)) )
        {
            if( length(Idx.AR.novelPMID)>0 & length(Idx.AR.novelPMIDsameCodon)>0 & !identical(Idx.AR.novelPMID, Idx.AR.novelPMIDsameCodon)  )
            {
                Idx.rec.onenovelPMIDonenovelPMIDsameCodon <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.onenovelPMIDonenovelPMIDsameCodon] <<- 36
            }
        }
        
        
        #### assign 37 for variants that are recessive, one PMID, one novel PMID related
        if( !identical(Idx.AR.PMID, integer(0)) & !identical(Idx.AR.novelPMIDsameCodon,integer(0)) )
        {
            if( length(Idx.AR.PMID)>0 & length(Idx.AR.novelPMIDsameCodon)>0 & !identical(Idx.AR.PMID, Idx.AR.novelPMIDsameCodon)  )
            {
                Idx.rec.onePMIDonenovelPMIDsameCodon <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.onePMIDonenovelPMIDsameCodon] <<- 37
            }
        }
        
        #### assign 38 for variants that are recessive, one novel PMID , one novel
        if( !identical(Idx.AR.novelPMID, integer(0)) & !identical(Idx.AR.novelnoPMID,integer(0)) )
        {
            if( length(Idx.AR.novelPMID)>0 & length(Idx.AR.novelnoPMID)>0 & !identical(Idx.AR.novelPMID,Idx.AR.novelnoPMID)  )
            {
                Idx.rec.onenovelPMIDonenovel <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.onePMIDonenovelPMIDsameCodon] <<- 38
            }
        }
        
        #### assign 39 for variants that are recessive, one PMID , one novel
        if( !identical(Idx.AR.PMID, integer(0)) & !identical(Idx.AR.novelnoPMID,integer(0)) ){
            if( length(Idx.AR.PMID)>0 & length(Idx.AR.novelnoPMID)>0 & !identical(Idx.AR.PMID,Idx.AR.novelnoPMID)  )
            {
                Idx.rec.onePMIDonenovel <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.onePMIDonenovel] <<- 39
            }
        }
        
        
        #### assign 40 for variants that are recessive, one novel PMID related, one novel
        if( !identical(Idx.AR.PMIDsameCodon, integer(0)) & !identical(Idx.AR.novelnoPMID,integer(0)) )
        {
            if( length(Idx.AR.PMIDsameCodon)>0 & length(Idx.AR.novelnoPMID)>0 & !identical(Idx.AR.PMIDsameCodon, Idx.AR.novelnoPMID ) )
            {
                Idx.rec.onePMIDsameCodononenovel <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.onePMIDsameCodononenovel] <<- 40
            }
        }
        
        #### assign 41 for variants that are recessive, two novel
        if( !identical(Idx.AR.novelnoPMID, integer(0)) ){
            if( length(Idx.AR.novelnoPMID)>1 | length(intersect(vrnt_df.filt$Zygosity[Idx.AR.novelnoPMID], c("Hom", "Hem"))) > 0 )
            {
                Idx.rec.twonovel <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.twonovel] <<- 41
            }
        }
        
        #### assign 42 for variants that are recessive, one novel one exac<6
        if( !identical(Idx.AR.exaclow, integer(0)) & !identical(Idx.AR.novelnoPMID,integer(0)) )
        {
            if( length(Idx.AR.exaclow)>0 & length(Idx.AR.novelnoPMID)>0 & !identical(Idx.AR.exaclow, Idx.AR.novelnoPMID ) )
            {
                Idx.rec.onenoveloneexaclow <- Idx.AR
                vrnt_df.filt$Inclusion.criteria[Idx.rec.onenoveloneexaclow] <<- 42
            }
        }
        
        
        
        
        c(Idx.rec.twotrunc, Idx.rec.twonovelPMID, Idx.rec.twoPMID, Idx.rec.novelPMIDsameCodon, Idx.rec.onetrunconenovelPMID,
        Idx.rec.onetrunconePMID, Idx.rec.onetrunconenovelPMIDsameCodon, Idx.rec.onetrunconenovel, Idx.rec.onenovelPMIDonePMID, Idx.rec.onenovelPMIDonenovelPMIDsameCodon,
        Idx.rec.onePMIDonenovelPMIDsameCodon, Idx.rec.onenovelPMIDonenovel,  Idx.rec.onePMIDonenovel, Idx.rec.onePMIDsameCodononenovel,
        Idx.rec.twonovel, Idx.rec.onenoveloneexaclow)
        
        
    }))
}


Prio.Idx <- unique(c(Idx.trunc.pLI,
Idx.Previous.P.LP, #Idx.dom.trunc.pLIsmallHGMDlarge,
Idx.domrec.nontrunc.PMID.novel,
Idx.domrec.nontrunc.noPMID.novel.PMIDsameCodon,
Idx.nontrunc.IDwithparents.novel.highpLI,
Idx.nontrunc.IDwithparents.lowpLI.highZ,
Idx.rec.mix,
Idx.SemanticSim.high))

Prio.Idx <- intersect(Prio.Idx, Idx.SemanticSim.pass)

Variants.individual.keep<- vrnt_df.filt[Prio.Idx,]

Variants.individual.keep
