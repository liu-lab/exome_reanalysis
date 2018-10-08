
library(shiny)

setwd("~/pengfeiProject_XS/shiny_app/")

ui <- fluidPage(
  h1("HPOmatcher", style="color:darkcyan"),
  p("Prioritizing variants based on the semantic similarity between patient's HPOs and disease gene-associated HPOs",
    style="color:dimgray;font-size: 18px"),
  p("contact:songxiaofei1002@gmail.com;cshaw@baylorgenetics.com;pengfeil@bcm.edu",style="color:dimgray;font-size: 15px"),
  hr(),
  sidebarLayout(
    
    sidebarPanel(
      # Input: Select a file ----
      fileInput("file1", "Choose CSV File",
                accept = c(".csv")),
      tags$hr(),
      p("Please input column names for",style="color:dimgray;font-size: 15px"),
      textInput('ID', 'sample id'),
      textInput('gene_name', 'gene symbol'), 
      textInput('ptHPO', 'patient HPOs'),
      actionButton('go', 'Go'),
      tags$hr(),
      p("Click to download",style="color:dimgray;font-size: 15px"),
      downloadButton("downloadData", "Download")
      
    ),
    
    
    mainPanel(
      tableOutput('table1')  
      
    )
  )
)




source("reanalysis_app_code.R")


server <- function(input, output){
  
  df <- eventReactive(input$go,{
    
    req(input$file1)
    
    ptdata<- read.csv(input$file1$datapath,
                   sep = ",", header = T)
    
    #ptdata <- read.csv("~/pengfeiProject_XS/shiny_app/data/test.input.csv",header = T)
    
    #colnames(ptdata)[colnames(ptdata)=="Pt_HPO"] <- 'annotated_terms'
    #colnames(ptdata)[colnames(ptdata)=="Gene_symbol"] <- 'Gene_name' #input$gene_name
    #colnames(ptdata)[colnames(ptdata)=="ID"] <- 'PtID' #input$ID
    
    colnames(ptdata)[colnames(ptdata)==input$ID] <- 'PtID' #
    colnames(ptdata)[colnames(ptdata)==input$gene_name] <- 'Gene_name' #
    colnames(ptdata)[colnames(ptdata)==input$ptHPO] <- 'annotated_terms'#
    
   variants_with_diseases = ptdata %>%
      left_join(entrez_to_hgnc ,by=c("Gene_name"= "Approved_Symbol")) %>%
      mutate(gene_diseases = map(Entrez_Gene_ID, ~filter(genes_to_diseases, entrez_gene_id == .x))) %>% 
      mutate(disease_terms = mclapply(gene_diseases, hpo_terms_from_diseases, mc.cores = 20)) %>%
      select(PtID, Gene_name, Entrez_Gene_ID, gene_diseases, disease_terms, annotated_terms)
    
    
    ####rank genes for each patient
    
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
            score_max = dplyr::first(annotation_disease_similarity)) %>%
      group_by(PtID, Gene_name, dz_ID_max, score_max ) %>%
      dplyr::summarise(
        dz_ID_all = paste(db_reference, collapse  = ";"),
        scores = paste(annotation_disease_similarity, collapse  = ";")) %>%
      right_join(ptdata, by = c("PtID", "Gene_name"))
    
    colnames(disease_rankings_by_patient_by_gene)[colnames(disease_rankings_by_patient_by_gene)=='PtID'] <- input$ID 
    colnames(disease_rankings_by_patient_by_gene)[colnames(disease_rankings_by_patient_by_gene)=='Gene_name'] <- input$gene_name #
    colnames(disease_rankings_by_patient_by_gene)[colnames(disease_rankings_by_patient_by_gene)=='annotated_terms'] <-input$ptHPO #
    
    return(disease_rankings_by_patient_by_gene)
    
  })
  
  output$table1 <- renderTable({
    df()
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$file1, "_HPOsim.csv", sep = "")
    },
    content = function(file) {
      write.csv(df(), file, row.names = FALSE)
    }
  )
  
  
}



shinyApp(ui, server)
