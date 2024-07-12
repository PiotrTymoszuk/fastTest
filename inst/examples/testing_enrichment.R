# Functionality of simulated enrichment tests in a GO term enrichment analysis

# tools -------

  library(tidyverse)
  library(fastTest)
  library(org.Hs.eg.db)
  library(biggrExtra)
  library(AnnotationDbi)

  library(microbenchmark)

  select <- dplyr::select
  reduce <- purrr::reduce

# analysis data ---------

  ## a list with mapping of the entrez IDs to GO terms

  all_go <- tcga_data$entrez_id %>%
    unique %>%
    mapIds(org.Hs.eg.db,
           keys = .,
           keytype = 'ENTREZID',
           column = 'GO')

  all_go <- tibble(entrez_id = names(all_go),
                   go_id = all_go) %>%
    filter(complete.cases(.))

  all_go <- split(all_go$entrez_id, all_go$go_id)

  ## vectors of IDs for all available genes, up- and downregulated genes

  all_genes <- tcga_data %>%
    filter(complete.cases(.)) %>%
    .$entrez_id %>%
    unique

  up_genes <- tcga_data %>%
    filter(complete.cases(.),
           regulation == 'upregulated') %>%
    .$entrez_id %>%
    unique

  down_genes <-  tcga_data %>%
    filter(complete.cases(.),
           regulation == 'downregulated') %>%
    .$entrez_id %>%
    unique

# Testing ---------

  testFisher <- f_enrichment(up_genes,
                             dict = all_go,
                             all = all_genes,
                             type = 'fisher',
                             as_data_frame = TRUE,
                             adj_method = 'BH')

  testEnrich <- f_enrichment(up_genes,
                             dict = all_go[1:1000],
                             all = all_genes,
                             type = 'random',
                             n_iter = 1000,
                             .parallel = TRUE,
                             .n_chunks = 100)

# END -----
