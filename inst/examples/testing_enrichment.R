# Functionality of simulated enrichment tests in a GO term enrichment analysis

# tools -------

  library(tidyverse)
  library(fastTest)
  library(org.Hs.eg.db)
  library(perich)
  library(trafo)
  library(AnnotationDbi)

  library(microbenchmark)

  select <- dplyr::select
  reduce <- purrr::reduce

# analysis data ---------

  ## a list with mapping of the entrez IDs to GO terms

  data("tcga_mutations")

  ## top mutations in the TCGA bladder cancer data and all
  ## recorded mutations

  mutations <- list()

  mutations$top <-
    sort(colMeans(tcga_mutations[, -2:-1]), decreasing = T)[1:200]

  mutations$all <- tcga_mutations[, -2:-1]

  mutations <- mutations %>%
    map(names)

  go_lexicon <-
    tibble(gene_symbol = mutations$all,
           go_id = mapIds(org.Hs.eg.db,
                          keys = mutations$all,
                          keytype = 'SYMBOL',
                          column = 'GO')) %>%
    blast(go_id) %>%
    map(~.x$gene_symbol)

# Testing ---------

  testFisher <- f_enrichment(mutations$top,
                             dict = go_lexicon,
                             all = mutations$all,
                             type = 'fisher',
                             as_data_frame = TRUE,
                             adj_method = 'BH')

  testEnrich <- f_enrichment(mutations$top,
                             dict = go_lexicon,
                             all = mutations$all,
                             type = 'random',
                             n_iter = 1000,
                             as_data_frame = TRUE,
                             .parallel = TRUE,
                             .n_chunks = 100)

# END -----
