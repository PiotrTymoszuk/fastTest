# Functionality of block-desing one-way ANOVA

# tools --------

  library(tidyverse)
  library(fastTest)
  library(microbenchmark)

# analysis data --------

  ## Box, Hunter and Hunter, book example, penicillin production

  box_hunter_penicillin <- rbind(c(89, 88, 97, 94),
                                 c(84, 77, 92, 76),
                                 c(81, 87, 87, 85),
                                 c(87, 92, 89, 84),
                                 c(79, 81, 80, 88))

  dimnames(box_hunter_penicillin) <- list(paste0('blend_', 1:5),
                                          LETTERS[1:4])

  box_hunter_penicillin <- box_hunter_penicillin %>%
    as.data.frame %>%
    rownames_to_column('block') %>%
    pivot_longer(cols = all_of(LETTERS[1:4]),
                 names_to = 'treatment',
                 values_to = 'yield') %>%
    mutate(treatment = factor(treatment, LETTERS[1:4]),
           block = factor(block, paste0('blend_', 1:5)))

  ## adding more variables to the table,
  ## for testing of the matrix functionality

  box_hunter_penicillin <- box_hunter_penicillin %>%
    mutate(yield2 = yield + rnorm(20, 5),
           yield3 = yield + runif(20))

# testing the functions --------

  ## single vectors

  aov(yield ~ treatment + Error(block),
      data = box_hunter_penicillin) %>%
    summary

  f_one_block_anova(box_hunter_penicillin$yield,
                    f = box_hunter_penicillin$treatment,
                    b = box_hunter_penicillin$block)

  ## data frames

  c('yield', 'yield2', 'yield3') %>%
    map(paste, '~ treatment + Error(block)') %>%
    map(as.formula) %>%
    map(aov, data = box_hunter_penicillin) %>%
    map(summary)

  f_one_block_anova(x = box_hunter_penicillin[, c("yield", "yield2", "yield3")],
                    f = box_hunter_penicillin$treatment,
                    b = box_hunter_penicillin$block,
                    as_data_frame = TRUE,
                    adj_method = 'bonferroni')

  ## benchmarking

  microbenchmark(c('yield', 'yield2', 'yield3') %>%
                   map(paste, '~ treatment + Error(block)') %>%
                   map(as.formula) %>%
                   map(aov, data = box_hunter_penicillin),

                 f_one_block_anova(x = box_hunter_penicillin[, c("yield", "yield2", "yield3")],
                                   f = box_hunter_penicillin$treatment,
                                   b = box_hunter_penicillin$block,
                                   as_data_frame = TRUE,
                                   adj_method = 'bonferroni'))

# END ------
