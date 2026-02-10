# Functionality of Friedman tests

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

  set.seed(12345)

  box_hunter_penicillin <- box_hunter_penicillin %>%
    mutate(yield2 = round(yield + rnorm(20, 5)),
           yield3 = round(yield + runif(20)))

# testing ---------

  ## vector variant

  f_friedman_test(x = box_hunter_penicillin$yield,
                  f = box_hunter_penicillin$treatment,
                  b = box_hunter_penicillin$block)

  friedman.test(box_hunter_penicillin$yield,
                box_hunter_penicillin$treatment,
                box_hunter_penicillin$block)

  ## data frame/matrix method

  f_friedman_test(box_hunter_penicillin[, c("yield", "yield2", "yield3")],
                  f = box_hunter_penicillin$treatment,
                  b = box_hunter_penicillin$block,
                  as_data_frame = TRUE,
                  adj_method = 'bonferroni')

  box_hunter_penicillin[, c("yield", "yield2", "yield3")] %>%
    map(friedman.test,
        group = box_hunter_penicillin$treatment,
        block = box_hunter_penicillin$block) %>%
    map(~.x[c('statistic', 'parameter', 'p.value')]) %>%
    map_dfr(unlist)

# END ------
