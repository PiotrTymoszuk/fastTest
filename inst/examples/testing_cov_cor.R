# functionality, benchmarking, and error handling by covariance and correlation
# tests

# tools --------

  library(tidyverse)
  library(fastTest)
  library(microbenchmark)

# analysis data -------

  ## MASS's biopsy data set
  ## with some NAs introduced in the first variable, the splitting factor,
  ## introducing also some unused levels as well

  num_data <- MASS::biopsy %>%
    mutate(class = factor(class, c('benign', 'malignant')))

  num_data$V1[100] <- NA
  num_data$class[500] <- NA

  num_data <- rbind(num_data,
                    data.frame(ID = 'dummy',
                               V1 = 1,
                               V2 = 2,
                               V3 = 3,
                               V4 = 4,
                               V5 = 5,
                               V6 = 6,
                               V7 = 7,
                               V8 = 8,
                               V9 = 9,
                               class = factor('benign',
                                              c('benign', 'malignant')))) %>%
    as_tibble

  ## introducing a variable used for paired tests

  num_data <- num_data %>%
    arrange(V1, V2, V3) %>%
    mutate(aggressivity = c(rep('low', 350), rep('high', 350)),
           aggressivity = factor(aggressivity, c('low', 'high')))

  test_tbl <- num_data[paste0('V', 1:9)]

# Testing of covariance --------

  ## there are some differences to be explained likely by numeric precision
  ## of R and Cpp

  cov(num_data$V1, num_data$V6, use = 'complete.obs', method = 'pearson')
  f_cov(num_data$V1, num_data$V6, method = 'pearson')

  cov(num_data$V1, num_data$V6, use = 'complete.obs', method = 'spearman')
  f_cov(num_data$V1, num_data$V6, method = 'spearman')

  cov(test_tbl, use = 'complete.obs', method = 'pearson')
  f_cov(test_tbl, method = 'pearson')

  cov(test_tbl, use = 'complete.obs', method = 'spearman')
  f_cov(test_tbl, method = 'spearman')

  ## our f_cov function is way slower than the R's genuine
  ## one that was written in C

  microbenchmark(cov(test_tbl, use = 'complete.obs', method = 'spearman'),
                 f_cov(test_tbl, method = 'spearman'))

# Testing of the correlations -------

  ## there are some differences to be explained likely by numeric precision
  ## of R and Cpp

  ## vector setting

  cor(num_data$V1, num_data$V6, use = 'complete.obs', method = 'pearson')
  f_cor(num_data$V1, num_data$V6, method = 'pearson')

  cor(num_data$V1, num_data$V6, use = 'complete.obs', method = 'spearman')
  f_cor(num_data$V1, num_data$V6, method = 'spearman')

  cor(num_data$V1, num_data$V6, use = 'complete.obs', method = 'kendall')
  f_cor(num_data$V1, num_data$V6, method = 'kendallB')

  f_cor(num_data$V1, num_data$V6, method = 'xiB')

  ## data frame or matrix

  cor(test_tbl, use = 'complete.obs', method = 'pearson')
  f_cor(test_tbl, method = 'pearson')

  cor(test_tbl, use = 'complete.obs', method = 'spearman')
  f_cor(test_tbl, method = 'spearman')

  cor(test_tbl, use = 'complete.obs', method = 'kendall')
  f_cor(test_tbl, method = 'kendallB')

  f_cor(test_tbl, method = 'xiA')
  f_cor(test_tbl, method = 'xiB')

  ## benchmarking in the matrix setting

  microbenchmark(cor(test_tbl, use = 'complete.obs', method = 'spearman'),
                 f_cor(test_tbl, method = 'spearman'))

  microbenchmark(cor(test_tbl, use = 'complete.obs', method = 'kendall'),
                 f_cor(test_tbl, method = 'kendallB'))

  microbenchmark(f_cor(test_tbl, method = 'xiA'),
                 f_cor(test_tbl, method = 'xiB'))

# END -------
