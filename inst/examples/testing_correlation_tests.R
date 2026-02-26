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

# Testing of functionality -----------

  f_cor_test(num_data$V1,
             num_data$V2,
             type = 'permutation',
             method = 'xiB',
             alternative = 'two.sided',
             n_iter = 1000)

  f_cor_test(num_data$V1,
             num_data$V2,
             type = 'bootstrap',
             method = 'kendallB',
             n_iter = 1000,
             as_data_frame = TRUE)

  f_cor_test(test_tbl,
             type = 'permutation',
             method = 'xiB',
             alternative = 'two.sided',
             n_iter = 100,
             as_data_frame = TRUE)

  f_cor_test(test_tbl,
             type = 'bootstrap',
             method = 'spearman',
             alternative = 'two.sided',
             n_iter = 100,
             as_data_frame = TRUE,
             adj_method = 'bonferroni')


# END -----
