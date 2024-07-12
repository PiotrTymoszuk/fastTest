# Testing of the package functionality: chi-square test

# tools --------

  library(tidyverse)
  library(fastTest)
  library(microbenchmark)

# analysis data -------

  ## MASS's biopsy data set
  ## with some NAs introduced in the first variable
  ## introducing also some unused levels as well

  num_data <- MASS::biopsy %>%
    mutate(class = ifelse(V1 >= 4 & V1 <= 7, 'uncertain', as.character(class)),
           class = factor(class,
                          c('benign', 'uncertain',
                            'malignant', 'not determined'))) %>%
    as_tibble

  num_data$V1[100] <- NA
  num_data$class[500] <- NA

  test_tbl <- num_data[paste0('V', 1:9)]

  ## a really big matrix for thorough benchmarking

  big_data <- as.matrix(test_tbl)

  for(i in 1:1000) {

    big_data <- cbind(big_data, as.matrix(test_tbl))

  }

  colnames(big_data) <- paste0('Var_', 1:ncol(big_data))

# Splitting -------

  ## functionality

  f_split(num_data$V6, droplevels(num_data$class))

  test_mtx <- as.matrix(num_data[2:5])

  dimnames(test_mtx) <- NULL

  f_split(test_mtx, num_data$class)

  rownames(test_mtx) <- paste0('obs_', 1:nrow(test_mtx))

  f_split(test_mtx, num_data$class)

  ## benchmarking

  microbenchmark(fastTest:::Split(num_data$V6, num_data$class),
                 f_split(num_data$V6, num_data$class),
                 split(num_data$V6, num_data$class),
                 times = 1000)

  microbenchmark(fastTest:::SplitMtx(test_mtx, num_data$class),
                 f_split(test_mtx, num_data$class),
                 split(test_mtx, num_data$class),
                 times = 1000)

# Contingency tables ---------

  ## functionality

  f_table(num_data$V1, num_data$class)
  table(num_data$V1, num_data$class)

  ## benchmarking

  microbenchmark(f_table(num_data$V1, num_data$class),
                 table(num_data$V1, num_data$class),
                 times = 1000)

# Chi-square test -------

  ## functionality

  test_cont_tbl <- f_table(num_data$V1, num_data$class)

  f_chisq_test(test_cont_tbl, correct = FALSE)
  chisq.test(test_cont_tbl)

  f_chisq_test(num_data$V1, num_data$class, correct = FALSE)
  chisq.test(num_data$V1, num_data$class)

  test_mtx <- as.matrix(test_tbl)
  dimnames(test_mtx) <- NULL


  f_chisq_test(test_mtx,
               num_data$class,
               correct = FALSE)

  f_chisq_test(as.matrix(test_tbl),
               num_data$class,
               correct = FALSE,
               adj_method = 'BH')

  f_chisq_test(as.matrix(test_tbl),
               num_data$class,
               correct = FALSE,
               adj_method = 'BH',
               as_data_frame = TRUE)

  map(test_tbl, ~chisq.test(.x, num_data$class)) %>%
    map_dbl(~.x[['p.value']])

  f_chisq_test(test_tbl,
               num_data$class,
               correct = FALSE)


  ## benchmarking

  microbenchmark(f_chisq_test(test_cont_tbl,
                              correct = FALSE),
                 f_chisq_test(test_cont_tbl,
                              correct = FALSE,
                              as_data_frame = TRUE),
                 chisq.test(test_cont_tbl))

  microbenchmark(f_chisq_test(num_data$V1,
                              num_data$class,
                              correct = FALSE),
                 f_chisq_test(num_data$V1,
                              num_data$class,
                              correct = FALSE,
                              as_data_frame = TRUE),
                 chisq.test(num_data$V1,
                            num_data$class))

  microbenchmark(f_chisq_test(as.matrix(test_tbl),
                              num_data$class,
                              correct = FALSE),

                 f_chisq_test(as.matrix(test_tbl),
                              num_data$class,
                              correct = FALSE,
                              adj_method = 'BH'),

                 f_chisq_test(as.matrix(test_tbl),
                              num_data$class,
                              correct = FALSE,
                              adj_method = 'BH',
                              as_data_frame = TRUE),

                 f_chisq_test(test_tbl,
                              num_data$class,
                              correct = FALSE,
                              adj_method = 'BH',
                              as_data_frame = TRUE),

                 map(test_tbl,
                     ~chisq.test(.x, num_data$class)))

  microbenchmark(f_chisq_test(as.matrix(test_tbl),
                              num_data$class,
                              correct = FALSE,
                              adj_method = 'BH'),

                 f_chisq_test(test_tbl,
                              num_data$class,
                              correct = FALSE,
                              adj_method = 'BH'))

  ## functionality with a big matrix

  f_chisq_test(big_data,
               num_data$class,
               correct = FALSE)

  f_chisq_test(big_data,
               num_data$class,
               correct = FALSE,
               adj_method = 'BH')

  f_chisq_test(big_data,
               num_data$class,
               correct = FALSE,
               adj_method = 'BH',
               as_data_frame = TRUE)

  map(as.data.frame(big_data),
      ~chisq.test(.x, num_data$class))

  #microbenchmark( f_chisq_test(big_data,
   #                            num_data$class,
    #                           correct = FALSE),
#
 #                 f_chisq_test(big_data,
  #                             num_data$class,
   #                            correct = FALSE,
    #                           adj_method = 'BH'),
#
 #                 f_chisq_test(big_data,
  #                             num_data$class,
   #                            correct = FALSE,
    #                           adj_method = 'BH',
     #                          as_data_frame = TRUE),
#
 #                 map(as.data.frame(big_data),
  #                    ~chisq.test(.x, num_data$class)),
   #               times = 10)

# Error handling ----------

  ## flawed data

  flaw_test_data <- test_tbl %>%
    mutate(V1 = NA, V5 = 1)

  flaw_ctg_table1 <- f_table(flaw_test_data$V1, num_data$class)
  flaw_ctg_table2 <- f_table(flaw_test_data$V5, num_data$class)

  ## do not execute

  # f_chisq_test(flaw_test_data$V1, num_data$class)
  # f_chisq_test(flaw_test_data$V2, num_data$class)
  # f_chisq_test(flaw_test_data$V5, num_data$class)

  ## safely mode

  f_chisq_test(flaw_test_data, num_data$class, safely = TRUE)

# END ------
