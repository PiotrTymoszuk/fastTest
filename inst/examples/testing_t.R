# functionality, benchmarking, and error handling by T tests

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

  ## a really big matrix for thorough benchmarking

  big_data <- as.matrix(test_tbl)

  for(i in 1:1000) {

    big_data <- cbind(big_data, as.matrix(test_tbl))

  }

  colnames(big_data) <- paste0('Var_', 1:ncol(big_data))

  big_data <- as.data.frame(big_data)

# testing the functionality -----

  class_split <- split(test_tbl$V1, num_data$class)
  aggres_split <- split(test_tbl$V6, num_data$aggressivity)

  ## independent T test

  f_t_test(class_split, type = 'welch')

  f_t_test(test_tbl$V1,
           num_data$class,
           type = 'welch')

  t.test(class_split[[2]],
         class_split[[1]],
         alternative = 'two.sided',
         var.equal = FALSE)

  test_mtx <- as.matrix(test_tbl)

  dimnames(test_mtx) <- NULL

  f_t_test(test_mtx, num_data$class)

  f_t_test(test_tbl,
           num_data$class,
           type = 'welch',
           alternative = 'two.sided')

  test_tbl %>%
    map(split, num_data$class) %>%
    map(~t.test(.x[[2]], .x[[1]],
                alternative = 'two.sided')) %>%
    map_dbl(~.x[['p.value']])

  ## paired T test

  f_t_test(aggres_split,
           type = 'paired',
           alternative = 'two.sided')

  f_t_test(test_tbl$V6,
           num_data$aggressivity,
           type = 'paired',
           alternative = 'two.sided')

  t.test(aggres_split[[2]],
         aggres_split[[1]],
         alternative = 'two.sided',
         var.equal = TRUE,
         paired = TRUE)

  f_t_test(test_tbl,
           num_data$aggressivity,
           type = 'paired',
           alternative = 'two.sided',
           safely = FALSE)

  test_tbl %>%
    map(split, num_data$aggressivity) %>%
    map(~t.test(.x[[2]], .x[[1]],
                paired = TRUE,
                alternative = 'two.sided'))


# benchmarking -------

  ## Welch's T test

  microbenchmark(f_t_test(test_tbl,
                          num_data$class,
                          type = 'welch',
                          alternative = 'two.sided'),

                 f_t_test(test_tbl,
                          num_data$class,
                          type = 'welch',
                          alternative = 'two.sided',
                          adj_method = 'BH',
                          as_data_frame = FALSE),

                 f_t_test(test_tbl,
                          num_data$class,
                          type = 'welch',
                          alternative = 'two.sided',
                          adj_method = 'BH',
                          as_data_frame = TRUE),

                 test_tbl %>%
                   map(split, num_data$class) %>%
                   map(~t.test(.x[[2]], .x[[1]],
                               alternative = 'two.sided')))

  ## paired T test

  microbenchmark(f_t_test(test_tbl,
                          num_data$aggressivity,
                          type = 'paired',
                          alternative = 'two.sided',
                          safely = FALSE),

                 f_t_test(test_tbl,
                          num_data$aggressivity,
                          type = 'paired',
                          alternative = 'two.sided',
                          safely = FALSE,
                          adj_method = 'BH',
                          as_data_frame = FALSE),

                 f_t_test(test_tbl,
                          num_data$aggressivity,
                          type = 'paired',
                          alternative = 'two.sided',
                          safely = FALSE,
                          adj_method = 'BH',
                          as_data_frame = TRUE),

                 test_tbl %>%
                   map(split, num_data$aggressivity) %>%
                   map(~t.test(.x[[2]], .x[[1]],
                               paired = TRUE,
                               alternative = 'two.sided')),

                 test_tbl %>%
                   map(split, num_data$aggressivity) %>%
                   map(~t.test(.x[[2]], .x[[1]],
                               paired = TRUE,
                               alternative = 'two.sided')) %>%
                   map(unclass) %>%
                   map_dfr(as.data.frame))

# exceptions --------

 # f_t_test(class_split,
  #         type = 'paired',
   #        alternative = 'two.sided')

# Large matrix handling --------

  f_t_test(big_data,
           num_data$class,
           type = 'welch',
           alternative = 'two.sided')

  big_data %>%
    map(split, num_data$class) %>%
    map(~t.test(.x[[2]], .x[[1]],
                alternative = 'two.sided'))

  microbenchmark(f_t_test(big_data,
                          num_data$class,
                          type = 'welch',
                          alternative = 'two.sided'),

                 big_data %>%
                   map(split, num_data$class) %>%
                   map(~t.test(.x[[2]], .x[[1]],
                               alternative = 'two.sided')),
                 times = 10)

# END ------
