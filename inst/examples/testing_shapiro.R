# Testing of functionality of Shapiro-wilk tests

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
                                              c('benign', 'malignant', 'not-determined')))) %>%
    as_tibble

  test_tbl <- num_data[paste0('V', 1:9)]

# Testing the functionality ------

  f_shapiro_test(num_data$V1)

  f_shapiro_test(num_data$V1,
                 num_data$class,
                 as_data_frame = FALSE)

  split(num_data$V1,
        num_data$class) %>%
    f_shapiro_test

  f_shapiro_test(test_tbl,
                 as_data_frame = TRUE)

  f_shapiro_test(test_tbl,
                 num_data$class,
                 as_data_frame = TRUE)

# END -----
