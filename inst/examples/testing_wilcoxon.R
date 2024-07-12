# Functionality, error handling, and banchmarking of Wilcoxon tests

# tools --------

  library(tidyverse)
  library(fastTest)
  library(effectsize)
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

# functionality: Mann-Whitney U test --------

  ## vectors

  class_split <- split(num_data$V4, num_data$class)

  f_wilcox_test(class_split,
                type = 'standard',
                alternative = 'two.sided',
                correct = TRUE)

  f_wilcox_test(num_data$V4,
                num_data$class,
                type = 'standard',
                alternative = 'two.sided',
                correct = TRUE)

  wilcox.test(class_split[[2]],
              class_split[[1]],
              paired = FALSE,
              alternative = 'two.sided',
              correct = TRUE,
              conf.int = TRUE) %>%
    unclass

  ## data frames

  f_wilcox_test(test_tbl,
                num_data$class,
                type = 'standard',
                alternative = 'two.sided',
                correct = TRUE)

  test_tbl %>%
    map(split, num_data$class) %>%
    map(~wilcox.test(.x[[1]], .x[[2]],
                     paired = FALSE,
                     alternative = 'two.sided',
                     correct = TRUE,
                     conf.int = TRUE)) %>%
    map(~.x[c('statistic', 'p.value', 'estimate')]) %>%
    map(unlist) %>%
    reduce(rbind)

  ## benchmarking of the data frame variant

  microbenchmark(f_wilcox_test(test_tbl,
                               num_data$class,
                               type = 'standard',
                               alternative = 'two.sided',
                               correct = TRUE),

                 test_tbl %>%
                   map(split, num_data$class) %>%
                   map(~wilcox.test(.x[[1]], .x[[2]],
                                    paired = FALSE,
                                    alternative = 'two.sided',
                                    correct = TRUE)))

# Paired test -------

  ## vectors

  aggressiv_split <- split(num_data$V6, num_data$aggressivity)

  f_wilcox_test(aggressiv_split,
                type = 'paired',
                alternative = 'two.sided',
                correct = TRUE)

  f_wilcox_test(num_data$V6,
                num_data$aggressivity,
                type = 'paired',
                alternative = 'two.sided',
                correct = TRUE)

  wilcox.test(aggressiv_split[[2]],
              aggressiv_split[[1]],
              paired = TRUE,
              alternative = 'two.sided',
              correct = TRUE,
              conf.int = TRUE) %>%
    .[c('p.value', 'estimate')]

  ## data frames

  f_wilcox_test(test_tbl,
                num_data$aggressivity,
                type = 'paired',
                alternative = 'two.sided',
                correct = TRUE,
                safely = TRUE,
                adj_method = 'BH')

  test_tbl %>%
    map(split, num_data$aggressivity) %>%
    map(~wilcox.test(.x[[1]], .x[[2]],
                     paired = TRUE,
                     alternative = 'two.sided',
                     correct = TRUE,
                     conf.int = TRUE)) %>%
    map(~.x[c('statistic', 'p.value', 'estimate')]) %>%
    map(unlist) %>%
    reduce(rbind)

  ## benchmarking the data frame variant

  microbenchmark(f_wilcox_test(test_tbl,
                               num_data$aggressivity,
                               type = 'paired',
                               alternative = 'two.sided',
                               correct = TRUE,
                               safely = TRUE),

                 test_tbl %>%
                   map(split, num_data$aggressivity) %>%
                   map(~wilcox.test(.x[[1]], .x[[2]],
                                    paired = TRUE,
                                    alternative = 'two.sided',
                                    correct = TRUE)))




# END -------
