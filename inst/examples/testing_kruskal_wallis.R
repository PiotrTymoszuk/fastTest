# Testing for Kruskal-Wallis tests

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

  ## introducing a new variable

  num_data <- num_data %>%
    arrange(V1, V2, V3) %>%
    mutate(aggressivity = c(rep('low', 150),
                            rep('intermediate', 200),
                            rep('high', 350)),
           aggressivity = factor(aggressivity,
                                 c('low', 'intermediate', 'high')))

  test_tbl <- num_data[paste0('V', 1:9)]

# Functionality and benchmarking --------

  kruskal.test(num_data$V6, num_data$aggressivity) %>%
    unclass

  f_kruskal_test(num_data$V6, num_data$aggressivity, as_data_frame = TRUE)

  kruskal.test(num_data$V6, num_data$class) %>%
    unclass

  split(num_data$V6, num_data$class) %>%
    f_kruskal_test

  f_kruskal_test(test_tbl,
                 num_data$aggressivity,
                 as_data_frame = TRUE,
                 adj_method = 'bonferroni')

  test_tbl %>%
    map(kruskal.test,
        g = num_data$aggressivity) %>%
    map(~.x[c('statistic', 'parameter', 'p.value')]) %>%
    map_dfr(unlist)

  microbenchmark(f_kruskal_test(test_tbl,
                                num_data$aggressivity,
                                as_data_frame = TRUE,
                                adj_method = 'bonferroni'),

                 test_tbl %>%
                   map(kruskal.test,
                       g = num_data$aggressivity))

# END -----
