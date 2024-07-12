# Testing for Levene test functions

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

# testing -------

  split(num_data$V1, num_data$aggressivity) %>%
    f_one_anova

  f_one_anova(num_data$V1, num_data$aggressivity)
  f_one_anova(num_data$V1, num_data$class)

  f_one_anova(test_tbl,
              num_data$aggressivity,
              as_data_frame = TRUE,
              adj_method = 'bonferroni')

  paste0('V', 1:9, ' ~ aggressivity') %>%
    map(as.formula) %>%
    map(aov, data = num_data) %>%
    map(summary) %>%
    map(~.x[[1]]) %>%
    map(as.data.frame) %>%
    map_dfr(~.x[1, ])


  microbenchmark(f_one_anova(test_tbl,
                             num_data$aggressivity,
                             as_data_frame = TRUE,
                             adj_method = 'bonferroni'),

                 paste0('V', 1:9, ' ~ aggressivity') %>%
                   map(as.formula) %>%
                   map(aov, data = num_data))


# END -----
