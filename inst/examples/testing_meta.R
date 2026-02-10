# Testing the very simple meta-analysis tools provided by the package

# tools ------

  library(tidyverse)
  library(metadat)
  library(meta)
  library(fastTest)

  library(microbenchmark)

  reduce <- purrr::reduce

# data: recidivism and mental health, Assink 2016 ------

  meta_data <- dat.assink2016 %>%
    mutate(semi = sqrt(vi)) %>% ## SEM
    as_tibble

  ## and a list of estimates and variances per study,
  ## used for analysis of the within study-meta-estimates

  meta_list <- split(meta_data, meta_data$study)

  estimate_list <- map(meta_list, ~.x$yi)
  sem_list <- map(meta_list, ~.x$semi)

  #names(estimate_list) <- NULL

  ## and lists of data frames with estimates and their errors

  estimate_df <- meta_list %>%
    map(~.x[c('id', 'yi')]) %>%
    map2(., names(.), ~set_names(.x, c('id', paste0('study_', .y)))) %>%
    reduce(full_join, by = 'id') %>%
    column_to_rownames('id')

  sem_df <- meta_list %>%
    map(~.x[c('id', 'semi')]) %>%
    map2(., names(.), ~set_names(.x, c('id', paste0('study_', .y)))) %>%
    reduce(full_join, by = 'id') %>%
    column_to_rownames('id')

# Vectors --------

  testMate_est <- metagen(meta_data$yi,
                          meta_data$semi,
                          method.tau = 'DL')

  fixed_vec_meta <- f_meta(meta_data$yi,
                           meta_data$semi,
                           type = 'fixed',
                           alternative = 'two.sided',
                           as_data_frame = TRUE)

  random_vec_meta <- f_meta(meta_data$yi,
                            meta_data$semi,
                            type = 'random',
                            alternative = 'two.sided',
                            as_data_frame = TRUE)

# Lists and data frames --------

  testMate_list <-
    map2(estimate_list,
         sem_list,
         metagen, method.tau = 'DL') %>%
    map(~.x[c('TE.fixed', 'seTE.fixed', 'lower.fixed', 'upper.fixed', 'zval.fixed',
              'TE.random', 'seTE.random', 'lower.random', 'upper.random', 'zval.random')]) %>%
    map_dfr(as.data.frame)

  ## list input

  fixed_list_meta <- f_meta(estimate_list,
                            sem_list,
                            type = 'fixed',
                            alternative = 'two.sided',
                            as_data_frame = TRUE,
                            safely = TRUE)

  random_list_meta <- f_meta(estimate_list,
                             sem_list,
                             type = 'random',
                             alternative = 'two.sided',
                             as_data_frame = TRUE,
                             safely = TRUE)

  ## data frame input

  fixed_df_meta <- f_meta(estimate_df,
                          sem_df,
                          type = 'fixed',
                          alternative = 'two.sided',
                          as_data_frame = TRUE,
                          safely = TRUE)

  random_df_meta <- f_meta(estimate_df,
                          sem_df,
                          type = 'random',
                          alternative = 'two.sided',
                          as_data_frame = TRUE,
                          safely = TRUE)

# benchmarking of the list and data frame variants --------

  microbenchmark(f_meta(estimate_list,
                        sem_list,
                        type = 'fixed',
                        alternative = 'two.sided',
                        as_data_frame = TRUE,
                        safely = TRUE),

                 f_meta(estimate_list,
                        sem_list,
                        type = 'random',
                        alternative = 'two.sided',
                        as_data_frame = TRUE,
                        safely = TRUE),

                 f_meta(estimate_df,
                        sem_df,
                        type = 'fixed',
                        alternative = 'two.sided',
                        as_data_frame = TRUE,
                        safely = TRUE),

                 f_meta(estimate_df,
                        sem_df,
                        type = 'random',
                        alternative = 'two.sided',
                        as_data_frame = TRUE,
                        safely = TRUE),

                 map2(estimate_list,
                      sem_list,
                      metagen, method.tau = 'DL'))



# END -------
