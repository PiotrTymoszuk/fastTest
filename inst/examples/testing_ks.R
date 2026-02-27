# Development and testing for functions computing two-sample
# Kolmogorov-Swmirnov tests

# packages --------

  library(tidyverse)
  library(fastTest)
  library(microbenchmark)

# testing data --------

  ## a matrix and a splitting factor

  mtx <- matrix(c(NA, rpois(n = 998, lambda = 2), NA), ncol = 10)

  colnames(mtx) <- paste0("Variable", 1:10)

  f <- c(1, NA, sample(c(0, 1), size = 98, replace = TRUE)) %>%
    factor(levels = 0:1)

  ## two matrices

  mtx1 <- matrix(c(NA, rt(n = 998, df = 20), NA), ncol = 10)
  mtx2 <- matrix(rnorm(n = 100), ncol = 10)

# Stats's KS test ----------

  ## a vector or a matrix and a splitting factor

  f_ks_test(mtx[, 1], f, alternative = "two.sided", as_data_frame = TRUE)

  ks.test(split(mtx[, 1], f)[[1]],
          split(mtx[, 1], f)[[2]],
          exact = FALSE)

  f_ks_test(mtx, f, alternative = "two.sided", as_data_frame = TRUE)
  f_ks_test(as.data.frame(mtx), f, alternative = "two.sided", as_data_frame = TRUE)

  ## pairwise KS tests for matrix columns

  f_ks_test(mtx, as_data_frame = TRUE)
  f_ks_test(as.data.frame(mtx), as_data_frame = TRUE)

  colnames(mtx) <- NULL
  f_ks_test(as.data.frame(mtx), as_data_frame = TRUE)

# END --------
