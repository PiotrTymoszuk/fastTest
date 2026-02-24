# Development and testing for functions computing two-sample
# Kolmogorov-Swmirnov tests

# packages --------

  library(tidyverse)
  library(fastTest)
  library(microbenchmark)

# testing data --------

  ## vectors

  tst_x <- c(rpois(n = 100, lambda = 2), NA)
  tst_y <- c(NA, rnorm(n = 150))

  ## a matrix and a splitting factor

  mtx <- matrix(c(NA, rpois(n = 998, lambda = 2), NA), ncol = 10)

  colnames(mtx) <- paste0("Variable", 1:10)

  f <- c(1, NA, sample(c(0, 1), size = 98, replace = TRUE)) %>%
    factor(levels = 0:1)

  ## two matrices

  mtx1 <- matrix(c(NA, rt(n = 998, df = 20), NA), ncol = 10)
  mtx2 <- matrix(rnorm(n = 100), ncol = 10)

# Stats's KS test ----------

  ks.test(tst_x, tst_y, alternative = "two.sided")
  fastTest:::ksTestCpp(tst_x, tst_y, alternative = "two.sided")

  microbenchmark(ks.test(tst_x, tst_y, alternative = "two.sided", exact = FALSE),
                 fastTest:::ksTestCpp(tst_x, tst_y, alternative = "two.sided"))

  ks.test(split(mtx[, 1], f)[[1]],
          split(mtx[, 1], f)[[2]],
          exact = FALSE)

  fastTest:::ksTestVec(mtx[, 1], f = f)

  fastTest:::ksTestMtx(mtx, f)

  fastTest:::ksTest2Mtx(mtx1, mtx2)

