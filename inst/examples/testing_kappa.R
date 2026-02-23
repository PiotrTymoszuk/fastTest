# Tests for functions computing Cohen's kappa

# tools --------

  library(tidyverse)
  library(fastTest)
  library(microbenchmark)

# test data -------

  set.seed(1232)

  tst_x <- factor(c(NA,
                    sample(1:4, size = 100, replace = TRUE)),
                  levels = 1:5)
  tst_y <- factor(c(sample(1:4, size = 100, replace = TRUE), 5),
                  levels = 1:5)

  mtx_x <- matrix(sample(1:4, size = 100, replace = TRUE), ncol = 10)
  colnames(mtx_x) <- paste0("variable_", 1:10)

  mtx_y <- matrix(c(NA, sample(1:4, size = 98, replace = TRUE), NA), ncol = 10)
  colnames(mtx_y) <- paste0("variable_", 1:10)


# VCD's kappa and fastTest kappa ---------

  ctg_tbl <- table(tst_x, tst_y)

  vcd::Kappa(ctg_tbl)
  vcd::Kappa(ctg_tbl, weights = "Fleiss-Cohen")

  f_kappa(tst_x, tst_y, "unweighted", as_data_frame = TRUE)

  f_kappa(tst_x, tst_y, "unweighted")
  f_kappa(tst_x, tst_y, "equal")
  f_kappa(tst_x, tst_y, "fleiss")

  f_kappa(mtx_x[, 1], mtx_y[, 1], "unweighted")

  f_kappa(mtx_x,
          mtx_y,
          "unweighted")

  f_kappa(mtx_x,
          mtx_y,
          "unweighted",
          as_data_frame = TRUE)

  f_kappa(as.data.frame(mtx_x),
          as.data.frame(mtx_y),
          "unweighted",
          as_data_frame = TRUE)

  f_kappa(map_dfc(as.data.frame(mtx_x), factor),
          map_dfc(as.data.frame(mtx_y), factor),
          "unweighted",
          as_data_frame = TRUE)

# Permutation tests ---------

  f_kappa_test(tst_x, tst_y)
  f_kappa_test(tst_x, tst_y,
               type = "bootstrap")

  f_kappa_test(mtx_x, mtx_y, alternative = "less")
  f_kappa_test(mtx_x, mtx_y, type = "bootstrap", as_data_frame = TRUE)


  f_kappa_test(map_dfc(as.data.frame(mtx_x), factor),
               map_dfc(as.data.frame(mtx_y), factor),
               method = "fleiss",
               as_data_frame = TRUE)

  f_kappa_test(map_dfc(as.data.frame(mtx_x), factor),
               map_dfc(as.data.frame(mtx_y), factor),
               type = "bootstrap",
               method = "fleiss",
               as_data_frame = TRUE)

# END ---------
