# Tests for functions computing Cohen's kappa

# tools --------

  library(tidyverse)
  library(fastTest)
  library(microbenchmark)

# test data -------

  #set.seed(1232)

  tst_x <- factor(sample(1:4, size = 100, replace = TRUE))
  tst_y <- factor(sample(1:4, size = 100, replace = TRUE))

# VCD's kappa ---------

  ctg_tbl <- table(tst_x, tst_y)

  vcd::Kappa(ctg_tbl)
  vcd::Kappa(ctg_tbl, weights = "Fleiss-Cohen")

  fastTest:::kappaCpp(tst_x, tst_y, "unweighted")
  fastTest:::kappaCpp(tst_x, tst_y, "equal")
  fastTest:::kappaCpp(tst_x, tst_y, "fleiss")




