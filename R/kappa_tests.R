# Unweighted and weighted Cohen's kappa for a pair of factors/integers
# or a pair of integer matrices: permutation and bootstrap tests.

#' Cohen's kappa inter-rater reliability statistic: permutation and bootstrap
#' tests
#'
#' @description
#' The function computes permuted and bootstrapped unweighted Cohen's kappa,
#' or Cohen's kappa weighted with equally spaced or Fleiss - Cohen weights.
#'
#' @details
#' If `x` and `y` are matrices or data frames, the kappas are computed for
#' the corresponding columns.
#' Hence, these matrices or data frames need to have equal dimensions.
#' If `x` is a matrix or a data frame and `y = NULL`, the kappas are tested
#' for all pairs of columns in `x`.
#' `NA` values are silently removed.
#'
#'
#' @references
#' Cohen J. A Coefficient of Agreement for Nominal Scales.
#' Educ Psychol Meas (1960) 20:37–46. doi:10.1177/001316446002000104
#'
#' @references
#' Fleiss JL, Cohen J, Everitt BS. Large sample standard errors of kappa and
#' weighted kappa. Psychol Bull (1969) 72:323–327. doi:10.1037/h0028106
#'
#' @return
#' A numeric vector, matrix, or a data frame with the following elements:
#' numbers of pairwise complete observations, the requested correlation
#' coefficients in the entire data set, bootstrapped means and confidence
#' interval bounds, numbers of successful iterations of the algorithm,
#' numbers of re-samples in favor of the null hypothesis H0, numbers of
#' re-samples in favor of the alternative hypothesis H1, p values, and,
#' optionally p values adjusted for multiple testing.
#'
#' @inheritParams f_cor_test
#' @param x a factor or integer vector, an integer matrix, or a data frame
#' with factor or integer columns.
#' @param y a factor or integer vector, an integer matrix, or a data frame
#' with factor or integer columns.
#' @param method weighting method: `"unweighted"` (default) returns unweighted kappa,
#' `"equal"` computes kappas with equally spaced weighting, and `"fleiss"`
#' calculates kappas with Fleiss - Cohen weights.
#' @param ... additional arguments passed to methods.
#'
#' @export

  f_kappa_test <- function(x, ...) UseMethod("f_kappa_test")

#' @rdname f_kappa_test
#' @export

  f_kappa_test.default <- function(x,
                                   y,
                                   type = c("permutation",
                                            "bootstrap"),
                                   method = c("unweighted", "equal", "fleiss"),
                                   alternative = c("two.sided",
                                                   "less",
                                                   "greater"),
                                   ci_type = c("bca", "percentile"),
                                   conf_level = 0.95,
                                   as_data_frame = FALSE,
                                   n_iter = 1000,
                                   ...) {

    ## input control ---------

    stopifnot(is.atomic(x))
    stopifnot(is.atomic(y))

    if(!is.factor(x) & !is.integer(x)) {

      stop("'x' has to be a factor or integer vector.", call. = FALSE)

    }

    if(!is.factor(y) & !is.integer(y)) {

      stop("'y' has to be a factor or integer vector.", call. = FALSE)

    }

    if((is.integer(x) & is.factor(y)) |
       (is.factor(x) & is.integer(y))) {

      stop("Incompatible data types: `x` and `y` must be both integer or factors.",
           call. = FALSE)

    }

    type <- match.arg(type[1], c("permutation", "bootstrap"))

    method <- match.arg(method[1], c("unweighted", "equal", "fleiss"))

    ci_type <- match.arg(ci_type[1], c("bca", "percentile"))

    stopifnot(is.numeric(conf_level))
    conf_level <- conf_level[1]

    stopifnot(is.logical(as_data_frame))

    stopifnot(is.numeric(n_iter))
    n_iter <- as.integer(n_iter[1])

    ## pre-processing --------

    if(is.factor(x)) {

      if(!identical(levels(x), levels(y))) {

        new_levs <- sort(union(levels(x), levels(y)))

        x <- factor(x, new_levs)
        y <- factor(y, new_levs)

      }

    }

    ## computation -------

    if(type == "permutation") {

      res <- permKappaVec(x, y, method, alternative, n_iter)

    } else {

      res <- bootKappaVec(x, y, method, ci_type, conf_level, n_iter)

    }

    res <- matrix(res,
                  nrow = 1,
                  dimnames = list(NULL,
                                  names(res)))

    if(!as_data_frame) return(res)

    return(as.data.frame(res))

  }

#' @rdname f_kappa_test
#' @export

  f_kappa_test.matrix <- function(x,
                                  y = NULL,
                                  type = c("permutation",
                                           "bootstrap"),
                                  method = c("unweighted", "equal", "fleiss"),
                                  alternative = c("two.sided",
                                                  "less",
                                                  "greater"),
                                  ci_type = c("bca", "percentile"),
                                  conf_level = 0.95,
                                  as_data_frame = FALSE,
                                  n_iter = 1000,
                                  adj_method = "none",
                                  ...) {

    ## entry control ----------

    stopifnot(is.matrix(x))
    x_colnames <- colnames(x)
    if(!is.integer(x)) stop("'x' has to be an integer matrix.", call. = FALSE)

    if(!is.null(y)) {

      stopifnot(is.matrix(y))

      if(ncol(x) != ncol(y)) {

        stop("'x' and 'y' have must have the same dimensions.", call. = FALSE)

      }

      if(!is.integer(y)) stop("'y' has to be an integer matrix.", call. = FALSE)

      y_colnames <- colnames(y)

      if(!is.null(x_colnames) & !is.null(x_colnames)) {

        if(!identical(x_colnames, y_colnames)) {

          stop("Column names of 'x' and 'y' differ: incompatible variables.",
               call. = FALSE)

        }

        y <- y[, x_colnames]

      }


    }

    type <- match.arg(type[1], c("permutation", "bootstrap"))

    method <- match.arg(method[1], c("unweighted", "equal", "fleiss"))

    ci_type <- match.arg(ci_type[1], c("bca", "percentile"))

    stopifnot(is.numeric(conf_level))
    conf_level <- conf_level[1]

    stopifnot(is.logical(as_data_frame))

    stopifnot(is.numeric(n_iter))
    n_iter <- as.integer(n_iter[1])

    stopifnot(is.character(adj_method))

    ## computation for a pair of matrices --------

    if(!is.null(y)) {

      if(type == "permutation") {

        res <- permKappa2Mtx(x, y, method, alternative, n_iter)

      } else {

        res <- bootKappa2Mtx(x, y, method, ci_type, conf_level, n_iter)

      }

      if(!is.null(x_colnames) & !is.null(x_colnames)) {

        rownames(res) <- x_colnames

      }

      ## p value adjustment and the output

      if(adj_method != "none") {

        p_adjusted <- NULL

        if(type == "permutation") {

          res <- cbind(res,
                       p_adjusted = p.adjust(res[, 6],
                                             method = adj_method))

        } else {

          res <- cbind(res,
                       p_adjusted = p.adjust(res[, 9],
                                             method = adj_method))

        }

      }

      if(!as_data_frame) return(res)

      return(rownames_to_column(as.data.frame(res), "variable"))

    }

    ## testing for single matrix ----------

    if(type == "permutation") {

      res <- permKappaMtx(x, method, alternative, n_iter)

    } else {

      res <- bootKappaMtx(x, method, ci_type, conf_level, n_iter)

    }

    if(adj_method != "none") {

      p_adjusted <- NULL

      if(type == "permutation") {

        res <- cbind(res,
                     p_adjusted = p.adjust(res[, 8],
                                           method = adj_method))

      } else {

        res <- cbind(res,
                     p_adjusted = p.adjust(res[, 11],
                                           method = adj_method))

      }

    }

    if(!as_data_frame) return(res)

    res <- as.data.frame(res)

    res[["variable1"]] <- x_colnames[res[["variable1"]]]
    res[["variable2"]] <- x_colnames[res[["variable2"]]]

    return(res)

  }

#' @rdname f_kappa_test
#' @export

  f_kappa_test.data.frame <- function(x,
                                      y = NULL,
                                      type = c("permutation",
                                               "bootstrap"),
                                      method = c("unweighted", "equal", "fleiss"),
                                      alternative = c("two.sided",
                                                      "less",
                                                      "greater"),
                                      ci_type = c("bca", "percentile"),
                                      conf_level = 0.95,
                                      as_data_frame = FALSE,
                                      n_iter = 1000,
                                      adj_method = "none",
                                      ...) {

    ## entry control ---------

    stopifnot(is.data.frame(x))
    x_colnames <- colnames(x)

    if(!is.null(y)) {

      stopifnot(is.data.frame(y))

      if(ncol(x) != ncol(y)) {

        stop("'x' and 'y' have must have the same dimensions.", call. = FALSE)

      }

      y_colnames <- colnames(y)

      if(!identical(x_colnames, y_colnames)) {

        stop("Column names of 'x' and 'y' differ: incompatible variables.",
             call. = FALSE)

      }

      y <- y[, x_colnames]

    }

    type <- match.arg(type[1], c("permutation", "bootstrap"))

    method <- match.arg(method[1], c("unweighted", "equal", "fleiss"))

    ci_type <- match.arg(ci_type[1], c("bca", "percentile"))

    stopifnot(is.numeric(conf_level))
    conf_level <- conf_level[1]

    stopifnot(is.logical(as_data_frame))

    stopifnot(is.numeric(n_iter))
    n_iter <- as.integer(n_iter[1])

    stopifnot(is.character(adj_method))

    ## type compatibility checks ----------

    x_int_check <- map_lgl(x, is.integer)
    x_factor_check <- map_lgl(x, is.factor)

    if(!(all(x_int_check) | all(x_factor_check))) {

      stop("All columns in 'x' have to be factors or integers.", call. = FALSE)

    }

    if(all(x_int_check)) x_type <- "integer" else x_type <- "factor"

    if(!is.null(y)) {

      y_int_check <- map_lgl(y, is.integer)
      y_factor_check <- map_lgl(y, is.factor)

      if(!(all(y_int_check) | all(y_factor_check))) {

        stop("All columns in 'y' have to be factors or integers.", call. = FALSE)

      }

      if(all(y_int_check)) y_type <- "integer" else y_type <- "factor"

      if(y_type != x_type) {

        stop(paste("Incompatible variable types: all columns",
                   "in 'x' and 'y' have to be either factors or integers."),
             call. = FALSE)

      }

    }

    ## factor processing --------

    if(x_type == "factor") {

      if(is.null(y)) {

        x_levs <- map(x, levels)

        all_levs <- reduce(x_levs, union)
        cmm_levs <- reduce(x_levs, intersect)

        if(length(setdiff(all_levs, cmm_levs)) > 0) {

          stop("Levels of all factor variables in data frame `x` must be identical.",
               call. = FALSE)

        }

        x <- map_dfc(x, ~as.integer(as.numeric(.x)))

      } else {

        for(i in 1:ncol(x)) {

          x_levs <- levels(x[[i]])
          y_levs <- levels(y[[i]])

          if(identical(x_levs, y_levs)) next

          new_levels <- sort(union(x_levs, y_levs))

          x[[i]] <- factor(x[[i]], new_levels)
          y[[i]] <- factor(y[[i]], new_levels)

        }

        x <- map_dfc(x, ~as.integer(as.numeric(.x)))
        y <- map_dfc(y, ~as.integer(as.numeric(.x)))

      }

    }

    ## calculation of the kappas -------

    f_kappa_test(x = as.matrix(x),
                 y = if(!is.null(y)) as.matrix(y),
                 type = type,
                 method = method,
                 alternative = alternative,
                 ci_type = ci_type,
                 conf_level = conf_level,
                 as_data_frame = as_data_frame,
                 n_iter = n_iter,
                 adj_method = adj_method)

  }

# END --------
