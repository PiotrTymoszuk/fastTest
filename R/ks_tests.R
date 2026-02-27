# R functions for computing of non-exact (asymptotic p values)
# Kolmogorov-Smirnov tests for a range of various inputs.

#' Two-sample Kolmogorov-Smirnov tests
#'
#' @description
#' Rcpp implementations of the two-sample Kolmogorov-Smirnov test.
#' Matrices and data frames: if `f = NULL`, pairwise Kolmogorov-Smirnov tests
#' for differences in distribution between the columns are performed.
#'
#' @return
#' A numeric vector, a matrix, or a data frame (if `as_data_frame = TRUE`) with
#' numbers of complete observations in the first and second sample, the testing
#' statistic d, and asymptotic (approximate) p value.
#' If `as_data_frame = TRUE`, the output data frame is appended with names of
#' variables.
#' In matrices and data frames, the testing is done column-wise.
#'
#' @inheritParams f_t_test
#'
#' @export

  f_ks_test <- function(x, ...) UseMethod("f_ks_test")

#' @rdname f_ks_test
#' @export

  f_ks_test.default <- function(x,
                                f,
                                alternative = c("two.sided", "greater", "less"),
                                as_data_frame = FALSE,
                                ...) {

    ## entry control ---------

    if(!is.numeric(x)) stop("'x' has to be a numeric vector.", call. = FALSE)

    if(!is.integer(f) & !is.factor(f)) {

      stop("'f' has to be a factor or integer vector.", call. = FALSE)

    }

    if(is.factor(f)) f <- as.integer(f)

    alternative = match.arg(alternative[1],
                            c('two.sided', 'less', 'greater'))

    stopifnot(is.logical(as_data_frame))

    ## testing ------

    result = ksTestVec(x, f, alternative)

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    result <- as.data.frame(result)

    result[["ties"]] <- result[["ties"]] == 1

    result

  }

#' @rdname f_ks_test
#' @export

  f_ks_test.matrix <- function(x,
                               f = NULL,
                               alternative = c("two.sided", "greater", "less"),
                               as_data_frame = FALSE,
                               adj_method = "none",
                               ...) {

    ## input control ------

    stopifnot(is.matrix(x))
    stopifnot(is.numeric(x))

    x_colnames <- colnames(x)

    if(!is.null(f)) {

      if(!is.integer(f) & !is.factor(f)) {

        stop("'f' has to be a factor or integer vector.", call. = FALSE)

      }

      if(is.factor(f)) f <- as.integer(f)

    }

    alternative = match.arg(alternative[1],
                            c('two.sided', 'less', 'greater'))

    stopifnot(is.character(adj_method))
    stopifnot(is.logical(as_data_frame))

    ## hypothesis testing -------

    if(is.null(f)) {

      result <- ksTestPairMtx(x, alternative)

    } else {

      result <- ksTestMtx(x, f, alternative)

    }

    if(adj_method != 'none') {

      p_adjusted <- NULL

      result <- cbind(result,
                      p_adjusted = p.adjust(result[, 5],
                                            method = adj_method))

    }

    if(!as_data_frame) return(result)

    result <- as.data.frame(result)

    if(is.null(f)) {

      if(!is.null(x_colnames)) {

        result[["variable1"]] <- x_colnames[result[["variable1"]]]
        result[["variable2"]] <- x_colnames[result[["variable2"]]]

      }

    } else {

      result <- rownames_to_column(result, 'variable')

    }

    result[["ties"]] <- result[["ties"]] == 1

    result

  }

#' @rdname f_ks_test
#' @export

  f_ks_test.data.frame <- function(x,
                                   f = NULL,
                                   alternative = c("two.sided", "greater", "less"),
                                   as_data_frame = FALSE,
                                   adj_method = "none",
                                   ...) {

    f_ks_test(as.matrix(x),
             f,
             alternative = alternative,
             as_data_frame = as_data_frame,
             adj_method = adj_method)

  }

# END -------
