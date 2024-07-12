# Kruskal-Wallis test with eta-square effect size statistic

#' Kruskal-Wallis test.
#'
#' @description
#' A faster and more versatile alternative to the base R's
#' \code{\link[stats]{kruskal.test}}.
#'
#' @details
#' If a numeric matrix or a numeric data frame is provided as the argument,
#' a series of Kruskal-Wallis tests is performed: one for each column.
#' The test statistic \eqn{H} is
#' \href{https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_test}{corrected for ties},
#' and the p-value is derived from chi-square distribution.
#' The eta-square effect size statistic of computed with the following formula:
#' \deqn{\eta^2 = \frac{H - k + 1}{N -k}}
#' where \eqn{k} represents the number
#' of groups, and \eqn{N} is the number of observations.
#'
#' @return
#' A numeric vector, numeric matrix, or a numeric data frame with the following
#' information: number of complete observations, number of groups,
#' tie-corrected test statistic H, tie correction factor, degrees of freedom,
#' p-value, eta-square effect size statistic, and, optionally, p-value
#' adjusted for multiple testing.
#'
#' @param x a numeric vector, a list of numeric vectors, numeric matrix, or a
#' numeric data frame.
#' @param f a factor or an integer vector used for splitting.
#' @param as_data_frame should the output be formatted as a data frame? This may
#' render the computation slower.
#' @param adj_method multiple testing adjustment method, as specified for
#' \code{\link[stats]{p.adjust}}.
#' @param safely if `TRUE`, most execution errors will be turned into warnings
#' and `NA` values will be returned as the testing results. This option is
#' particularly useful for large analysis pipelines with possible data quality
#' issues.
#' @param ... extra arguments passed to methods.
#'
#' @export

  f_kruskal_test <- function(x, ...) UseMethod('f_kruskal_test')

#' @rdname f_kruskal_test
#' @export

  f_kruskal_test.default <- function(x, f, as_data_frame = FALSE, ...) {

    ## entry control ------

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric vector.", call. = FALSE)

    }

    if(!is.factor(f) & !is.integer(f)) {

      stop("'f' has to be a factor or an integer vector.", call. = FALSE)

    }

    if(length(x) != length(f)) {

      stop("Incompatible lengths of 'x' and 'f'.", call. = FALSE)

    }

    if(is.factor(f)) f <- as.integer(f)

    stopifnot(is.logical(as_data_frame))

    ## testing --------

    result = kruskalVec(x, f, TRUE)

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    as.data.frame(result)

  }

#' @rdname f_kruskal_test
#' @export

  f_kruskal_test.list <- function(x, as_data_frame = FALSE, ...) {

    ## entry control --------

    stopifnot(is.list(x))

    lens <- map_dbl(x, length)

    lenCheck <- map_lgl(lens, ~.x < 2)

    if(any(lenCheck)) stop("Not enough observations.", call. = FALSE)

    stopifnot(is.logical(as_data_frame))

    ## construction of the x and f vectors ------

    x_vec <- unname(unlist(x))

    f_vec <- map2(1:length(x), lens, rep)
    f_vec <- as.integer(unname(unlist(f_vec)))

    ## testing --------

    result <- kruskalVec(x_vec, f_vec, TRUE)

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    as.data.frame(result)

  }

#' @rdname f_kruskal_test
#' @export

  f_kruskal_test.matrix <- function(x,
                                    f,
                                    as_data_frame = FALSE,
                                    adj_method = 'none',
                                    safely = FALSE, ...) {

    ## entry control --------

    stopifnot(is.matrix(x))

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric matrix or a numeric data frame.",
           call. = FALSE)

    }

    if(!is.factor(f) & !is.integer(f)) {

      stop("'f' has to be a factor or an integer vector.", call. = FALSE)

    }

    if(nrow(x) != length(f)) {

      stop("Incompatible sizes of 'x' and 'f'.", call. = FALSE)

    }

    if(is.factor(f)) f <- as.integer(f)

    ## testing --------

    result <- kruskalMtx(x, f, !safely)

    if(adj_method != 'none') {

      result <- cbind(result,
                      p_adjusted = p.adjust(result[, 6],
                                            method = adj_method))

    }

    if(!as_data_frame) return(result)

    result <- as.data.frame(result)

    rownames_to_column(result, 'variable')

  }

#' @rdname f_kruskal_test
#' @export

  f_kruskal_test.data.frame <- function(x,
                                        f,
                                        as_data_frame = FALSE,
                                        adj_method = 'none',
                                        safely = FALSE, ...) {

    f_kruskal_test(as.matrix(x),
                   f = f,
                   as_data_frame = as_data_frame,
                   adj_method = adj_method,
                   safely = safely)

  }

# END -------
