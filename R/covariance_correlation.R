# Covariance and correlation coefficients

# Covariance coefficients ---------

#' Covariance and correlation coefficients.
#'
#' @description
#' Alternatives to base R \code{\link[stats]{cov}}
#' and \code{\link[stats]{cor}} provided by `f_cov()` and `f_cor()` function
#' family.
#'
#' @details
#' The functions always use pairwise complete observations. Except for
#' Kednall's tau, they are not faster than the usual R functions which are
#' implemented in C.
#' If a matrix or a data frame is provided as the input, the columns are
#' handled as variables. The correlation-calculating functions offer also the
#' novel correlation coefficient Xi proposed by Chatterjee. Note that the later
#' correlation coefficients are generated with random resolution of rank ties
#' and may hence vary from run to run. Furthermore, those correlation
#' coefficients are not symmetric.
#'
#' @references
#' Chatterjee S. A New Coefficient of Correlation. J Am Stat Assoc (2021)
#' 116:2009–2022. doi:10.1080/01621459.2020.1758115
#'
#' @return
#' Methods for vectors and lists return always two-element numeric vector with
#' the number of complete observations and the requested coefficient.
#'
#' @param x a numeric vector, a list of numeric vectors, matrix or a data frame.
#' @param y a numeric vector.
#' @param method type of covariance or correlation to be analyzed:
#' Pearson's r (default), Spearman's rho, Kendall's TauA, Kendall's TauB,
#' Chatterjee's Xi (xiA: no assumption on ties, xiB: tie correction included).
#' @param ... extra arguments passed to methods.
#'
#' @export

  f_cov <- function(x, ...) UseMethod('f_cov')

#' @rdname f_cov
#' @export

  f_cor <- function(x, ...) UseMethod('f_cor')

#' @rdname f_cov
#' @export

  f_cov.default <- function(x,
                            y,
                            method = c('pearson', 'spearman'), ...) {

    stopifnot(is.atomic(x))
    stopifnot(is.atomic(y))

    if(!is.numeric(x) | !is.numeric(y)) {

      stop("'x' and 'y' have to be numeric vectors.", call. = FALSE)

    }

    method = match.arg(method[1], c('pearson', 'spearman'))

    Cov(x, y, method)

  }

#' @rdname f_cov
#' @export

  f_cov.list <- function(x,
                         method = c('pearson', 'spearman'), ...) {

    if(length(x) < 2) {

      stop("At least two elements are required.", call. = FALSE)

    }

    if(!is.numeric(x[[1]]) | !is.numeric(x[[2]])) {

      stop("A list of numeric vector is required.", call. = FALSE)

    }

    method = match.arg(method[1], c('pearson', 'spearman'))

    Cov(x[[1]], x[[2]], method)

  }

#' @rdname f_cov
#' @export

  f_cov.matrix <- function(x,
                           method = c('pearson', 'spearman'), ...) {

    stopifnot(is.matrix(x))

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric matrix or data frame.",
           call. = FALSE)

    }

    method = match.arg(method[1], c('pearson', 'spearman'))

    CovMtx(x, method)

  }

#' @rdname f_cov
#' @export

  f_cov.data.frame <- function(x,
                               method = c('pearson', 'spearman'), ...) {

    f_cov(as.matrix(x), method)

  }

#' @rdname f_cov
#' @export

  f_cor.default <- function(x, y,
                            method = c('pearson',
                                       'spearman',
                                       'kendallA',
                                       'kendallB',
                                       'xiA',
                                       'xiB'), ...) {

    stopifnot(is.atomic(x))
    stopifnot(is.atomic(y))

    if(!is.numeric(x) | !is.numeric(y)) {

      stop("'x' and 'y' have to be numeric vectors.", call. = FALSE)

    }

    method = match.arg(method[1],
                       c('pearson',
                         'spearman',
                         'kendallA',
                         'kendallB',
                         'xiA',
                         'xiB'))

    Cor(x, y, method)

  }

#' @rdname f_cov
#' @export

  f_cor.list <- function(x, y,
                         method = c('pearson',
                                    'spearman',
                                    'kendallA',
                                    'kendallB',
                                    'xiA',
                                    'xiB'), ...) {

    if(length(x) < 2) {

      stop("At least two elements are required.", call. = FALSE)

    }

    if(!is.numeric(x[[1]]) | !is.numeric(x[[2]])) {

      stop("A list of numeric vector is required.", call. = FALSE)

    }

    method = match.arg(method[1],
                       c('pearson',
                         'spearman',
                         'kendallA',
                         'kendallB',
                         'xiA',
                         'xiB'))

    Cor(x[[1]], x[[2]], method)

  }

#' @rdname f_cov
#' @export

  f_cor.matrix <- function(x,
                           method = c('pearson',
                                      'spearman',
                                      'kendallA',
                                      'kendallB',
                                      'xiA',
                                      'xiB'), ...) {

    stopifnot(is.matrix(x))

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric matrix or data frame.",
           call. = FALSE)

    }

    method = match.arg(method[1],
                       c('pearson',
                         'spearman',
                         'kendallA',
                         'kendallB',
                         'xiA',
                         'xiB'))

    CorMtx(x, method)

  }

#' @rdname f_cov
#' @export

  f_cor.data.frame <- function(x,
                               method = c('pearson',
                                          'spearman',
                                          'kendallA',
                                          'kendallB',
                                          'xiA',
                                          'xiB'), ...) {

    f_cor(as.matrix(x), method)

  }

# END ------
