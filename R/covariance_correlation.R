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
#' If `as_data_frame = TRUE`, the output is a data frame with variable indexes or
#' names (if present), numbers of complete observations used to compute the
#' coefficients of interest, and the covariance/correlation coefficients.
#' If `as_square_matrix = TRUE`, the output is a numeric matrix with pairwise
#' correlation coefficients, like for base R's \code{\link[stats]{cov}} and
#' \code{\link[stats]{cor}}.
#'
#' @param x a numeric vector, a list of numeric vectors, matrix or a data frame.
#' @param y a numeric vector.
#' @param method type of covariance or correlation to be analyzed:
#' Pearson's r (default), Spearman's rho, Kendall's TauA, Kendall's TauB,
#' Chatterjee's Xi (xiA: no assumption on ties, xiB: tie correction included).
#' @param as_data_frame logical: should the output be coerced to a data frame?
#' @param as_square_matrix logical: should the output be a square matrix? Ignored
#' if `as_data_frame = TRUE`.
#' @param ... additional arguments passed to methods.
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
                           method = c('pearson', 'spearman'),
                           as_data_frame = FALSE,
                           as_square_matrix = FALSE, ...) {

    ## input control ---------

    stopifnot(is.matrix(x))

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric matrix or data frame.",
           call. = FALSE)

    }

    x_colnames <- colnames(x)

    method = match.arg(method[1], c('pearson', 'spearman'))

    stopifnot(is.logical(as_data_frame))
    stopifnot(is.logical(as_square_matrix))

    ## covariance ------

    if(!as_data_frame & as_square_matrix) return(CovMtxSquare(x, method))

    res <- CovMtx(x, method)

    if(!as_data_frame) return(res)

    res <- as.data.frame(res)

    if(is.null(x_colnames)) return(res)

    res[["variable1"]] <- x_colnames[res[["variable1"]]]
    res[["variable2"]] <- x_colnames[res[["variable2"]]]

    return(res)

  }

#' @rdname f_cov
#' @export

  f_cov.data.frame <- function(x,
                               method = c('pearson', 'spearman'),
                               as_data_frame = FALSE,
                               as_square_matrix = FALSE, ...) {

    f_cov(as.matrix(x),
          method,
          as_data_frame,
          as_square_matrix)

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
                                      'xiB'),
                           as_data_frame = FALSE,
                           as_square_matrix = FALSE, ...) {

    ## the input control -------

    stopifnot(is.matrix(x))

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric matrix or data frame.",
           call. = FALSE)

    }

    x_colnames <- colnames(x)

    method = match.arg(method[1],
                       c('pearson',
                         'spearman',
                         'kendallA',
                         'kendallB',
                         'xiA',
                         'xiB'))

    stopifnot(is.logical(as_data_frame))
    stopifnot(is.logical(as_square_matrix))

    ## correlation coefficients ------

    if(!as_data_frame & as_square_matrix) return(CorMtxSquare(x, method))

    if(method %in% c("xiA", "xiB")) {

      warning(paste("Xi correlation coefficients are not symmetric.",
                    "To compute the whole correlation matrix, please",
                    "choose `as_data_frame = FALSE` and `as_square_matrix = TRUE`."),
              call. = FALSE)

    }

    res <- CorMtx(x, method)

    if(!as_data_frame) return(res)

    res <- as.data.frame(res)

    if(is.null(x_colnames)) return(res)

    res[["variable1"]] <- x_colnames[res[["variable1"]]]
    res[["variable2"]] <- x_colnames[res[["variable2"]]]

    return(res)

  }

#' @rdname f_cov
#' @export

  f_cor.data.frame <- function(x,
                               method = c('pearson',
                                          'spearman',
                                          'kendallA',
                                          'kendallB',
                                          'xiA',
                                          'xiB'),
                               as_data_frame = FALSE,
                               as_square_matrix = FALSE, ...) {

    f_cor(as.matrix(x),
          method,
          as_data_frame,
          as_square_matrix)

  }

# END ------
