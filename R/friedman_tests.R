# Friedman tests with Kendall's W effect size statistic

#' Friedman test.
#'
#' @description
#' Faster and more flexible counterparts of base R's
#' \code{\link[stats]{friedman.test}}.
#'
#' @details
#' The matrix and data frame methods perform a series of Friedman tests: one for
#' each column. Kednall's W serves as an effect size statistic.
#'
#' @return
#' A numeric vector, matrix, or data frame with the following information:
#' total observation number (`n`), number of treatment groups (`k`), number of
#' blocks (`b`), test statistic Q, tiee correction factor, degrees of freedom,
#' p-value, Kendall's W effect size statistic, and, optionally, p value adjusted
#' for multiple testing.
#'
#' @param x a vector, list of numeric vectors, matrix or a data frame.
#' Vectors, matrices, and data frames have to contain numeric values.
#' @param f a factor or integer vector used for splitting by the treatment.
#' @param b a factor of integer vector which defines the blocks.
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

  f_friedman_test <- function(x, ...) UseMethod('f_friedman_test')

#' @rdname f_friedman_test
#' @export

  f_friedman_test.default <- function(x, f, b,
                                      as_data_frame = FALSE, ...) {

    ## input control --------


    if(!is.numeric(x)) {

      stop("'x' has to be a numeric vector.", call. = FALSE)

    }

    if(!is.factor(f) & !is.integer(f)) {

      stop("'f' has to be a factor or integer vector.", call. = FALSE)

    }

    if(!is.factor(b) & !is.integer(b)) {

      stop("'b' has to be a factor or integer vector.", call. = FALSE)

    }

    if((length(x) != length(f)) |
       (length(x) != length(b)) |
       (length(b) != length(f))) {

      stop("Incompatible lengths of 'x', 'f' and 'b'.", call. = FALSE)

    }

    if(is.factor(f)) f <- as.integer(f)
    if(is.factor(b)) b <- as.integer(b)

    stopifnot(is.logical(as_data_frame))

    ## hypothesis testing -------

    result <- friedmanVec(x, f, b, TRUE)

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    as.data.frame(result)

  }

#' @rdname f_friedman_test
#' @export

  f_friedman_test.matrix <- function(x, f, b,
                                     as_data_frame = FALSE,
                                     adj_method = 'none',
                                     safely = FALSE, ...) {

    ## entry control -------

    stopifnot(is.matrix(x))

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric matrix or a numeric data frame.",
           call. = FALSE)

    }

    if(!is.factor(f) & !is.integer(f)) {

      stop("'f' has to be a factor or integer vector.", call. = FALSE)

    }

    if(!is.factor(b) & !is.integer(b)) {

      stop("'b' has to be a factor or integer vector.", call. = FALSE)

    }

    if((nrow(x) != length(f)) |
       (nrow(x) != length(b)) |
       length(f) != length(b)) {

      stop("Incompatible sizes of 'x', 'f' and 'b'.", call. = FALSE)

    }

    if(is.factor(f)) f <- as.integer(f)
    if(is.factor(b)) b <- as.integer(b)

    stopifnot(is.logical(as_data_frame))
    stopifnot(is.character(adj_method))
    stopifnot(is.logical(safely))

    ## testing -------

    result <- friedmanMtx(x = x,
                          f = f,
                          b = b,
                          crash = !safely)

    if(adj_method != 'none') {

      result <- cbind(result,
                      p_adjusted = p.adjust(result[, 7],
                                            adj_method))

    }

    if(!as_data_frame) return(result)

    result <- as.data.frame(result)

    rownames_to_column(result, 'variable')

  }

#' @rdname f_friedman_test
#' @export

  f_friedman_test.data.frame <- function(x, f, b,
                                         as_data_frame = FALSE,
                                         adj_method = 'none',
                                         safely = FALSE, ...) {

    f_friedman_test(as.matrix(x),
                    f = f,
                    b = b,
                    as_data_frame = as_data_frame,
                    adj_method = adj_method,
                    safely = safely)

  }

# END ------
