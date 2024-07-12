# Levene tests for equality of variances

#' Levene test for equality of variances.
#'
#' @description
#' Levene test checks for equality of variances of a set of samples (H0: all
#' variances equal, H1: at least one variance not equal as compared with the
#' rest). The family of `f_levene_test` functions provide faster and more
#' versatile alternatives to \code{\link[car]{leveneTest}}.
#'
#' @details
#' In a numeric matrix or data frame is provided as an argument, a series of
#' Levene tests with the columns is performed. The argument `type` allow for
#' choice of the group's centrality statistic: mean (`type = 'standard'`) or
#' median (`type = 'bf'`). This later option results in a Brown - Forsythe test
#' and may be a better option for strongly skewed distributions.
#'
#' @return a numeric vector, matrix, or a data frame with the following
#' information: number of complete observations, number of groups, the test
#' statistic F, number of degrees of freedom (`df1`: nominator, `df2`
#' denominator) and p value.
#'
#' @param x a vector, list of numeric vectors, matrix or a data frame.
#' Vectors, matrices, and data frames have to contain numeric values.
#' @param f a factor or integer vector used for splitting.
#' @param type type of the test: standard which employs group mean as a
#' centrality measure, or Brown-Forsythe, which operates with group medians.
#' @param as_data_frame should the output be formatted as a data frame? This may
#' render the computation slower.
#' @param safely if `TRUE`, most execution errors will be turned into warnings
#' and `NA` values will be returned as the testing results. This option is
#' particularly useful for large analysis pipelines with possible data quality
#' issues.
#' @param ... extra arguments passed to methods.
#'
#' @export

  f_levene_test <- function(x, ...) UseMethod('f_levene_test')

#' @rdname f_levene_test
#' @export

  f_levene_test.list <- function(x,
                                 type = c('standard', 'bf'),
                                 as_data_frame = FALSE, ...) {

    ## entry control

    stopifnot(is.list(x))

    type = match.arg(type[1], c('standard', 'bf'))

    stopifnot(is.logical(as_data_frame))

    ## computation

    result <- leveneBase(x, type)

    result <- result[-1]

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    as.data.frame(result)

  }

#' @rdname f_levene_test
#' @export

  f_levene_test.default <- function(x, f,
                                    type = c('standard', 'bf'),
                                    as_data_frame = FALSE, ...) {

    ## entry control -------

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric vector.", call. = FALSE)

    }

    if(!is.factor(f) & !is.integer(f)) {

      stop("'f' has to be factor or an integer vector,", call. = FALSE)

    }

    if(length(x) != length(f)) {

      stop("Incompatible lengths of 'x' and 'f'", call. = FALSE)

    }

    if(is.factor(f)) f <- as.integer(f);

    type <- match.arg(type[1], c('standard', 'bf'))

    stopifnot(is.logical(as_data_frame))

    ## testing --------

    result <- leveneVec(x, f, type, TRUE)

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    as.data.frame(result)

  }

#' @rdname f_levene_test
#' @export

  f_levene_test.matrix <- function(x, f,
                                   type = c('standard', 'bf'),
                                   as_data_frame = FALSE,
                                   safely = FALSE, ...) {

    ## entry control --------

    stopifnot(is.matrix(x))

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric matrix or data frame.", call. = FALSE)

    }

    if(!is.factor(f) & !is.integer(f)) {

      stop("'f' has to be factor or an integer vector,", call. = FALSE)

    }

    if(length(f) != nrow(x)) {

      stop("Incompatible sizes of 'x' and 'f'.", call. = FALSE)

    }

    if(is.factor(f)) f <- as.integer(f)

    type <- match.arg(type[1], c('standard', 'bf'))

    stopifnot(is.logical(as_data_frame))
    stopifnot(is.logical(safely))

    ## testing -------

    result <- leveneMtx(x, f, type, !safely)

    if(!as_data_frame) return(result)

    result <- as.data.frame(result)

    rownames_to_column(result, 'variable')

  }

#' @rdname f_levene_test
#' @export

  f_levene_test.data.frame <-  function(x, f,
                                        type = c('standard', 'bf'),
                                        as_data_frame = FALSE,
                                        safely = FALSE, ...) {

    f_levene_test(as.matrix(x),
                  f = f,
                  type = type,
                  as_data_frame = as_data_frame,
                  safely = safely)

  }

# END ------
