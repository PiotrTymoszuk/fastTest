# classical one-way ANOVA

#' One-way analysis of variance (ANOVA).
#'
#' @description
#' A faster and more versatile version of one-way ANOVA genuinely implemented
#' by \code{\link[stats]{aov}}.
#'
#' @details
#' If a numeric matrix or a numeric data frame is provided as the first
#' argument, a series of one-way ANOVAs is performed: one for each column.
#'
#' @return
#' A numeric vector, matrix, or a data frame containing the following
#' information: number of complete observations, number of groups,
#' sum-of-squares: between the groups, within the groups and total, test
#' statistic F, degrees of freedom (nominator: `df1`, denominator `df2`),
#' p value, eta-square effect size statistic, and, optionally, multiple
#' testing-adjusted p value.
#'
#' @param x a vector, list of numeric vectors, matrix or a data frame.
#' Vectors, matrices, and data frames have to contain numeric values.
#' @param f a factor or integer vector used for splitting.
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

  f_one_anova <- function(x, ...) UseMethod('f_one_anova')

#' @rdname f_one_anova
#' @export

  f_one_anova.list <- function(x, as_data_frame = FALSE, ...) {

    ## entry control

    stopifnot(is.list(x))

    stopifnot(is.logical(as_data_frame))

    ## computation

    result <- oneAnovaBase(x)

    result <- result[-1]

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    as.data.frame(result)

  }

#' @rdname f_one_anova
#' @export

  f_one_anova.default <- function(x,
                                  f,
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

    stopifnot(is.logical(as_data_frame))

    ## testing --------

    result <- oneAnovaVec(x, f, TRUE)

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    as.data.frame(result)

  }

#' @rdname f_one_anova
#' @export

  f_one_anova.matrix <- function(x,
                                 f,
                                 as_data_frame = FALSE,
                                 adj_method = 'none',
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

    stopifnot(is.logical(as_data_frame))
    stopifnot(is.logical(safely))

    ## testing -------

    result <- oneAnovaMtx(x = x, f = f, crash = !safely)

    if(adj_method != 'none') {

      result <- cbind(result,
                      p_adjusted = p.adjust(result[, 9],
                                            method = adj_method))

    }

    if(!as_data_frame) return(result)

    result <- as.data.frame(result)

    rownames_to_column(result, 'variable')

  }

#' @rdname f_one_anova
#' @export

  f_one_anova.data.frame <- function(x,
                                     f,
                                     as_data_frame = FALSE,
                                     adj_method = 'none',
                                     safely = FALSE, ...) {

    f_one_anova(as.matrix(x),
                f = f,
                as_data_frame = as_data_frame,
                adj_method = adj_method,
                safely = safely)

  }

# END ------
