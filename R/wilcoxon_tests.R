# Wilcoxon tests

#' Wilcoxon or Mann-Whitney test.
#'
#' @description
#' A faster and more versatile alternative to R
#' \code{\link[stats]{wilcox.test}} performing Wilcoxon rank tests
#' (also known as Mann-Whitney test) and Wilcoxon signed rank test for paired
#' data.
#'
#' @details
#' In case of matrices or data frames, a series of Wilcoxon tests is performed
#' in a column-wise manner, i.e. each column is treated as a separate variable.
#'
#' @return a vector, matrix, or a data frame with the following information;
#' numbers of samples in the first and the second group, Wilcoxon/Mann-Whitney
#' test statistic U (also known as V for the signed rank test), p value,
#' total number of observation pairs, fractions of pairs in favor of the NULL
#' hypothesis H0, fraction of pairs in favor of the alternative hypothesis H1,
#' estimate of difference between the groups (Hodges-Lehman location parameter
#' or difference of medians, group 2 against group 1, as defined by levels of
#' the splitting factor `f`),
#' biserial rank correlation coefficient r as an effect size metric, and,
#' optionally, multiple testing-adjusted p value.
#'
#' @param x a list of two numeric vectors, vector, matrix or a data frame.
#' Vectors, matrices, and data frames have to contain numeric values.
#' @param f a factor or integer vector used for splitting. If it contains more
#' than two non-empty levels, the first two will be used in the analysis.
#' @param type type of the test: standard Wilcoxon rank test or paired Wilcoxon
#' signed rank test.
#' @param alternative type of the alternative hypothesis concerning difference
#' between the groups (group 2 against group 1).
#' @param exact logical, should p values be derived from Wilcoxon
#' distributions? Otherwise, p values will be computed with a normal
#' distribution approximation. Note that calculation of the exact p value is
#' possible only if there are no ties in ranks of the groups. If there are any
#' ties, normal distribution approximation will be used.
#' @param correct logical, should continuity correction be applied to the test
#' statistic?
#' @param hl_estimate logical, should Hodges-Lehman estimate of location be
#' computed? This option is computationally intensive, as default
#' (`hl_estimate = FALSE`), difference of medians of the analysis groups is
#' returned.
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

  f_wilcox_test <- function(x, ...) UseMethod('f_wilcox_test')

#' @rdname f_wilcox_test
#' @export

  f_wilcox_test.list <- function(x,
                                 type = c('standard', 'paired'),
                                 alternative = c('two.sided', 'less', 'greater'),
                                 exact = TRUE,
                                 correct = TRUE,
                                 hl_estimate = FALSE,
                                 as_data_frame = FALSE, ...) {

    ## input control --------

    stopifnot(is.list(x))
    stopifnot(length(x) >= 2)

    type = match.arg(type[1], c('standard', 'paired'))

    alternative = match.arg(alternative[1],
                            c('two.sided', 'less', 'greater'))

    stopifnot(is.logical(exact))
    stopifnot(is.logical(correct))
    stopifnot(is.logical(as_data_frame))

    ## hypothesis testing ------

    result <- testWilcoxonStd(x[[1]], x[[2]],
                              type = type,
                              alternative = alternative,
                              exact = exact,
                              correct = correct,
                              hl_estimate = hl_estimate,
                              crash = TRUE)

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    as.data.frame(result)

  }

#' @rdname f_wilcox_test
#' @export

  f_wilcox_test.default <- function(x, f,
                                    type = c('standard', 'paired'),
                                    alternative = c('two.sided', 'less', 'greater'),
                                    exact = TRUE,
                                    correct = TRUE,
                                    hl_estimate = FALSE,
                                    as_data_frame = FALSE, ...) {

    ## input control -------

    if(!is.numeric(x)) stop("'x' has to be a numeric vector.", call. = FALSE)

    if(!is.integer(f) & !is.factor(f)) {

      stop("'f' has to be a factor or integer vector.", call. = FALSE)

    }

    if(is.factor(f)) f <- as.integer(f)

    type = match.arg(type[1], c('standard', 'paired'))

    alternative = match.arg(alternative[1],
                            c('two.sided', 'less', 'greater'))

    stopifnot(is.logical(exact))
    stopifnot(is.logical(correct))
    stopifnot(is.logical(as_data_frame))

    ## hypothesis testing --------

    result <- testWilcoxonVec(x, f,
                              type = type,
                              alternative = alternative,
                              exact = exact,
                              correct = correct,
                              hl_estimate = hl_estimate,
                              crash = TRUE)

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    as.data.frame(result)

  }

#' @rdname f_wilcox_test
#' @export

  f_wilcox_test.matrix <- function(x, f,
                                   type = c('standard', 'paired'),
                                   alternative = c('two.sided', 'less', 'greater'),
                                   exact = TRUE,
                                   correct = TRUE,
                                   hl_estimate = FALSE,
                                   adj_method = 'none',
                                   safely = FALSE,
                                   as_data_frame = FALSE, ...) {

    ## input control ------

    stopifnot(is.matrix(x))
    stopifnot(is.numeric(x))

    if(!is.integer(f) & !is.factor(f)) {

      stop("'f' has to be a factor or integer vector.", call. = FALSE)

    }

    if(is.factor(f)) f <- as.integer(f)

    type = match.arg(type[1], c('standard', 'paired'))

    alternative = match.arg(alternative[1],
                            c('two.sided', 'less', 'greater'))

    stopifnot(is.logical(exact))
    stopifnot(is.logical(correct))

    stopifnot(is.character(adj_method))
    stopifnot(is.logical(safely))
    stopifnot(is.logical(as_data_frame))

    ## hypothesis testing -------

    result <- testWilcoxonMtx(x, f,
                              type = type,
                              alternative = alternative,
                              exact = exact,
                              correct = correct,
                              hl_estimate = hl_estimate,
                              crash = !safely)

    p_adjusted <- NULL

    if(adj_method != 'none') {

      result <- cbind(result,
                      p_adjusted = p.adjust(result[, 4],
                                            method = adj_method))

    }

    if(!as_data_frame) return(result)

    result <- as.data.frame(result)

    rownames_to_column(result, 'variable')

  }

#' @rdname f_wilcox_test
#' @export

  f_wilcox_test.data.frame <- function(x, f,
                                       type = c('standard', 'paired'),
                                       alternative = c('two.sided', 'less', 'greater'),
                                       exact = TRUE,
                                       correct = TRUE,
                                       hl_estimate = FALSE,
                                       adj_method = 'none',
                                       safely = FALSE,
                                       as_data_frame = FALSE, ...) {

    f_wilcox_test(as.matrix(x),
                  f,
                  type = type,
                  alternative = alternative,
                  exact = exact,
                  correct = correct,
                  hl_estimate = hl_estimate,
                  adj_method = adj_method,
                  safely = safely,
                  as_data_frame = as_data_frame)

  }

# END ---------
