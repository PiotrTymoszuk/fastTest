# Pooled variance, T tests, and Cohen's d effect size statistics

#' Student's, Welch's, and paired T test.
#'
#' @description
#' A faster and more versatile alternative to the standard R
#' \code{\link[stats]{t.test}}.
#'
#' @details
#' In case of matrices or data frames, a series of T tests is performed
#' in a column-wise manner, i.e. each column is treated as a separate variable.
#'
#' @return a vector/matrix or a data frame with numbers of complete
#' observations in the investigated groups, T statistics, degrees of freedom,
#' p values, means in the first group, means in the second group, estimate aka
#' difference of means (the second groups versus the first one, or, in case of
#' paired tests, mean of differences), confidence intervals of the difference,
#' pooled variance used for effect size estimation, and Cohen's d effect size
#' statistic.
#' Optionally, multiple-testing-adjusted p values are returned in the last
#' column.
#'
#' @param x a list of two numeric vectors, vector, matrix or a data frame.
#' Vectors, matrices, and data frames have to contain numeric values.
#' @param f a factor or integer vector used for splitting. If it contains more
#' than two non-empty levels, the first two will be used in the analysis.
#' @param type type of the test: standard/Student's, Welch's, or paired
#' T test.
#' @param alternative type of the alternative hypothesis concerning difference
#' of means.
#' @param conf_level confidence level used for computation of the confidence
#' intervals.
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

  f_t_test <- function(x, ...) UseMethod('f_t_test')

#' @rdname f_t_test
#' @export

  f_t_test.list <- function(x,
                            type = c('standard', 'welch', 'paired'),
                            alternative = c('two.sided', 'less', 'greater'),
                            conf_level = 0.95,
                            as_data_frame = FALSE, ...) {

    ## input control ------

    stopifnot(is.list(x))
    stopifnot(length(x) >= 2)

    type = match.arg(type[1],
                     c('standard', 'welch', 'paired'))

    alternative = match.arg(alternative[1],
                            c('two.sided', 'less', 'greater'))

    stopifnot(is.numeric(conf_level))

    if(conf_level < 0 | conf_level > 1) {

      stop("'conf_level' must be within the [0, 1] range.", call. = FALSE)

    }

    stopifnot(is.logical(as_data_frame))

    ## hypothesis testing ---------

    result <- tTestStd(x[[1]], x[[2]],
                       type = type,
                       alternative = alternative,
                       conf_level = conf_level,
                       crash = TRUE)

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    as.data.frame(result)

  }

#' @rdname f_t_test
#' @export

  f_t_test.default <- function(x,
                               f,
                               type = c('standard', 'welch', 'paired'),
                               alternative = c('two.sided', 'less', 'greater'),
                               conf_level = 0.95,
                               as_data_frame = FALSE, ...) {

    ## input control --------

    if(!is.numeric(x)) stop("'x' has to be a numeric vector.", call. = FALSE)

    if(!is.integer(f) & !is.factor(f)) {

      stop("'f' has to be a factor or integer vector.", call. = FALSE)

    }

    if(is.factor(f)) f <- as.integer(f)

    type = match.arg(type[1],
                     c('standard', 'welch', 'paired'))

    alternative = match.arg(alternative[1],
                            c('two.sided', 'less', 'greater'))

    stopifnot(is.numeric(conf_level))

    if(conf_level < 0 | conf_level > 1) {

      stop("'conf_level' must be within the [0, 1] range.", call. = FALSE)

    }

    stopifnot(is.logical(as_data_frame))

    ## hypothesis testing --------

    result <- tTestVec(x, f,
                       type = type,
                       alternative = alternative,
                       conf_level = conf_level,
                       crash = TRUE)

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    as.data.frame(result)

  }

#' @rdname f_t_test
#' @export

  f_t_test.matrix <- function(x,
                              f,
                              type = c('standard', 'welch', 'paired'),
                              alternative = c('two.sided', 'less', 'greater'),
                              conf_level = 0.95,
                              as_data_frame = FALSE,
                              adj_method = 'none',
                              safely = FALSE, ...) {

    ## input control ------

    stopifnot(is.matrix(x))
    stopifnot(is.numeric(x))

    if(!is.integer(f) & !is.factor(f)) {

      stop("'f' has to be a factor or integer vector.", call. = FALSE)

    }

    if(is.factor(f)) f <- as.integer(f)

    type = match.arg(type[1],
                     c('standard', 'welch', 'paired'))

    alternative = match.arg(alternative[1],
                            c('two.sided', 'less', 'greater'))

    stopifnot(is.numeric(conf_level))

    if(conf_level < 0 | conf_level > 1) {

      stop("'conf_level' must be within the [0, 1] range.", call. = FALSE)

    }

    stopifnot(is.character(adj_method))
    stopifnot(is.logical(safely))
    stopifnot(is.logical(as_data_frame))

    ## hypothesis testing -------

    result <- tTestMtx(x, f,
                       type = type,
                       alternative = alternative,
                       conf_level = conf_level,
                       crash = !safely)

    p_adjusted <- NULL

    if(adj_method != 'none') {

      result <- cbind(result,
                      p_adjusted = p.adjust(result[, 5],
                                            method = adj_method))

    }

    if(!as_data_frame) return(result)

    result <- as.data.frame(result)

    rownames_to_column(result, 'variable')

  }

#' @rdname f_t_test
#' @export

  f_t_test.data.frame <- function(x,
                                  f,
                                  type = c('standard', 'welch', 'paired'),
                                  alternative = c('two.sided', 'less', 'greater'),
                                  conf_level = 0.95,
                                  as_data_frame = FALSE,
                                  adj_method = 'none',
                                  safely = FALSE, ...) {

    f_t_test(as.matrix(x),
             f,
             type = type,
             alternative = alternative,
             conf_level = conf_level,
             as_data_frame = as_data_frame,
             adj_method = adj_method,
             safely = safely)

  }

# END -------
