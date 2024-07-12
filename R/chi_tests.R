# Chi-square test with contingency tables

#' Perason's chi-square test.
#'
#' @description
#' A faster and more versatile alternative to base R's
#' \code{\link[stats]{chisq.test}}.
#'
#' @details
#' In case of matrices or data frames, a series of chi-square tests is performed
#' in a column-wise manner, i.e. each column is treated as a separate variable.
#'
#' @return a vector/matrix or a data frame with numbers of complete observations,
#' chi-square statistic, degrees of freedom, p values, Cramer's V
#' effect size statistics, and, optionally, multiple testing-adjusted p value.
#'
#' @param x a table object, vector, matrix or a data frame. Vectors, matrices,
#' and data frames have to contain numeric values, factors, or data compatible
#' with numeric conversion.
#' @param f a factor or integer vector used for splitting.
#' @param correct should Yates correction be applied?
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

  f_chisq_test <- function(x, ...) UseMethod('f_chisq_test')

#' @rdname f_chisq_test
#' @export

  f_chisq_test.table <- function(x,
                                 correct = TRUE,
                                 as_data_frame = FALSE, ...) {

    stopifnot(is.table(x))
    stopifnot(is.logical(correct))
    stopifnot(is.logical(as_data_frame))

    res <- chiSqTestTbl(x, correct)

    if(!as_data_frame) return(res)

    res <- matrix(res,
                  nrow = 1,
                  dimnames = list(NULL,
                                  names(res)))

    as.data.frame(res)

  }

#' @rdname f_chisq_test
#' @export

  f_chisq_test.default <- function(x,
                                   f,
                                   correct = TRUE,
                                   as_data_frame = FALSE, ...) {

    ## entry control ------

    if(!is.numeric(x) & is.factor(x)) {

      stop("'x' has to be a numeric vector or a factor.", call. = FALSE)

    }

    if(!is.integer(f) & !is.factor(f)) {

      stop("'f' has to be a factor or integer vector.", call. = FALSE)

    }

    if(is.factor(x)) x <- as.numeric(x)
    if(!is.numeric(x)) x <- as.numeric(factor(x))

    if(is.factor(f)) f <- as.integer(f)

    stopifnot(is.logical(correct))
    stopifnot(is.logical(as_data_frame))

    ## hypothesis testing -------

    res <- chiSqTestVec(x, f, correct)

    if(!as_data_frame) return(res)

    res <- matrix(res,
                  nrow = 1,
                  dimnames = list(NULL,
                                  names(res)))

    as.data.frame(res)

  }

#' @rdname f_chisq_test
#' @export

  f_chisq_test.matrix <- function(x,
                                  f,
                                  correct = TRUE,
                                  as_data_frame = FALSE,
                                  adj_method = 'none',
                                  safely = FALSE, ...) {

    ## entry control -------

    stopifnot(is.matrix(x))
    stopifnot(is.numeric(x))

    if(!is.integer(f) & !is.factor(f)) {

      stop("'f' has to be a factor or integer vector.", call. = FALSE)

    }

    if(is.factor(f)) f <- as.integer(f)

    stopifnot(is.logical(correct))
    stopifnot(is.logical(as_data_frame))

    ## hypothesis testing -------

    res <- chiSqTestMtx(x, f, correct, !safely)

    p_adjusted <- NULL

    if(adj_method != 'none') {

      res <- cbind(res,
                   p_adjusted = p.adjust(res[, 4],
                                         method = adj_method))

    }

    if(!as_data_frame) return(res)

    res <- as.data.frame(res)

    rownames_to_column(res, 'variable')

  }

#' @rdname f_chisq_test
#' @export

  f_chisq_test.data.frame <- function(x,
                                      f,
                                      correct = TRUE,
                                      as_data_frame = FALSE,
                                      adj_method = 'none',
                                      safely = FALSE, ...) {

    ## input control ----------

    var_conv <- function(x) {

      if(is.numeric(x)) return(x)

      if(is.factor(x)) return(as.numeric(x))

      return(as.numeric(factor(as.character(x))))

    }

    input_data <- map_dfc(x, var_conv)

    ## calculation ---------

    f_chisq_test(as.matrix(input_data),
                 f,
                 correct = correct,
                 as_data_frame = as_data_frame,
                 adj_method = adj_method,
                 safely = safely)

  }

# END ------
