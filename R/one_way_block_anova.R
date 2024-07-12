# One-way ANOVA for block designs

#' One-way ANOVA for block designs.
#'
#' @description
#' The `f_one_block_anova` family allows for computation of onw-way analysis
#' of variance (ANOVA) with effects of a single treatment factor `f` and
#' partitioning of error between blocks defined by a single grouping factor `b`.
#'
#' @details
#' As such, the results are equivalent to the traditional R's ANOVA build with
#' the `x ~ f  + Error(b)` formula. If your data has a more sophisticated block
#' design, you are better served with the base R's function
#' \code{\link[stats]{aov}}.
#' If a numeric matrix or a numeric data frame is provided ast the function
#' argument, a series of ANOVAs is computed in a column-wise manner.
#'
#' @return
#' A numeric vector, numeric matrix, or a numeric data frame with the following
#' information:
#'
#' * `n`, `k`, `b`: numbers of complete observations, treatment groups, and
#' blocks, respectively
#'
#' * `ss_treatment`, `ss_block`, `ss_within`, `ss_total`: sum of squares, for,
#' respectively treatment effect, blocks, errors (i.e. within
#' block and treatment groups), and all observations
#'
#' * `f`, `df1`, `df2`, `p`, `etasq`: test statistic F, nominator and
#' denominator degrees of freedom, p-value, and eta-square effect size statistic
#' for the treatment and block effect
#'
#' * `p_adjusted`: optional, multiple testing-adjusted p values for the
#' treatment and block effect. The adjustment is done for matrix and data frame
#' inputs in a column-wise manner
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

  f_one_block_anova <- function(x, ...) UseMethod('f_one_block_anova')

#' @rdname f_one_block_anova
#' @export

  f_one_block_anova.default <- function(x,
                                        f,
                                        b,
                                        as_data_frame = FALSE, ...) {
    ## entry control -------

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

    ## testing --------

    result <- oneBlockAnovaVec(x, f, b)

    if(result[1] == 1.0) stop('Computation failed.', call. = FALSE)

    result <- result[-1]

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    as.data.frame(result)

  }

#' @rdname f_one_block_anova
#' @export

  f_one_block_anova.matrix <- function(x,
                                       f,
                                       b,
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

    result <- oneBlockAnovaMtx(x = x,
                               f = f,
                               b = b,
                               crash = !safely)

    if(adj_method != 'none') {

      result <- cbind(result,
                      p_adjusted_treatment = p.adjust(result[, 11],
                                                      adj_method),
                      p_adjusted_block = p.adjust(result[, 16],
                                                  adj_method))

    }

    if(!as_data_frame) return(result)

    result <- as.data.frame(result)

    rownames_to_column(result, 'variable')


  }

#' @rdname f_one_block_anova
#' @export

  f_one_block_anova.data.frame <- function(x,
                                           f,
                                           b,
                                           as_data_frame = FALSE,
                                           adj_method = 'none',
                                           safely = FALSE, ...) {

    f_one_block_anova(as.matrix(x),
                      f = f,
                      b = b,
                      as_data_frame = as_data_frame,
                      adj_method = adj_method,
                      safely = safely)

  }

# END --------
