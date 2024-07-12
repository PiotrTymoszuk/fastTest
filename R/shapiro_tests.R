# Shapiro-Wilk tests for various inputs.
# Those are simple R wrappers around the canonical `shapiro.test`, whose
# core is implemented in C - a C++ translation is not likely to be faster!

#' Shapiro-Wilk tests.
#'
#' @description
#' The functions compute Shapiro-Wilk tests for various inputs: single vectors,
#' a list of vectors, numeric matrices and data frames.
#'
#' @details
#' Because the base R's \code{\link[stats]{shapiro.test}} is already a fast C
#' implementation, we decided not to include an extra C++ version: functions
#' of the `f_shapiro_test` are just simple R wrappers designed for various input
#' objects.
#' If a data frame or a matric is provided as the `x` argument, a series of
#' Shapiro-Wilk tests is performed, where columns are treated as single
#' variables.
#'
#' @return
#' If the splitting vector `f = NULL`, the functions return numeric vectors,
#' matrices or data frames. Otherwise, a list of those objects are returned,
#' where each element represents a single non-empty level of `f`. The output
#' contains information on the number of complete cases, test statistic W, and p
#' values.
#'
#' @param x a numeric vector, list of numeric vectors, matrix or a data frame.
#' @param f an optional factor or integer vector used for splitting.
#' @param as_data_frame should the output be formatted as a data frame?
#' @param ... extra arguments passed to methods.
#'
#' @export

  f_shapiro_test <- function(x, ...) UseMethod('f_shapiro_test')

#' @rdname f_shapiro_test
#' @export

  f_shapiro_test.default <- function(x,
                                     f = NULL,
                                     as_data_frame = FALSE, ...) {

    ## entry control --------

    stopifnot(is.atomic(x))

    if(!is.numeric(x)) stop("'x' has to be a numeric vector.", call. = FALSE)

    if(!is.null(f)) {

      if(!is.factor(f) & !is.integer(f)) {

        stop("'f' has to be a factor or an integer vector.", call. = FALSE)

      }

    }

    stopifnot(is.logical(as_data_frame))

    ## no splitting factor --------

    if(is.null(f)) {

      x <- x[!is.na(x)]

      tst <- unclass(shapiro.test(x))

      result <- c('n' = length(x),
                  'w' = unname(tst[[1]]),
                  'p_value' = tst[[2]])

      if(!as_data_frame) return(result)

      result <- matrix(result, nrow = 1)

      return(as.data.frame(result))

    } else {

      if(is.factor(f)) f <- droplevels(f)

      complete_pairs <- complete.cases(x, f)

      x <- x[complete_pairs]
      f <- f[complete_pairs]

      x_splits <- split(x, f)

      result <- map(x_splits,
                    f_shapiro_test.default,
                    f = NULL,
                    as_data_frame = as_data_frame)

      return(map(result,
                 set_names,
                 c('n', 'w', 'p_value')))

    }

  }

#' @rdname f_shapiro_test
#' @export

  f_shapiro_test.list <- function(x, as_data_frame = FALSE, ...) {

    stopifnot(is.list(x))
    stopifnot(is.logical(as_data_frame))

    x <- compact(x)

    err_txt <- "'x' has to be a list of numeric vectors."

    atoms <- map_lgl(x, is.atomic)
    numbers <- map_lgl(x, is.numeric)

    if(any(!atoms) | any(!numbers)) stop(err_txt, call. = FALSE)

    map(x,
        f_shapiro_test.default,
        f = NULL,
        as_data_frame = as_data_frame)

  }

#' @rdname f_shapiro_test
#' @export

  f_shapiro_test.matrix <- function(x,
                                    f = NULL,
                                    as_data_frame = FALSE, ...) {

    ## entry control --------

    stopifnot(is.matrix(x))

    if(!is.numeric(x)) {

      stop("'x' has to be numeric matrix or data frame.", call = FALSE)

    }

    if(!is.null(f)) {

      if(!is.factor(f) & !is.integer(f)) {

        stop("'f' has to be a factor or an integer vector.", call. = FALSE)

      }

    }

    stopifnot(is.logical(as_data_frame))

    ## testing -------

    x_splits <- map(1:ncol(x), ~x[, .x])

    if(is.null(f)) {

      result <- map(x_splits,
                    f_shapiro_test.default,
                    f = NULL,
                    as_data_frame = FALSE)

      result <- do.call('rbind', result)

      rownames(result) <- colnames(x)

      if(!as_data_frame) return(result)

      result <- as.data.frame(result)

      return(rownames_to_column(result, 'variable'))

    } else {

      result <- map(x_splits,
                    f_shapiro_test.default,
                    f = f,
                    as_data_frame = FALSE)

      if(!is.null(colnames))

      result <- map(result, ~do.call('rbind', .x))

      if(!is.null(colnames(x))) names(result) <- colnames(x)

      if(!as_data_frame) return(result)

      result <- map(result, as.data.frame)

      return(map(result,
                 rownames_to_column,
                 'f'))

    }

  }

#' @rdname f_shapiro_test
#' @export

  f_shapiro_test.data.frame <- function(x,
                                        f = NULL,
                                        as_data_frame = FALSE, ...) {

    f_shapiro_test(as.matrix(x),
                   f = f,
                   as_data_frame = as_data_frame)

  }

# END --------
