# General purpose utilities

# Subsetting with a factor --------

#' Sub-setting of a vector or a matrix by a splitting factor.
#'
#' @description
#' Alternatives to base R `split()` for numeric vectors, matrices,
#' and data frames.
#'
#' @return a list of vectors, matrices, or data frames.
#'
#' @param x a numeric vector, matrix, or data frame.
#' @param f a factor or integer vector used for splitting.
#' @param ... extra arguments passed to methods.
#'
#' @export

  f_split <- function(x, f, ...) UseMethod('f_split')

#' @rdname f_split
#' @export

  f_split.default <- function(x, f, ...) {

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric vector.", call. = FALSE)

    }

    if(!is.integer(f) & !is.factor(f)) {

      stop("'f' has to be a factor or integer vector.", call. = FALSE)

    }

    #comp_indexes <- complete.cases(x, f)

    #x <- x[comp_indexes]
    #f <- f[comp_indexes]

    if(is.factor(f)) {

      f <- droplevels(f)

      new_names <- levels(f)
      f <- as.integer(f)

    } else {

      new_names <- sort(unique(f))

    }

    set_names(compact(Split(x, f)), new_names)

  }

#' @rdname f_split
#' @export

  f_split.matrix <- function(x, f, ...) {

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric matrix or data frame.", call. = FALSE)

    }

    if(!is.integer(f) & !is.factor(f)) {

      stop("'f' has to be a factor or integer vector.", call. = FALSE)

    }

    #comp_indexes <- complete.cases(x, f)

    #x <- x[comp_indexes, ]
    #f <- f[comp_indexes]

    if(is.factor(f)) {

      f <- droplevels(f)

      new_names <- levels(f)
      f <- as.integer(f)

    } else {

      new_names <- sort(unique(f))

    }

    set_names(compact(SplitMtx(x, f)), new_names)

  }

#' @rdname f_split
#' @export

  f_split.data.frame <- function(x, f, ...) {

    f_split.matrix(as.matrix(x), f)

  }

# Contingency table --------

#' Contingency tables.
#'
#' @description
#' Alternative to base R's `table()`.
#'
#' @return a numeric matrix with counts.
#'
#' @param x a vector, matrix, or data frame: numeric or factor,
#' @param f a factor or integer vector used for splitting, if `NULL`, single
#' dimension table is returned.
#' @param as_table should an instance of `table` class be returned?
#' @param ... extra arguments passed to methods.
#'
#' @export

  f_table <- function(x, f, ...) UseMethod('f_table')

#' @rdname f_table
#' @export

  f_table.default <- function(x, f = NULL, as_table = TRUE, ...) {

    stopifnot(is.logical(as_table))

    if(!is.numeric(x) & is.factor(x)) {

      stop("'x' has to be a numeric data frame.", call. = FALSE)

    }

    if(is.null(f)) return(table(x))

    if(!is.integer(f) & !is.factor(f)) {

      stop("'f' has to be a factor or integer vector.", call. = FALSE)

    }

    comp_indexes <- complete.cases(x, f)

    x <- x[comp_indexes]
    f <- f[comp_indexes]

    if(is.factor(f)) {

      f <- droplevels(f)

      new_names <- levels(f)
      f <- as.integer(f)

    } else {

      new_names <- sort(unique(f))

    }

    if(is.factor(x)) {

      levs <- levels(x)

      x <- as.numeric(x)

    } else {

      levs <- sort(unique(x))

    }

    res <- xTable(x, f)

    rownames(res) <- levs
    colnames(res) <- new_names

    if(as_table) {

      res <- structure(res, class = 'table')

    }

    res

  }

# Asymptotic p values in Smirnov distribution ---------

#' Asymptotic p values in Smirnov distribution.
#'
#' @description
#' An internal, simplified wrapper function around \code{\link{psmirnov}}.
#'
#' @return a single numeric value: p value for the given quantile öf
#' Smirnov distribution.
#'
#' @param q quantile of the distribution.
#' @param n1 number of complete observations in the first sample.
#' @param n2 number of complete observations in the second sample.
#' @param alternative type of alternative hypothesis
#'
#' @export

  asympt_smirnov <- function(q,
                             n1,
                             n2,
                             alternative = c("two.sided", "less", "greater")) {

    ## the input control is done by the wrapped function

    psmirnov(q = q,
             sizes = c(n1, n2),
             alternative = alternative,
             exact = FALSE,
             lower.tail = FALSE)

  }

# END ------
