# Utilities for R-side paralellization

# Splitting of the input into chunks by columns ------

#' Split an object into chunks.
#'
#' @description
#' Splits an object (vector, list, matrix or a data frame) into `n` chunks.
#' If the object length or its column number is not divisible by `n`, the
#' output accomplished by `cut` is intended to comprise elements of a similar
#' lengths.
#'
#' @return a list of objects split by elements or columns.
#'
#' @param x an object.
#' @param n an integer defining the chunk number.
#' @param ... extra arguments passed to methods.
#'
#' @export

  chunk <- function(x, ...) UseMethod('chunk')

#' @rdname chunk
#' @export

  chunk.default <- function(x, n, ...) {

    stopifnot(is.atomic(x))
    stopifnot(is.numeric(n))
    stopifnot(n >= 1)

    split_idx <- cut(seq_along(x), n)

    new_x <- split(x, split_idx)

    set_names(new_x, paste0('chunk_', seq_along(new_x)))

  }

#' @rdname chunk
#' @export

  chunk.list <- function(x, n, ...) {

    stopifnot(is.list(x))

    idx <- chunk(seq_along(x), n)

    map(idx, ~x[.x])

  }

#' @rdname chunk
#' @export

  chunk.matrix <- function(x, n, ...) {

    stopifnot(is.matrix(x))

    idx <- chunk(1:ncol(x), n)

    map(idx, ~x[, .x, drop = FALSE])

  }

#' @rdname chunk
#' @export

  chunk.data.frame <- function(x, n, ...) {

    stopifnot(is.data.frame(x))

    idx <- chunk(1:ncol(x), n)

    map(idx, ~x[, .x, drop = FALSE])

  }

# END ------
