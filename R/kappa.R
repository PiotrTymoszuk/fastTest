# Unweighted and weighted Cohen's kappa for a pair of factors/integers
# or a pair of integer matrices.

#' Cohen's kappa inter-rater reliability statistic.
#'
#' @description
#' The function computes unweighted Cohen's kappa, or Cohen's kappa weighted
#' with equally spaced or Fleiss - Cohen weights.
#'
#' @details
#' If `x` and `y` are matrices, the kappas are computed for the corresponding
#' columns.
#' Hence, these matrices or data frames need to have equal dimensions.
#' `NA` values are silently removed.
#'
#' @references
#' Cohen J. A Coefficient of Agreement for Nominal Scales.
#' Educ Psychol Meas (1960) 20:37–46. doi:10.1177/001316446002000104
#'
#' @references
#' Fleiss JL, Cohen J, Everitt BS. Large sample standard errors of kappa and
#' weighted kappa. Psychol Bull (1969) 72:323–327. doi:10.1037/h0028106
#'
#' @return
#' The default `f_kappa()` method for a pair of vectors returns a numeric vector
#' with the number of complete observations and kappa value.
#' The method for a pair of matrices or data frames returns
#' a numeric matrix or a data frame with the numbers of complete observations and
#' kappa values.
#' If only one matrix or data frame is provided, i.e. `y = NULL`, kappa
#' coefficients are computed for all pairs of columns in `x`.
#' If `as_data_frame = TRUE`, the output is coerced to a data frame, and,
#' optionally appended with names of variables in matrices or data frames.
#' Note: because of checks of level compatibility, the function is way faster
#' for integers than factors.
#'
#' @param x a factor or integer vector, an integer matrix, or a data frame
#' with factor or integer columns.
#' @param y `NULL` or a factor or integer vector, an integer matrix, or
#' a data frame with factor or integer columns.
#' @param method weighting method: `"unweighted"` (default) returns unweighted kappa,
#' `"equal"` computes kappas with equally spaced weighting, and `"fleiss"`
#' calculates kappas with Fleiss - Cohen weights.
#' @param as_data_frame logical: should the output be coerced to a data frame?
#' @param ... additional arguments passed to methods.
#'
#' @export

  f_kappa <- function(x, ...) UseMethod("f_kappa")

#' @rdname f_kappa
#' @export

  f_kappa.default <- function(x,
                              y,
                              method = c("unweighted", "equal", "fleiss"),
                              as_data_frame = FALSE,
                              ...) {

    ## input control ---------

    stopifnot(is.atomic(x))
    stopifnot(is.atomic(y))

    if(!is.factor(x) & !is.integer(x)) {

      stop("'x' has to be a factor or integer vector.", call. = FALSE)

    }

    if(!is.factor(y) & !is.integer(y)) {

      stop("'y' has to be a factor or integer vector.", call. = FALSE)

    }

    if((is.integer(x) & is.factor(y)) |
       (is.factor(x) & is.integer(y))) {

      stop("Incompatible data types: `x` and `y` must be both integer or factors.",
           call. = FALSE)

    }

    method <- match.arg(method[1], c("unweighted", "equal", "fleiss"))

    stopifnot(is.logical(as_data_frame))

    ## pre-processing --------

    if(is.factor(x)) {

      if(!identical(levels(x), levels(y))) {

        new_levs <- sort(union(levels(x), levels(y)))

        x <- factor(x, new_levs)
        y <- factor(y, new_levs)

      }

    }

    ## computation -------

    res <- kappaCpp(x, y, method)

    if(!as_data_frame) return(res)

    return(as.data.frame(as.list(res)))

  }

#' @rdname f_kappa
#' @export

  f_kappa.matrix <- function(x,
                             y = NULL,
                             method = c("unweighted", "equal", "fleiss"),
                             as_data_frame = FALSE,
                             ...) {

    ## entry control ----------

    stopifnot(is.matrix(x))
    if(!is.integer(x)) stop("'x' has to be an integer matrix.", call. = FALSE)

    x_colnames <- colnames(x)

    if(!is.null(y)) {

      stopifnot(is.matrix(y))
      if(!is.integer(y)) stop("'y' has to be an integer matrix.", call. = FALSE)

      if(ncol(x) != ncol(y)) {

        stop("'x' and 'y' have must have the same dimensions.", call. = FALSE)

      }

      y_colnames <- colnames(y)

      if(!is.null(x_colnames) & !is.null(x_colnames)) {

        if(!identical(x_colnames, y_colnames)) {

          stop("Column names of 'x' and 'y' differ: incompatible variables.",
               call. = FALSE)

        }

        y <- y[, x_colnames]

      }

    }

    method <- match.arg(method[1], c("unweighted", "equal", "fleiss"))

    stopifnot(is.logical(as_data_frame))

    ## computation for a pair of matrices --------

    if(!is.null(y)) {

      res <- kappa2Mtx(x, y, method)

      if(!is.null(x_colnames) & !is.null(x_colnames)) {

        rownames(res) <- x_colnames

      }

      if(!as_data_frame) return(res)

      return(rownames_to_column(as.data.frame(res), "variable"))

    }

    ## computation of kappas for a single matrix --------

    res <- kappaMtx(x, method)

    if(!as_data_frame) return(res)

    res <- as.data.frame(res)

    if(is.null(x_colnames)) return(res)

    res[["variable1"]] <- x_colnames[res[["variable1"]]]
    res[["variable2"]] <- x_colnames[res[["variable2"]]]

    return(res)

  }

#' @rdname f_kappa
#' @export

  f_kappa.data.frame <- function(x,
                                 y = NULL,
                                 method = c("unweighted", "equal", "fleiss"),
                                 as_data_frame = FALSE, ...) {

    ## entry control ---------

    stopifnot(is.data.frame(x))
    x_colnames <- colnames(x)

    if(!is.null(y)) {

      stopifnot(is.data.frame(y))

      if(ncol(x) != ncol(y)) {

        stop("'x' and 'y' have must have the same dimensions.", call. = FALSE)

      }

      y_colnames <- colnames(y)

      if(!identical(x_colnames, y_colnames)) {

        stop("Column names of 'x' and 'y' differ: incompatible variables.",
             call. = FALSE)

      }

      y <- y[, x_colnames]

    }

    method <- match.arg(method[1], c("unweighted", "equal", "fleiss"))

    stopifnot(is.logical(as_data_frame))

    ## type compatibility checks ----------

    x_int_check <- map_lgl(x, is.integer)
    x_factor_check <- map_lgl(x, is.factor)

    if(!(all(x_int_check) | all(x_factor_check))) {

      stop("All columns in 'x' have to be factors or integers.", call. = FALSE)

    }

    if(all(x_int_check)) x_type <- "integer" else x_type <- "factor"

    if(!is.null(y)) {

      y_int_check <- map_lgl(y, is.integer)
      y_factor_check <- map_lgl(y, is.factor)

      if(!(all(y_int_check) | all(y_factor_check))) {

        stop("All columns in 'y' have to be factors or integers.", call. = FALSE)

      }

      if(all(y_int_check)) y_type <- "integer" else y_type <- "factor"

      if(y_type != x_type) {

        stop(paste("Incompatible variable types: all columns",
                   "in 'x' and 'y' have to be either factors or integers."),
             call. = FALSE)

      }

    }

    ## factor processing --------

    if(x_type == "factor") {

      if(is.null(y)) {

        x_levs <- map(x, levels)

        all_levs <- reduce(x_levs, union)
        cmm_levs <- reduce(x_levs, intersect)

        if(length(setdiff(all_levs, cmm_levs)) > 0) {

          stop("Levels of all factor variables in data frame `x` must be identical.",
               call. = FALSE)

        }

        x <- map_dfc(x, ~as.integer(as.numeric(.x)))

      } else {

        for(i in 1:ncol(x)) {

          x_levs <- levels(x[[i]])
          y_levs <- levels(y[[i]])

          if(identical(x_levs, y_levs)) next

          new_levels <- sort(union(x_levs, y_levs))

          x[[i]] <- factor(x[[i]], new_levels)
          y[[i]] <- factor(y[[i]], new_levels)

        }

        x <- map_dfc(x, ~as.integer(as.numeric(.x)))
        y <- map_dfc(y, ~as.integer(as.numeric(.x)))


      }

    }

    ## computation of kappas --------

    f_kappa(x = as.matrix(x),
            y = if(!is.null(y)) as.matrix(y),
            method = method,
            as_data_frame = as_data_frame)

  }

# END --------
