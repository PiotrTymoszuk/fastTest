# Meta-estimates calculated with the inverse-variance weighting

#' Inverse variance-weighted meta-estimates.
#'
#' @description
#' The functions compute inverse variance-weighted meta-estimates of an effect
#' with the fixed and random between-study variance (DerSimonian-Lair algorithm).
#' While those tools are well suited for high throughput analyses of e.g.
#' transcriptomic data, the excellent R package `meta` with richer algorithm
#' and quality control options is recommended for  meta-analyses of clinical
#' studies.
#'
#' @references
#' Schwarzer G, Carpenter JR, Rücker G. Meta-Analysis with R. (2015)
#' doi:10.1007/978-3-319-21416-0
#'
#' @details
#' In case of lists, matrices or data frames provided as the function arguments,
#' the meta-estimates are computed in an element- and column-wise manner.
#' This means that each element or each column is treated as a a vector of
#' estimates or errors obtained in a single study.
#'
#' @return
#' A numeric vector, matrix, or data frame with the following statistics:
#' number of studies, Cochran's Q statistic of between-study variance,
#' p value for significant between-study variance determined from chi-square
#' distribution of Q, tau-square estimate of between-study variance,
#' meta-estimate calculated by the inverse-variance method, standard error of
#' the meta-estimate (SEM), lower and upper bound of the confidence interval,
#' test statistic Z, and p-value for significant meta-effect.
#'
#' @param y a numeric vector, list, matrix, or a data frame with estimates
#' obtained in single studies.
#' @param e a numeric vector, list, matrix, or a data frame with errors, e.g.
#' standard errors or other errors at the same scale as `y`.
#' @param type algorithm used for calculation of the meta-estimates. `fixed`
#' (default) assumes that there's no significant between-study variance.
#' `random` estimates the between-study variance (i.e. tau-square) with the
#' DerSimonian-Lair algorithm.
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

  f_meta <- function(y, e, ...) UseMethod('f_meta')

#' @rdname f_meta
#' @export

  f_meta.default <- function(y, e,
                             type = c('fixed', 'random'),
                             alternative = c('two.sided',
                                             'less',
                                             'greater'),
                             conf_level = 0.95,
                             as_data_frame = FALSE, ...) {

    ## entry control ------

    stopifnot(is.numeric(y))
    stopifnot(is.numeric(e))

    type <- match.arg(type[1], c('fixed', 'random'))

    alternative <- match.arg(alternative[1],
                             c('two.sided', 'less', 'greater'))

    stopifnot(is.numeric(conf_level))
    stopifnot(is.logical(as_data_frame))

    ## meta estimates --------

    result <- metaVec(y, e,
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

#' @rdname f_meta
#' @export

  f_meta.list <- function(y, e,
                          type = c('fixed', 'random'),
                          alternative = c('two.sided',
                                          'less',
                                          'greater'),
                          conf_level = 0.95,
                          as_data_frame = FALSE,
                          adj_method = 'none',
                          safely = FALSE, ...) {

    ## input control -----

    stopifnot(is.list(y))
    stopifnot(is.list(e))

    type <- match.arg(type[1], c('fixed', 'random'))

    alternative <- match.arg(alternative[1],
                             c('two.sided', 'less', 'greater'))

    stopifnot(is.numeric(conf_level))
    stopifnot(is.logical(as_data_frame))
    stopifnot(is.logical(safely))

    ## meta-estimates -------

    result <- metaList(y, e,
                       type = type,
                       alternative = alternative,
                       conf_level = conf_level,
                       crash = !safely)

    p_adjusted <- NULL

    if(adj_method != 'none') {

      result <- cbind(result,
                      p_adjusted = p.adjust(result[, 10],
                                            method = adj_method))

    }

    if(!as_data_frame) return(result)

    result <- as.data.frame(result)

    rownames_to_column(result, 'variable')

  }

#' @rdname f_meta
#' @export

  f_meta.matrix <- function(y, e,
                            type = c('fixed', 'random'),
                            alternative = c('two.sided',
                                            'less',
                                            'greater'),
                            conf_level = 0.95,
                            as_data_frame = FALSE,
                            adj_method = 'none',
                            safely = FALSE, ...) {

    ## input control ------

    stopifnot(is.matrix(y))
    stopifnot(is.matrix(e))

    if(!is.numeric(y) | !is.numeric(e)) {

      stop("'y' and 'e' have to be numeric matrices or data frames.")

    }

    type <- match.arg(type[1], c('fixed', 'random'))

    alternative <- match.arg(alternative[1],
                             c('two.sided', 'less', 'greater'))

    stopifnot(is.numeric(conf_level))
    stopifnot(is.logical(as_data_frame))
    stopifnot(is.character(adj_method))
    stopifnot(is.logical(safely))

    ## meta estimates --------

    result <- metaMtx(y, e,
                      type = type,
                      alternative = alternative,
                      conf_level = conf_level,
                      crash = !safely)


    p_adjusted <- NULL

    if(adj_method != 'none') {

      result <- cbind(result,
                      p_adjusted = p.adjust(result[, 10],
                                            method = adj_method))

    }

    if(!as_data_frame) return(result)

    result <- as.data.frame(result)

    rownames_to_column(result, 'variable')

  }

#' @rdname f_meta
#' @export

  f_meta.data.frame <- function(y, e,
                                type = c('fixed', 'random'),
                                alternative = c('two.sided',
                                                'less',
                                                'greater'),
                                conf_level = 0.95,
                                as_data_frame = FALSE,
                                adj_method = 'none',
                                safely = FALSE, ...) {

    f_meta(as.matrix(y),
           as.matrix(e),
           type = type,
           alternative = alternative,
           conf_level = conf_level,
           as_data_frame = as_data_frame,
           adj_method = adj_method,
           safely = safely)

  }

# END -------
