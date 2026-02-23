# Correlation tests via permutation and bootstrap

#' Permutation and bootstrap tests for correlations.
#'
#' @description
#' Permutation and bootstrap tests for correlation which are not covered by
#' base R.
#'
#' @details
#' The permutation tests are done by permuting one of the correlation partners
#' (i.e. re-sampling without replacement). The bootstrap tests employ
#' re-sampling of data points (i.e. X and Y pairs) with replacement.
#' Please note that bootstrapping is likely going to break trends in data
#' investigated by the xi family of correlation coefficients: the bootstrapped
#' mean is likely to be substantially lower than the correlation coefficient
#' for the entire data set. Nevertheless, bootstrapped xi coefficients with
#' the respective confidence intervals are good estimates of interpolation of the
#' data trends.
#' Please note, that bootstrapped p values refer to the two-sided setting.
#' If a numeric matrix or data frame is provided as the `x` argument, the
#' correlations are performed in pairwise manner between columns of the input
#' object.
#'
#' @references
#' Chatterjee S. A New Coefficient of Correlation. J Am Stat Assoc (2021)
#' 116:2009–2022. doi:10.1080/01621459.2020.1758115
#'
#' @return
#' A numeric vector, matrix, or a data frame with the following elements:
#' numbers of pairwise complete observations, the requested correlation
#' coefficients in the entire data set, bootstrapped means and confidence
#' interval bounds, numbers of successful iterations of the algorithm,
#' numbers of re-samples in favor of the null hypothesis H0, numbers of
#' re-samples in favor of the alternative hypothesis H1, p values, and,
#' optionally p values adjusted for multiple testing.
#'
#' @param x a numeric vector, list, or data frame.
#' @param y a numeric vector.
#' @param type type of the test: permutation (default) or bootstrap.
#' @param method type of covariance or correlation to be analyzed:
#' Pearson's r (default), Spearman's rho, Kendall's TauA, Kendall's TauB,
#' Chatterjee's Xi (xiA: no assumption on ties, xiB: tie correction included).
#' @param alternative type of the alternative hypothesis concerning correlation
#' coefficient sign. Ignored if `type = 'bootstrap'`.
#' @param ci_type type of confidence intervals: BCA (default) or percentile.
#' @param conf_level confidence level used for computation of the confidence
#' intervals.
#' @param as_data_frame should the output be formatted as a data frame? This may
#' render the computation slower.
#' @param adj_method multiple testing adjustment method, as specified for
#' \code{\link[stats]{p.adjust}}.
#' @param n_iter number of the algorithm's iterations.
#' @param ... extra arguments passed to methods.
#'
#' @export

  f_cor_test <- function(x, ...) UseMethod('f_cor_test')

#' @rdname f_cor_test
#' @export

  f_cor_test.default <- function(x, y,
                                 type = c('permutation',
                                          'bootstrap'),
                                 method = c('pearson',
                                            'spearman',
                                            'kendallA',
                                            'kendallB',
                                            'xiA',
                                            'xiB'),
                                 alternative = c("two.sided",
                                                 "less",
                                                 "greater"),
                                 ci_type = c('bca', 'percentile'),
                                 conf_level = 0.95,
                                 as_data_frame = FALSE,
                                 n_iter = 1000, ...) {

    ## entry control -------

    stopifnot(is.atomic(x))
    stopifnot(is.atomic(y))

    if(!is.numeric(x) | !is.numeric(y)) {

      stop("'x' and 'y' have to be numeric vectors.", call. = FALSE)

    }

    type <- match.arg(type[1], c('permutation',  'bootstrap'))

    method <- match.arg(method[1],
                        c('pearson',
                          'spearman',
                          'kendallA',
                          'kendallB',
                          'xiA',
                          'xiB'))

    alternative <- match.arg(alternative[1],
                             c("two.sided", "less", "greater"))

    ci_type <- match.arg(ci_type[1], c('bca', 'percentile'))

    stopifnot(is.numeric(conf_level))

    if(conf_level < 0 | conf_level > 1) {

      stop("'conf_level' has to be in the [0, 1] range.", call. = FALSE)

    }

    stopifnot(is.logical(as_data_frame))
    stopifnot(is.numeric(n_iter))
    stopifnot(n_iter >= 1)

    n_iter <- as.integer(n_iter)

    ## correlation testing  --------

    if(type == 'permutation') {

      result <- permCorVec(x, y,
                           method = method,
                           alternative = alternative,
                           n_iter = n_iter)

    } else {

      result <- bootCorVec(x, y,
                           method = method,
                           ci_type = ci_type,
                           conf_level = conf_level,
                           n_iter = n_iter)

    }

    if(!as_data_frame) return(result)

    result <- matrix(result,
                     nrow = 1,
                     dimnames = list(NULL,
                                     names(result)))

    as.data.frame(result)

  }

#' @rdname f_cor_test
#' @export

  f_cor_test.list <- function(x,
                              type = c('permutation',
                                       'bootstrap'),
                              method = c('pearson',
                                         'spearman',
                                         'kendallA',
                                         'kendallB',
                                         'xiA',
                                         'xiB'),
                              alternative = c("two.sided",
                                              "less",
                                              "greater"),
                              ci_type = c('bca', 'percentile'),
                              conf_level = 0.95,
                              as_data_frame = FALSE,
                              n_iter = 1000, ...) {

    ## input control -------

    stopifnot(is.list(x))

    if(length(x) < 2) {

      stop("'x' has to have at least two elements.", call. = FALSE)

    }

    if(!is.numeric(x[[1]]) | !is.numeric(x[[2]])) {

      stop("A list of numeric elements is required.", call. = FALSE)

    }

    ## correlations ---------

    f_cor_test(x[[1]], x[[2]],
               type = type,
               method = method,
               alternative = alternative,
               ci_type = ci_type,
               conf_level = conf_level,
               as_data_frame = as_data_frame,
               n_iter = n_iter)

  }

#' @rdname f_cor_test
#' @export

  f_cor_test.matrix <- function(x,
                                type = c('permutation',
                                         'bootstrap'),
                                method = c('pearson',
                                           'spearman',
                                           'kendallA',
                                           'kendallB',
                                           'xiA',
                                           'xiB'),
                                alternative = c("two.sided",
                                                "less",
                                                "greater"),
                                ci_type = c('bca', 'percentile'),
                                conf_level = 0.95,
                                as_data_frame = FALSE,
                                n_iter = 1000,
                                adj_method = 'none', ...) {

    stopifnot(is.matrix(x))

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric matrix.", call. = FALSE)

    }

    type <- match.arg(type[1], c('permutation',  'bootstrap'))

    method <- match.arg(method[1],
                        c('pearson',
                          'spearman',
                          'kendallA',
                          'kendallB',
                          'xiA',
                          'xiB'))

    alternative <- match.arg(alternative[1],
                             c("two.sided", "less", "greater"))

    ci_type <- match.arg(ci_type[1], c('bca', 'percentile'))

    stopifnot(is.numeric(conf_level))

    if(conf_level < 0 | conf_level > 1) {

      stop("'conf_level' has to be in the [0, 1] range.", call. = FALSE)

    }

    stopifnot(is.logical(as_data_frame))
    stopifnot(is.numeric(n_iter))
    stopifnot(n_iter >= 1)

    n_iter <- as.integer(n_iter)

    stopifnot(is.character(adj_method))

    ## correlations -------

    if(type == 'permutation') {

      result <- permCorMtx(x, method, alternative, n_iter)

    } else {

      result <- bootCorMtx(x, method, ci_type, conf_level, n_iter)

    }

    if(adj_method != 'none') {

      result <- cbind(result,
                      p_adjusted = p.adjust(result[, 11],
                                            method = adj_method))

    }

    if(!as_data_frame) return(result)

    result <- as.data.frame(result)

    if(!is.null(colnames(x))) {

      result[['variable1']] <- colnames(x)[result[['variable1']]]
      result[['variable2']] <- colnames(x)[result[['variable2']]]

    }

    result

  }

#' @rdname f_cor_test
#' @export

  f_cor_test.data.frame <- function(x,
                                     type = c('permutation',
                                              'bootstrap'),
                                     method = c('pearson',
                                                'spearman',
                                                'kendallA',
                                                'kendallB',
                                                'xiA',
                                                'xiB'),
                                     alternative = c("two.sided",
                                                     "less",
                                                     "greater"),
                                     ci_type = c('bca', 'percentile'),
                                     conf_level = 0.95,
                                     as_data_frame = FALSE,
                                     n_iter = 1000,
                                     adj_method = 'none', ...) {

    f_cor_test(as.matrix(x),
               type = type,
               method = method,
               alternative = alternative,
               ci_type = ci_type,
               conf_level = conf_level,
               as_data_frame = as_data_frame,
               n_iter = n_iter,
               adj_method = adj_method)

  }

# END ------
