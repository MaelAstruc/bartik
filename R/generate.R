#' Generate shift-share data
#'
#' `generate` produces a data set for a shift-share analysis. The data is
#' constructed following the Bartik model where we are interested in the impact
#' of \eqn{x_l} the employment growth in location \eqn{l} on \eqn{y_l} the wage growth in
#' the same region. \eqn{x_l} is decomposed as \eqn{z_{lk}} the share of industry \eqn{k} in
#' the location \eqn{l} and \eqn{g_{lk}} the employment growth for the industry \eqn{k} in
#' location \eqn{l}. The last step is to decompose this local industry growth
#' between the global industry growth \eqn{g_k} and the local specific growth
#' \eqn{\tilde{g}_{lk}}.
#'
#' Given this framework, we have the model:
#' \deqn{
#' y_l = \beta x_l + e_l
#' x_l = \sum_k z_l{lk} g_{lk}
#' g_{lk} = \gamma g_k + (1-\gamma) \tilde{g}_{lk}
#' }
#'
#' With \eqn{\beta} the parameter of interest and \eqn{e_l} the structural error.
#'
#' The endogeneity of the model comes from the correlation between
#' \eqn{\tilde{g}_{lk}} and \eqn{e_l}.
#'
#' The Bartik instrument is constructed as :
#' \deqn{
#' B_l = \sum_k z_{lk} g_{k}
#' }
#'
#' This instrument is exogenous if \eqn{z_{lk}} is not correlated with \eqn{e_l}.
#'
#' @param L an integer. The number of locations.
#' @param K an integer. The number of groups.
#' @param alpha a numeric. The model intercept.
#' @param beta a numeric. The parameter of interest.
#' @param g_global a numeric between 0 and 1. The parameter \eqn{\gamma}.
#' @param g_endo a numeric between 0 and 1. The share of \eqn{\tilde{g}_{lk}}
#'   determined by \eqn{e_l}.
#'
#' @return A list with the generated data, the share and the shift.
#' @export
#'
#' @examples
#' # 1. No endogeneity and the local specific growth is not important
#'
#' dataset <- generate(L = 10000, K = 50, alpha = 1, beta = 3,
#'                     g_global = 1, g_endo = 0)
#' lm(y ~ x + e, data = dataset$data)
#' lm(x ~ e, data = dataset$data)
#' lm(y ~ x, data = dataset$data)
#' lm(x ~ bartik, data = dataset$data)
#' ivreg::ivreg(y ~ 1 | x | bartik, data = dataset$data)
#'
#' # Here the OLS estimation of beta is close to the true parameter
#' # The bartik gives an equivalent result
#'
#' # 2. Strong endogeneity and the local specific growth is not important
#'
#' dataset <- generate(L = 10000, K = 50, alpha = 1, beta = 3,
#'                     g_global = 1, g_endo = 1)
#' lm(y ~ x + e, data = dataset$data)
#' lm(x ~ e, data = dataset$data)
#' lm(y ~ x, data = dataset$data)
#' lm(x ~ bartik, data = dataset$data)
#' ivreg::ivreg(y ~ 1 | x | bartik, data = dataset$data)
#'
#' # The estimation of beta is still accurate
#' # The Bartik instrument gives the same results
#'
#' # 3. Strong endogeneity and the local specific growth is important
#'
#' dataset <- generate(L = 10000, K = 50, alpha = 1, beta = 3,
#'                     g_global = 0, g_endo = 1)
#' lm(y ~ x + e, data = dataset$data)
#' lm(x ~ e, data = dataset$data)
#' lm(y ~ x, data = dataset$data)
#' lm(x ~ bartik, data = dataset$data)
#' ivreg::ivreg(y ~ 1 | x | bartik, data = dataset$data)
#'
#' # The estimation of beta is biased
#' # The Bartik instrument gives the same biased results
#' # The instrument is too weak
#'
#' # 4. No endogeneity and the local specific growth is important
#'
#' dataset <- generate(L = 10000, K = 50, alpha = 1, beta = 3,
#'                     g_global = 0, g_endo = 0)
#' lm(y ~ x + e, data = dataset$data)
#' lm(x ~ e, data = dataset$data)
#' lm(y ~ x, data = dataset$data)
#' lm(x ~ bartik, data = dataset$data)
#' ivreg::ivreg(y ~ 1 | x | bartik, data = dataset$data)
#'
#' # The estimation of beta is accurate
#' # The Bartik instrument gives inaccurate results
#' # The instrument is too weak
#'
#' # 5. Strong endogeneity and the local specific growth is moderately important
#'
#' dataset <- generate(L = 10000, K = 50, alpha = 1, beta = 3,
#'                     g_global = 0.5, g_endo = 1)
#' lm(y ~ x + e, data = dataset$data)
#' lm(y ~ x, data = dataset$data)
#' lm(x ~ e, data = dataset$data)
#' lm(x ~ bartik, data = dataset$data)
#' ivreg::ivreg(y ~ 1 | x | bartik, data = dataset$data)
#'
#' # The estimation of beta is biased
#' # The Bartik instrument gives more accurate result
generate <- function(L, K, alpha = 0, beta = 1,
                     g_global = 0.5, g_endo = 0.5) {
  # Check that the dimensions are correct
  stopifnot(is.numeric(L),
            is.numeric(K))

  # Check that the parameters are well defined
  stopifnot(g_global >= 0 && g_global <= 1,
            g_endo >= 0 && g_endo <= 1)

  # First we have the error
  e_l <- stats::rnorm(L)

  # Then we construct the growth with two components
  g_k <- stats::rnorm(K)
  g_k <- matrix(rep(g_k, L), nrow = L, byrow = TRUE)

  g_lk <- stats::rnorm(L*K)
  g_lk <- matrix(g_lk, nrow = L, byrow = TRUE)
  g_lk <- (1 - g_endo) * g_lk + g_endo * e_l

  shift <- list(global = g_k,
                local = g_lk)

  g_lk <- g_global * g_k + (1 - g_global) * g_lk

  # Then we construct the share
  z_lk <- stats::runif(L*K)
  z_lk <- matrix(z_lk, nrow = L, byrow = TRUE)
  z_lk <- z_lk / rowSums(z_lk)
  share <- z_lk

  # We construct the x variable
  x_l <- share * g_lk
  x_l <- rowSums(x_l)

  # We can also construct the Bartik
  bartik <- share * g_k
  bartik <- rowSums(bartik)

  # Finally we construct the y variable
  y_l = alpha + beta * x_l + e_l

  data <- data.frame(y = y_l, x = x_l, e = e_l, bartik = bartik)

  list(data = data,
       share = share,
       shift = shift)
}
