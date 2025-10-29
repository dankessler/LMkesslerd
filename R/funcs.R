##' Construct a confidence interval for a linear combination of regression coefficients
##'
##' @param y Response vector
##' @param X Design matrix
##' @param alpha Significance level (e.g., 0.05 for 95% CI)
##' @param gamma Vector specifying the linear combination
##' @param sigma2 Optional known variance. If not provided, uses t-distribution
##' @return Numeric vector of length 2: lower and upper confidence interval bounds
##' @examples
##' y <- c(1,2,3,4)
##' X <- cbind(1, c(1,2,3,4))
##' gamma <- c(0,1)
##' make_CI(y, X, 0.05, gamma)
##' make_CI(y, X, 0.05, gamma, sigma2=1)
##' @export
make_CI <- function(y, X, alpha, gamma, sigma2 = NULL) {
  n <- nrow(X)
  p <- ncol(X)

  # Fit the linear model
  fit <- lm(y ~ X - 1)  # -1 to exclude intercept if X already includes it
  beta_hat <- coef(fit)
  residuals <- resid(fit)

  if (is.null(sigma2)) {
    # Estimate variance
    sigma2_hat <- sum(residuals^2) / (n - p)

    # Compute standard error for gamma^T * beta_hat
    se_gamma_beta <- sqrt(sigma2_hat * t(gamma) %*% solve(t(X) %*% X) %*% gamma)

    # Compute the critical value from the t-distribution
    crit_value <- qt(1 - alpha / 2, df = n - p)
  } else {
    # Compute standard error for gamma^T * beta_hat using known variance
    se_gamma_beta <- sqrt(sigma2 * t(gamma) %*% solve(t(X) %*% X) %*% gamma)

    # Compute the critical value from the normal distribution
    crit_value <- qnorm(1 - alpha / 2)
  }

  # Compute the confidence interval
  lower_bound <- t(gamma) %*% beta_hat - crit_value * se_gamma_beta
  upper_bound <- t(gamma) %*% beta_hat + crit_value * se_gamma_beta

  return(c(lower_bound, upper_bound))
}



##' Compute the OLS estimate of regression coefficients
##'
##' @param y Response vector
##' @param X Design matrix
##' @return Numeric vector of OLS estimates
##' @examples
##' y <- c(1,2,3,4)
##' X <- cbind(1, c(1,2,3,4))
##' find_beta(y, X)
##' @export
find_beta <- function(y, X) {
  fit <- lm(y ~ X - 1)  # -1 to exclude intercept if X already includes it
  beta_hat <- coef(fit)
  return(beta_hat)
}

##' Compute constrained OLS estimate subject to linear constraints
##'
##' @param y Response vector
##' @param X Design matrix
##' @param C Constraint matrix
##' @param d Constraint vector
##' @return Numeric vector of constrained OLS estimates
##' @examples
##' y <- c(1,2,3,4)
##' X <- cbind(1, c(1,2,3,4))
##' C <- matrix(c(0,1), nrow=1)
##' d <- 0
##' # find_betaH(y, X, C, d)
##' @export
find_betaH <- function(y, X, C, d) {
  beta_hat <- find_beta(y, X)
  beta_hat_H <- beta_hat - solve(t(X) %*% X) %*% t(C) %*% solve(C %*% solve(t(X) %*% X) %*% t(C)) %*% (C %*% beta_hat - d)
  return(beta_hat_H)
}

##' Compute F-statistic and p-value for linear restrictions (model comparison)
##'
##' @param y Response vector
##' @param X Design matrix
##' @param C Constraint matrix
##' @param d Constraint vector
##' @return List with F_statistic and p_value
##' @examples
##' y <- c(1,2,3,4)
##' X <- cbind(1, c(1,2,3,4))
##' C <- matrix(c(0,1), nrow=1)
##' d <- 0
##' # compute_F_stat(y, X, C, d)
##' @export
f_test_RSS <- function(y, X, C, d) {
  n <- nrow(X)
  p <- ncol(X)
  q <- nrow(C)

  # Unrestricted model
  beta_hat_unrestricted <- find_beta(y, X)
  residuals_unrestricted <- y - X %*% beta_hat_unrestricted
  ssr_unrestricted <- sum(residuals_unrestricted^2)

  # Restricted model
  beta_hat_restricted <- find_betaH(y, X, C, d)
  residuals_restricted <- y - X %*% beta_hat_restricted
  ssr_restricted <- sum(residuals_restricted^2)

  # Compute F-statistic
  F_stat <- ((ssr_restricted - ssr_unrestricted) / q) / (ssr_unrestricted / (n - p))

  # Compute p-value
  p_value <- 1 - pf(F_stat, df1 = q, df2 = n - p)

  return(list(F_statistic = F_stat, p_value = p_value))
}

##' Compute F-statistic and p-value using quadratic form of constraints
##'
##' @param y Response vector
##' @param X Design matrix
##' @param C Constraint matrix
##' @param d Constraint vector
##' @return List with F_statistic and p_value
##' @examples
##' y <- c(1,2,3,4)
##' X <- cbind(1, c(1,2,3,4))
##' C <- matrix(c(0,1), nrow=1)
##' d <- 0
##' # compute_F_stat_quadratic(y, X, C, d)
##' @export
f_test_quad <- function(y, X, C, d) {
  n <- nrow(X)
  p <- ncol(X)
  q <- nrow(C)

  beta_hat_unrestricted <- find_beta(y, X)
  s2 = sum((y - X %*% beta_hat_unrestricted)^2) / (n - p)
  diff <- C %*% beta_hat_unrestricted - d
  cov_matrix <- C %*% solve(t(X) %*% X) %*% t(C)
  F_stat <- as.numeric(t(diff) %*% solve(cov_matrix) %*% diff) / (q * s2)
  p_value <- 1 - pf(F_stat, df1 = q, df2 = n - p)
  return(list(F_statistic = F_stat, p_value = p_value))
}
