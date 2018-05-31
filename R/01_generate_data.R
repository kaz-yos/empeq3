################################################################################
### Data generation functions for three-group empirical equipoise study
##
## Created on: 2018-05-24
## Author: Kazuki Yoshida
################################################################################

##' Generate bivariate standard normal distribution with correlation rho
##'
##' .. content for details ..
##'
##' @param n  Sample size
##' @param rho Correlation coefficient between X1 and X2
##'
##' @return data_frame containing two covariates X1 and X2. Both are marginally N(0,1). Their correlation is rho.
##'
##' @export
generate_bivariate_standard_normal_covariate <- function(n, rho) {
    X <- datagen3::generate_mvn_covariates(n = n,
                                           mu = c(0,0),
                                           Sigma = matrix(c(1,rho,
                                                            rho,1),
                                                          nrow = 2, byrow = TRUE))
    ## Fix covariate names
    names(X) <- gsub("Z", "X", names(X))
    X
}

##' Generate p-variate standard normal distribution with correlation rho
##'
##' .. content for details ..
##'
##' @param n Sample size
##' @param p Dimension
##' @param rho Correlation coefficient. corr(Xi, Xj) = rho^abs(i - j).
##'
##' @return data_frame containing p covariates X1 through Xp. Each one is marginally N(0,1). Their correlation structure is compound symmetry.
##'
##' @export
generate_p_dimensional_standard_normal_covariates <- function(n, p, rho) {
    mu <- rep(0, p)
    ## Create a compound symmetry type correlation matrix
    Sigma <- matrix(rep(NA, p*p), nrow = p)
    for (i in seq_len(p)) {
        for (j in seq_len(p)) {
            Sigma[i,j] <- rho^abs(i - j)
        }
    }

    X <- datagen3::generate_mvn_covariates(n = n,
                                           mu = mu,
                                           Sigma = Sigma)
    ## Fix covariate names
    names(X) <- gsub("Z", "X", names(X))
    X
}

##' Generate p covariate, one continuous, one count, and binary
##'
##' .. content for details ..
##'
##' @param n Sample size
##' @param p Dimension
##' @param rho Correlation coefficient. corr(Xi, Xj) = rho^abs(i - j) for the latent p-variate normal distribution.
##' @param lambda mean parameter for Poisson X2 variable
##' @param prev prevalence parameter vector for remaining binary variables. This must be p - 2.
##'
##' @return data_frame containing p covariates X1 through Xp. Each one is marginally N(0,1). Their correlation structure is compound symmetry.
##'
##' @export
generate_cont_count_bin_covariates <- function(n, p, rho, lambda, prev) {
    assertthat::assert_that(length(n) == 1)
    assertthat::assert_that(length(p) == 1)
    assertthat::assert_that(length(rho) == 1)
    assertthat::assert_that(length(lambda) == 1)
    assertthat::assert_that(length(prev) == (p - 2))

    X <- generate_p_dimensional_standard_normal_covariates(n, p, rho)
    ## Generate a Poisson variable with mean lambda
    X[[2]] <- qpois(pnorm(X[[2]]), lambda = lambda)
    ## Generate Bernoulli variables
    for (i in seq_along(prev)) {
        X[[2 + i]] <- qbinom(pnorm(X[[2 + i]]), size = 1, prob = prev[i])
    }

    X
}


##' Generate three-valued treatment with constraint on
##'
##' .. content for details ..
##'
##' @param X data_frame containing two covariates
##' @param alpha01 True intercept for the first linear predictor.
##' @param alpha02 True intercept for the second linear predictor.
##' @param alphaXm1 True coefficients for the first variable for the first linear predictor.
##' @param alphaXm2 True coefficients for the first variable for the second linear predictor.
##' @param gamma Ratio of the true coefficient for the second variable and the first variable.
##' @param sigma scaling of all covariate effects. A higher value means stronger covariate effects, i.e., less clinical equipoise.
##'
##' @return df with treatment (A) added.
##'
##' @export
generate_tri_treatment_from_two_covariates <- function(X, alpha01, alpha02, alphaXm1, alphaXm2, gamma, sigma) {
    assertthat::assert_that(length(alpha01) == 1)
    assertthat::assert_that(length(alpha02) == 1)
    assertthat::assert_that(length(alphaXm1) == 1)
    assertthat::assert_that(length(alphaXm2) == 1)

    datagen3::generate_tri_treatment(X,
                                     alphas1 = c(alpha01, sigma * alphaXm1, sigma * gamma * alphaXm1),
                                     alphas2 = c(alpha02, sigma * alphaXm2, sigma * gamma * alphaXm2))
}


##' Generate data from bivariate normal covariates (count outcome)
##'
##' .. content for details ..
##'
##' @param n Sample size
##'
##' @param rho correlation coefficient
##'
##' @param alphas True coefficients for the first and second treatment linear predictors. This vector should contain the intercept. alphas = c(alpha01, alphaXm1, alpha02, alphaXm2)
##' @param gamma Ratio of the true coefficient for the second variable and the first variable.
##' @param sigma scaling of all covariate effects. A higher value means stronger covariate effects, i.e., less clinical equipoise.
##'
##' @param beta0 Outcome model intercept coefficient
##' @param betaA Outcome model coefficient for I(A_i = 1) and I(A_i = 2)
##' @param betaX Outcome model coefficient vector for covariates X_i
##' @param betaXA Outcome model interaction coefficients for covariates. betaXA = c(betaXA1, betaXA2)
##'
##' @return a complete simulated data_frame
##'
##' @author Kazuki Yoshida
##'
##' @export
generate_bivariate_normal_data_count <- function(n,
                                                 ## Covariate generation
                                                 rho,
                                                 ## Treatment assignment
                                                 alphas,
                                                 gamma,
                                                 sigma,
                                                 ## Outcome assignment
                                                 beta0,
                                                 betaA,
                                                 betaX,
                                                 betaXA) {

    n_covariates <- 2
    assertthat::assert_that(length(n) == 1)
    assertthat::assert_that(length(rho) == 1)
    ## These alphas
    assertthat::assert_that(length(alphas) == 4)
    assertthat::assert_that(length(gamma) == 1)
    assertthat::assert_that(length(sigma) == 1)
    ## betaA = c(betaA1, betaA2)
    assertthat::assert_that(length(betaA) == 2)
    ## Only two covariates
    assertthat::assert_that(length(betaX) == n_covariates)
    ## Four interaction coefficients
    assertthat::assert_that(length(betaXA) == 2 * n_covariates)


    ## Extract parameters for use
    alpha01 <- alphas[1]
    alphaXm1 <- alphas[2]
    alpha02 <- alphas[3]
    alphaXm2 <- alphas[4]
    betaA1 <- betaA[1]
    betaA2 <- betaA[2]
    betaXA1 <- betaXA[seq_len(n_covariates)]
    betaXA2 <- betaXA[n_covariates + seq_len(n_covariates)]

    n %>%
        generate_bivariate_standard_normal_covariate(rho = rho) %>%
        generate_tri_treatment_from_two_covariates(alpha01 = alpha01,
                                                   alpha02 = alpha02,
                                                   alphaXm1 = alphaXm1,
                                                   alphaXm2 = alphaXm2,
                                                   gamma = gamma,
                                                   sigma = sigma) %>%
        datagen3::generate_count_outcome_log_tri_treatment(beta0 = beta0,
                                                           betaA1 = betaA1,
                                                           betaA2 = betaA2,
                                                           betaX = betaX,
                                                           betaXA1 = betaXA1,
                                                           betaXA2 = betaXA2)
}


##' Generate data from p-variate normal covariates (count outcome)
##'
##' .. content for details ..
##'
##' @param n Sample size
##'
##' @param p number of covariates. Must be 3 or more.
##' @param rho correlation coefficient
##' @param lambda mean parameter for X2
##' @param prev prevalence vector for X3 through Xp
##'
##' @param alphas True coefficients for the first and second treatment linear predictors. This vector should contain the intercept. alphas = c(alpha01, alphaXm1, alpha02, alphaXm2)
##' @param gamma Ratio of the true coefficient for the second variable and the first variable.
##' @param sigma scaling of all covariate effects. A higher value means stronger covariate effects, i.e., less clinical equipoise.
##'
##' @param beta0 Outcome model intercept coefficient
##' @param betaA Outcome model coefficient for I(A_i = 1) and I(A_i = 2)
##' @param betaX Outcome model coefficient vector for covariates X_i
##' @param betaXA Outcome model interaction coefficients for covariates. betaXA = c(betaXA1, betaXA2)
##'
##' @return a complete simulated data_frame
##'
##' @author Kazuki Yoshida
##'
##' @export
generate_p_norm_count_bin_data_count <- function(n,
                                                 ## Covariate generation
                                                 p,
                                                 rho,
                                                 lambda,
                                                 prev,
                                                 ## Treatment assignment
                                                 alphas,
                                                 gamma,
                                                 sigma,
                                                 ## Outcome assignment
                                                 beta0,
                                                 betaA,
                                                 betaX,
                                                 betaXA) {

    assertthat::assert_that(length(p) == 1)
    n_covariates <- p
    assertthat::assert_that(n_covariates >= 3)
    assertthat::assert_that(length(n) == 1)
    assertthat::assert_that(length(rho) == 1)
    assertthat::assert_that(length(lambda) == 1)
    assertthat::assert_that(length(prev) == (n_covariates - 2))
    ## These alphas
    assertthat::assert_that(length(alphas) == (n_covariates + 1) * 2)
    assertthat::assert_that(length(gamma) == 1)
    assertthat::assert_that(length(sigma) == 1)
    ## betaA = c(betaA1, betaA2)
    assertthat::assert_that(length(betaA) == 2)
    ## Only two covariates
    assertthat::assert_that(length(betaX) == n_covariates)
    ## Four interaction coefficients
    assertthat::assert_that(length(betaXA) == 2 * n_covariates)


    ## Extract parameters for use
    alpha01 <- alphas[1]
    alphaX1 <- alphas[1 + seq_len(n_covariates)]
    alpha02 <- alphas[n_covariates + 2]
    alphaX2 <- alphas[n_covariates + 2 + seq_len(n_covariates)]
    betaA1 <- betaA[1]
    betaA2 <- betaA[2]
    betaXA1 <- betaXA[seq_len(n_covariates)]
    betaXA2 <- betaXA[n_covariates + seq_len(n_covariates)]

    n %>%
        generate_p_dimensional_standard_normal_covariates(p = p, rho = rho) %>%
        datagen3::generate_tri_treatment(alphas1 = c(alpha01, sigma * alphaX1),
                                         alphas2 = c(alpha02, sigma * alphaX2)) %>%
        datagen3::generate_count_outcome_log_tri_treatment(beta0 = beta0,
                                                           betaA1 = betaA1,
                                                           betaA2 = betaA2,
                                                           betaX = betaX,
                                                           betaXA1 = betaXA1,
                                                           betaXA2 = betaXA2)
}
