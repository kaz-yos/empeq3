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
                                                 ## Covariate geenration
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
