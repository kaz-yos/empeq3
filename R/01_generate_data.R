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
##'
##' @return df with treatment (A) added.
##'
##' @export
generate_tri_treatment_from_two_covariates <- function(X, alpha01, alpha02, alphaXm1, alphaXm2, gamma) {
    assertthat::assert_that(length(alpha01) == 1)
    assertthat::assert_that(length(alpha02) == 1)
    assertthat::assert_that(length(alphaXm1) == 1)
    assertthat::assert_that(length(alphaXm2) == 1)

    datagen3::generate_tri_treatment(X,
                                     alphas1 = c(alpha01, alphaXm1, gamma * alphaXm1),
                                     alphas2 = c(alpha02, alphaXm2, gamma * alphaXm2))
}


##' Generate data from bivariate normal covariates (count outcome)
##'
##' .. content for details ..
##'
##' @param n Sample size
##'
##' @param alpha01 True intercept for the first linear predictor.
##' @param alpha02 True intercept for the second linear predictor.
##' @param alphaXm1 True coefficients for the first variable for the first linear predictor.
##' @param alphaXm2 True coefficients for the first variable for the second linear predictor.
##' @param gamma Ratio of the true coefficient for the second variable and the first variable.
##'
##' @param beta0 Outcome model intercept coefficient
##' @param betaA Outcome model coefficient for I(A_i = 1) and I(A_i = 2)
##' @param betaX Outcome model coefficient vector for covariates X_i
##' @param betaXA1 Outcome model interaction coefficients for covariates. betaXA = c(betaXA1, betaXA2)
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
                                                 alpha01,
                                                 alpha02,
                                                 alphaXm1,
                                                 alphaXm2,
                                                 gamma,
                                                 ## Outcome assignment
                                                 beta0,
                                                 betaA1,
                                                 betaA2,
                                                 betaX,
                                                 betaXA1,
                                                 betaXA2) {

    n %>%
        generate_bivariate_standard_normal_covariate(rho = rho) %>%
        generate_tri_treatment_from_two_covariates(alpha01 = alpha01,
                                                   alpha02 = alpha02,
                                                   alphaXm1 = alphaXm1,
                                                   alphaXm2 = alphaXm2,
                                                   gamma = gamma) %>%
        datagen3::generate_count_outcome_log_tri_treatment(beta0 = beta0,
                                                           betaA1 = betaA1,
                                                           betaA2 = betaA2,
                                                           betaX = betaX,
                                                           betaXA1 = betaXA1,
                                                           betaXA2 = betaXA2)
}
