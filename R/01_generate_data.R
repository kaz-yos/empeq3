################################################################################
### Data generation functions for three-group empirical equipoise study
##
## Created on: 2018-05-24
## Author: Kazuki Yoshida
################################################################################

##' Generate bivariate standard normal distribution with correlation rho
##'
##' .. content for \details{} ..
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
##' .. content for \details{} ..
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
    assert_that(length(alpha01) == 1)
    assert_that(length(alpha02) == 1)
    assert_that(length(alphaXm1) == 1)
    assert_that(length(alphaXm2) == 1)

    datagen3::generate_tri_treatment(X,
                                     alphas1 = c(alpha01, alphaXm1, gamma * alphaXm1),
                                     alphas2 = c(alpha02, alphaXm2, gamma * alphaXm2))
}
