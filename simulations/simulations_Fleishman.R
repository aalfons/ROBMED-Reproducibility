# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## load required packages
library("SimMultiCorrData")
library("car")
library("robmed")


## additional functions

# compute Fleishman's polynomial transformation of normally distributed values
get_Fleishman <- function(x, c0, c1, c2, c3) c0 + c1 * x + c2 * x^2 + c3 * x^3

# compute mean of distribution after Fleishman's polynomial transformation
get_mean <- function(c0, c1, c2, c3, mu = 0, sigma = 1) {
  # compute moments
  mu_1 <- mu
  mu_2 <- mu^2 + sigma^2
  mu_3 <- mu^3 + 3 * mu * sigma^2
  # compute mean
  c0 + c1 * mu_1 + c2 * mu_2 + c3 * mu_3
}


## control parameters for simulation
seed <- 20200309             # seed for the random number generator
sample_sizes <- c(100, 250)  # number of observations
K <- 1000                    # number of simulation runs
R <- 5000                    # number of bootstrap samples
# coefficients in mediation model
effect <- 0.4
a <- rep(effect, 2)
b <- c(effect, 0)
c <- rep(effect, 2)
# error scales in mediation model
sigma_m <- sqrt(1 - a^2)
sigma_y <- sqrt(1 - b^2 - c^2 - 2*a*b*c)

## skewness and kurtosis
# define grid of target values for those parameters
skewness <- seq(-1, 1, by = 0.5)
kurtosis <- seq(-1, 4, by = 1)
grid <- expand.grid(Skewness = skewness, Kurtosis = kurtosis)
# compute coefficients of Fleishman's polynomial transformation
parameter_list <- mapply(function(skewness, kurtosis) {
  tryCatch({
    coefficients <- find_constants(method = "Fleishman", skews = skewness,
                                   skurts = kurtosis)
    data.frame(Skewness = skewness, Kurtosis = kurtosis,
               t(coefficients$constants), Valid = coefficients$valid)
  }, warning = function(condition) NULL, error = function(condition) NULL)
}, skewness = grid$Skewness, kurtosis = grid$Kurtosis, SIMPLIFY = FALSE)
parameters <- do.call(rbind, parameter_list)

## control parameters for methods
level <- 0.95                                        # confidence level
lmrob_control <- reg_control(max_iterations = 5000)  # MM-estimator


## run simulation
cat(paste(Sys.time(), ": starting ...\n"))
results_list <- lapply(sample_sizes, function(n) {

  ## use the same seed for each sample size
  ## (then simulations can be run separately for different sample sizes)
  set.seed(seed)

  ## print sample size
  cat(paste(Sys.time(), sprintf(": n = %d\n", n)))

  ## perform simulation for current sample size
  results_n <- lapply(1:K, function(k) {

    ## print simulation run
    cat(paste(Sys.time(), sprintf(":   run = %d\n", k)))

    ## generate independent variable and error terms such that results are
    ## comparable across different values of the parameters
    x <- rnorm(n)
    # First standard normal values are drawn, which are later transformed
    # to the given skewness and kurtorsis by a polynomial transformation.
    e_normal_m <- rnorm(n)
    e_normal_y <- rnorm(n)

    ## ensure that the same bootstrap samples are used for all parameter
    ## settings and all bootstrap tests for maximum comparability
    indices <- boot_samples(n, R = R)

    ## loop over different error distributions
    results_nk <- mapply(function(skewness, kurtosis, c0, c1, c2, c3) {

      ## transform errors
      # Fleishman's polynomial transformation
      e_m <- get_Fleishman(e_normal_m, c0 = c0, c1 = c1, c2 = c2, c3 = c3)
      e_y <- get_Fleishman(e_normal_y, c0 = c0, c1 = c1, c2 = c2, c3 = c3)
      # center error terms
      e_m <- e_m - get_mean(c0 = c0, c1 = c1, c2 = c2, c3 = c3)
      e_y <- e_y - get_mean(c0 = c0, c1 = c1, c2 = c2, c3 = c3)

      ## loop over different parameter settings
      results_nke <- mapply(function(a, b, c, sigma_m, sigma_y) {

        # generate hypothesized mediator and dependent variable
        m <- a * x + sigma_m * e_m
        y <- b * m + c * x + sigma_y * e_y
        simulated_data <- data.frame(x, y, m)

        # OLS regression
        df_ols <- tryCatch({
          # fit mediation model
          ols_fit <- fit_mediation(simulated_data, "x", "y", "m",
                                   method = "regression", robust = FALSE,
                                   family = "gaussian")
          # perform test for indirect effect
          ols_boot <- test_mediation(ols_fit, test = "boot", level = level,
                                     indices = indices)
          ols_sobel <- test_mediation(ols_fit, test = "sobel")
          # extract indirect effect
          ab_boot <- c(coef(ols_boot, parm = "ab", type = "boot"),
                       coef(ols_boot, parm = "ab", type = "data"),
                       use.names = FALSE)
          ab_sobel <- c(NA_real_, coef(ols_sobel, parm = "ab"),
                        use.names = FALSE)
          # check for significant indirect effect
          reject_boot <- prod(ols_boot$ci) > 0
          reject_sobel <- ols_sobel$p_value < 1 - level
          # combine into data frame
          data.frame(n = n, Run = k, Skewness = skewness, Kurtosis = kurtosis,
                     a = a, b = b, Method = c("ols_boot", "ols_sobel"),
                     ab_boot = c(ab_boot[1], ab_sobel[1]),
                     ab_data = c(ab_boot[2], ab_sobel[2]),
                     reject = c(reject_boot, reject_sobel),
                     stringsAsFactors = FALSE)
        }, error = function(condition) NULL)

        # Box-Cox transformations (with negative values allowed) followed by OLS
        df_bcn <- tryCatch({
          # transform variables
          transformed_list <- lapply(simulated_data, function(variable) {
            bcn <- powerTransform(variable, family = "bcnPower")
            bcnPower(variable, lambda = bcn$lambda, gamma = bcn$gamma)
          })
          transformed_data <- as.data.frame(transformed_list)
          # fit mediation model
          bcn_fit <- fit_mediation(transformed_data, "x", "y", "m",
                                   method = "regression", robust = FALSE,
                                   family = "gaussian")
          # perform test for indirect effect
          bcn_boot <- test_mediation(bcn_fit, test = "boot", level = level,
                                     indices = indices)
          bcn_sobel <- test_mediation(bcn_fit, test = "sobel")
          # extract indirect effect
          ab_boot <- c(coef(bcn_boot, parm = "ab", type = "boot"),
                       coef(bcn_boot, parm = "ab", type = "data"),
                       use.names = FALSE)
          ab_sobel <- c(NA_real_, coef(bcn_sobel, parm = "ab"),
                        use.names = FALSE)
          # check for significant indirect effect
          reject_boot <- prod(bcn_boot$ci) > 0
          reject_sobel <- bcn_sobel$p_value < 1 - level
          # combine into data frame
          data.frame(n = n, Run = k, Skewness = skewness, Kurtosis = kurtosis,
                     a = a, b = b, Method = c("bcn_boot", "bcn_sobel"),
                     ab_boot = c(ab_boot[1], ab_sobel[1]),
                     ab_data = c(ab_boot[2], ab_sobel[2]),
                     reject = c(reject_boot, reject_sobel),
                     stringsAsFactors = FALSE)
        }, error = function(condition) NULL)

        # Huberized covariance matrix
        df_huber <- tryCatch({
          # fit mediation model
          huber_fit <- fit_mediation(simulated_data, "x", "y", "m",
                                     method = "covariance", robust = TRUE)
          # perform test for indirect effect
          huber_boot <- test_mediation(huber_fit, test = "boot", level = level,
                                       indices = indices)
          huber_sobel <- test_mediation(huber_fit, test = "sobel")
          # extract indirect effect
          ab_boot <- c(coef(huber_boot, parm = "ab", type = "boot"),
                       coef(huber_boot, parm = "ab", type = "data"),
                       use.names = FALSE)
          ab_sobel <- c(NA_real_, coef(huber_sobel, parm = "ab"),
                        use.names = FALSE)
          # check for significant indirect effect
          reject_boot <- prod(huber_boot$ci) > 0
          reject_sobel <- huber_sobel$p_value < 1 - level
          # combine into data frame
          data.frame(n = n, Run = k, Skewness = skewness, Kurtosis = kurtosis,
                     a = a, b = b, Method = c("huber_boot", "huber_sobel"),
                     ab_boot = c(ab_boot[1], ab_sobel[1]),
                     ab_data = c(ab_boot[2], ab_sobel[2]),
                     reject = c(reject_boot, reject_sobel),
                     stringsAsFactors = FALSE)
        }, error = function(condition) NULL)

        # median regression
        df_median <- tryCatch({
          # fit mediation model
          median_fit <- fit_mediation(simulated_data, "x", "y", "m",
                                      method = "regression", robust = "median")
          # perform test for indirect effect
          median_boot <- test_mediation(median_fit, test = "boot", level = level,
                                        indices = indices)
          median_sobel <- test_mediation(median_fit, test = "sobel")
          # extract indirect effect
          ab_boot <- c(coef(median_boot, parm = "ab", type = "boot"),
                       coef(median_boot, parm = "ab", type = "data"),
                       use.names = FALSE)
          ab_sobel <- c(NA_real_, coef(median_sobel, parm = "ab"),
                        use.names = FALSE)
          # check for significant indirect effect
          reject_boot <- prod(median_boot$ci) > 0
          reject_sobel <- median_sobel$p_value < 1 - level
          # combine into data frame
          data.frame(n = n, Run = k, Skewness = skewness, Kurtosis = kurtosis,
                     a = a, b = b, Method = c("median_boot", "median_sobel"),
                     ab_boot = c(ab_boot[1], ab_sobel[1]),
                     ab_data = c(ab_boot[2], ab_sobel[2]),
                     reject = c(reject_boot, reject_sobel),
                     stringsAsFactors = FALSE)
        }, error = function(condition) NULL)

        # regression with selection of error distribution
        df_select <- NULL
        # df_select <- tryCatch({
        #   # fit mediation model
        #   select_fit <- fit_mediation(simulated_data, "x", "y", "m",
        #                               method = "regression", robust = FALSE,
        #                               family = "select", total_model = FALSE)
        #   # perform test for indirect effect
        #   select_boot <- test_mediation(select_fit, test = "boot", level = level,
        #                                 indices = indices)
        #   select_sobel <- test_mediation(select_fit, test = "sobel")
        #   # extract indirect effect
        #   ab_boot <- c(coef(select_boot, parm = "ab", type = "boot"),
        #                coef(select_boot, parm = "ab", type = "data"),
        #                use.names = FALSE)
        #   ab_sobel <- c(NA_real_, coef(select_sobel, parm = "ab"),
        #                 use.names = FALSE)
        #   # check for significant indirect effect
        #   reject_boot <- prod(select_boot$ci) > 0
        #   reject_sobel <- select_sobel$p_value < 1 - level
        #   # combine into data frame
        #   data.frame(n = n, Run = k, Skewness = skewness, Kurtosis = kurtosis,
        #              a = a, b = b, Method = c("select_boot", "select_sobel"),
        #              ab_boot = c(ab_boot[1], ab_sobel[1]),
        #              ab_data = c(ab_boot[2], ab_sobel[2]),
        #              reject = c(reject_boot, reject_sobel),
        #              stringsAsFactors = FALSE)
        # }, error = function(condition) NULL)

        # robust regression
        df_robust <- tryCatch({
          # fit mediation model
          robust_fit <- fit_mediation(simulated_data, "x", "y", "m",
                                      method = "regression", robust = "MM",
                                      control = lmrob_control)
          # perform test for indirect effect
          robust_boot <- test_mediation(robust_fit, test = "boot", level = level,
                                        indices = indices)
          robust_sobel <- test_mediation(robust_fit, test = "sobel")
          # extract indirect effect
          ab_boot <- c(coef(robust_boot, parm = "ab", type = "boot"),
                       coef(robust_boot, parm = "ab", type = "data"),
                       use.names = FALSE)
          ab_sobel <- c(NA_real_, coef(robust_sobel, parm = "ab"),
                        use.names = FALSE)
          # check for significant indirect effect
          reject_boot <- prod(robust_boot$ci) > 0
          reject_sobel <- robust_sobel$p_value < 1 - level
          # combine into data frame
          data.frame(n = n, Run = k, Skewness = skewness, Kurtosis = kurtosis,
                     a = a, b = b, Method = c("robust_boot", "robust_sobel"),
                     ab_boot = c(ab_boot[1], ab_sobel[1]),
                     ab_data = c(ab_boot[2], ab_sobel[2]),
                     reject = c(reject_boot, reject_sobel),
                     stringsAsFactors = FALSE)
        }, error = function(condition) NULL)

        ## results for current parameter settings
        rbind(df_ols, df_bcn, df_huber, df_median, df_select, df_robust)

      }, a = a, b = b, c = c, sigma_m = sigma_m, sigma_y = sigma_y,
      SIMPLIFY = FALSE)

      ## combine results for current error distribution into data frame
      do.call(rbind, results_nke)

    }, skewness = parameters$Skewness, kurtosis = parameters$Kurtosis,
    c0 = parameters$c0, c1 = parameters$c1, c2 = parameters$c2,
    c3 = parameters$c3, SIMPLIFY = FALSE)

    ## combine results for current simulation run into data frame
    do.call(rbind, results_nk)

  })

  ## combine results for current sample size into data frame
  do.call(rbind, results_n)

})

## combine results into data frame
results <- do.call(rbind, results_list)
cat(paste(Sys.time(), ": finished.\n"))

## store results
file <- "simulations/results/results_Fleishman.RData"
save(results, sample_sizes, a, b, c, level, file = file)
