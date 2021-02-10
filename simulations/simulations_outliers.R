# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## load required packages
library("car")
library("robmed")


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

## control parameters for outliers
shift <- c(x = 0, y = 1, m = -1)          # base outlier shift configuration
multiplier <- c(x = 1, y = 0.1, m = 0.1)  # multiplication factor for outliers
# full range of outlier settings
outlier_settings <- rbind(
  data.frame(epsilon = 0, d = 0),
  expand.grid(epsilon = seq(0.01, 0.04, by = 0.01),
              d = seq(1, 6, by = 1))
)

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
    e_m <- rnorm(n)
    e_y <- rnorm(n)

    ## probability of being an outlier for each observation
    ## (only used in setting with outliers)
    p_out <- runif(n)

    ## ensure that the same bootstrap samples are used for all parameter
    ## settings and all bootstrap tests for maximum comparability
    indices <- boot_samples(n, R = R)

    ## loop over different parameter settings
    results_nk <- mapply(function(a, b, c, sigma_m, sigma_y) {

      # generate hypothesized mediator and dependent variable
      m <- a * x + sigma_m * e_m
      y <- b * m + c * x + sigma_y * e_y
      clean_data <- data.frame(x, y, m)

      # loop over outlier settings
      results_nkl <- mapply(function(epsilon, d) {

        # replace selected observations by outliers except in first iteration
        simulated_data <- clean_data
        if (epsilon > 0 && d > 0) {
          replace <- p_out < epsilon
          simulated_data[replace, ] <- mapply(function(v, m, s) m * v + s,
                                              v = simulated_data[replace, ],
                                              m = multiplier, s = d * shift,
                                              SIMPLIFY = FALSE)
        }

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
          data.frame(n = n, Run = k, a = a, b = b, epsilon = epsilon,
                     d = d, Method = c("ols_boot", "ols_sobel"),
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
          data.frame(n = n, Run = k, a = a, b = b, epsilon = epsilon,
                     d = d, Method = c("bcn_boot", "bcn_sobel"),
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
          data.frame(n = n, Run = k, a = a, b = b, epsilon = epsilon,
                     d = d, Method = c("huber_boot", "huber_sobel"),
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
          data.frame(n = n, Run = k, a = a, b = b, epsilon = epsilon,
                     d = d, Method = c("median_boot", "median_sobel"),
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
        #   data.frame(n = n, Run = k, a = a, b = b, epsilon = epsilon,
        #              d = d, Method = c("select_boot", "select_sobel"),
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
          data.frame(n = n, Run = k, a = a, b = b, epsilon = epsilon,
                     d = d, Method = c("robust_boot", "robust_sobel"),
                     ab_boot = c(ab_boot[1], ab_sobel[1]),
                     ab_data = c(ab_boot[2], ab_sobel[2]),
                     reject = c(reject_boot, reject_sobel),
                     stringsAsFactors = FALSE)
        }, error = function(condition) NULL)

        ## results for current parameter settings
        rbind(df_ols, df_bcn, df_huber, df_median, df_select, df_robust)

      }, epsilon = outlier_settings$epsilon, d = outlier_settings$d,
      SIMPLIFY = FALSE)

      # combine results for different outlier settings into data frame
      do.call(rbind, results_nkl)

    }, a = a, b = b, c = c, sigma_m = sigma_m, sigma_y = sigma_y,
    SIMPLIFY = FALSE)

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
file <- "simulations/results/results_outliers.RData"
save(results, sample_sizes, a, b, c, level, file = file)
