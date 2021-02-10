# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## load required packages
library("car")
library("parallel")
library("robmed")
library("sn")


## control parameters for simulation
seed <- 20200309                            # seed for random number generator
sample_sizes <- c(50, 100, 250, 500, 1000)  # number of observations
K <- 1000                                   # number of simulation runs
R <- 5000                                   # number of bootstrap samples
settings <- c("normal", "outliers", "skewnormal", "t")
# coefficients in mediation model
effect_sizes <- seq(0.2, 0.8, by = 0.1)
n_effects <- length(effect_sizes)
a <- rep(effect_sizes, times = 2)
b <- c(effect_sizes, rep(0, n_effects))
c <- rep(effect_sizes, times = 2)
# error scales in mediation model
sigma_m <- sqrt(1 - a^2)
sigma_y <- rep(1, times = 2 * n_effects)

## control parameters for outliers
epsilon <- 0.02
multiplier <- c(x = 1, y = 0.1, m = 0.1)
shift <- c(x = 0, y = 3, m = -3)

## control parameters for methods
level <- 0.95                                        # confidence level
lmrob_control <- reg_control(max_iterations = 5000)  # MM-estimator


## function to generate different error terms
## (for maximum comparability of the results, the different error terms
## are transformed from the normal errors that are generated first)
generate_errors <- function(n, alpha = Inf, df = 2) {
  # generate standard normal errors
  e_normal <- rnorm(n)
  # first apply the standard normal CDF to transform to uniformly distributed
  # values, then apply the inverse CDF of other error distributions
  u <- pnorm(e_normal)
  # gnerate skew-normal errors
  e_skewnormal <- qsn(u, alpha = alpha)
  # generate t distributed errors
  e_t <- qt(u, df = df)
  # combine the error terms as matrix
  cbind(normal = e_normal, outliers = e_normal,
        skewnormal = e_skewnormal, t = e_t)
}


## run simulation
cat(paste(Sys.time(), ": starting ...\n"))
results_list <- lapply(sample_sizes, function(n) {

  ## use the same seed for each sample size
  ## (then simulations can be run separately for different dimensions)
  set.seed(seed)

  ## print sample size
  cat(paste(Sys.time(), sprintf(": n = %d\n", n)))

  ## loop over simulation runs
  results_n <- lapply(1:K, function(k) {

    ## print simulation run
    cat(paste(Sys.time(), sprintf(":   run = %d\n", k)))

    ## generate independent variable and error terms such that results are
    ## comparable across different values of the parameters
    x <- rnorm(n)
    e_m <- generate_errors(x)
    e_y <- generate_errors(x)

    ## probability of being an outlier for each observation
    ## (only used in setting with outliers)
    p_out <- runif(n)

    ## ensure that the same bootstrap samples are used for all parameter
    ## settings and all bootstrap tests for maximum comparability
    indices <- boot_samples(n, R = R)

    ## loop over different settings
    results_nk <- lapply(settings, function(setting) {

      ## print simulation run
      cat(paste(Sys.time(), sprintf(":    setting = \"%s\"\n", setting)))

      ## loop over different parameter values
      results_nks <- mcmapply(function(a, b, c, sigma_m, sigma_y) {

        # ## print simulation run
        # cat(paste(Sys.time(), sprintf(":      ab = %.02f\n", a * b)))

        # generate hypothesized mediator and dependent variable
        m <- a * x + sigma_m * e_m[, setting]
        y <- b * m + c * x + sigma_y * e_y[, setting]
        simulated_data <- data.frame(x, y, m)

        # for setting with outliers, replace selected observations
        if (setting == "outliers") {
          # observations to be replaced
          replace <- p_out < epsilon
          # outlier shift needs to scale with standard deviations
          std_dev <- c(1, sqrt(b^2 + c^2 + 2*a*b*c + 1), 1)
          # generate outliers
          simulated_data[replace, ] <- mapply(function(v, m, s) m * v + s,
                                              v = simulated_data[replace, ],
                                              m = multiplier,
                                              s = shift * std_dev,
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
          data.frame(n = n, Run = k, Setting = setting, a = a, b = b,
                     Method = c("ols_boot", "ols_sobel"),
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
          data.frame(n = n, Run = k, Setting = setting, a = a, b = b,
                     Method = c("bcn_boot", "bcn_sobel"),
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
          data.frame(n = n, Run = k, Setting = setting, a = a, b = b,
                     Method = c("huber_boot", "huber_sobel"),
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
          data.frame(n = n, Run = k, Setting = setting, a = a, b = b,
                     Method = c("median_boot", "median_sobel"),
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
        #   data.frame(n = n, Run = k, Setting = setting, a = a, b = b,
        #              Method = c("select_boot", "select_sobel"),
        #              ab_boot = c(ab_boot[1], ab_sobel[1]),
        #              ab_data = c(ab_boot[2], ab_sobel[2]),
        #              reject = c(reject_boot, reject_sobel),
        #              stringsAsFactors = FALSE)
        # }, error = function(condition) NULL)

        # robust regression
        df_robust <- tryCatch({
          # fit mediation model
          robust_fit <- fit_mediation(simulated_data, "x", "y", "m",
                                      method = "regression", robust = "MM")
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
          data.frame(n = n, Run = k, Setting = setting, a = a, b = b,
                     Method = c("robust_boot", "robust_sobel"),
                     ab_boot = c(ab_boot[1], ab_sobel[1]),
                     ab_data = c(ab_boot[2], ab_sobel[2]),
                     reject = c(reject_boot, reject_sobel),
                     stringsAsFactors = FALSE)
        }, error = function(condition) NULL)

        ## results for current parameter settings
        rbind(df_ols, df_bcn, df_huber, df_median, df_select, df_robust)

      }, a = a, b = b, c = c, sigma_m = sigma_m, sigma_y = sigma_y,
      SIMPLIFY = FALSE, mc.preschedule = FALSE, mc.set.seed = FALSE,
      mc.cores = 2)

      ## 'mc.preschedule = FALSE' and 'mc.set.seed = FALSE' are set such that
      ## the seed of the random number generator is passed on from main thread
      ## to the worker threads.  The only place inside the parallelized loop
      ## where random numbers are needed is the subsampling algorithm of the MM
      ## estimator.  By setting these options, the same subsamples are used
      ## across parameter settings.  Note that the variability in the samples
      ## is guaranteed by performing the random draws outside the parallelized
      ## loop.

      ## combine results for current simulation run into data frame
      do.call(rbind, results_nks)

    })

    ## combine results for current error distribution into data frame
    do.call(rbind, results_nk)

  })

  ## combine results for current sample size into data frame
  do.call(rbind, results_n)

})

## combine results into data frame
results <- do.call(rbind, results_list)
cat(paste(Sys.time(), ": finished.\n"))

## store results
file <- "simulations/results/results_extended.RData"
save(results, sample_sizes, a, b, c, level, file = file)
