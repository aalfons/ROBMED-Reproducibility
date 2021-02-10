# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## load required packages
library("dplyr")
library("stringr")
library("tidyr")

## information and methods to include in tables
parameters <- c("n", "epsilon", "d")
methods <- c("ols_boot", "ols_sobel", "bcn_boot",
             "huber_boot", "median_boot", "robust_boot")

## nice labels to be used in tables
parameter_labels <- c(n = "n", epsilon = "Probability of outliers",
                      d = "Outlier shift d")
method_labels <- c(ols_boot = "OLS bootstrap",
                   ols_sobel = "OLS Sobel",
                   bcn_boot = "Box-Cox bootstrap",
                   huber_boot = "Winsorized bootstap",
                   median_boot = "Median bootstrap",
                   robust_boot = "ROBMED")

## file containing results
file_results <- "simulations/results/results_outliers.RData"
load(file_results)

## filter results and prepare point estimates for indirect effect (bootstrap
## estimates for bootstrap tests, estimates on original sample otherwise)
results <- results %>%
  filter(Method %in% methods) %>%
  mutate(ab = if_else(str_detect(Method, "boot"), ab_boot, ab_data),
         bias = ab - a*b)


## loop over different mediation scenarios
scenarios <- results %>% distinct(a, b)
for (i in 1:nrow(scenarios)) {

  # filter results for current mediation scenario
  results_ab <- results %>% filter(a == scenarios$a[i], b == scenarios$b[i])

  # evaluate estimates of the indirect effect
  df_ab <- results_ab %>%
    group_by(n, epsilon, d, Method) %>%
    summarize(Bias = mean(bias), SD = sd(ab))

  # different evaluation of tests in case of mediation and nonmediation
  if (scenarios$a[i] * scenarios$b[i] != 0) {
    # for mediation, compute rate of rejection with correct sign
    results_ab <- results_ab %>%
      mutate(correct = reject & sign(ab) == sign(a * b))
    # aggregate results
    df_reject <- results_ab %>%
      group_by(n, epsilon, d, Method) %>%
      summarize(Value = mean(correct))
  } else {
    # for nonmediation, compute rejection rate
    df_reject <- results_ab %>%
      group_by(n, epsilon, d, Method) %>%
      summarize(Value = mean(reject))
  }

  # further preparation of tables of results
  # evaluation of estimates
  df_ab <- df_ab %>%
    pivot_wider(names_from = Method, values_from = Bias:SD) %>%
    arrange(n, epsilon, d)
  method_vars <- sapply(methods,
                        function(m, prefix) paste(prefix, m, sep = "_"),
                        prefix = c("Bias", "SD"))
  df_ab <- df_ab[, c(parameters, method_vars)]
  names(df_ab) <- c(parameter_labels[parameters],
                    sapply(method_labels[methods], paste, c("Bias", "Std. dev.")))
  df_ab[, -(1:3)] <- lapply(df_ab[, -(1:3)], function(x) sprintf("%.3f", x))
  # evaluation of tests
  df_reject <- df_reject %>%
    pivot_wider(names_from = Method, values_from = Value) %>%
    arrange(n, epsilon, d)
  df_reject <- df_reject[, c(parameters, methods)]
  names(df_reject) <- c(parameter_labels, method_labels)[names(df_reject)]
  df_reject[, -(1:3)] <- lapply(df_reject[, -(1:3)],
                                function(x) sprintf("%.3f", x))

  # file name for .csv file with table
  if (scenarios$b[i] == 0) {
    file_table <- "simulations/tables/table_outliers_a=%.1f_b=%d_%s.csv"
  } else {
    file_table <- "simulations/tables/table_outliers_a=%.1f_b=%.1f_%s.csv"
  }
  # write tables to .csv file
  write.csv(df_ab,
            file = sprintf(file_table, scenarios$a[i], scenarios$b[i], "ab"),
            quote = TRUE, na = "", row.names = FALSE)
  write.csv(df_reject,
            file = sprintf(file_table, scenarios$a[i], scenarios$b[i], "reject"),
            quote = FALSE, na = "", row.names = FALSE)

}
