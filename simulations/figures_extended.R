# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## load required packages
library("dplyr")
library("ggplot2")
library("grid")
library("scales")
library("stringr")
library("tidyr")

## colors and line types
colors <- c(rep("#F8766D", 2), "#B79F00", "#F564E3", "#619CFF", "#00BFC4")
line_types <- c("dashed", "dotted", "longdash", "twodash", "dotdash", "solid")
line_sizes <- c(2/3, 0.75, rep(2/3, 4))  # dotted line ever so slightly thicker

## methods to include in plot
methods <- c("ols_boot", "ols_sobel", "bcn_boot",
             "huber_boot", "median_boot", "robust_boot")

## nice labels to be used for plots
method_labels <- c(ols_boot = "OLS bootstrap",
                   ols_sobel = "OLS Sobel",
                   bcn_boot = "Box-Cox bootstrap",
                   huber_boot = "Winsorized bootstap",
                   median_boot = "Median bootstrap",
                   robust_boot = "ROBMED")
facet_labels <- c(bias = "\"Bias\"",
                  relative = "\"Relative bias\"",
                  correct = "\"Rate of rejection with correct sign\"",
                  reject = "\"Rejection rate\"")
setting_labels <- c(normal = "\"Normal\"",
                    outliers = "\"Outliers\"",
                    skewnormal = "\"Skewness\"",
                    t = "\"Heavy tails\"")

## file containing results
file_results <- "simulations/results/results_extended.RData"
load(file_results)

## filter results and prepare point estimates for indirect effect (bootstrap
## estimates for bootstrap tests, estimates on original sample otherwise)
results <- results %>%
  filter(Method %in% methods) %>%
  mutate(Scenario = if_else(b == 0, "nonmediation", "mediation"),
         ab = if_else(str_detect(Method, "boot"), ab_boot, ab_data),
         bias = ab - a*b)


## loop over different mediation scenarios
scenarios <- results %>% distinct(n, Scenario)
for (i in 1:nrow(scenarios)) {

  # slightly different plots for mediation and nonmediation
  # rejection rate (with correct sign) is plotted in the bottom row
  if (scenarios$Scenario[i] == "mediation") {
    bias <- "bias"
    reject <- "correct"
    # for mediation, compute rate of rejection with correct sign
    results_ab <- results %>%
      filter(n == scenarios$n[i], b != 0) %>%
      mutate(correct = reject & sign(ab) == sign(a * b))
    # aggregate results
    df_results <- results_ab %>%
      group_by(a, b, Setting, Method) %>%
      summarize(bias = mean(bias, na.rm = TRUE),
                correct = mean(correct, na.rm = TRUE)) %>%
      gather(key = "Parameter", value = "Value", bias, correct)
    # make the longer facet label a little smaller
    size_reject <- 11
  } else {
    bias <- "relative"
    reject <- "reject"
    # for nonmediation, compute relative bias with respect to varying
    # coefficient a
    results_ab <- results %>%
      filter(n == scenarios$n[i], b == 0) %>%
      mutate(relative = bias / a)
    # aggregate results
    df_results <- results_ab %>%
      group_by(a, b, Setting, Method) %>%
      summarize(relative = mean(relative, na.rm = TRUE),
                reject = mean(reject, na.rm = TRUE)) %>%
      gather(key = "Parameter", value = "Value", relative, reject)
    # the short facet label can be the usual size
    size_reject <- 12
  }

  # nice labels for plotting
  df_results <- df_results %>% ungroup %>%
    mutate(Method = factor(method_labels[Method], levels = method_labels),
           Parameter = factor(facet_labels[Parameter], levels = facet_labels),
           Setting = factor(setting_labels[Setting], levels = setting_labels))

  # data frame for adding true value
  df_true <- data.frame(Parameter = facet_labels[bias], Value = 0)
  if (scenarios$Scenario[i] == "nonmediation") {
    df_alpha <- data.frame(Parameter = facet_labels[reject], Value = 1 - level)
    df_true <- rbind(df_true, df_alpha)
  }

  # expand axis limits of power facet with blank layer
  df_ylim <- data.frame(Parameter = unname(facet_labels[reject]),
                        a = 0.4, Value = 0:1)

  # x-axis label
  if (scenarios$Scenario[i] == "mediation") {
    x_lab <- expression(paste("Effect size of ", italic(a), " and ", italic(b)))
  } else {
    x_lab <- expression(paste("Effect size of ", italic(a), " (",
                              italic(b == 0), ")"))
  }

  # plot results
  mapping <- aes(x = a, y = Value, color = Method,
                 linetype = Method, size = Method)
  p <- ggplot() +
    geom_hline(aes(yintercept = Value), data = df_true) +
    geom_blank(mapping = aes(x = a, y = Value), data = df_ylim) +
    geom_line(mapping = mapping, data = df_results) +
    facet_grid(Parameter ~ Setting, scales = "free_y",
               labeller = "label_parsed") +
    labs(x = x_lab, y = NULL) +
    scale_color_manual("", values = colors) +
    scale_linetype_manual("", values = line_types) +
    scale_size_manual("", values = line_sizes) +
    theme_bw() + theme(axis.text = element_text(size = 12),
                       axis.title = element_text(size = 13),
                       legend.position = "top",
                       legend.direction = "horizontal",
                       legend.key = element_rect(color = NA),
                       legend.key.width = unit(4.9, "line"),
                       legend.text = element_text(size = 12),
                       panel.spacing.x = unit(0.6, "line"),
                       panel.spacing.y = unit(0.5, "line"),
                       strip.text.x = element_text(size = 12),
                       strip.text.y = element_text(size = size_reject))

  # save plot to file (in color)
  file_plot <- "simulations/figures/figure_extended_%s_n=%d.pdf"
  pdf(file = sprintf(file_plot, scenarios$Scenario[i], scenarios$n[i]),
      width = 8.5, height = 6.75)
  print(p)
  dev.off()

}


## closer look for outliers under increasing sample size
results_out <- results %>%
  filter(Setting == "outliers")

## loop over different mediation scenarios
scenarios <- results_out %>% distinct(Scenario)
for (i in 1:nrow(scenarios)) {

  # slightly different plots for mediation and nonmediation
  # rejection rate (with correct sign) is plotted in the bottom row
  if (scenarios$Scenario[i] == "mediation") {
    bias <- "bias"
    reject <- "correct"
    # for mediation, compute rate of rejection with correct sign
    results_ab <- results_out %>%
      filter(b != 0) %>%
      mutate(correct = reject & sign(ab) == sign(a * b))
    # aggregate results
    df_results <- results_ab %>%
      group_by(n, a, b, Method) %>%
      summarize(bias = mean(bias, na.rm = TRUE),
                correct = mean(correct, na.rm = TRUE)) %>%
      gather(key = "Parameter", value = "Value", bias, correct)
  } else {
    bias <- "relative"
    reject <- "reject"
    # for nonmediation, compute relative bias with respect to varying
    # coefficient a
    results_ab <- results_out %>%
      filter(b == 0) %>%
      mutate(relative = bias / a)
    # aggregate results
    df_results <- results_ab %>%
      group_by(n, a, b, Method) %>%
      summarize(relative = mean(relative, na.rm = TRUE),
                reject = mean(reject, na.rm = TRUE)) %>%
      gather(key = "Parameter", value = "Value", relative, reject)
  }

  # nice labels for plotting
  sample_sizes <- unique(df_results$n)
  if (scenarios$Scenario[i] == "mediation") {
    df_results <- df_results %>% ungroup %>%
      mutate(Effect = paste0("\"a = b = ", a, "\""))
  } else {
    df_results <- df_results %>% ungroup %>%
      mutate(Effect = paste0("\"a = ", a, ", b = ", b, "\""))
  }
  df_results <- df_results %>%
    mutate(Method = factor(method_labels[Method], levels = method_labels),
           Parameter = factor(facet_labels[Parameter], levels = facet_labels))

  # data frame for adding true value
  df_true <- data.frame(Parameter = facet_labels[bias], Value = 0)
  if (scenarios$Scenario[i] == "nonmediation") {
    df_alpha <- data.frame(Parameter = facet_labels[reject], Value = 1 - level)
    df_true <- rbind(df_true, df_alpha)
  }

  # expand axis limits of power facet with blank layer
  df_ylim <- data.frame(Parameter = unname(facet_labels[reject]),
                        n = 100, Value = 0:1)

  # plot results
  mapping <- aes(x = n, y = Value, color = Method,
                 linetype = Method, size = Method)
  text_pos <- c(1, 0, 0, 0.5, 0.5)
  p_out <- ggplot() +
    geom_hline(aes(yintercept = Value), data = df_true) +
    geom_blank(mapping = aes(x = n, y = Value), data = df_ylim) +
    geom_line(mapping = mapping, data = df_results) +
    facet_grid(Parameter ~ Effect, scales = "free_y",
               labeller = "label_parsed") +
    labs(x = "Sample size", y = NULL) +
    scale_x_continuous(breaks = sample_sizes) +
    scale_color_manual("", values = colors) +
    scale_linetype_manual("", values = line_types) +
    scale_size_manual("", values = line_sizes) +
    theme_bw() + theme(axis.text.x = element_text(size = 11, angle = 270,
                                                  hjust = 0, vjust = text_pos),
                       axis.text.y = element_text(size = 11),
                       axis.title = element_text(size = 13),
                       legend.position = "top",
                       legend.direction = "horizontal",
                       legend.key = element_rect(color = NA),
                       legend.key.width = unit(4.9, "line"),
                       legend.text = element_text(size = 12),
                       panel.grid.minor.x = element_blank(),
                       panel.spacing.x = unit(0.6, "line"),
                       panel.spacing.y = unit(0.5, "line"),
                       strip.text.x = element_text(size = 11),
                       strip.text.y = element_text(size = 11))

  # save plot to file (in color)
  file_plot <- "simulations/figures/figure_extended_outliers_%s.pdf"
  pdf(file = sprintf(file_plot, scenarios$Scenario[i]),
      width = 8.5, height = 6.75)
  print(p_out)
  dev.off()

}
