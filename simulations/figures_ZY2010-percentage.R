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
facet_labels <- c(ab = "paste(\"Indirect effect \", italic(ab))",
                  correct = "\"Rate of rejection with correct sign\"",
                  reject = "\"Rejection rate\"")

## file containing results
file_results <- "simulations/results/results_ZY2010.RData"
load(file_results)

## filter results and prepare point estimates for indirect effect (bootstrap
## estimates for bootstrap tests, estimates on original sample otherwise)
results <- results %>%
  filter(Method %in% methods, NrOutliers == 0 | Shift == 6) %>%
  mutate(Epsilon = NrOutliers / n,
         ab = if_else(str_detect(Method, "boot"), ab_boot, ab_data))

## loop over different mediation scenarios
scenarios <- results %>% distinct(a, b)
for (i in 1:nrow(scenarios)) {

  # filter results for current mediation scenario
  results_ab <- results %>% filter(a == scenarios$a[i], b == scenarios$b[i])

  # slightly different plots for mediation and nonmediation
  if (scenarios$a[i] * scenarios$b[i] != 0) {
    reject <- "correct"
    # for mediation, compute rate of rejection with correct sign
    results_ab <- results_ab %>%
      mutate(correct = reject & sign(ab) == sign(a * b))
    # aggregate results
    df_results <- results_ab %>%
      group_by(n, a, b, Epsilon, Method) %>%
      summarize(ab = mean(ab, na.rm = TRUE),
                correct = mean(correct, na.rm = TRUE)) %>%
      gather(key = "Parameter", value = "Value", ab, correct)
  } else {
    reject <- "reject"
    # aggregate results
    df_results <- results_ab %>%
      group_by(n, a, b, Epsilon, Method) %>%
      summarize(ab = mean(ab, na.rm = TRUE),
                reject = mean(reject, na.rm = TRUE)) %>%
      gather(key = "Parameter", value = "Value", ab, reject)
  }

  # nice labels for plotting
  df_results <- ungroup(df_results) %>%
    mutate(n = paste("paste(italic(n),", sprintf("\" = %d\")", n)),
           Method = factor(method_labels[Method], levels = method_labels),
           Parameter = factor(facet_labels[Parameter], levels = facet_labels))

  # data frame for adding true value
  df_true <- data.frame(Parameter = facet_labels["ab"],
                        Value = scenarios$a[i] * scenarios$b[i])
  if (scenarios$a[i] * scenarios$b[i] == 0) {
    df_alpha <- data.frame(Parameter = facet_labels["reject"],
                           Value = 1 - level)
    df_true <- rbind(df_true, df_alpha)
  }

  # expand axis limits if there is little variation with blank layer
  df_ab <- df_results %>% filter(Parameter == facet_labels["ab"])
  ab_ylim <- mean(df_ab$Value) + c(-0.01, 0.01)
  reject_ylim <- c(0, 0.1)
  parameters <- rep(c("ab", reject), each = 2)
  df_ylim <- data.frame(Parameter = factor(facet_labels[parameters],
                                           levels = facet_labels),
                        Epsilon = 0, Value = c(ab_ylim, reject_ylim))

  # plot results
  mapping <- aes(x = Epsilon, y = Value, color = Method,
                 linetype = Method, size = Method)
  p <- ggplot() +
    geom_blank(mapping = aes(x = Epsilon, y = Value), data = df_ylim) +
    geom_hline(aes(yintercept = Value), data = df_true) +
    geom_line(mapping = mapping, data = df_results) +
    facet_grid(Parameter ~ n, scales = "free_y", labeller = "label_parsed") +
    labs(x = "Percentage of outliers", y = NULL) +
    scale_x_continuous(labels = percent) +
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
                       # panel.spacing.x = unit(0.8, "line"),
                       # panel.spacing.y = unit(0.7, "line"),
                       strip.text = element_text(size = 12),
                       plot.margin = margin(0, 1, 0.5, 0.5, unit = "line"))

  # save plot to file
  if (scenarios$b[i] == 0) {
    file_plot <- "simulations/figures/figure_ZY2010-percentage_a=%.1f_b=%d.pdf"
  } else {
    file_plot <- "simulations/figures/figure_ZY2010-percentage_a=%.1f_b=%.1f.pdf"
  }
  pdf(file = sprintf(file_plot, scenarios$a[i], scenarios$b[i]),
      width = 8.5, height = 6.75)
  print(p)
  dev.off()

}
