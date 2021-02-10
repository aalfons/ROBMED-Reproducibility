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

## colors and plot symbols
colors <- c(rep("#F8766D", 2), "#B79F00", "#F564E3", "#619CFF", "#00BFC4")
symbols <- c(21:22, rep(21, 5))

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
file_results <- "simulations/results/results_lognormal.RData"
load(file_results)

## filter results and prepare point estimates for indirect effect (bootstrap
## estimates for bootstrap tests, estimates on original sample otherwise)
results <- results %>%
  filter(Method %in% methods) %>%
  mutate(ab = if_else(str_detect(Method, "boot"), ab_boot, ab_data))


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
    df_reject <- results_ab %>%
      group_by(n, a, b, Method) %>%
      summarize(Value = mean(correct)) %>%
      mutate(Parameter = reject)
    # make the longer facet label a little smaller
    size_reject <- 11
  } else {
    reject <- "reject"
    # aggregate results
    df_reject <- results_ab %>%
      group_by(n, a, b, Method) %>%
      summarize(Value = mean(reject)) %>%
      mutate(Parameter = reject)
    # the short facet label can be the usual size
    size_reject <- 12
  }

  # boxplots of the point estimates are plotted in the top row
  df_ab <- results_ab %>%
    dplyr::select(n, a, b, Method, Value = ab) %>%
    mutate(Parameter = "ab") %>%
    group_by(n, a, b, Method, Parameter) %>%
    summarize(Min = boxplot.stats(Value)$stats[1],
              Lower = boxplot.stats(Value)$stats[2],
              Middle = boxplot.stats(Value)$stats[3],
              Upper = boxplot.stats(Value)$stats[4],
              Max = boxplot.stats(Value)$stats[5])

  # nice labels for plotting
  df_ab <- df_ab %>% ungroup %>%
    mutate(n = paste("paste(italic(n),", sprintf("\" = %d\")", n)),
           Method = factor(method_labels[Method], levels = method_labels),
           Parameter = factor(facet_labels[Parameter], levels = facet_labels))
  df_reject <- df_reject %>% ungroup %>%
    mutate(n = paste("paste(italic(n),", sprintf("\" = %d\")", n)),
           Method = factor(method_labels[Method], levels = method_labels),
           Parameter = factor(facet_labels[Parameter], levels = facet_labels))

  # data frame for adding true value
  df_true <- data.frame(Parameter = facet_labels["ab"],
                        Value = scenarios$a[i] * scenarios$b[i])
  if (scenarios$a[i] * scenarios$b[i] == 0) {
    df_alpha <- data.frame(Parameter = facet_labels[reject],
                           Value = 1 - level)
    df_true <- rbind(df_true, df_alpha)
  }

  # expand axis limits of power facet with blank layer
  df_ylim <- data.frame(Parameter = unname(facet_labels[reject]),
                        Method = unname(method_labels[1]), Value = 0:1)

  # plot results
  p <- ggplot() +
    geom_boxplot(mapping = aes(x = Method, ymin = Min, lower = Lower,
                               middle = Middle, upper = Upper, ymax = Max,
                               fill = Method),
                 data = df_ab, stat = "identity", show.legend = FALSE) +
    geom_hline(aes(yintercept = Value), data = df_true) +
    geom_blank(mapping = aes(x = Method, y = Value), data = df_ylim) +
    geom_point(mapping = aes(x = Method, y = Value,
                             shape = Method, fill = Method),
               data = df_reject, size = 3, show.legend = FALSE) +
    facet_grid(Parameter ~ n, scales = "free_y", labeller = "label_parsed") +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = symbols) +
    labs(x = NULL, y = NULL) + theme_bw() +
    theme(axis.text.x = element_text(size = 12, angle = 270,
                                     hjust = 0, vjust = 0.5),
          axis.text.y = element_text(size = 12),
          # panel.spacing.x = unit(0.8, "line"),
          # panel.spacing.y = unit(0.7, "line"),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = size_reject))

  # save plot to file
  if (scenarios$b[i] == 0) {
    file_plot <- "simulations/figures/figure_lognormal_a=%.1f_b=%d.pdf"
  } else {
    file_plot <- "simulations/figures/figure_lognormal_a=%.1f_b=%.1f.pdf"
  }
  # in color
  pdf(file = sprintf(file_plot, scenarios$a[i], scenarios$b[i]),
      width = 8.5, height = 6.75)
  print(p)
  dev.off()

}
