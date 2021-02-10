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

## function to convert colors to grayscale
col2gray <- function (color, method = c("luminosity", "average")) {
  # initializations
  method <- match.arg(method)
  # convert to RGB values
  rgb <- col2rgb(color)
  # compute weighted combination of RGB values
  weights <- switch(method, luminosity = c(0.2126, 0.7152, 0.0722),
                    average = rep.int(1/3, 3))
  gray <- crossprod(weights, rgb)
  # convert to text string
  rgb(gray, gray, gray, maxColorValue = 255)
}

## colors and plot symbols
colors <- c(rep("#F8766D", 2), "#B79F00", "#00BA38",
            "#F564E3", "#619CFF", "#00BFC4")
symbols <- c(21:22, rep(21, 5))

## methods to include in paper
methods <- c("ols_boot", "ols_sobel", "bcn_boot", "select_boot",
             "huber_boot", "median_boot", "robust_boot")

## nice labels to be used for plots
method_labels <- c(ols_boot = "OLS bootstrap",
                   ols_sobel = "OLS Sobel",
                   bcn_boot = "Box-Cox bootstrap",
                   select_boot = "SNT bootstrap",
                   huber_boot = "Winsorized bootstap",
                   median_boot = "Median bootstrap",
                   robust_boot = "ROBMED")
facet_labels <- c(ab = "paste(\"Indirect effect \", italic(ab))",
                  correct = "\"Rate of rejection with correct sign\"",
                  reject = "\"Rejection rate\"")
setting_labels <- c(normal = "\"Normal\"",
                    outliers = "\"Outliers\"",
                    skewnormal = "\"Skewness\"",
                    t = "\"Heavy tails\"")

## file containing results
file_results <- "simulations/results/results_ORM.RData"
load(file_results)

## filter results and prepare point estimates for indirect effect (bootstrap
## estimates for bootstrap tests, estimates on original sample otherwise)
results <- results %>%
  filter(Method %in% methods) %>%
  mutate(ab = if_else(str_detect(Method, "boot"), ab_boot, ab_data))


## loop over different mediation scenarios
scenarios <- results %>% distinct(a, b)
for (i in 1:nrow(scenarios)) {

  # filter results for current mediation setting
  results_ab <- results %>% filter(a == scenarios$a[i], b == scenarios$b[i])

  # slightly different plots for mediation and nonmediation
  # rejection rate (with correct sign) is plotted in the bottom row
  if (scenarios$a[i] * scenarios$b[i] != 0) {
    reject <- "correct"
    # for mediation, compute rate of rejection with correct sign
    results_ab <- results_ab %>%
      mutate(correct = reject & sign(ab) == sign(a * b))
    # aggregate results
    df_reject <- results_ab %>%
      group_by(a, b, Setting, Method) %>%
      summarize(Value = mean(correct)) %>%
      mutate(Parameter = reject)
    # make the longer facet label a little smaller
    size_reject <- 11
  } else {
    reject <- "reject"
    # aggregate results
    df_reject <- results_ab %>%
      group_by(a, b, Setting, Method) %>%
      summarize(Value = mean(reject)) %>%
      mutate(Parameter = reject)
    # the short facet label can be the usual size
    size_reject <- 12
  }

  # boxplots of the point estimates are plotted in the top row
  df_ab <- results_ab %>%
    dplyr::select(a, b, Setting, Method, Value = ab) %>%
    mutate(Parameter = "ab") %>%
    group_by(a, b, Setting, Method, Parameter) %>%
    summarize(Min = boxplot.stats(Value)$stats[1],
              Lower = boxplot.stats(Value)$stats[2],
              Middle = boxplot.stats(Value)$stats[3],
              Upper = boxplot.stats(Value)$stats[4],
              Max = boxplot.stats(Value)$stats[5])

  # nice labels for plotting
  df_ab <- df_ab %>% ungroup %>%
    mutate(Method = factor(method_labels[Method], levels = method_labels),
           Parameter = factor(facet_labels[Parameter], levels = facet_labels),
           Setting = factor(setting_labels[Setting], levels = setting_labels))
  df_reject <- df_reject %>% ungroup %>%
    mutate(Method = factor(method_labels[Method], levels = method_labels),
           Parameter = factor(facet_labels[Parameter], levels = facet_labels),
           Setting = factor(setting_labels[Setting], levels = setting_labels))

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
    facet_grid(Parameter ~ Setting, scales = "free_y",
               labeller = "label_parsed") +
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
    file_plot_color <- "simulations/figures/figure_ORM_a=%.1f_b=%d_color.pdf"
    file_plot_bw <- "simulations/figures/figure_ORM_a=%.1f_b=%d_bw.pdf"
  } else {
    file_plot_color <- "simulations/figures/figure_ORM_a=%.1f_b=%.1f_color.pdf"
    file_plot_bw <- "simulations/figures/figure_ORM_a=%.1f_b=%.1f_bw.pdf"
  }
  # in color
  pdf(file = sprintf(file_plot_color, scenarios$a[i], scenarios$b[i]),
      width = 8.5, height = 6.75)
  print(p + scale_fill_manual(values = colors))
  dev.off()
  # in grayscale
  pdf(file = sprintf(file_plot_bw, scenarios$a[i], scenarios$b[i]),
      width = 8.5, height = 6.75)
  print(p + scale_fill_manual(values = col2gray(colors)))
  dev.off()

}
