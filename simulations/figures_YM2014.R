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

# additional utility function
as_index0 <- function(x, levels) {
  if (missing(levels)) f <- as.factor(x)
  else f <- factor(x, levels = levels)
  as.numeric(f) - 1
}


## colors and plot symbols
colors <- c(rep("#F8766D", 2), "#B79F00", "#F564E3", "#619CFF", "#00BFC4")
symbols <- c(21:22, rep(21, 4))

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
setting_labels <- c(normal = "\"Standard normal\"",
                    # t = "paste(italic(t), \"(2)\")",
                    t = "\"Heavy tails\"",
                    mixture = "\"Contaminated normal\"")

## file containing results
file_results <- "simulations/results/results_YM2014.RData"
load(file_results)

## filter results and prepare point estimates for indirect effect (bootstrap
## estimates for bootstrap tests, estimates on original sample otherwise)
results <- results %>%
  filter(Method %in% methods) %>%
  mutate(Scenario = if_else(a == 0, "nonmediation", "mediation"),
         ab = if_else(str_detect(Method, "boot"), ab_boot, ab_data))

## determine positions on the x-axis and x-axis labels
n_methods <- length(methods)
box_width <- 1
box_margin <- 0.3
offset <- 2
width <- n_methods * box_width + (n_methods - 1) * box_margin + offset
step <- box_width + box_margin
groups <- results %>%
  distinct(b, Method) %>%
  mutate(x = as_index0(b) * width + as_index0(Method, levels = methods) * step + box_width / 2)
axis_labels <- groups %>%
  group_by(b) %>%
  summarize(x = mean(x)) %>%
  rename(Label = b)


## loop over different mediation scenarios
scenarios <- results %>% distinct(n, Scenario)
for (i in 1:nrow(scenarios)) {

  # slightly different plots for mediation and nonmediation
  # rejection rate (with correct sign) is plotted in the bottom row
  if (scenarios$Scenario[i] == "mediation") {
    reject <- "correct"
    # for mediation, compute rate of rejection with correct sign
    results_ab <- results %>%
      filter(n == scenarios$n[i], a != 0) %>%
      mutate(correct = reject & sign(ab) == sign(a * b))
    # aggregate results
    df_reject <- results_ab %>%
      group_by(a, b, Setting, Method) %>%
      summarize(Value = mean(correct)) %>%
      mutate(Parameter = reject)
  } else {
    reject <- "reject"
    # filter results
    results_ab <- results %>%
      filter(n == scenarios$n[i], a == 0)
    # aggregate results
    df_reject <- results %>%
      filter(n == scenarios$n[i], a == 0) %>%
      group_by(a, b, Setting, Method) %>%
      summarize(Value = mean(reject)) %>%
      mutate(Parameter = reject)
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

  # determine positions on the x-axis
  df_ab <- df_ab %>% left_join(groups, by = c("b", "Method"))
  df_reject <- df_reject %>% left_join(groups, by = c("b", "Method"))

  # nice labels for plotting
  df_ab <- df_ab %>% ungroup %>%
    mutate(Method = factor(method_labels[Method], levels = method_labels),
           Parameter = factor(facet_labels[Parameter], levels = facet_labels),
           Setting = factor(setting_labels[Setting], levels = setting_labels))
  df_reject <- df_reject %>% ungroup %>%
    mutate(Method = factor(method_labels[Method], levels = method_labels),
           Parameter = factor(facet_labels[Parameter], levels = facet_labels),
           Setting = factor(setting_labels[Setting], levels = setting_labels))

  # add lines for subpanels given by different effect sizes
  n_effects <- nrow(axis_labels)
  df_subpanel <- data.frame(x = (axis_labels$x[-1]+axis_labels$x[-n_effects])/2)

  # add reference lines for true values
  if (scenarios$Scenario[i] == "mediation") {
    start <- c(-offset/2, df_subpanel$x)
    end <- c(df_subpanel$x, max(groups$x) + box_width/2 + offset / 2)
    true_ab <- df_ab %>% ungroup() %>% distinct(a, b)
    df_true <- data.frame(Parameter = unname(facet_labels["ab"]),
                          x = start, xend = end, y = true_ab$a * true_ab$b)
    expand <- 0
  } else {
    df_true <- data.frame(Parameter = facet_labels["ab"], Value = 0)
    df_size <- data.frame(Parameter = facet_labels[reject], Value = 1 - level)
    expand <- offset / 2
  }

  # expand axis limits of power facet with blank layer
  df_ylim <- data.frame(Parameter = unname(facet_labels[reject]),
                        Method = unname(method_labels[1]),
                        x = axis_labels$x[1], Value = 0:1)

  # x-axis label
  if (scenarios$Scenario[i] == "mediation") {
    x_lab <- expression(paste("Effect size of ", italic(a), " and ", italic(b)))
  } else {
    x_lab <- expression(paste("Effect size of ", italic(b), " (",
                              italic(a == 0), ")"))
  }

  # plot results
  p <- ggplot() +
    geom_boxplot(mapping = aes(x = x, ymin = Min, lower = Lower,
                               middle = Middle, upper = Upper, ymax = Max,
                               fill = Method, group = x),
                 data = df_ab, stat = "identity", width = box_width,
                 show.legend = FALSE)
  if (scenarios$Scenario[i] == "mediation") {
    p <- p +
      geom_segment(aes(x = x, y = y, xend = xend, yend = y), data = df_true)
  } else {
    p <- p +
      geom_hline(aes(yintercept = Value), data = df_true) +
      geom_hline(aes(yintercept = Value), data = df_size)
  }
  p <- p +
    geom_vline(aes(xintercept = x), data = df_subpanel, color = "darkgray") +
    geom_blank(mapping = aes(x = x, y = Value), data = df_ylim) +
    geom_point(mapping = aes(x = x, y = Value,
                             shape = Method, fill = Method),
               data = df_reject, size = 2.75) +
    facet_grid(Parameter ~ Setting, scales = "free_y",
               labeller = "label_parsed") +
    scale_x_continuous(breaks = axis_labels$x,
                       minor_breaks = groups$x,
                       labels = axis_labels$Label,
                       expand = expansion(add = expand)) +
    # coord_cartesian(xlim = c(11, 38)) +
    scale_shape_manual("", values = symbols) +
    scale_fill_manual("", values = colors) +
    labs(x = x_lab, y = NULL) + theme_bw() +
    theme(axis.text = element_text(size = 12),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.key = element_rect(color = NA),
          legend.key.width = unit(4.9, "line"),
          legend.text = element_text(size = 12),
          panel.grid.major.x = element_blank(),
          # panel.spacing.x = unit(0.8, "line"),
          # panel.spacing.y = unit(0.7, "line"),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12))

  # save plot to file
  file_plot <- "simulations/figures/figure_YM2014_%s_n=%d.pdf"
  pdf(file = sprintf(file_plot, scenarios$Scenario[i], scenarios$n[i]),
      width = 8.5, height = 6.75)
  print(p)
  dev.off()

}
