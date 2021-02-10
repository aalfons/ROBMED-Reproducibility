# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# This script reproduces the illustration of (robust) mediation analysis on
# simulated data from Figure 1 of the paper:
#
# Alfons, A., Ates, N.Y., & Groenen, P.J.F. (2021). A Robust Bootstrap Test for
# Mediation Analysis. Organizational Research Methods, accepted for publication.
#
# Before running this script, please make sure that all necessary packages are
# installed.


## load packages
library("robmed")
library("colorspace")
library("viridis")

## function to let color fade
fade_color <- function(color, alpha = 1) {
  # 'alpha; behaves similarly to transparancy in ggplot2, but colors are faded
  # to white and not transparant (the default is to use fully opaque colors)
  rgb <- col2rgb(color) / 255
  faded <- mixcolor(1 - alpha, RGB(rgb[1, ], rgb[2, ], rgb[3, ]), RGB(1, 1, 1))
  rgb(faded@coords[, 1], faded@coords[, 2], faded@coords[, 3])
}

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

## control parameters for data generation
n <- 100            # number of observations
a <- b <- c <- 0.8  # true effects
sigma <- 0.4        # standard deviation of error terms
level <- 0.95       # confidence level
seed <- 20170627    # seed for the random number generator

## generate data
RNGversion("3.5.3")
set.seed(seed)
X <- sample(0:1, n, replace = TRUE)
M <- a * X + sigma * rnorm(n)
Y <- b * M + c * X + sigma * rnorm(n)

## different situations for outlier
mu <- c(1, 3.2, -3.2)    # means for contaminated data
affected <- arrow(angle = 0, length = unit(0, "cm"))

## graphical parameters for pointers to outlier
offset <- c(NA, -0.275, 1)
label <- "Outlier"
col_pointer <- "darkgray"
pointer <- arrow(length = unit(0.15, "cm"), ends = "last", type = "closed")

## control parameters for plots
lower_X <- -1.25
lower_Y <- mu[3]
max_Y <- max(Y)
col_lines <- viridis_pal(begin = 0.3, end = 0.7, direction = -1)(2)
# col_points <- fade_color(col_lines, alpha = 0.8)  # faded colors for points
col_points <- col_lines                           # no color fade for points
line_types <- c("dotdash", "solid", "dashed")
line_size <- 2/3
effect <- arrow(length = unit(0.2, "cm"), ends = "both", type = "open")
arrow_size <- 2/3
scenarios <- c("No deviations from normality", "Including outlier")
methods <- c("OLS", "ROBMED")


## loop over scenario and method
df_points <- NULL
df_MX <- df_YX <- df_YMX <- NULL
df_a <- df_c <- df_c_prime <- NULL
df_ab <- list()
df_label_top <- df_label_right <- df_label_arrow <- NULL
df_outlier_arrow <- df_outlier_label <- NULL
for (scenario in scenarios) {

  ## add outlier in second iteration
  if (scenario == scenarios[2]) {
    X <- c(X, mu[1])
    M <- c(M, mu[2])
    Y <- c(Y, mu[3])
  }

  for (method in methods) {

    ## generate data frame and fit mediation
    data <- data.frame(X, Y, M, Scenario = scenario, Method = method)
    if (method == methods[1]) {
      fit <- fit_mediation(data, "X", "Y", "M", robust = FALSE)
    } else if (method == methods[2]) {
      fit <- fit_mediation(data, "X", "Y", "M", robust = TRUE)
    }

    ## combine data frames for plotting
    data$X <- as.factor(data$X)
    df_points <- rbind(df_points, data)

    ## extract coefficients
    coef_MX <- coef(fit$fit_mx)
    coef_YMX <- coef(fit$fit_ymx)
    if (method == methods[1]) {
      coef_YX <- coef(fit$fit_yx)
    } else if (method == methods[2]) {
      i2 <- coef_YMX[1] + coef_YMX[2] * coef_MX[1]
      coef_YX <- c(i2, i2 + fit$total)
    }

    ## equations to be plotted
    tmp_MX <- data.frame(Equation = c("i1", "i1 + a"),
                         Value = c(unname(coef_MX)[1], sum(coef_MX)),
                         Scenario = scenario, Method = method)
    df_MX <- rbind(df_MX, tmp_MX)
    tmp_YMX <- data.frame(Equation = c("i2 + b*M", "i2 + b*M + c"),
                          Intercept = c(unname(coef_YMX)[1], sum(coef_YMX[-2])),
                          Slope = rep.int(unname(coef_YMX)[2], 2),
                          Scenario = scenario, Method = method)
    df_YMX <- rbind(df_YMX, tmp_YMX)
    tmp_YX <- data.frame(Equation = c("i3", "i3 + c'"),
                         Value = c(unname(coef_YX)[1], sum(coef_YX)),
                         Scenario = scenario, Method = method)
    df_YX <- rbind(df_YX, tmp_YX)

    ## effects to be plotted
    tmp_a <- data.frame(x = tmp_MX$Value, y = lower_Y, Scenario = scenario,
                        Method = method)
    df_a <- rbind(df_a, tmp_a)
    tmp_c_prime <- data.frame(x = lower_X, y = tmp_YX$Value,
                              Scenario = scenario, Method = method)
    df_c_prime <- rbind(df_c_prime, tmp_c_prime)
    intersection_Y <- tmp_YMX[1, "Intercept"] + tmp_YMX[1, "Slope"] * tmp_MX[2, "Value"]
    tmp_c <- data.frame(x = tmp_MX[2, "Value"],
                        y = c(intersection_Y, tmp_YX[2, "Value"]),
                        Scenario = scenario, Method = method)
    df_c <- rbind(df_c, tmp_c)
    tmp_ab <- data.frame(x = tmp_MX[2, "Value"],
                         y = c(tmp_YX[1, "Value"], intersection_Y),
                         Scenario = scenario, Method = method)
    df_ab[[paste(scenario, method, sep = "_")]] <- tmp_ab

    ## effect labels to be plotted
    # labels plotted on top of arrows
    tmp <- data.frame(x = mean(tmp_a$x), y = mean(tmp_a$y),
                      Label = "italic(hat(a))", Scenario = scenario,
                      Method = method)
    df_label_top <- rbind(df_label_top, tmp)
    # labels plotted to the right of arrows
    labels_right <- c("italic(hat(c))", "italic(paste(hat(c), \"'\"))",
                      "italic(widehat(ab))")
    if(scenario == scenarios[2] && method == methods[1]) {
      # no space for arrow of standard method when the outlier is included,
      # so label should be offset a bit
      offset_ab <- c(x = 0.4, y = -1/3)
      tmp <- data.frame(x = c(mean(tmp_c$x), mean(tmp_c_prime$x),
                              mean(tmp_ab$x) + offset_ab["x"]),
                        y = c(mean(tmp_c$y), mean(tmp_c_prime$y),
                              mean(tmp_ab$y) + offset_ab["y"]),
                        Label = labels_right, Scenario = scenario,
                        Method = method)
      df_label_right <- rbind(df_label_right, tmp)
      # add an arrow to connect the label to the effect
      start <- c(mean(tmp_ab$x) + offset_ab["x"],
                 mean(tmp_ab$y) + offset_ab["y"])
      end <- c(x = mean(tmp_ab$x), y = mean(tmp_ab$y))
      cut_start <- start - 0.1 * (end - start)
      cut_end <- start + 0.95 * (end - start)
      tmp <- data.frame(x = cut_start["x"], y = cut_start["y"],
                        x_end = cut_end["x"], y_end = cut_end["y"],
                        Scenario = scenario, Method = method)
      df_label_arrow <- rbind(df_label_arrow, tmp)
    } else {
      # enough space for all the arrows and labels
      tmp <- data.frame(x = c(mean(tmp_c$x), mean(tmp_c_prime$x),
                              mean(tmp_ab$x)),
                        y = c(mean(tmp_c$y), mean(tmp_c_prime$y),
                              mean(tmp_ab$y)),
                        Label = labels_right, Scenario = scenario,
                        Method = method)
      df_label_right <- rbind(df_label_right, tmp)
    }

    ## extra information on outlier
    if(scenario == scenarios[2]) {
      ## arrow to outlier to be plotted
      start <- c(x = mu[2] + offset[2], y = mu[3] + offset[3])
      end <- c(x = mu[2], y = mu[3])
      cut_start <- start + 0.1 * (end - start)
      cut_end <- start + 0.8 * (end - start)
      tmp <- data.frame(x = cut_start["x"], y = cut_start["y"],
                        x_end = cut_end["x"], y_end = cut_end["y"],
                        Scenario = scenario, Method = method)
      df_outlier_arrow <- rbind(df_outlier_arrow, tmp)

      ## outlier label to be plotted
      tmp <- data.frame(x = mu[2] + offset[2], y = mu[3] + offset[3],
                        Label = label, Scenario = scenario, Method = method)
      df_outlier_label <- rbind(df_outlier_label, tmp)
    }
  }
}

# fix data frame for arrows for indirect effect
df_ab <- list(arrow = rbind(df_ab[[1]], df_ab[[2]], df_ab[[4]]),
              line = df_ab[[3]])

# equation labels for legend
equations <- expression(italic(hat(M) == hat(i)[1]),
                        italic(paste(hat(M) == hat(i)[1] + hat(a), "  ")),
                        italic(hat(Y) == hat(i)[2] + hat(b) %.% M),
                        italic(paste(hat(Y) == hat(i)[2] + hat(b) %.% M + hat(c), "'  ")),
                        italic(hat(Y) == hat(i)[3]),
                        italic(hat(Y) == hat(i)[3] + hat(c)))


## plot data and mediation effects
# create plot object
p <- ggplot() +
  geom_point(aes(x = M, y = Y, fill = X), data = df_points,
             color = "transparent", shape = 21, size = 2.5) +
  geom_vline(aes(xintercept = Value, linetype = Equation, color = Equation),
             data = df_MX, size = line_size, show.legend = FALSE) +
  geom_hline(aes(yintercept = Value, linetype = Equation, color = Equation),
             data = df_YX, size = line_size) +
  geom_abline(aes(intercept = Intercept, slope = Slope, linetype = Equation,
                  color = Equation), data = df_YMX, size = line_size,
              show.legend = FALSE) +
  geom_segment(aes(x = x, y = y, xend = x_end, yend = y_end),
               data = df_label_arrow, color = col_pointer,
               arrow = pointer) +
  geom_line(aes(x = x, y = y), data = df_a, arrow = effect,
            size = arrow_size) +
  geom_line(aes(x = x, y = y), data = df_c, arrow = effect,
            size = arrow_size) +
  geom_line(aes(x = x, y = y), data = df_c_prime, arrow = effect,
            size = arrow_size) +
  geom_line(aes(x = x, y = y), data = df_ab[[1]], arrow = effect,
            size = arrow_size) +
  geom_line(aes(x = x, y = y), data = df_ab[[2]], arrow = affected,
            size = arrow_size) +
  geom_text(aes(x = x, y = y, label = Label), data = df_label_top, size = 4,
            hjust = "center", vjust = "bottom", parse = TRUE, nudge_y = 0.1) +
  geom_text(aes(x = x, y = y, label = Label), data = df_label_right, size = 4,
            hjust = "left", vjust = "middle", parse = TRUE, nudge_x = 0.1) +
  geom_segment(aes(x = x, y = y, xend = x_end, yend = y_end),
               data = df_outlier_arrow, color = col_pointer,
               arrow = pointer) +
  geom_text(aes(x = x, y = y, label = Label), data = df_outlier_label, size = 4,
            hjust = "center", vjust = "bottom") +
  scale_fill_manual(values = col_points) +
  scale_color_manual(values = rep.int(col_lines, 3), labels = equations,
                     guide = guide_legend(label.hjust = 0)) +
  scale_linetype_manual(values = rep(line_types, each = 2), labels = equations,
                        guide = guide_legend(override.aes = aes(xintercept = NA))) +
  facet_grid(Method~Scenario, scales = "free") + theme_bw() +
  xlim(lower_X, mu[2]) + ylim(lower_Y, max_Y) +
  theme(axis.text = element_blank(),
        axis.title = element_text(size = 13, face = "italic"),
        legend.position = "top",
        legend.direction = "vertical",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        plot.margin = unit(c(0, 0.25, 0.25, 0.25), "line"),
        strip.text = element_text(size = 12)) +
  guides(color = guide_legend(label.hjust = 0, keywidth = unit(1.5, "line"),
                              nrow = 2),
         fill = guide_legend(title.theme = element_text(face = "italic",
                                                        angle = 0)))
# plot to file
pdf("illustration/mediation_color.pdf", width = 7.5, height = 6)
print(p)
dev.off()

## plot in grayscale
p <- ggplot() +
  geom_point(aes(x = M, y = Y, fill = X), data = df_points,
             color = "transparent", shape = 21, size = 2.5) +
  geom_vline(aes(xintercept = Value, linetype = Equation, color = Equation),
             data = df_MX, size = line_size, show.legend = FALSE) +
  geom_hline(aes(yintercept = Value, linetype = Equation, color = Equation),
             data = df_YX, size = line_size) +
  geom_abline(aes(intercept = Intercept, slope = Slope, linetype = Equation,
                  color = Equation), data = df_YMX, size = line_size,
              show.legend = FALSE) +
  geom_segment(aes(x = x, y = y, xend = x_end, yend = y_end),
                        data = df_label_arrow, color = col_pointer,
                        arrow = pointer) +
  geom_line(aes(x = x, y = y), data = df_a, arrow = effect,
            size = arrow_size) +
  geom_line(aes(x = x, y = y), data = df_c, arrow = effect,
            size = arrow_size) +
  geom_line(aes(x = x, y = y), data = df_c_prime, arrow = effect,
            size = arrow_size) +
  geom_line(aes(x = x, y = y), data = df_ab[[1]], arrow = effect,
            size = arrow_size) +
  geom_line(aes(x = x, y = y), data = df_ab[[2]], arrow = affected,
            size = arrow_size) +
  geom_text(aes(x = x, y = y, label = Label), data = df_label_top, size = 4,
            hjust = "center", vjust = "bottom", parse = TRUE, nudge_y = 0.1) +
  geom_text(aes(x = x, y = y, label = Label), data = df_label_right, size = 4,
            hjust = "left", vjust = "middle", parse = TRUE, nudge_x = 0.1) +
  geom_segment(aes(x = x, y = y, xend = x_end, yend = y_end),
               data = df_outlier_arrow, color = col_pointer,
               arrow = pointer) +
  geom_text(aes(x = x, y = y, label = Label), data = df_outlier_label, size = 4,
            hjust = "center", vjust = "bottom") +
  scale_fill_manual(values = col2gray(col_points)) +
  scale_color_manual(values = rep.int(col2gray(col_points), 3), labels = equations,
                     guide = guide_legend(label.hjust = 0)) +
  scale_linetype_manual(values = rep(line_types, each = 2), labels = equations,
                        guide = guide_legend(override.aes = aes(xintercept = NA))) +
  facet_grid(Method~Scenario, scales = "free") + theme_bw() +
  xlim(lower_X, mu[2]) + ylim(lower_Y, max_Y) +
  theme(axis.text = element_blank(),
        axis.title = element_text(size = 13, face = "italic"),
        legend.position = "top",
        legend.direction = "vertical",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        plot.margin = unit(c(0, 0.25, 0.25, 0.25), "line"),
        strip.text = element_text(size = 12)) +
  guides(color = guide_legend(label.hjust = 0, keywidth = unit(1.5, "line"),
                              nrow = 2),
         fill = guide_legend(title.theme = element_text(face = "italic",
                                                        angle = 0)))
pdf("illustration/mediation_bw.pdf", width = 7.5, height = 6)
print(p)
dev.off()
