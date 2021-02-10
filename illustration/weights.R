# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# This script reproduces the plot of the loss functions and the corresponding
# weighting functions of (robust) mediation analysis from Figure 2 of the paper:
#
# Alfons, A., Ates, N.Y., & Groenen, P.J.F. (2021). A Robust Bootstrap Test for
# Mediation Analysis. Organizational Research Methods, accepted for publication.
#
# Before running this script, please make sure that all necessary packages are
# installed.


# load packages
library("dplyr")
library("robmed")
library("scales")

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

# control parameters
tuning.psi <- reg_control()$tuning.psi
col <- hue_pal()(4)[c(1, 3)]
lty <- c(2, 1)
line_size <- 2/3

# grid of x-values
x <- seq(-4, 4, length.out = 100)

# loss functions
OLS_loss <- data.frame(x = x,
                       y = x^2,
                       Method = "OLS",
                       Function = "Loss")
robust_loss <- data.frame(x = x,
                          y = Mpsi(x, tuning.psi, "bisquare", deriv = -1),
                          Method = "Robust",
                          Function = "Loss")

# weight functions
OLS_weight <- data.frame(x = x,
                         y = 1,
                         Method = "OLS",
                         Function = "Weight")
robust_weight <- data.frame(x = x,
                            y = Mwgt(x, tuning.psi, "bisquare"),
                            Method = "Robust",
                            Function = "Weight")

# data frame for plot
df <- rbind(OLS_loss, robust_loss, OLS_weight, robust_weight) %>% filter(y <= 4)

# create plot
p <- ggplot() +
  geom_line(aes(x = x, y = y, color = Method, linetype = Method),
            data = df, size = line_size) +
  labs(x = NULL, y = NULL) +
  scale_linetype_manual(values = lty) +
  facet_wrap(~ Function, scales = "free_y") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 13),
                     legend.key = element_rect(color = NA),
                     legend.key.width = unit(1.5, "line"),
                     legend.text = element_text(size = 12),
                     legend.title = element_text(size = 13),
                     strip.text = element_text(size = 12))

# plot to file (in color)
pdf(file = "illustration/weights_color.pdf", width = 8.5, height = 3.75)
print(p + scale_color_manual(values = col))
dev.off()

# plot to file (in grayscale)
pdf(file = "illustration/weights_bw.pdf", width = 8.5, height = 3.75)
print(p + scale_color_manual(values = col2gray(col)))
dev.off()
