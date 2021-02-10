# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# This script reproduces the results of the illustrative empirical case,
# including Table 4 and Figures 3 and 4, from the paper:
#
# Alfons, A., Ates, N.Y., & Groenen, P.J.F. (2021). A Robust Bootstrap Test for
# Mediation Analysis. Organizational Research Methods, accepted for publication.
#
# Before running this script, please make sure that all necessary packages are
# installed.


# load packages and data
library("robmed")
library("car")
library("dplyr")
library("scales")
data("BSG2014")

# function to convert colors to grayscale
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


# illustrative empirical case
x <- "ValueDiversity"
y <- "TeamCommitment"
m <- "TaskConflict"

# plot data
labs <- c("Value diversity", "Team commitment", "Task conflict")
names(labs) <- c(x, y, m)
# pdf(file = "example/illustrative_case_data.pdf", width = 7, height = 7)
# plot(BSG2014[, c(x, y, m)], labels = gsub(" ", "\n", labs),
#      pch = 21, bg = "black", cex = 1.9, las = 1, cex.axis = 1.4,
#      oma = rep.int(2.5, 4))
# dev.off()

# seed of random number generator
RNGversion("3.5.3")
seed <- 20150601

# OLS bootstrap and OLS Sobel tests
set.seed(seed)
ols_fit <- fit_mediation(BSG2014, x = x, y = y, m = m, method = "regression",
                         robust = FALSE, family = "gaussian")
ols_boot <- test_mediation(ols_fit, test = "boot")
ols_sobel <- test_mediation(ols_fit, test = "sobel")
# Box-Cox transformations followed by OLS bootstrap
BSG2014_bc <- BSG2014[, c(x, y, m)] %>%
  lapply(function(variable) {
    bc <- powerTransform(variable, family = "bcPower")
    bcPower(variable, lambda = bc$lambda)
  }) %>%
  as.data.frame()
set.seed(seed)
bc_boot <- test_mediation(BSG2014_bc, x = x, y = y, m = m, test = "boot",
                          method = "regression", robust = FALSE,
                          family = "gaussian")
# regression with normal, skew-normal, t, or skew-t error distribution
# within bootstrap, selection of error distribution via BIC
set.seed(seed)
snt_boot <- test_mediation(BSG2014, x = x, y = y, m = m, test = "boot",
                           method = "regression", robust = FALSE,
                           family = "select")
# bootstrap following winsorization of the data
set.seed(seed)
winsorized_boot <- test_mediation(BSG2014, x = x, y = y, m = m, test = "boot",
                                  method = "covariance", robust = TRUE)
# bootstrap using median regression
set.seed(seed)
median_boot <- test_mediation(BSG2014, x = x, y = y, m = m, test = "boot",
                              method = "regression", robust = "median")
# ROBMED
set.seed(seed)
robust_boot <- test_mediation(BSG2014, x = x, y = y, m = m, test = "boot",
                              method = "regression", robust = TRUE)

# show results
summary(ols_boot)
summary(ols_sobel)
summary(bc_boot)
summary(snt_boot)
summary(winsorized_boot)
summary(median_boot)
summary(robust_boot, plot = FALSE)

# compute or extract p-values of indirect effect
p_value(ols_boot, parm = "ab")
p_value(ols_sobel, parm = "ab")
p_value(bc_boot, parm = "ab")
p_value(snt_boot, parm = "ab")
p_value(winsorized_boot, parm = "ab")
p_value(median_boot, parm = "ab")
p_value(robust_boot, parm = "ab")


# let's have a closer look at the effect of x on m

# extract information for diagnostic plot with tolerance ellipse
boot_list <- list("OLS bootstrap" = ols_boot,
                  "ROBMED" = robust_boot)
ellipses <- setup_ellipse_plot(boot_list, horizontal = x, vertical = m)

# control parameters for plots
colors <- c("#F8766D", "#00BFC4")
linetypes <- c("dashed", "solid")
line_size <- 2/3

# create plot
p_ellipse <- ggplot() +
  geom_path(aes(x = x, y = y, color = Method, linetype = Method),
            data = ellipses$ellipse, size = line_size) +
  geom_point(aes(x = x, y = y, fill = Weight), data = ellipses$data,
             shape = 21, size = 3) +
  geom_abline(aes(intercept = intercept, slope = slope,
                  color = Method, linetype = Method),
              data = ellipses$line, size = line_size,
              show.legend = FALSE) +
  scale_linetype_manual(values = linetypes) +
  scale_fill_gradient(limits = 0:1, low = "white", high = "black") +
  labs(x = labs[x], y = labs[m]) + theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.key.width = unit(1.5, "line"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13))

# plot to file (in color)
pdf(file = "example/illustrative_case_mx_color.pdf", width = 6.5, height = 5)
print(p_ellipse + scale_color_manual(values = colors))
dev.off()

# plot to file (in grayscale)
pdf(file = "example/illustrative_case_mx_bw.pdf", width = 6.5, height = 5)
print(p_ellipse + scale_color_manual(values = col2gray(colors)))
dev.off()

# diagnostic plot of robust regression weights
p_weights <- weight_plot(robust_boot, outcome = m) +
  labs(title = NULL) + scale_y_continuous(labels = percent) +
  theme_bw() + theme(legend.position = "top")

# plot to file (in color)
pdf(file = "example/illustrative_case_weights_color.pdf",
    width = 6.5, height = 5)
print(p_weights + scale_color_manual("", values = c("black", colors[2])))
dev.off()

# plot to file (in grayscale)
pdf(file = "example/illustrative_case_weights_bw.pdf",
    width = 6.5, height = 5)
print(p_weights + scale_color_manual(values = col2gray(c("black", colors[2]))))
dev.off()
