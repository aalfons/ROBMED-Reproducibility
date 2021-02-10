# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# This script reproduces the descriptive statistics (Tables 1 and 2) of the
# data used in the illustrative empirical case of the paper:
#
# Alfons, A., Ates, N.Y., & Groenen, P.J.F. (2021). A Robust Bootstrap Test for
# Mediation Analysis. Organizational Research Methods, accepted for publication.
#
# Before running this script, please make sure that all necessary packages are
# installed.


# load packages and data
library("dplyr")
library("tidyr")
data("BSG2014", package = "robmed")

# select the relevant variables (aggregated from multiple Likert items)
BSG2014 <- BSG2014 %>%
  select(TaskConflict, TeamCommitment, ValueDiversity)

# compute descriptive statistics
descriptives <- BSG2014 %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  summarize(Mean = mean(Value),
            SD = sd(Value),
            Median = median(Value),
            MAD = mad(Value),
            Min = min(Value),
            Max = max(Value))

# round to 3 digits and print descriptive statistics
descriptives %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  as.data.frame()

# Spearman's rank correlation is more robust than the Pearson correlation,
# but it needs to be transformed to be consistent for the Pearson correlation
# (see Croux & Dehon 2010)
corSpearman <- function(...) 2 * sin(pi/6 * cor(..., method = "spearman"))

# compute correlation table and round to 3 digits
R <- corSpearman(BSG2014)
round(R, digits = 3)
