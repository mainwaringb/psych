library(psych)

# Generate four variables - two on a 6-pt scale, and two on a 4-pt scale
set.seed(6171729)
mydf <- data.frame(x1 = sample(1:6, 1000, replace = TRUE), x2 = sample(1:6, 1000, replace = TRUE), x3 = sample(1:4, 1000, replace = TRUE), x4 = sample(1:4, 1000, replace = TRUE))

# ISSUE A: try running with "global" as true - we should get a warning message and force global = FALSE, but don't
polychoric(mydf, global = TRUE)

# ISSUE B: Turn set "global" to false as we have unbalanced scales, and run "polychoric" with a continuity correction
polychoric(mydf, global = FALSE, correct = TRUE)
polychoric(mydf, global = FALSE, correct = .5)
polychoric(mydf, global = FALSE)

#When the correction is disabled, the code runs
polychoric(mydf, global = FALSE, correct = FALSE)
