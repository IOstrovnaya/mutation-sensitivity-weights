# ------------------------------------------------
# Example: mutation ~ stage with VAF-based weights
# ------------------------------------------------

# Simulate a small example dataset
set.seed(123)

n <- 200
stage <- rbinom(n, size = 1, prob = 0.5)           # y: group (e.g. 0 = early, 1 = late stage)
mut   <- rbinom(n, size = 1, prob = plogis(-1 + 0.8 * stage))  # x: mutation indicator

# VAF: positive for mutation carriers, 0 otherwise
vaf <- ifelse(mut == 1,
              runif(n, min = 0.2, max = 0.8),      # carriers
              0)                                   # non-carriers

# Compute weights
wt <- estimate.weights(
  y        = stage,
  x        = mut,
  vaf      = vaf,
  type     = "uncond",
  var_type = "separate"
)

# Check a few rows
head(data.frame(stage, mut, vaf, wt))

# Fit unweighted logistic regression: mutation ~ stage
fit_unweighted <- glm(mut ~ stage, family = binomial())
summary(fit_unweighted)

# Fit weighted logistic regression using VAF-based weights
fit_weighted <- glm(mut ~ stage, family = binomial(), weights = wt)
summary(fit_weighted)
