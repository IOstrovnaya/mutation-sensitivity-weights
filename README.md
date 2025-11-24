# VAF-based weights for mutation–group association models

This repository contains an R function to compute patient-level weights
based on variant allele frequencies (VAFs) and a simple example of how to
use these weights in logistic regression models of mutation status on a
binary group (e.g. stage).

## Files

- `estimate_weights.R`  
  Defines the function `estimate.weights()`.

- `example_usage.R`  
  Simulated example showing how to:
  - compute weights from `y` (group), `x` (mutation), and `vaf`;
  - fit unweighted and weighted logistic regression models.

## Usage

```r
# Download the function
source("https://raw.githubusercontent.com/<your-user>/mutation-vaf-weights/main/estimate_weights.R")

# y: group (0/1), x: mutation (0/1), vaf: variant allele frequency
wt <- estimate.weights(
  y        = stage,
  x        = mut,
  vaf      = vaf,
  type     = "uncond",
  var_type = "separate"
)

fit <- glm(mut ~ stage, family = binomial(), weights = wt)
summary(fit)
