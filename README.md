# VAF-based weights for mutation–group association models

This repository contains an R function to compute patient-level weights
based on variant allele frequencies (VAFs) and a simple example of how to
use these weights in logistic regression models of mutation status on a
binary group (e.g. stage).
The preprint is located at https://www.biorxiv.org/content/10.64898/2026.10.02.756288v1




## Usage

```r
# Download the function
source("https://raw.githubusercontent.com/<your-user>/mutation-vaf-weights/main/Data examples code.R")

# y: group (0/1), x: mutation (0/1), vaf: variant allele frequency

wt<-estimate.weights.geom(y,x,pmin(1,2*vaf))
fit <- glm(mut ~ stage, family = binomial(), weights = wt)
summary(fit)
