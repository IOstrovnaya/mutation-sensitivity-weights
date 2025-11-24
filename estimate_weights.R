#' Estimate observation-level weights using VAFs
#'
#' @param y        Binary grouping variable (e.g. clinical group, stage; 0/1).
#' @param x        Binary mutation indicator (0 = no mutation, 1 = mutation).
#' @param vaf      Numeric vector of variant allele frequencies (VAFs), same length as x.
#' @param type     "uncond" (default) uses unconditional mean reconstruction;
#'                 "cond" uses observed conditional means directly.
#' @param var_type "separate" (default) uses separate variances by y group;
#'                 "common" uses a common variance across groups.
#'
#' @return Numeric vector of weights of length(x).
#'         For x == 1 (mutated), weight = 1.
#'         For x == 0, weight depends on y and VAF distribution.
estimate.weights <- function(y, x, vaf,
                             type     = c("uncond", "cond"),
                             var_type = c("separate", "common")) {

  ## ---- argument handling ----
  type     <- match.arg(type)
  var_type <- match.arg(var_type)

  if (length(y) != length(x) || length(x) != length(vaf)) {
    stop("y, x, and vaf must have the same length.")
  }
  if (!all(x %in% c(0, 1))) {
    stop("x must be a binary (0/1) mutation indicator.")
  }
  if (!all(y %in% c(0, 1))) {
    stop("y must be a binary (0/1) group indicator.")
  }

  ## ---- helper: solve quadratic equation a * z^2 + b * z + c = 0 ----
  quad <- function(a, b, c) {
    a <- as.complex(a)

    roots <- c(
      (-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
      (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
    )

    # If both roots are real, drop imaginary part
    if (all(Im(roots) == 0)) {
      roots <- Re(roots)
    }

    # If both roots are identical, just return one
    if (roots[1] == roots[2]) {
      return(roots[1])
    }

    roots
  }

  ## ---- helper: unconditional mean given conditional mean/variance ----
  # Given conditional mean (condmean) and (assumed) equal unconditional variance condvar,
  # solve for an unconditional mean parameter.
  uncondmean <- function(condmean, condvar) {
    quad(1, -condmean, condvar)
  }

  ## ---- conditional means of VAF among mutation carriers by group ----
  # x == 1 : mutation carriers
  # y == 0/1 : group indicator (e.g. stage 0 vs 1)

  m.theta0 <- mean(vaf[x == 1 & y == 0], na.rm = TRUE)
  m.theta1 <- mean(vaf[x == 1 & y == 1], na.rm = TRUE)

  if (type == "cond") {
    # Use conditional means directly
    m.theta0.uncond <- m.theta0
    m.theta1.uncond <- m.theta1

  } else { # type == "uncond"

    ## ---- variances ----
    if (var_type == "common") {
      # Common variance across groups, estimated among carriers
      v.common <- stats::var(vaf[x == 1], na.rm = TRUE)
      v.theta0 <- v.common
      v.theta1 <- v.common
    } else {
      # Separate variances by group among carriers
      v.theta0 <- stats::var(vaf[x == 1 & y == 0], na.rm = TRUE)
      v.theta1 <- stats::var(vaf[x == 1 & y == 1], na.rm = TRUE)
    }

    ## ---- reconstruct unconditional means from conditional moments ----
    # For each group, solve quadratic and pick the root closest to the
    # observed conditional mean. If the solution is complex, fall back
    # to the conditional mean.

    a0 <- uncondmean(m.theta0, v.theta0)
    m.theta0.uncond <- if (is.complex(a0[1])) {
      m.theta0
    } else {
      a0[which.min(abs(m.theta0 - a0))]
    }

    a1 <- uncondmean(m.theta1, v.theta1)
    m.theta1.uncond <- if (is.complex(a1[1])) {
      m.theta1
    } else {
      a1[which.min(abs(m.theta1 - a1))]
    }
  }

  ## ---- construct weights ----
  wt <- rep(NA_real_, length(x))

  # No mutation observed, group y = 0
  wt[x == 0 & y == 0] <- m.theta0.uncond

  # No mutation observed, group y = 1
  wt[x == 0 & y == 1] <- m.theta1.uncond

  # Mutation observed: full weight
  wt[x == 1] <- 1

  wt
}
