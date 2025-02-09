# ------------------------------
# Example R code
# ------------------------------

set.seed(123)  # for reproducibility

# --- Parameters ---
nGroups       <- 100     # total number of groups
nNoEffect    <- 50      # number of groups with no effect
nEffect       <- 50      # number of groups with effect
baselineProb  <- 0.05    # baseline response rate (control)
ratio   <- 10      # treatment effect ratio for "with effect" groups
nControl      <- 1000    # sample size in control arm
nTreatment    <- 1000    # sample size in treatment arm

# --- Create data frame structure ---
df <- data.frame(
  group_id = 1:nGroups,
  # For simplicity, let's assign cluster=1 to the first 50 groups, cluster=2 to the last 50
  cluster = rep(c(1, 2), each = 50),
  # Label which groups have effect vs. no effect
  effect_type = rep(c("no_effect", "effect"), c(nNoEffect, nEffect)),
  stringsAsFactors = FALSE
)

# --- Generate responses ---
# All control arms draw from Binomial(1000, 0.05)
df$control_response <- rbinom(nGroups, size = nControl, prob = baselineProb)

# Treatment arms:
# - If no effect: same prob = 0.05
# - If effect: prob = 0.05 * 10 = 0.5
df$treatment_response <- ifelse(
  df$effect_type == "no_effect",
  rbinom(nGroups, size = nTreatment, prob = baselineProb),
  rbinom(nGroups, size = nTreatment, prob = baselineProb * ratio)
)

# --- View the first few rows ---
# head(df)
