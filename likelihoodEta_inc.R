#-------------------------------------
# Imcomplete Beta Version
#-------------------------------------
library(zipfR)

# # Example parameters
# n1 <- 30
# y1 <- 20
# y0 <- 2
# n0 <- 30
# a_pi <- 1
# b_pi <- 1
# eta <--2.1
# likelihoodEta_inc(n0, y0, n1, y1, a_pi, b_pi, eta)

# Define the function to calculate the summation
likelihoodEta_inc <- function(n0, y0, n1, y1, a_pi, b_pi, eta) {
  # Initialize the summation result
  if(eta <= 0){
    positives <- c()
    negatives <- c()
    for (k in 0:(n1 - y1)) {
      log_term_k = lchoose(n1 - y1, k) + k*eta + lbeta(y0 + y1 + a_pi + k, n0 - y0 + b_pi)  # your p, q, etc.
      if (k %% 2 == 0) {
        # even k => positive sign
        positives <- c(positives, log_term_k)
      } else {
        # odd k => negative sign
        negatives <- c(negatives, log_term_k)
      }
    }

    # log-sum-exp for positives
    maxPos <- max(positives)
    sumPos <- 0
    for (Lp in positives) {
      sumPos <- sumPos + exp(Lp - maxPos)
    }
    logSumPos <- maxPos + log(sumPos)   # log of sum of positive terms

    # log-sum-exp for negatives
    maxNeg <- max(negatives)
    sumNeg <- 0
    for (Ln in negatives) {
      sumNeg <- sumNeg + exp(Ln - maxNeg)
    }
    logSumNeg <- maxNeg + log(sumNeg)   # log of sum of negative terms

    log_diff <- log(1 - exp(logSumNeg - logSumPos))
    # If '1 - exp(logSumNeg - logSumPos)' is slightly negative due to rounding:
    if (is.nan(log_diff) || log_diff == -Inf) {
      # They basically canceled out to a negative or zero
      # Theoretically, that shouldn't happen, so clamp:
      log_diff <- -Inf
    }
    log_final <- logSumPos + log_diff

    lresult <- lchoose(n0, y0) + lchoose(n1,y1) + eta * y1 - lbeta(a_pi,b_pi)  + log_final
    result <- exp(lresult)

  }else{

    positives <- c()
    negatives <- c()
    for (k in 0:(n1 - y1)) {
      log_term_k <- lchoose(n1 - y1, k) + k * eta +
        Ibeta(exp(-eta), y0 + y1 + a_pi + k, n0 - y0 + b_pi, lower = TRUE, log = TRUE)
      if (k %% 2 == 0) {
        # even k => positive sign
        positives <- c(positives, log_term_k)
      } else {
        # odd k => negative sign
        negatives <- c(negatives, log_term_k)
      }
    }

    # log-sum-exp for positives
    maxPos <- max(positives)
    sumPos <- 0
    for (Lp in positives) {
      sumPos <- sumPos + exp(Lp - maxPos)
    }
    logSumPos <- maxPos + log(sumPos)   # log of sum of positive terms

    # log-sum-exp for negatives
    maxNeg <- max(negatives)
    sumNeg <- 0
    for (Ln in negatives) {
      sumNeg <- sumNeg + exp(Ln - maxNeg)
    }
    logSumNeg <- maxNeg + log(sumNeg)   # log of sum of negative terms

    suppressWarnings(log_diff <- log(1 - exp(logSumNeg - logSumPos)))
    # If '1 - exp(logSumNeg - logSumPos)' is slightly negative due to rounding:
    if (is.nan(log_diff) || log_diff == -Inf) {
      # They basically canceled out to a negative or zero
      # Theoretically, that shouldn't happen, so clamp:
      log_diff <- -Inf
    }
    log_final <- logSumPos + log_diff

    lresult <- lchoose(n0, y0) + lchoose(n1,y1) + eta * y1 - lbeta(a_pi,b_pi)  + log_final
    result <- exp(lresult)
  }

  return(result)
}


# likelihoodEta <- function(n0, y0, n1, y1, a_pi, b_pi, eta) {
# 
#   if(eta <= 0){
# 
#     integral <- 0
#     # Loop over k from 0 to (n1j - y1j)
#     for (k in 0:(n1 - y1)) {
#       # Compute each term in the summation
#       # log_term <- lchoose(n1 - y1, k) + k * eta +
#       #   lgamma(y0 + y1 + a_pi + k) + lgamma(n0 - y0 + b_pi) -
#       #   lgamma(n0 + y1 + a_pi + b_pi + k)
# 
#       log_term <- lchoose(n1 - y1, k) + k * eta +
#         lbeta(y0 + y1 + a_pi + k, n0 - y0 + b_pi)
# 
#       term <- (-1)^k * exp(log_term)
# 
#       # Add the term to the result
#       integral <- integral + term
#     }
# 
#     lresult <- lchoose(n0, y0) + lchoose(n1,y1) + eta * y1 - lbeta(a_pi,b_pi)  + log(integral)
#     result <- exp(lresult)
#   }else{
# 
#     integral <- 0
#     # Loop over k from 0 to (n1j - y1j)
#     for (k in 0:(n1 - y1)) {
#       # Compute each term in the summation
#       # log_term <- lchoose(n1 - y1, k) + k * eta +
#       #   lgamma(y0 + y1 + a_pi + k) + lgamma(n0 - y0 + b_pi) -
#       #   lgamma(n0 + y1 + a_pi + b_pi + k)
# 
#       log_term <- lchoose(n1 - y1, k) + k * eta +
#         Ibeta(exp(-eta), y0 + y1 + a_pi + k, n0 - y0 + b_pi, lower = TRUE, log = TRUE)
# 
#       term <- (-1)^k * exp(log_term)
# 
#       # Add the term to the result
#       integral <- integral + term
#     }
# 
#     lresult <- lchoose(n0, y0) + lchoose(n1,y1) + eta * y1 - lbeta(a_pi,b_pi)  + log(integral)
#     result <- exp(lresult)
#   }
# 
#   return(result)
# }
# 
# 
# 
# likelihoodEta <- function(n0, y0, n1, y1, a_pi, b_pi, eta) {
#   
#   if(eta <= 0){
#     
#     integral <- 0
#     log_term <- numeric(n1 - y1 + 1)
#     signs <- numeric(n1 - y1 + 1)
#     # Loop over k from 0 to (n1j - y1j)
#     for (k in 0:(n1 - y1)) {
#       # Compute each term in the summation
#       # log_term <- lchoose(n1 - y1, k) + k * eta +
#       #   lgamma(y0 + y1 + a_pi + k) + lgamma(n0 - y0 + b_pi) -
#       #   lgamma(n0 + y1 + a_pi + b_pi + k)
#       
#       log_term[k+1] <- lchoose(n1 - y1, k) + k * eta +
#         lbeta(y0 + y1 + a_pi + k, n0 - y0 + b_pi)
#       
#       signs[k+1] <- (-1)^k
#       # term <- (-1)^k * exp(log_term)
#       
#       # # Add the term to the result
#       # integral <- integral + term
#     }
#     log_integral <- signedLogSumExp(signs, log_term)
#     
#     lresult <- lchoose(n0, y0) + lchoose(n1,y1) + eta * y1 - lbeta(a_pi,b_pi)  +  log_integral
#     result <- exp(lresult)
#   }else{
#     
#     integral <- 0
#     # Loop over k from 0 to (n1j - y1j)
#     for (k in 0:(n1 - y1)) {
#       # Compute each term in the summation
#       # log_term <- lchoose(n1 - y1, k) + k * eta +
#       #   lgamma(y0 + y1 + a_pi + k) + lgamma(n0 - y0 + b_pi) -
#       #   lgamma(n0 + y1 + a_pi + b_pi + k)
#       
#       log_term <- lchoose(n1 - y1, k) + k * eta +
#         Ibeta(exp(-eta), y0 + y1 + a_pi + k, n0 - y0 + b_pi, lower = TRUE, log = TRUE)
#       
#       term <- (-1)^k * exp(log_term)
#       
#       # Add the term to the result
#       integral <- integral + term
#     }
#     
#     lresult <- lchoose(n0, y0) + lchoose(n1,y1) + eta * y1 - lbeta(a_pi,b_pi)  + log(integral)
#     result <- exp(lresult)
#   }
#   
#   return(result)
# }

# # Generate a range of eta values to plot
# eta_values <- seq(-5, 5, by = 0.01)
# 
# results <- sapply(eta_values, function(eta) likelihoodEta_inc(n0, y0, n1, y1, a_pi, b_pi, eta))
# # results <- NA
# # for (k in seq_len(length(eta_values))){
# #   results[k] <- likelihoodEta(n0, y0, n1, y1, a_pi, b_pi, eta_values[k])
# # }
# # eta_values[72]
# # likelihoodEta(n0, y0, n1, y1, a_pi, b_pi, eta)
# # eta <- eta_values[72]
# #
# #
# # Plot the likelihood function
# plot(eta_values, results, type = "l", col = "blue", lwd = 2,
#      xlab = expression(eta), ylab = "Likelihood",
#      main = "Likelihood Function with treatment(2/100), control(50/100)",
#      xaxt = "n")
# axis(1, at = seq(min(eta_values), max(eta_values), by = 0.1))