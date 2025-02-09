library(hypergeo)
library(gsl)
source("log-hypergeo.R")
# # Example parameters
# n1 <- 30
# y1 <- 2
# y0 <- 20
# n0 <- 30
# a <- 1
# b <- 1
# a_pi <- 1
# b_pi <- 1
# eta <-2
# likelihoodEta_hyp(n0, y0, n1, y1, a_pi, b_pi, eta)


likelihoodEta_hyp <- function(n0, y0, n1, y1, a, b, eta) {
  if (eta < 0) {
    log_hype_term <- logHyperg2F1(n0 + n1 + a + b, n0 - y0 + b, n0 + y1 + a + b, exp(eta))
    log_likelihood <- 
      y1 * eta + 
      (n0 + n1 - y0 - y1 + b) * log(1 - exp(eta)) +
      lchoose(n0, y0) + 
      lchoose(n1, y1) +
      (lbeta(y0 + y1 + a, n0 - y0 + b) - lbeta(a, b)) +
      log_hype_term
    
  } else if (eta > 0) {
    log_hype_term <- logHyperg2F1(n0 + n1 + a + b, n1 - y1, n1 + y0 + a + 1, exp(-eta))
    log_likelihood <- 
      -eta * (y0 + a) + 
      (n0 + n1 - y0 - y1 + b) * log(1 - exp(-eta)) +
      lchoose(n0, y0) + 
      lchoose(n1, y1) +
      (lbeta(y0 + y1 + a, n1 - y1) - lbeta(a, b)) +
      log_hype_term
  } else {
    log_likelihood <- y1 * eta + lchoose(n0, y0) + 
      lchoose(n1, y1) + lbeta(y0 + y1 + a, n0 + n1 - y0 - y1 + b)
  }
  
  # Return the likelihood in log space, or exponentiate to return in normal space
  return(exp(log_likelihood))  # For log-likelihood
  # return(exp(log_likelihood))  # For regular likelihood
}


# ### Original
# likelihoodEta<-function(n0,y0,n1,y1,a,b,eta){
#   if(eta<0){
#     exp(y1*eta) *
#       choose(n0,y0) * choose(n1,y1) *
#       (beta(y0+y1+a,n0-y0+b)/beta(a,b)) *
#       logHyperg2F1(y1 - n1, y0 + y1 + a, n0 + y1 + a + b, exp(eta))
#     
#   }else if(eta>0){
#     exp(-eta * (y0 + a)) *
#       choose(n0,y0) * choose(n1,y1) *
#       (beta(y0 + y1 + a, n1 - y1)/beta(a,b))*
#       logHyperg2F1(y0 - n0 - b + 1, y0 + y1 + a, n1 + y0 + a + 1, exp(-eta))
#   } else{
#     beta(y0 + y1 + a, n0 + n1 - y0 - y1 + b)
#   }
# }
# 
# 
# hyperg_2F1(y0 - n0 - b + 1, y0 + y1 + a, n1 + y0 + a + 1, exp(-eta))
# ### Kummer Trans
# likelihoodEta<-function(n0,y0,n1,y1,a,b,eta){
#   if(eta<0){
#     exp(y1*eta) * (1 - exp(eta))^(n0 + n1 - y0 - y1 + b) *
#       choose(n0,y0) * choose(n1,y1) *
#       (beta(y0+y1+a,n0-y0+b)/beta(a,b)) *
#       logHyperg2F1(n0 + n1 + a + b, n0 - y0 + b, n0 + y1 + a + b, exp(eta))
# 
#   }else if(eta>0){
#     exp(-eta * (y0 + a)) * (1 - exp(-eta))^(n0 + n1 - y0 - y1 + b) *
#       choose(n0,y0) * choose(n1,y1) *
#       (beta(y0 + y1 + a, n1 - y1)/beta(a,b))*
#       logHyperg2F1(n0 + n1 + a + b, n1 - y1, n1 + y0 + a + 1, exp(-eta))
#   } else{
#     beta(y0 + y1 + a, n0 + n1 - y0 - y1 + b)
#   }
# }
# 
# ### Log Kummer Trans
# likelihoodEta <- function(n0, y0, n1, y1, a, b, eta) {
#   if (eta < 0) {
#     hype_term <- Re(hypergeo::hypergeo(n0 + n1 + a + b, n0 - y0 + b, n0 + y1 + a + b, exp(eta)))
#     hype_term <- ifelse(hype_term < 0, -hype_term, hype_term)
#     hype_term <- ifelse(hype_term == Inf, .Machine$double.xmax, hype_term)
#     log_likelihood <- 
#       y1 * eta + 
#       (n0 + n1 - y0 - y1 + b) * log(1 - exp(eta)) +
#       lchoose(n0, y0) + 
#       lchoose(n1, y1) +
#       (lbeta(y0 + y1 + a, n0 - y0 + b) - lbeta(a, b)) +
#       log(hype_term)
#     
#   } else if (eta > 0) {
#     hype_term <- Re(hypergeo::hypergeo(n0 + n1 + a + b, n1 - y1, n1 + y0 + a + 1, exp(-eta)))
#     hyperg_2F1(n0 + n1 + a + b,n1 - y1,n1 + y0 + a + 1,exp(-eta))
#     hype_term <- ifelse(hype_term < 0, -hype_term, hype_term)
#     hype_term <- ifelse(hype_term == Inf, .Machine$double.xmax, hype_term)
#     log_likelihood <- 
#       -eta * (y0 + a) + 
#       (n0 + n1 - y0 - y1 + b) * log(1 - exp(-eta)) +
#       lchoose(n0, y0) + 
#       lchoose(n1, y1) +
#       (lbeta(y0 + y1 + a, n1 - y1) - lbeta(a, b)) +
#       log(hype_term)
#     
#   } else {
#     log_likelihood <- lbeta(y0 + y1 + a, n0 + n1 - y0 - y1 + b)
#   }
#   
#   # Return the likelihood in log space, or exponentiate to return in normal space
#   return(log_likelihood)  # For log-likelihood
#   # return(exp(log_likelihood))  # For regular likelihood
# }

# logHyperg2F1(2002, 950, 1052, 0.11)
  
  



# a <- 2002
# b <- 1000
# c <- 1050
# z <- 0.4
# log_hypergeo <- function(a, b, c, z) {
#   # Check if z is near 1 to avoid divergence
#   if (z == 1) {
#     stop("Hypergeometric function diverges for z = 1.")
#   }
#   
#   # Compute the hypergeometric series in log-space
#   max_terms <- 500  # Limit the series to prevent infinite summation
#   log_sum <- 0
#   for (k in 0:max_terms) {
#     term_log <- lfactorial(a + k - 1) - lfactorial(a - 1) +
#       lfactorial(b + k - 1) - lfactorial(b - 1) -
#       lfactorial(c + k - 1) + lfactorial(c - 1) -
#       lfactorial(k) +
#       k * log(z)
#     
#     # Add the terms in log-space
#     if (k == 0) {
#       log_sum <- term_log
#     } else {
#       log_sum <- log(exp(log_sum) + exp(term_log))
#     }
#     
#     # Stop if terms are vanishingly small
#     if (abs(exp(term_log)) < 1e-12) break
#   }
#   
#   return(log_sum)
# }
# 
# ?hypergeo
# hyperg_2F1
# # 
# # likelihoodEta(n0,y0,n1,y1,a,b,4)
# # likelihoodEta(n0,y0,n1,y1,a,b,0.1)
# # 
# # 
# # hypergeo(y1-n1,y0+y1+a,n0+y1+a+b,-4)
# # 
# # 
# # n0 <- 100
# # y0 <- 50
# # n1 <- 100
# # y1 <- 2
# # a <- 1
# # b <- 1
# # 
# Generate a range of eta values
# eta_values <- seq(-5, 5, by = 0.01)
# # eta_values <- c(seq(-5, -0.2, by = 0.10), seq(0.2, 5, by = 0.10))
# results <- sapply(eta_values, function(eta) likelihoodEta_hyp(n0, y0, n1, y1, a_pi, b_pi, eta))
# #
# #
# # Plot the likelihood function
# plot(eta_values, results, type = "l", col = "blue", lwd = 2,
#      xlab = expression(eta), ylab = "Likelihood",
#      main = "Likelihood Function with treatment(2/100), control(50/100)",
#      xaxt = "n")
# axis(1, at = seq(min(eta_values), max(eta_values), by = 0.1))
# 
# likelihoodEta(n0, y0, n1, y1, a, b, 2)
# 
# log(1.25)
# 
# log(6)
