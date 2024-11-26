update_mu_tau <- function(theta, mu, tau, alpha, J,
                          mu_prior_mean, mu_prior_sd,
                          tau_prior_loc, tau_prior_scale,
                          m) {
  # Initialize proposal vectors
  mu_pro <- numeric(m)
  tau_pro <- numeric(m)
  prob_c <- numeric(J + m)
  
  # Step 1: Update class labels of object j
  for (j in 1:J) {
    is_unique_j <- 1
    max_logprob_c <- -Inf
    sum_prob_c <- 0
    
    mu_j <- mu[j]
    tau_j <- tau[j]
    
    mu_pro[] <- 0
    tau_pro[] <- 0
    prob_c[] <- 0
    
    # Check if j is a singleton (unique cluster)
    for (j1 in 1:J) {
      if (j1 != j && mu[j1] == mu_j && tau[j1] == tau_j) {
        is_unique_j <- 0
        break
      }
    }
    
    # Compute log-probabilities for existing clusters
    for (j1 in 1:J) {
      if (j1 != j) {
        prob_c[j1] <- dnorm(theta[j], mu[j1], tau[j1], log = TRUE)
        if (max_logprob_c < prob_c[j1]) max_logprob_c <- prob_c[j1]
      }
    }
    
    # If j is a singleton, keep current value and adjust probabilities
    if (is_unique_j == 1) {
      mu_pro[1] <- mu[j]
      tau_pro[1] <- tau[j]
      prob_c[J + 1] <- dnorm(theta[j], mu_pro[1], tau_pro[1], log = TRUE)
      prob_c[J + 1] <- prob_c[J + 1] + log(alpha) - log(m)
      if (max_logprob_c < prob_c[J + 1]) max_logprob_c <- prob_c[J + 1]
    }
    
    # Draw proposals from priors for new clusters
    for (m1 in (is_unique_j + 1):m) {
      mu_pro[m1] <- rnorm(1, mu_prior_mean, mu_prior_sd)
      # Draw tau from half-Cauchy distribution
      repeat {
        tau_pro_temp <- rcauchy(1, tau_prior_loc, tau_prior_scale)
        if (tau_pro_temp > 0) break
      }
      tau_pro[m1] <- tau_pro_temp
      
      prob_c[J + m1] <- dnorm(theta[j], mu_pro[m1], tau_pro[m1], log = TRUE)
      prob_c[J + m1] <- prob_c[J + m1] + log(alpha) - log(m)
      if (max_logprob_c < prob_c[J + m1]) max_logprob_c <- prob_c[J + m1]
    }
    
    # Normalize probabilities
    for (j1 in 1:(J + m)) {
      if (j1 != j) {
        prob_c[j1] <- prob_c[j1] - max_logprob_c
        prob_c[j1] <- exp(prob_c[j1])
        sum_prob_c <- sum_prob_c + prob_c[j1]
      }
    }
    prob_c[j] <- 0  # Exclude current index
    
    # Sample new class label
    u <- runif(1)
    cumulative_prob <- cumsum(prob_c / sum_prob_c)
    j1 <- which(u <= cumulative_prob)[1]
    
    # Assign new mu and tau based on sampled class label
    if (j1 <= J) {
      mu[j] <- mu[j1]
      tau[j] <- tau[j1]
    } else {
      mu[j] <- mu_pro[j1 - J]
      tau[j] <- tau_pro[j1 - J]
    }
  }
  
  # Step 2: Update distinct mu and tau values
  mu_tau_updated <- rep(0, J)
  sd_prop <- 0.5
  
  for (j in 1:J) {
    if (mu_tau_updated[j] == 0) {
      mu_j <- mu[j]
      tau_j <- tau[j]
      
      # Propose new tau using log-normal proposal distribution
      norm_prop <- rnorm(1, 0.0, sd_prop)
      tau_j_pro <- tau_j * exp(norm_prop)
      
      # Find indices of data points in the same cluster
      same_cluster_indices <- which(mu == mu_j & tau == tau_j)
      n_c <- length(same_cluster_indices)
      theta_c_sum <- sum(theta[same_cluster_indices])
      
      # Sample new mu from conjugate normal posterior
      mu_post_var <- 1 / ((1 / mu_prior_sd^2) + (n_c / tau_j^2))
      mu_post_mean <- mu_post_var * ((mu_prior_mean / mu_prior_sd^2) + (theta_c_sum / tau_j^2))
      mu_j_new <- rnorm(1, mu_post_mean, sqrt(mu_post_var))
      
      # Calculate log-posterior probabilities for current and proposed tau
      logpost_cur <- sum(dnorm(theta[same_cluster_indices], mu_j_new, tau_j, log = TRUE))
      logpost_pro <- sum(dnorm(theta[same_cluster_indices], mu_j_new, tau_j_pro, log = TRUE))
      
      # Add prior distributions and Jacobian adjustments
      logpost_cur <- logpost_cur + dcauchy(tau_j, tau_prior_loc, tau_prior_scale, log = TRUE) + log(tau_j)
      logpost_pro <- logpost_pro + dcauchy(tau_j_pro, tau_prior_loc, tau_prior_scale, log = TRUE) + log(tau_j_pro)
      
      # Metropolis-Hastings acceptance step
      u <- runif(1)
      if (log(u) < (logpost_pro - logpost_cur)) {
        tau_j_new <- tau_j_pro
      } else {
        tau_j_new <- tau_j
      }
      
      # Update mu and tau for all points in the cluster
      mu[same_cluster_indices] <- mu_j_new
      tau[same_cluster_indices] <- tau_j_new
      mu_tau_updated[same_cluster_indices] <- 1
    }
  }
  
  # Return the updated mu and tau
  return(list(mu = mu, tau = tau))
}


# Example data
set.seed(123)
J <- 100
theta <- rnorm(J, mean = 5, sd = 2)
mu <- rnorm(J, mean = 5, sd = 1)
tau <- runif(J, min = 1, max = 3)
alpha <- 1
mu_prior_mean <- 0
mu_prior_sd <- 10
tau_prior_loc <- 0
tau_prior_scale <- 1
a_tau <- 
m <- 5

# Update mu and tau
result <- update_mu_tau(theta, mu, tau, alpha, J,
                        mu_prior_mean, mu_prior_sd,
                        tau_prior_loc, tau_prior_scale,
                        m)

# Access updated mu and tau
updated_mu <- result$mu
updated_tau <- result$tau
