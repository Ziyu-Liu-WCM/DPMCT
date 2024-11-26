
update_mu_tau_alpha <- function(theta, mu, tau, alpha, J,
                          mu_prior_mean, mu_prior_sd,
                          a_tau, b_tau, a_alpha, b_alpha,
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
    
    # Identify clusters and calculate cluster sizes
    unique_clusters <- unique(data.frame(mu = mu[-j], tau = tau[-j]))
    num_unique_clusters <- nrow(unique_clusters)
    
    # Calculate log-probabilities for existing clusters
    for (k in 1:num_unique_clusters) {
      mu_c <- unique_clusters$mu[k]
      tau_c <- unique_clusters$tau[k]
      
      # Number of data points in cluster c excluding theta[j]
      n_c_minus_j <- sum(mu[-j] == mu_c & tau[-j] == tau_c)
      
      # Only consider clusters with at least one data point
      if (n_c_minus_j > 0) {
        # Find an index of a data point in cluster c (excluding j)
        j1 <- which((mu == mu_c & tau == tau_c) & (1:J != j))[1]
        
        # Compute the log-probability
        prob_c[j1] <- log(n_c_minus_j) + dnorm(theta[j], mu_c, tau_c, log = TRUE)
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
      # Draw tau from Gamma distribution
      repeat {
        tau_pro_temp <- rgamma(1, a_tau, b_tau)
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
      logpost_cur <- logpost_cur + dgamma(tau_j, a_tau, b_tau, log = TRUE) + log(tau_j)
      logpost_pro <- logpost_pro + dgamma(tau_j_pro, a_tau, b_tau, log = TRUE) + log(tau_j_pro)
      
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
  
  # Step 3: Update alpha
  # Compute the number of clusters K
  cluster_labels <- unique(data.frame(mu = mu, tau = tau))
  K <- nrow(cluster_labels)
  
  # Sample eta from Beta distribution
  eta <- rbeta(1, alpha + 1, J)
  
  # Compute s based on a threshold
  epsilon <- 1e-6  # Small constant to prevent division by zero
  s <- (a_alpha + K - 1) / (J * (b_alpha - log(eta) + epsilon))
  pi_eta <- s / (1 + s)
  
  # Sample from Bernoulli to decide on new shape parameter
  u <- runif(1)
  if (u < pi_eta) {
    a_alpha_new <- a_alpha + K
  } else {
    a_alpha_new <- a_alpha + K - 1
  }
  
  b_alpha_new <- b_alpha - log(eta)
  
  # Sample new alpha from Gamma distribution
  alpha <- rgamma(1, a_alpha_new, b_alpha_new)
  
  # Return the updated mu and tau
  return(list(mu = mu, tau = tau, alpha = alpha))
}



# Example data
set.seed(123)
J <- 100
theta <- rnorm(J, mean = 5, sd = 2)

# Initial values
mu <- rnorm(J, mean = 5, sd = 1)
tau <- rgamma(J, shape = 2, rate = 1)
alpha <- rgamma(1, shape = 1, rate = 1)
mu_prior_mean <- 0
mu_prior_sd <- 10
a_tau <- 2
b_tau <- 1
a_alpha <- 1
b_alpha <- 1
m <- 5
sd_prop <- 0.5


result <- update_mu_tau_alpha(theta, mu, tau, alpha, J, mu_prior_mean, mu_prior_sd,
                              a_tau, b_tau, a_alpha, b_alpha, m)
