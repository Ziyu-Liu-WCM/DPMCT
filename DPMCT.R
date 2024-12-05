source("likelihoodEta.R")

J <- 5
y0 <- c(30,30,40,30,50)
y1 <- c(30,40,50,60,70)
n0 <- rep(100, J)
n1 <- rep(100, J)
### Check when eta = 0
# MCMC start value
eta <- rep(-1, J)
mu <- rep(0, J)
tau <- rep(1, J)
alpha <- 2
m <- 5

# Prior specification
a_pi <- 1
b_pi <- 1
mu_prior_mean <- 0
mu_prior_sd <- 2
a_tau <- 1
b_tau <- 1
a_alpha <- 1
b_alpha <- 1
# j <- 1

update_eta(n0, y0, n1, y1, eta_current, a_pi, b_pi, mu_current, tau_current)


update_eta <- function(n0, y0, n1, y1, eta, 
                       a_pi, b_pi, mu, tau){
  J <- length(y0)
  if(length(y1) != J | length(n0) != J | length(n1) != J | length(mu) != J | length(tau) != J) stop("y0, y1, n0 and n1 must be of equal length!")
  
  for (j in 1:J) {
    # Current value of eta at position j
    eta_j_cur <- eta[j]
    # Propose a new eta value from N(mu[j], tau[j]^2)
    eta_j_pro <- rnorm(1, mean = mu[j], sd = tau[j])
    
    ## Prior Probability
    logpost_cur <- dnorm(eta_j_cur, mean = mu[j], sd = tau[j], log = TRUE)
    logpost_pro <- dnorm(eta_j_pro, mean = mu[j], sd = tau[j], log = TRUE)
    
    ## Likelihood
    logpost_cur <- logpost_cur + log(Re(likelihoodEta(n0[j], y0[j], n1[j], y1[j], a_pi, b_pi, eta_j_cur)))
    logpost_pro <- logpost_pro + log(Re(likelihoodEta(n0[j], y0[j], n1[j], y1[j], a_pi, b_pi, eta_j_pro)))
    
    # Metropolis-Hastings acceptance step
    u <- runif(1)
    if (log(u) < (logpost_pro - logpost_cur)) {
      eta_j_new <- eta_j_pro
    } else {
      eta_j_new <- eta_j_cur
    }
    eta[j] <- eta_j_new
    print(j)
  }
  
  return(eta)
}





update_mu_tau <- function(theta, mu, tau, alpha,
                          mu_prior_mean, mu_prior_sd,
                          a_tau, b_tau, m) {
  
  J <- length(mu)
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
  
  # Return the updated mu and tau
  return(list(mu = mu, tau = tau))
}




update_alpha <- function(alpha, mu, tau, a_alpha, b_alpha){
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
  
  return(alpha)
}


## Single iteration
# eta_updated <- update_eta(n0, y0, n1, y1, eta, a_pi, b_pi, mu, tau)
# mu_tau_updated <- update_mu_tau(eta_updated, mu, tau, alpha,
#                             mu_prior_mean, mu_prior_sd,
#                             a_tau, b_tau, m)
# mu_updated <- mu_tau_updated$mu
# tau_updated <- mu_tau_updated$tau
# alpha_updated <- update_alpha(alpha, mu_updated, tau_updated, a_alpha, b_alpha)


eta_current <- eta
mu_current <- mu
tau_current <- tau
alpha_current <- alpha
num_iter <- 10000

eta_samples <- matrix(NA, nrow = num_iter, ncol = J)
mu_samples <- matrix(NA, nrow = num_iter, ncol = J)
tau_samples <- matrix(NA, nrow = num_iter, ncol = J)
alpha_samples <- numeric(num_iter)

# for (iter in 1:num_iter) {
#   eta_current <- update_eta(n0, y0, n1, y1, eta_current, a_pi, b_pi, mu_current, tau_current)
#   eta_samples[iter, ] <- eta_current
# }

for (iter in 1:num_iter) {
  # Update eta
  eta_current <- update_eta(n0, y0, n1, y1, eta_current, a_pi, b_pi, mu_current, tau_current)
  
  # Update mu and tau
  mu_tau_updated <- update_mu_tau(eta_current, mu_current, tau_current, alpha_current,
                                  mu_prior_mean, mu_prior_sd, a_tau, b_tau, m)
  mu_current <- mu_tau_updated$mu
  tau_current <- mu_tau_updated$tau
  
  # Update alpha
  alpha_current <- update_alpha(alpha_current, mu_current, tau_current, a_alpha, b_alpha)
  
  # Store samples
  eta_samples[iter, ] <- eta_current
  mu_samples[iter, ] <- mu_current
  tau_samples[iter, ] <- tau_current
  alpha_samples[iter] <- alpha_current
}

burn_in <- 9000

posterior_samples <- eta_samples[(burn_in : num_iter),1]

posterior_mean <- mean(posterior_samples)
posterior_mode <- as.numeric(density(posterior_samples)$x[which.max(density(posterior_samples)$y)])


hist(posterior_samples, breaks = 50, probability = TRUE, 
     main = "Posterior Distribution with treatment(6/10), control(3/10)", xlab = "Parameter Value", 
     col = "lightblue", border = "white")

# Overlay density plot
lines(density(posterior_samples), col = "blue", lwd = 2)

# Add vertical lines for mean and mode
abline(v = posterior_mean, col = "red", lwd = 2, lty = 2)  # Mean
abline(v = posterior_mode, col = "green", lwd = 2, lty = 2)  # Mode

# Add legend
legend("topright", legend = c("Mean", "Mode"), 
       col = c("red", "green"), lty = 2, lwd = 2, bty = "n")

# Annotate the mean and mode
text(posterior_mean, 0.05, labels = paste0("Mean = ", round(posterior_mean, 2)), col = "red", pos = 4)
text(posterior_mode, 0.04, labels = paste0("Mode = ", round(posterior_mode, 2)), col = "green", pos = 4)
