update_eta <- function(n0, y0, n1, y1, eta, 
                       a_pi, b_pi, mu, tau, step_size, version = "hyp"){
  if(version == "hyp"){
    likelihoodEta <- likelihoodEta_hyp
  } else{
    likelihoodEta <- likelihoodEta_inc
  }
  J <- length(y0)
  if(length(y1) != J | length(n0) != J | length(n1) != J | length(mu) != J | length(tau) != J) stop("y0, y1, n0 and n1 must be of equal length!")
  
  for (j in 1:J) {
    # Current value of eta at position j
    eta_j_cur <- eta[j]
    # Propose a new eta value from N(mu[j], tau[j]^2)
    eta_j_pro <- rnorm(1, mean = eta_j_cur, sd = step_size)
    
    ## Prior Probability
    logpost_cur <- dnorm(eta_j_cur, mean = mu[j], sd = tau[j], log = TRUE)
    logpost_pro <- dnorm(eta_j_pro, mean = mu[j], sd = tau[j], log = TRUE)
    
    ## Likelihood
    logpost_cur <- logpost_cur + log(likelihoodEta(n0[j], y0[j], n1[j], y1[j], a_pi, b_pi, eta_j_cur))
    logpost_pro <- logpost_pro + log(likelihoodEta(n0[j], y0[j], n1[j], y1[j], a_pi, b_pi, eta_j_pro))
    
    if (is.na(logpost_pro - logpost_cur)){
      next
    }
    # Metropolis-Hastings acceptance step
    u <- runif(1)
    if (log(u) < (logpost_pro - logpost_cur)) {
      eta_j_new <- eta_j_pro
    } else {
      eta_j_new <- eta_j_cur
    }
    eta[j] <- eta_j_new
  }
  
  return(eta)
}





update_mu_tau <- function(theta, mu, tau, alpha,
                          mu_prior_mean, mu_prior_sd,
                          a_tau, b_tau, m) {
  
  n <- length(theta)
  
  for (i in sample(seq_len(n))) {
    
    old_mu  <- mu[i]
    old_tau <- tau[i]
    
    # Remove i from its old cluster by "temporarily" clearing out mu[i], tau[i].
    # We'll reassign it below.  For checking if the old cluster was a singleton,
    # see how many other j != i share the same (old_mu, old_tau).
    mu[i]  <- NA
    tau[i] <- NA
    
    # Check if the old cluster is now empty
    # i.e., no other j has (mu[j], tau[j]) = (old_mu, old_tau).
    is_singleton <- !any((mu == old_mu & tau == old_tau), na.rm=TRUE)
    
    # Gather unique clusters among the other data points
    df_others <- data.frame(mu = mu, tau = tau)
    df_others <- df_others[complete.cases(df_others), ]  
    # Above excludes the row i, which is now NA.
    unique_clusters <- unique(df_others)
    K <- nrow(unique_clusters)
    
    #-----------------------------------
    # Build new proposals for i
    #-----------------------------------
    mu_pro  <- numeric(m)
    tau_pro <- numeric(m)
    
    if (is_singleton) {
      # Recycle i's old parameter as first new proposal:
      mu_pro[1]  <- old_mu
      tau_pro[1] <- old_tau
      
      # Draw remaining m-1 from base distribution G0
      for (mm in 2:m) {
        mu_pro[mm]  <- rnorm(1, mu_prior_mean, mu_prior_sd)
        tau_pro[mm] <- rgamma(1, shape=a_tau, rate=b_tau)
      }
    } else {
      # i's old cluster is not a singleton => standard approach
      for (mm in 1:m) {
        mu_pro[mm]  <- rnorm(1, mu_prior_mean, mu_prior_sd)
        tau_pro[mm] <- rgamma(1, shape=a_tau, rate=b_tau)
      }
    }
    
    #-------------------------------------
    # Compute unnormalized log probabilities
    #-------------------------------------
    # 1) for existing clusters
    log_prob <- numeric(K + m)
    for (k in seq_len(K)) {
      cluster_mu  <- unique_clusters$mu[k]
      cluster_tau <- unique_clusters$tau[k]
      # how many points are in that cluster among the "others"?
      nk <- sum(df_others$mu == cluster_mu & df_others$tau == cluster_tau)
      # unnormalized log probability = log(n_k) + log likelihood
      log_prob[k] <- log(nk) + dnorm(theta[i], mean=cluster_mu, sd=cluster_tau, log=TRUE)
    }
    # 2) for new proposals
    for (mm in seq_len(m)) {
      log_prob[K + mm] <- log(alpha/m) + 
        dnorm(theta[i], mean=mu_pro[mm], sd=tau_pro[mm], log=TRUE)
    }
    
    #-------------------------------------
    # Sample a new cluster for data i
    #-------------------------------------
    lp_max <- max(log_prob)
    w <- exp(log_prob - lp_max)
    w <- w / sum(w)
    pick <- which(cumsum(w) > runif(1))[1]
    
    # Assign i
    if (pick <= K) {
      # existing cluster
      mu[i]  <- unique_clusters$mu[pick]
      tau[i] <- unique_clusters$tau[pick]
    } else {
      # new cluster
      mu[i]  <- mu_pro[pick - K]
      tau[i] <- tau_pro[pick - K]
    }
  }  # end loop over i
  
  # Step 2: Update distinct mu and tau values
  mu_tau_updated <- rep(0, J)
  sd_prop <- 0.5
  
  for (j in sample(seq_len(n))) {
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


































# eta_current <- eta
# mu_current <- mu
# tau_current <- tau
# alpha_current <- alpha
# num_iter <- 10000
# 
# eta_samples <- matrix(NA, nrow = num_iter, ncol = J)
# mu_samples <- matrix(NA, nrow = num_iter, ncol = J)
# tau_samples <- matrix(NA, nrow = num_iter, ncol = J)
# alpha_samples <- numeric(num_iter)
# 
# # for (iter in 1:num_iter) {
# #   eta_current <- update_eta(n0, y0, n1, y1, eta_current, a_pi, b_pi, mu_current, tau_current)
# #   eta_samples[iter, ] <- eta_current
# # }
# 
# for (iter in 1:num_iter) {
#   # Update eta
#   eta_current <- update_eta(n0, y0, n1, y1, eta_current, a_pi, b_pi, mu_current, tau_current)
#   
#   # Update mu and tau
#   mu_tau_updated <- update_mu_tau(eta_current, mu_current, tau_current, alpha_current,
#                                   mu_prior_mean, mu_prior_sd, a_tau, b_tau, m)
#   mu_current <- mu_tau_updated$mu
#   tau_current <- mu_tau_updated$tau
#   
#   # Update alpha
#   alpha_current <- update_alpha(alpha_current, mu_current, tau_current, a_alpha, b_alpha)
#   
#   # Store samples
#   eta_samples[iter, ] <- eta_current
#   mu_samples[iter, ] <- mu_current
#   tau_samples[iter, ] <- tau_current
#   alpha_samples[iter] <- alpha_current
# }
# 
# burn_in <- 9000
# 
# posterior_samples <- eta_samples[(burn_in : num_iter),1]
# 
# posterior_mean <- mean(posterior_samples)
# posterior_mode <- as.numeric(density(posterior_samples)$x[which.max(density(posterior_samples)$y)])
# 
# 
# hist(posterior_samples, breaks = 50, probability = TRUE, 
#      main = "Posterior Distribution with treatment(6/10), control(3/10)", xlab = "Parameter Value", 
#      col = "lightblue", border = "white")
# 
# # Overlay density plot
# lines(density(posterior_samples), col = "blue", lwd = 2)
# 
# # Add vertical lines for mean and mode
# abline(v = posterior_mean, col = "red", lwd = 2, lty = 2)  # Mean
# abline(v = posterior_mode, col = "green", lwd = 2, lty = 2)  # Mode
# 
# # Add legend
# legend("topright", legend = c("Mean", "Mode"), 
#        col = c("red", "green"), lty = 2, lwd = 2, bty = "n")
# 
# # Annotate the mean and mode
# text(posterior_mean, 0.05, labels = paste0("Mean = ", round(posterior_mean, 2)), col = "red", pos = 4)
# text(posterior_mode, 0.04, labels = paste0("Mode = ", round(posterior_mode, 2)), col = "green", pos = 4)

