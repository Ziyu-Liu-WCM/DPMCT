source("likelihoodEta_hyp.R")
source("likelihoodEta_inc.R")
source("DPMCT.R")
source("generate_data.R")

J <- nrow(df)
y0 <- df$control_response
n0 <- rep(nControl,J)
y1 <- df$treatment_response
n1 <- rep(nTreatment,J)
### Check when eta = 0
# MCMC start value
eta <- rep(1, J)
mu <- rep(0, J)
tau <- rep(1, J)
alpha <- 100
m <- 5
likelihood_version <- "hyp"

# Prior specification
a_pi <- 1
b_pi <- 1
step_size <- 1
mu_prior_mean <- 0
mu_prior_sd <- 2
a_tau <- 1
b_tau <- 1
a_alpha <- 1
b_alpha <- 1
# j <- 5

nIter  <- 2000
burnIn <- 500

chain_alpha <- numeric(nIter)
chain_mu    <- matrix(NA, nrow = nIter, ncol = J)
chain_tau   <- matrix(NA, nrow = nIter, ncol = J)
chain_eta   <- matrix(NA, nrow = nIter, ncol = J)



for (iter in 1:nIter) {
  # A) Update alpha
  alpha <- update_alpha(alpha, mu, tau, a_alpha, b_alpha)
  
  # B) Update mu, tau (DP clustering) with the current eta
  mu_tau_res  <- update_mu_tau(
    theta         = eta, 
    mu            = mu,
    tau           = tau,
    alpha         = alpha,
    mu_prior_mean = mu_prior_mean,
    mu_prior_sd   = mu_prior_sd,
    a_tau         = a_tau,
    b_tau         = b_tau,
    m             = m
  )
  mu  <- mu_tau_res$mu
  tau <- mu_tau_res$tau
  
  # C) Update eta using binomial data
  eta <- update_eta(n0, y0, n1, y1, eta, a_pi, b_pi, mu, tau, step_sizem, version = likelihood_version)
  
  # D) Store samples
  chain_alpha[iter]    <- alpha
  chain_mu[iter, ]     <- mu
  chain_tau[iter, ]    <- tau
  chain_eta[iter, ]    <- eta
}


idx_keep   <- (burnIn + 1):nIter
alpha_post <- chain_alpha[idx_keep]
mu_post    <- chain_mu[idx_keep, ]
tau_post   <- chain_tau[idx_keep, ]
eta_post   <- chain_eta[idx_keep, ]

mu_post[1500,]
tau_post[1500,]
eta_post[1500,]
alpha_post[1500]