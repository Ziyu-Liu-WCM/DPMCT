library(hypergeo)

n0<-100
n1<-100
y0<-50
y1<-55

a<-1
b<-1

eta<--1
exp(-log(1.1))



# Define the function to calculate the summation
likelihoodEta <- function(n1j, y1j, y0j, n0j, a_pi_j, b_pi_j, e_pi_j) {
  # Initialize the summation result
  result <- 0
  
  # Loop over k from 0 to (n1j - y1j)
  for (k in 0:(n1j - y1j)) {
    # Compute each term in the summation
    term <- choose(n1j - y1j, k) * 
      (-e_pi_j)^k * 
      (gamma(y0j + y1j + a_pi_j + k) * gamma(n0j - y0j + b_pi_j) /
         gamma(y1j + a_pi_j + k + n0j + b_pi_j))
    
    # Add the term to the result
    result <- result + term
  }
  
  return(result)
}

# Example parameters
n1j <- 10
y1j <- 3
y0j <- 4
n0j <- 12
a_pi_j <- 2
b_pi_j <- 3
e_pi_j <- 0.5

# Call the function with example inputs
result <- calculate_formula(n1j, y1j, y0j, n0j, a_pi_j, b_pi_j, e_pi_j)

# Print the result
cat("The result of the summation is:", result, "\n")



































# likelihoodEta<-function(n0,y0,n1,y1,a,b,eta){
#   if(eta<=0){
#     # choose(n0,y0)*choose(n1,y1)*exp(y1*eta)*
#     #   beta(y0+y1+a,n0-y0+b)/beta(a,b)*
#     #   hypergeo((y1-n1),y0+y1+a,n0+y1+a+b,exp(eta))
#     exp(y1*eta) * (1 - exp(eta))^(n0 + n1 - y0 - y1 + b) *
#       choose(n0,y0) * choose(n1,y1) * 
#       (beta(y0+y1+a,n0-y0+b)/beta(a,b)) *
#       hypergeo(n0 + n1 + a + b, n0 - y0 + b, n0 + y1 + a + b, exp(eta))
#     
#   }else{
#     # choose(n0,y0)*choose(n1,y1)*
#     #   beta(y0+y1+a,n1-y1+b)/exp(y1*eta)/beta(a,b)*
#     #   hypergeo(-(y0-n0),y0+y1+a,n1+y0+a+b,exp(-eta))
#     exp(-y1*eta) * (1 - exp(-eta))^(n0 + n1 - y0 - y1 + b) *
#       choose(n0,y0) * choose(n1,y1) * 
#       (beta(y0+y1+a,n1-y1+b)/beta(a,b))*
#       hypergeo(n0 + n1 + a + b, n1 - y1 + b, n1 + y0 + a + b, exp(-eta))
#   }
# }
# 
# likelihoodEta(n0,y0,n1,y1,a,b,4)
# likelihoodEta(n0,y0,n1,y1,a,b,0.1)
# 
# 
# hypergeo(y1-n1,y0+y1+a,n0+y1+a+b,-4)
# 
# 
# n0 <- 100
# y0 <- 50
# n1 <- 100
# y1 <- 2
# a <- 1
# b <- 1
# 
# # Generate a range of eta values
# eta_values <- c(seq(-21.1, -0.1, length.out = 1000), seq(0.1, 21.1, length.out = 1000))
# 
# # Calculate likelihood for each eta
# likelihood_values <- sapply(eta_values, function(eta) {
#   tryCatch(likelihoodEta(n0, y0, n1, y1, a, b, eta), error = function(e) NA)
# })
# 
# 
# # Plot the likelihood function
# plot(eta_values, likelihood_values, type = "l", col = "blue", lwd = 2,
#      xlab = expression(eta), ylab = "Likelihood",
#      main = "Likelihood Function with treatment(2/100), control(50/100)")
# 
# 
# likelihoodEta(n0, y0, n1, y1, a, b, 2)


