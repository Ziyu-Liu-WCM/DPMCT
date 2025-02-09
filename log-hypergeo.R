# -----------------------------------------------
# logHyperg2F1(a, b, c, z, ...)
#
# Computes log( 2F1(a, b; c; z) ) by summing the series
#   sum_{n=0 to Inf} [ (a)_n (b)_n / (c)_n ] * (z^n / n!)
# in log form to avoid overflow when terms get large.
#
# The partial sums are accumulated in log-space.
# Convergence is checked via the difference of consecutive partial sums.
#
# Arguments:
#   a, b, c    : (Possibly real) parameters of the hypergeometric function.
#   z          : Complex or real argument with |z| < 1 for series convergence.
#   maxit      : Maximum number of terms in the series.
#   tol        : Tolerance for convergence on the *log-scale* partial sums.
#
# Returns:
#   A numeric value (double) giving log( 2F1(a,b;c;z) ), if convergent.
#   Possibly NA if the iteration does not converge (or if c is a negative integer).
# -----------------------------------------------


logHyperg2F1 <- function(a, b, c, z, maxit = 1e+6, tol = 1e-14) {
  # If c is a non-positive integer, the series definition breaks down
  # (poles in the gamma function). We'll do a basic check:
  if (c <= 0 && abs(c - round(c)) < 1e-15) {
    warning("Parameter 'c' is a non-positive integer; 2F1 may be undefined or require a limiting process.")
    return(NA_real_)
  }
  
  # For |z| >= 1, the naive series expansion does not converge in general.
  if (Mod(z) > 1) {
    warning("The power series for 2F1 diverges or is not valid for |z| >= 1.")
    return(NA_real_)
  }
  
  # A small helper for log-sum-exp:
  logspace_add <- function(x, y) {
    # log( e^x + e^y ) = max(x,y) + log(1 + exp(-|x-y|))
    m <- max(x, y)
    return(m + log1p(exp(-abs(x - y))))
  }
  
  # We'll build the sum term-by-term in log space.
  # 2F1(a,b;c;z) = sum_{n=0 to inf} T(n),
  # where T(n) = [ (a)_n (b)_n / (c)_n ] * (z^n / n! ).
  #
  # We accumulate partial sums in log form: partial_log = log( sum of terms ).
  # The 0th term is T(0) = 1, so its log is 0.
  partial_log <- 0.0  # log(1) = 0
  old_partial_log <- -Inf
  
  # We'll track the log of the nth term. Let:
  #   T(0) = 1 => logTerm = log(1) = 0
  logTerm <- 0.0
  
  for (n in seq_len(maxit)) {
    # Add the current term's log to the partial sum in log-space
    partial_log <- logspace_add(partial_log, logTerm)
    
    # Check for convergence in log-space:
    #   if partial_log hasn't changed by more than `tol`,
    #   we consider we’ve converged. More precisely,
    #   we check if difference is < tol in linear scale:
    #   difference in linear scale = exp(partial_log) - exp(old_partial_log).
    #
    # But we only have logs. A simple rule of thumb is:
    #   if abs(partial_log - old_partial_log) < tol
    #   then we assume it's converged enough on a log scale.
    if (abs(partial_log - old_partial_log) < tol/2) {
      break
    }
    old_partial_log <- partial_log
    
    # Now compute logTerm(n+1) from logTerm(n) in a *recursive* way:
    #
    # (a)_n+1 = (a)_n * (a + n)
    # so in log: log((a)_n+1) = log((a)_n) + log(a + n).
    #
    # T(n+1) / T(n) = [ (a + n)(b + n) / ((c + n)*(n+1)) ] * z
    # => logTerm_{n+1} = logTerm_n + log(a + n) + log(b + n)
    #                           - log(c + n) - log(n+1) + log(z).
    logTerm <- logTerm +
      log(a + (n - 1)) +
      log(b + (n - 1)) -
      log(c + (n - 1)) -
      log(n) +
      log(z)
  }
  
  # After the loop, we must add the final term as well
  # (the loop as written adds T(n) before computing T(n+1)).
  partial_log <- logspace_add(partial_log, logTerm)
  
  if (abs(partial_log - old_partial_log) > tol) {
    warning("logHyperg2F1 did not fully converge within maxit = ", maxit,
            " iterations (log difference ~ ", 
            format(abs(partial_log - old_partial_log)), ").")
  }
  
  return(partial_log)
}