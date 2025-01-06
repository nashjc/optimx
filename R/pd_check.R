pd_check <- function(A, pivot=TRUE, tol=1.e-07) {
# checks whether the input (symmetric) matrix is positive-definite
# Checks for errors or warnings
  result <- tryCatch({
    chol(A, pivot=pivot, tol=tol)   # Attempt Cholesky decomposition
    TRUE      # If successful, return TRUE
  }, warning = function(w) {
    FALSE      # Return FALSE if there's a warning
  }, error = function(e) {
    FALSE      # Return FALSE if there's an error
  })
  return(result)
}
