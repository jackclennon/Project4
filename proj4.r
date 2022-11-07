newt <- function(theta, func, grad, hess=NULLL, ..., tol=1e-8) {
  g <- grad(theta)
  while (length(g[g <= tol]) > 0) {
    print(theta)
    #Need to find a way to get rid of the solve
    theta <- theta - solve(hess(theta)) %*% grad(theta)
    g <- grad(theta)
  }
  theta
}

#Example function, with gradient and hessian

func <- function(theta) {
  log(theta[1]) + theta[2] - theta[1]*exp(theta[2])
}

grad <- function(theta) {
  c(1/theta[1] - exp(theta[1]), 1-theta[1]*exp(theta[2]))
}

hess <- function(theta) {
  v <- c(-theta[1]^-2, -exp(theta[2]), -exp(theta[2]), -theta[1]*exp(theta[2]))
  matrix(v, 2, 2)
}