hessian <- function(theta, grad, eps) {
  g<-grad
  n<-length(theta)
  hess<-matrix(nrow=n,ncol=n)
  h<-rep(0, 2)
  h[1] <- eps
  h[2] <- eps
  for (i in 1:n){
    for (j in 1:n){
      e_i=rep(0,n)
      e_i[i]=1
      e_j=rep(0,n)
      e_j[j]=1
      hess[i, j] = (g(theta + h[j]*e_j)[i] - g(theta-h[j]*e_j)[i])/4*eps +  (g(theta+h[i]*e_i)[j] - g(theta-h[i]*e_i)[j])/4*eps
    }
  }
  return(hess)
}

newt <- function(theta,func,grad,hess,...,tol=1e-8,max.half=20,eps=1e-6) {#fscale=1,maxit=100,max.half=20){
  g <- grad(theta)
  if (is.null(hess)){
    hess<-function(theta){hessian(theta,grad,eps)}
  }
  
  #check if the hessian is positive definite. If it isn't, the we perturb it
  if (inherits(try(chol(hess), silent=TRUE), "try-error")) {
    #perturb the hessian if it isn't positive definite
    
    #check that this value for the dk works
    h0 <- hess
    dk <- eps/sqrt(eps)
    hess <- function(theta) {h0(theta) + diag(length(theta))*dk}
  }
  
  
  while (length(g[g > tol]) > 0) {
    #Need to find a way to get rid of the solve
    delta <- solve(hess(theta)) %*% grad(theta)
    
    #step-halving
    #we need to check that the function's value is actually being minimised
    counter <- 0
    while(func(theta + delta) <= func(theta)) {
      if (counter >= max.half) {
        #error goes here
        break
      }
      delta <- delta/2
      counter <- counter + 1
    }
    
    theta <- theta - delta
    
    g <- grad(theta)
  }
  theta
}
rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}
gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}

rc<-function(theta){
  return(theta[1]^2)
}

gc<-function(theta){
  return(c(2*theta[1],0))
}
hc<-function(theta){
  h<-matrix(0,2,2)
  h[1,1]<-2
  return(h)
}

newt(c(1, 1), rc, gc, hc)


