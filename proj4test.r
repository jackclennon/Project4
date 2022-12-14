hessian1 <- function(theta, grad, eps,...) {
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
      hess[i, j] = ((g(theta + eps*e_j,...)[i] - g(theta-eps*e_j,...)[i])/(4*eps) +  ((g(theta+eps*e_i,...)[j] - g(theta-eps*e_i,...)[j])/(4*eps)))
    }
  }
  return(hess)
}

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,max.half=20,eps=1e-6) {#fscale=1,maxit=100,max.half=20){
  g <- grad(theta,...)
  if (is.null(hess)){
    hess<-function(theta,...){hessian1(theta,grad,eps,...)}
  }
  if(is.infinite(func) || all(is.finite(grad))==FALSE)
    stop(' The objective is not finite at the initial theta')
  
  counter2<-0
  while (length(g[abs(g) > tol]) > 0) { #breaks when there is convergence
    counter2<-counter2+1
    #check if the hessian is positive definite. If it isn't, the we perturb it
    multiple<-1
    H<-hess(theta,...)
    
    while (inherits(try(chol(H), silent=TRUE), "try-error")) {
      #perturb the hessian if it isn't positive definite
      print('converting non-pos to pos')
      print(norm(H, type = 'F'))
      H <- H + diag(length(theta))*(abs(min(eigen(H)$values))+1)
      multiple<-multiple+10
    }
    #Need to find a way to get rid of the solve
    delta <- -(solve(H) %*% grad(theta,...))
    
    #step-halving
    #we need to check that the function's value is actually being minimised
    counter <- 0
    while(func(theta + delta,...) > func(theta,...)) {
      print(counter)
      if (counter >= max.half) {
        warning("The step fails to reduce the objective despite trying max.half step halvings")
      }
      delta <- delta/2
      counter <- counter + 1
    }
    
    theta <- theta + delta
    g <- grad(theta,...)
  }
  
  return(theta)
}

#test case 1
rb <- function(th,k) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}
gb <- function(th,k) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}
hb <- function(th,k) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}

#Test case 2
dummy<-function(theta){
  x<-theta[1];y<-theta[2]
  return((1.5-x+x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625-x+x*y^3)^2)
}


dummy_grad<-function(theta){
  return(grad(dummy,theta))
}
dummy_hess<-function(theta){
  return(hessian_csd(dummy,theta))
}

findgradient<-function(theta){
  trig.exp <- expression((1.5-x+x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625-x+x*y^3)^2)
  ( dxy <- deriv(trig.exp, c("x", "y")) )
  y <- theta[2]
  x<- theta[1]
  eval(dxy)
}

t80 <- 1:13 ## years since 1980
y1 <- c(12,14,33,50,67,74,123,141,165,204,253,246,240) ## AIDS cases
nll <- function(theta) {
  ## -ve log likelihood for AIDS model y_i ~ Poi(alpha*exp(beta*t_i))
  ## theta = (alpha,beta)
  mu <- theta[1] * exp(theta[2] * t80) ## mu = E(y)
  -sum(dpois(y1,mu,log=TRUE)) ## the negative log likelihood
} ## nll
gll <- function(theta) {
  ## grad of -ve log lik of Poisson AIDS early epidemic model
  alpha <- theta[1];beta <- theta[2] ## enhances readability
  ebt <- exp(beta*t80) ## avoid computing twice
  -c(sum(y1)/alpha - sum(ebt), ## -dl/dalpha
     sum(y1*t80) - alpha*sum(t80*ebt)) ## -dl/dbeta
} ## gll
newt(c(400,1),dummy,dummy_grad)


