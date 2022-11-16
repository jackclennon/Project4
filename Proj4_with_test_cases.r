hessian <- function(theta, grad, eps, ...) {
  #function for calculating the Hessian
  #takes theta, the variables of the function to find the Hessian of
  #grad, the gradient of the function
  #eps, a small number
  #... extra parameters to pass to grad
  
  #start by putting grad into g to make it easier to see what's going on
  g<-grad
  
  #get the number of elements in theta
  n<-length(theta)
  
  #create an empty nxn matrix where we'll put the Hessian
  hess<-matrix(nrow=n,ncol=n)
  
  #loop over each row in the Hessian
  for (i in 1:n){
    
    #loop over each column up to the current row
    #this makes the lower triangle of the Hessian
    for (j in 1:i){
      #make the vector that will add a slight change to the ith element of theta
      e_i=rep(0,n)
      e_i[i]=1
      
      #make the vector that will add a slight change to the jth element of theta
      e_j=rep(0,n)
      e_j[j]=1
      
      #this is the equation for the element H_ij
      hess[i, j] = ((g(theta + eps*e_j,...)[i] - g(theta-eps*e_j,...)[i])/(4*eps) +  ((g(theta+eps*e_i,...)[j] - g(theta-eps*e_i,...)[j])/(4*eps)))
      
      #the Hessian is symmetric, so we can skip calculating the upper triangle
      #and just use the transpose of the lower triangle for the upper triangle
      hess[j, i] = hess[i, j]
    }
  }
  hess
}

newt <- function(theta, func, grad, hess=NULL, ..., tol=1e-8, fscale=1, maxit=100, max.half=20,eps=1e-6) {
  #function for carrying out Newton's method on a function
  
  g <- grad(theta,...)
  if (is.null(hess)){
    hess<-function(theta,...){hessian(theta,grad,eps,...)}
  }
  
  if(is.infinite(func(theta, ...)))
    stop('The objective function is not finite at the initial theta')
  
  if (all(is.finite(g))==FALSE)
    stop('The gradient vector contains non-finite values at initial theta')
  
  # for counting how many times we iterate Newton's method
  iterations<-0
  
  # breaks when there is convergence
  while (length(g[abs(g) >= tol*(abs(func(theta, ...))+fscale)]) > 0) { 
    
    # ends the algorithm if we exceed the maximum number of iterations maxit
    if (iterations > maxit)
      stop("Iterations exceeded maxit ", maxit)
    
    #check if the hessian is positive definite. If it isn't, the we perturb it
    H<-hess(theta,...)
    
    while (inherits(try(chol(H), silent=TRUE), "try-error")) {
      if (all(is.finite(H))==FALSE)
        stop('The Hessian contains non-finite values')
      #perturb the hessian if it isn't positive definite
      H <- H + diag(length(theta))*(abs(min(eigen(H)$values))+1)
    }
    # chol2inv(chol(H)) <-> solve(H) <-> H^-1 because H is positive definite
    Hi <- chol2inv(chol(H))
    delta <- -Hi %*% grad(theta,...)
    
    # step-halving
    # we need to check that the function's value is actually being minimised
    counter <- 0
    while(func(theta + delta,...) > func(theta,...) || is.infinite(func(theta, ...))) {
      if (counter >= max.half) {
        stop("Step halving fails to reduce the objective function despite trying max.half step halvings")
      }
      
      #halves the step size
      delta <- delta/2
      
      #counts the number of step halvings
      counter <- counter + 1
    }
    
    #update the value for theta
    theta <- theta + delta
    g <- grad(theta,...)
    
    iterations<-iterations+1
  }
  
  H<-hess(theta,...)
  
  if (inherits(try(chol(H), silent=TRUE), "try-error")) {
    warning("Hessian is not positive definite at minimum and thus the inverse cannot be computed")
    Hi <- NULL
  } else {
    Hi <- chol2inv(chol(H))
  }
  
  return(list(f=func(theta, ...),
              theta=theta,
              iter=iterations-1,
              g=grad(theta, ...),
              Hi=Hi))
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
dummy0<-function(theta){ # (3,0.5)
  x<-theta[1];y<-theta[2]
  return((1.5-x+x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625-x+x*y^3)^2)
}

dummy_grad0<-function(theta){
  return(grad(dummy,theta))
}
dummy_hess0<-function(theta){
  return(hessian_csd(dummy,theta))
}

#Test case 3
dummy1<-function(theta){  # x=(-1), y=1.5
  x<-theta[1];y<-theta[2]
  return(x-y+2*x^2+2*x*y+y^2)
}
dummy_grad1<-function(theta){
  return(grad(dummy1,theta))
}
dummy_hess1<-function(theta){
  return(hessian_csd(dummy1,theta)) # Hessian for complex step derivatives (csd)
}

#Test Case 4
dummy2<-function(theta){  # The Hessian contains non-finite values
  x<-theta[1];y<-theta[2]
  return((exp(x*y)-exp(-x*y))/2)
}
dummy_grad2<-function(theta){
  return(grad(dummy2,theta))
}
dummy_hess2<-function(theta){
  return(hessian_csd(dummy2,theta)) # Hessian for complex step derivatives (csd)
}

# Test Case 5
dummy3<-function(theta){ # This is negative definite
  x<-theta[1];y<-theta[2]
  return(-y-x^2+2*x-2)
}
dummy_grad3<-function(theta){
  return(grad(dummy3,theta))
}
dummy_hess3<-function(theta){
  return(hessian_csd(dummy3,theta)) # Hessian for complex step derivatives (csd)
}

# Test Case 6
dummy4<-function(theta){ # This is positive in-definite
  x<-theta[1];y<-theta[2]
  return(-y+2*x^2-x-5)
}
dummy_grad4<-function(theta){
  return(grad(dummy4,theta))
}
dummy_hess4<-function(theta){
  return(hessian_csd(dummy4,theta)) # Hessian for complex step derivatives (csd)
}

# Test Case 7
dummy5<-function(theta){ # Hess
  x<-theta[1];y<-theta[2];z<-theta[3]
  return(x^2*y^2*z^2+1)
}
dummy_grad5<-function(theta){
  return(grad(dummy5,theta))
}
dummy_hess5<-function(theta){
  return(hessian_csd(dummy5,theta)) # Hessian for complex step derivatives (csd)
}

#Test Case 8
dummy6<-function(theta){ # This is positive definite
  return((theta[1]+theta[2]+theta[3])^2-1)
}
dummy_grad6<-function(theta){
  return(grad(dummy6,theta))
}
dummy_hess6<-function(theta){
  return(hessian_csd(dummy6,theta)) # Hessian for complex step derivatives (csd)
}

#Test Case 9
dummy7<-function(theta){ # This is positive definite
  x<-theta[1];y<-theta[2];z<-theta[3]
  return(100*(x^2+y^2+z^2)+2*x*z+4*y*z)
}
dummy_grad7<-function(theta){
  return(grad(dummy7,theta))
}
dummy_hess7<-function(theta){
  return(hessian_csd(dummy7,theta)) # Hessian for complex step derivatives (csd)
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

# 2 Variable Functions
newt(c(400,1),dummy0,dummy_grad0, maxit=500)

newt(c(400,1),dummy1,dummy_grad1, maxit=500)

newt(c(400,1),dummy2,dummy_grad2, maxit=500)

newt(c(400,1),dummy3,dummy_grad3, maxit=500)

newt(c(400,1),dummy4,dummy_grad4, maxit=500)

# 3 Variable Functions
newt(c(10,10,10),dummy5,dummy_grad5, maxit=500)

newt(c(10,10,10),dummy6,dummy_grad6, maxit=500)

newt(c(10,10,10),dummy7,dummy_grad7, maxit=500)

f<-function(x){(x[1]+x[2]+x[3])^2-1}
p=(c(10,10,10))
nlm(f,p)
