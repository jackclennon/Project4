hessian<-function(theta,func,eps){
  n<-length(theta)
  hess<-matrix(nrow=n,ncol=n)
  for (i in 1:n){
    for (j in 1:n){
      e_i=rep(0,n)
      e_i[i]=1
      e_j=rep(0,n)
      e_j[j]=1
      plus_plus<-func(theta+eps*e_i + eps*e_j)
      plus_minus<-func(theta+eps*e_i - eps*e_j)
      minus_plus<-func(theta-eps*e_i + eps*e_j)
      minus_minus<-func(theta-eps*e_i - eps*e_j)
      finitedifferencing<-(plus_plus-plus_minus-minus_plus+minus_minus)/(4*eps^2)
      hess[i,j]<-finitedifferencing
      
    }
  }
  return(hess)
}

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  g <- grad(theta)
  if (hess==NULL){
    hess<-hessian(theta,func,eps)
  }
  while (length(g[g <= tol]) > 0) {
    #Need to find a way to get rid of the solve
    theta <- theta - solve(hess(theta)) %*% grad(theta)
    g <- grad(theta)
  }
  theta
}
