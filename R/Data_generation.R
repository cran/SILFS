Cov_fac <- function(n,p,r,beta,alpha){
  B <- matrix(0,p,r)
  for (i in 1:p){
    for (j in 1:r){
      B[i,j]=rnorm(1,0,1)
    }
  }

  #generate F
  Phi <- matrix(0,r,r)
  for (i in 1:r){
    for (j in 1:r){
      Phi[i,j]=(0.3)^abs(i-j)*0.5^(i==j)
    }
  }
  F0 <- matrix(0,2*n,r)
  f0 <- c(0,0,0)
  F0[1,] <- f0
  for (i in 1:(2*n-1)){
    F0[(i+1),]<-t(Phi%*%F0[i,])+mvrnorm(1,c(0,0,0),diag(0.1,r))
  }
  F_t <- F0[(n+1):(2*n),]

  #generate U
  U <- matrix(0,n,p)
  for (i in 1:n){
    U[i,]=mvrnorm(1,c(rep(0,p)),diag(0.1,p))
  }

  #generate X
  X <- F_t%*%t(B)+U
  #Estimate
  Fhat <- FacF(X,r)
  Uhat <- FacU(X,r)
  Xaug <- cbind(Fhat,Uhat)
  #regression data checked
  Y <- X%*%beta+alpha+mvrnorm(1,c(rep(0,n)),diag(0.1,n))
  res <- list(Fhat,F_t,Uhat,X,Xaug,Y,B)
  return(res)
}


Intercept_alpha <- function(type,n){
  alpha <- switch(type,
                  "1"<- sample(c(-3,3),n,replace=T,prob=c(1/2,1/2)),
                  "2"<- sample(c(-5,5),n,replace=T,prob=c(1/2,1/2)),
                  "3"<- sample(c(-3,0,3),n,replace=T,prob=c(1/3,1/3,1/3)),
                  "4"<- sample(c(-5,0,5),n,replace=T,prob=c(1/3,1/3,1/3)),
  )
}


FacF = function(X,k){
  n = nrow(X)
  p = ncol(X)
  XXT = X%*%t(X)
  eigenvectors = eigen(XXT)$vectors
  Fhat = sqrt(n)*(eigenvectors[,1:k])
  return(Fhat)
}


FacU = function(X,k){
n = nrow(X)
p = ncol(X)
XXT = X%*%t(X)
eigenvectors = eigen(XXT)$vectors
Fhat = sqrt(n)*(eigenvectors[,1:k])
Bhat = 1/n*t(X)%*%Fhat
Uhat = X - Fhat%*%t(Bhat)
return(Uhat)
}
