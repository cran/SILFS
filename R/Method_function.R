#' Initialization Function for the Intercept Parameter
#'
#' This function computes initial values for intercept parameter by solving a ridge regression problem.
#' 
#' @param Y The response vector of length \eqn{n}.
#' @param X The design matrix of size \eqn{n\times p}.
#' @param lam_ridge The tuning parameter for ridge regression.
#' @return A numeric vector of length \eqn{n}, representing the initial estimation for intercept parameter.
#' @export
#' @examples
#' n <- 100
#' p <- 100
#' beta <- rep(1,p)
#' X <- matrix(rnorm(100*100), n, p)
#' Y <- sample(c(-3,3),n,replace=TRUE,prob=c(1/2,1/2)) +  X%*%beta
#' lam_ridge <- 0.1
#' alpha_init <- INIT(Y, X, lam_ridge)
INIT = function(Y,X,lam_ridge){
  n <- nrow(X)
  cov_aug <- cbind(diag(1,n),X)
  dim <- ncol(cov_aug)
  theta <- solve(t(cov_aug)%*%cov_aug+lam_ridge*diag(1,dim))%*%t(cov_aug)%*%Y
  alpha0 <-  theta[1:n]
  return(alpha0)
}


#' SILFS-Based Subgroup Identification and Variable Selection Optimized by Coordinate Descent under the L2 Distance
#' 
#' This function employs SILFS method under L2 distance and uses the Coordinate Descent Algorithm for optimization to effectively identify subgroup structures and perform variable selection.
#'
#' @import MASS
#' @importFrom MASS mvrnorm
#' @param Y The response vector of length \eqn{n}.
#' @param X_aug The augmented design matrix created by row concatenation of common and idiosyncratic factor matrices, with a size of \eqn{n \times (r+p)}.
#' @param r The user supplied number of common factors.
#' @param lam_CAR The tuning parameter for Center-Augmented Regularization.
#' @param lam_lasso The tuning parameter for LASSO.
#' @param alpha_init The initialization of intercept parameter.
#' @param K The user-supplied group number.
#' @param epsilon The user-supplied stopping tolerance.
#' @return A vector containing the following components:
#' \item{alpha_m}{The estimated intercept parameter vector of length \eqn{n}.}
#' \item{gamma}{The estimated vector of subgroup centers of length \eqn{K}.}
#' \item{theta_m}{The estimated regression coefficient vector, matched with common factor terms, with a dimension of \eqn{r}.}
#' \item{beta_m}{The estimated regression coefficients matched with idiosyncratic factors, with a dimension of \eqn{p}.}
#' @author Yong He, Liu Dong, Fuxin Wang, Mingjuan Zhang, Wenxin Zhou.
#' @references He, Y., Liu, D., Wang, F., Zhang, M., Zhou, W., 2024. High-Dimensional Subgroup Identification under Latent Factor Structures.
#' @export
#' @examples
#' n <- 50
#' p <- 50
#' r <- 3
#' K <- 2
#' alpha <- sample(c(-3,3),n,replace=TRUE,prob=c(1/2,1/2))
#' beta <- c(rep(1,2),rep(0,48))
#' B <- matrix((rnorm(p*r,1,1)),p,r)
#' F_1 <- matrix((rnorm(n*r,0,1)),n,r)
#' U <- matrix(rnorm(p*n,0,0.1),n,p)
#' X <- F_1%*%t(B)+U
#' Y <- alpha + X%*%beta + rnorm(n,0,0.5)
#' alpha_init <- INIT(Y,F_1,0.1)
#' SILFS(Y,cbind(F_1,U),3,0.01,0.05,alpha_init,K,0.3)
SILFS <- function(Y,X_aug,r,lam_CAR,lam_lasso,alpha_init,K,epsilon){
  m <- 1
#  if(is.null(Tng_init)){
#    Tng.init =
#  }
  n <- nrow(X_aug)

  d <- ncol(X_aug)

  F_hat <- X_aug[,1:r]

  U_hat <- X_aug[,(r+1):d]

  Res <- list(rep(1,n),rep(0,n))

  alpha_m <- alpha_init

  theta_m <- itertheta(alpha_m,F_hat,Y)

  beta_m <- iterbeta(r,alpha_m,theta_m,X_aug,Y,lam_lasso)

  gamma_m <- as.matrix(kmeans(alpha_m,K)$centers)

  while(sum(abs(Res[[m+1]]-Res[[m]]))/sqrt(n)>epsilon){

    alpha_m <- iteralpha(Y,X_aug,r,alpha_m,beta_m,theta_m,gamma_m,lam_CAR)

    theta_m <- itertheta(alpha_m,F_hat,Y)

    beta_m <- iterbeta(r,alpha_m,theta_m,X_aug,Y,lam_lasso)

    kmeans_m <- kmeans(alpha_m,K,nstart = 50)

    gamma_m <- as.matrix(kmeans_m$centers)

    Res[[m+2]] <- alpha_m

    m <- m+1
  }
  center <- unlist(kmeans_m$centers)
  cluster<- unlist(kmeans_m$cluster)
  gamma <- center[cluster]
  result <- list(c(alpha_m),c(gamma),c(theta_m),c(beta_m))
  return(result)
}


GradS2 <- function(alpha0,gamma0,lambda){
  n = length(alpha0)
  K = length(gamma0)
  sortedgamma = sort(gamma0)
  e_temp = 0
  ekj = matrix(c(rep(0,n)),n)
  for (i in 1:n){
    for (k in 2:K){
      e_temp <- max(abs(alpha0[i]-sortedgamma[k-1]),
                    abs(alpha0[i]-sortedgamma[k]))
      e_temp <- e_temp*sign(alpha0[i] - 0.5*(sortedgamma[k-1]+sortedgamma[k]))
      ekj[i] <- ekj[i] + e_temp
    }
  }
  return(2*lambda*ekj)
}

iteralpha <- function(Y,X,r,alpha0,beta0,theta0,gamma0,lambda){
  n <- nrow(X)
  dim <- ncol(X)
  K <- nrow(gamma0)
  Fm1 <- X[,1:r]
  U <- X[,(r+1):dim]
  Ytil <- Y-Fm1%*%theta0-U%*%beta0
  alphaupdate <- 1/(1/n+2*lambda*K)*((1/n)*Ytil+2*lambda*matrix(c(rep(sum(gamma0),n)),n)+GradS2(alpha0,gamma0,lambda))
  return(alphaupdate)
}

itertheta <- function(alpha,Fhat,Y){
  theta01 <-  matrix(solve(t(Fhat)%*%Fhat)%*%t(Fhat)%*%(Y-alpha))
  return(theta01)
}


#' @import glmnet
iterbeta <- function(r,alpha,theta,X,Y,lamlasso){
  dim_X <- ncol(X)
  Fm2 <- X[,1:r]
  U <- X[,(r+1):dim_X]
  Yc = Y-alpha-Fm2%*%theta
  #lambda_lasso = cv.glmnet(U,Yc,Intercept = "F")$lambda.min
  updatebeta = matrix(glmnet(U,Yc,Intercept = "F",lambda = lamlasso)$beta)
  return(updatebeta)
}



#' Factor Adjusted-Pairwise Fusion Penalty (FA-PFP) Method for Subgroup Identification and Variable Selection
#'
#' This function utilizes the FA-PFP method implemented via the Alternating Direction Method of Multipliers (ADMM) algorithm to identify subgroup structures and conduct variable selection.
#'
#' @import MASS
#' @importFrom MASS mvrnorm
#' @param Fhat The estimated common factors matrix of size \eqn{n \times r}.
#' @param Uhat The estimated idiosyncratic factors matrix of size \eqn{n \times p}.
#' @param Y The response vector of length \eqn{n}.
#' @param vartheta The Lagrangian augmentation parameter for intercepts.
#' @param lam The tuning parameter for Pairwise Fusion Penalty.
#' @param gam The user-supplied parameter for Alternating Direction Method of Multipliers (ADMM) algorithm.
#' @param alpha_init The initialization of intercept parameter.
#' @param lam_lasso The tuning parameter for LASSO.
#' @param epsilon The user-supplied stopping tolerance.
#' @return A list with the following components:
#' \item{alpha_m}{The estimated intercept parameter vector of length \eqn{n}.}
#' \item{theta_m}{The estimated regression coefficient vector, matched with common factor terms, with a dimension of \eqn{r}.}
#' \item{beta_m}{The estimated regression coefficients matched with idiosyncratic factors, with a dimension of \eqn{p}.}
#' \item{eta_m}{A numeric matrix storing the pairwise differences of the estimated intercepts, with size of \eqn{n \times (n\times(n-1)/2)}.}
#' @author Yong He, Liu Dong, Fuxin Wang, Mingjuan Zhang, Wenxin Zhou.
#' @references Ma, S., Huang, J., 2017. A concave pairwise fusion approach to subgroup analysis.
#' @export
#' @examples
#' n <- 50
#' p <- 50
#' r <- 3
#' alpha <- sample(c(-3,3),n,replace=TRUE,prob=c(1/2,1/2))
#' beta <- c(rep(1,2),rep(0,48))
#' B <- matrix((rnorm(p*r,1,1)),p,r)
#' F_1 <- matrix((rnorm(n*r,0,1)),n,r)
#' U <- matrix(rnorm(p*n,0,0.1),n,p)
#' X <- F_1%*%t(B)+U
#' Y <- alpha + X%*%beta + rnorm(n,0,0.5)
#' alpha_init <- INIT(Y,F_1,0.1)
#' FA_PFP(Y,F_1,U,1,0.67,3,alpha_init,0.05,0.3)
FA_PFP<-function(Y,Fhat,Uhat,vartheta,lam,gam,alpha_init,lam_lasso,epsilon){

  alpha_m <- alpha_init

  n <- nrow(Fhat)

  eta_m <- Del(n)%*%alpha_m

  v_m <- matrix(0,n,n)

  m <- 1

  Res <- list(rep(1,n),rep(0,n))

  while((sum(abs(Res[[m+1]]-Res[[m]])))/sqrt(n)>epsilon){

    delta_m <- delu(n,alpha_m,v_m,vartheta)

    eta_m <- etu(n,delta_m,lam,gam,vartheta)

    theta_m <- thetau(n,alpha_m,Fhat,Y)

    beta_m <- betau(alpha_m,theta_m,Fhat,Uhat,Y,lam_lasso)

    alpha_m <- alphau(n,Fhat,Uhat,Y,theta_m,beta_m,vartheta,MtoV(n,eta_m),MtoV(n,v_m))

    v_m <- vu(vartheta,n,alpha_m,eta_m,v_m)

    Res[[m+2]] <- alpha_m

    m <- m+1
  }
  return(list(c(alpha_m),c(theta_m),c(beta_m),eta_m))
}



ST<-function(x,l){
  return(sign(x)*(abs((abs(x)-l))+(abs(x)-l))/2)
}



Del<-function(n){
  D<-matrix(0,1,n)
  for (i in 1:(n-1)){
    Db <- matrix(0,(n-i),n)
    for (j in (i+1):n){
      v<-c(rep(0,n))
      v[i]<- 1
      v[j]<- (-1)
      Db[(j-i),]<-v
    }
    D <- rbind(D,Db)
  }
  return(D[-1,])
}



MtoV <- function(n,M){
  vec<-c()
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      vec<- append(vec,M[i,j])
    }
  }
  return(vec)
}



Qproj<- function(X){
  Q<-X%*%solve(t(X)%*%X)%*%t(X)
  return(Q)
}



delu <- function(n,alpha,v,vartheta){
  d<-matrix(0,n,n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      d[i,j] = alpha[i]-alpha[j]+(vartheta)^(-1)*v[i,j]
    }
  }
  return(d)
}



etu <- function(n,delta,lam,gam,vartheta){
  eta <- matrix(0,n,n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if (abs(delta[i,j])<=lam+(lam/vartheta)){
        eta[i,j]<-ST(delta[i,j],lam/vartheta)
      }
      if (abs((lam+lam/vartheta)<delta[i,j])<=gam*lam){
        eta[i,j]<-(ST(delta[i,j],gam*lam/((gam-1)*vartheta)))/(1-1/((gam-1)*vartheta))
      }
      if (abs(delta[i,j])>gam*lam){
        eta[i,j]<-delta[i,j]
      }
    }
  }
  return(eta)
}


vu<-function(vartheta,n,alpha,eta,v){
  vn<-matrix(0,n,n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      vn[i,j]<-v[i,j]+vartheta*(alpha[i]-alpha[j]-eta[i,j])
    }
  }
  return(vn)
}



alphau <- function(n,Fhat,Uhat,Y,theta,beta,vartheta,eta,v){
  alpha <- solve((diag(n)+vartheta*t(Del(n))%*%Del(n)))%*%(Y-Fhat%*%theta-Uhat%*%beta+vartheta*t(Del(n))%*%(eta-1/(vartheta)*v))
  return(alpha)
}




thetau <- function(n,alpha,Fhat,Y){
  theta<-1/n*t(Fhat)%*%(Y-alpha)
  return(theta)
}


#' @import glmnet
betau <- function(alpha,theta,Fhat,Uhat,Y,lasso){
  Yc <- Y-alpha-Fhat%*%theta
  beta <- matrix(glmnet(Uhat,Yc,Intercept = "F",lambda = lasso)$beta)
  return(beta)
}


#' @import glmnet
Oracle_fac <- function(Y,alpha,Fhat,Uhat,lam,K){
  n <- nrow(Fhat)
  gamma_m <- rep(0,K)
  O <- Omega(alpha)
  m <- 1
  Res <- list(rep(0,K),rep(1,K))
  while((sum(abs(Res[[m+1]]-Res[[m]]))/sqrt(2)>0.01)){
    theta_m <- (1/n)*t(Fhat)%*%(Y-O%*%gamma_m)
    beta_m <- matrix(glmnet(Uhat,(Y-O%*%gamma_m-Fhat%*%theta_m),Intercept = "F", lambda = lam)$beta)
    gamma_m <- solve(t(O)%*%O)%*%t(O)%*%(Y-Fhat%*%theta_m-Uhat%*%beta_m)
    Res[[m+2]] <- gamma_m
    m <- m+1
  }
  return(list(c(gamma_m),c(O%*%gamma_m),c(theta_m),c(beta_m)))
}


#' Standard Center Augmented Regularization (S-CAR) Method for Subgroup Identification and Variable Selection
#'
#' This function employs the S-CAR method under L2 distance and uses the Coordinate Descent Algorithm for optimization to identify subgroup structures and execute variable selection.
#'
#' @import MASS
#' @import stats
#' @importFrom MASS mvrnorm
#' @param Y The response vector of length \eqn{n}.
#' @param X The design matrix of size \eqn{n \times p}.
#' @param lam_CAR The tuning parameter for Center-Augmented Regularization.
#' @param lam_lasso The tuning parameter for lasso.
#' @param alpha_init The initialization of intercept parameter.
#' @param K The estimated group number.
#' @param epsilon The user-supplied stopping tolerance.
#' @return A list with the following components:
#' \item{alpha_m}{The estimated intercept parameter vector of length \eqn{n}.}
#' \item{gamma}{The estimated vector of subgroup centers of length \eqn{K}.}
#' \item{beta_m}{The estimated regression coefficient vector of dimension \eqn{p}.}
#' @author Yong He, Liu Dong, Fuxin Wang, Mingjuan Zhang, Wenxin Zhou.
#' @export
#' @examples
#' n <- 50
#' p <- 50
#' r <- 3
#' K <- 2
#' alpha <- sample(c(-3,3),n,replace=TRUE,prob=c(1/2,1/2))
#' beta <- c(rep(1,2),rep(0,48))
#' B <- matrix((rnorm(p*r,1,1)),p,r)
#' F_1 <- matrix((rnorm(n*r,0,1)),n,r)
#' U <- matrix(rnorm(p*n,0,0.1),n,p)
#' X <- F_1%*%t(B)+U
#' Y <- alpha + X%*%beta + rnorm(n,0,0.5)
#' alpha_init <- INIT(Y,X,0.1)
#' SCAR(Y,X,0.01,0.05,alpha_init,K,0.3)
SCAR <- function(Y,X,lam_CAR,lam_lasso,alpha_init,K,epsilon){
  n <- nrow(X)
  p <- ncol(X)
  Res <- list(rep(1,n),rep(0,n))
  m <- 1
  alpha_m <- alpha_init
  beta_m <- iterbeta_SCAR(alpha_m,X,Y,lam_lasso)
  gamma_m <- as.matrix(kmeans(alpha_m,centers = K)$centers)

  while(sum(abs(Res[[m+1]]-Res[[m]]))/sqrt(n)>epsilon && m<101){
    alpha_m <- iteralpha_SCAR(Y,X,alpha_m,beta_m,gamma_m,lam_CAR)
    beta_m <- iterbeta_SCAR(alpha_m,X,Y,lam_lasso)
    kmeans_m <- kmeans(alpha_m,K,nstart = 50)
    gamma_m <- as.matrix(kmeans_m$centers)
    Res[[m+2]] <- alpha_m
    m <- m+1
  }
  center=unlist(kmeans_m$centers)
  cluster=unlist(kmeans_m$cluster)
  gamma <- center[cluster]
  result = list(c(alpha_m),c(gamma),c(beta_m))
  return(result)
}



GradS2_SCAR = function(alpha0,gamma0,lambda){
  n = length(alpha0)
  K = length(gamma0)
  sortedgamma = sort(gamma0)
  e_temp = 0
  ekj = matrix(c(rep(0,n)),n)
  for (i in 1:n){
    for (k in 2:K){
      e_temp <- max(abs(alpha0[i]-sortedgamma[k-1]),
                    abs(alpha0[i]-sortedgamma[k]))
      e_temp <- e_temp*sign(alpha0[i] - 0.5*(sortedgamma[k-1]+sortedgamma[k]))
      ekj[i] <- ekj[i] + e_temp
    }
  }
  return(2*lambda*ekj)
}



iteralpha_SCAR = function(Y,X,alpha0,beta0,gamma0,lambda){
  n <- nrow(X)
  dim <- ncol(X)
  K <- nrow(gamma0)
  Ytil = Y-X%*%beta0
  alphaupdate = 1/(1/n+2*lambda*K)*((1/n)*Ytil+2*lambda*matrix(c(rep(sum(gamma0),n)),n)+GradS2_SCAR(alpha0,gamma0,lambda))
  return(alphaupdate)
}

#' @import glmnet
iterbeta_SCAR = function(alpha,X,Y,lamlasso){
  Yc = Y-alpha
  updatebeta = matrix(glmnet(X,Yc,Intercept = "F",lambda = lamlasso)$beta)
  return(updatebeta)
}


grad_Z2 <- function(delta){
  n <- nrow(delta)
  K <- ncol(delta)

  Grad <- matrix(0,n,K)


  for(i in 1:n){
    Grad[i,K] <- sign(delta[i,K])*ifelse(abs(delta[i,K]) > abs(delta[i,(K-1)]), 1, 0)
    Grad[i,1] <- sign(delta[i,1])*ifelse(abs(delta[i,1]) > abs(delta[i,2]), 1, 0)
  }

  if(K>2){
    for(j in 1:n){
      for(l in 2:(K-1)){
        Grad[j,l] <- sign(delta[j,l])*ifelse(abs(delta[j,l]) > abs(delta[j,(l-1)]), 1, 0)+sign(delta[j,l])*ifelse(abs(delta[j,l]) > abs(delta[j,(l+1)]), 1, 0)
      }
    }
  }


  return(Grad)
}

C1_mat <- function(n,K){
  block <- rep(1, K)

  C1 <- matrix(0, n * K, n)
  for (i in 1:n) {
    C1[((i - 1) * K + 1):(i * K), i] <- block
  }
  return(C1)
}


C2_mat <- function(n,K){
  C2 <- matrix(0,n*K,K)
  for(i in 1:n){
    C2[((i-1)*K+1):(i*K),] <- -diag(1,K)
  }
  return(C2)
}



D_mat <- function(K){
  D_mat <- matrix(0, K - 1, K)
  for (i in 1:(K - 1)) {
    D_mat[i, i] <- 1
    D_mat[i, i + 1] <- -1
  }
  return(D_mat)

}





Theta_update <- function(Y,U_hat,r_1,r_2,r_3,u,v,w,eta,delta,y,C1,C2,n,K,P,D){
  n <- nrow(delta)
  K <- ncol(delta)
  p <- length(eta)

  Hess <- matrix(0,(n+p+K),(n+p+K))
  Hess[(1:n),(1:n)] <- P/n + r_1*t(C1)%*%C1
  Hess[(1:n),((n+1):(n+p))] <- U_hat/n
  Hess[(1:n),((n+p+1):(n+p+K))] <- r_1*t(C1)%*%C2
  Hess[((n+1):(n+p)),(1:n)] <- t(U_hat)/n
  Hess[((n+1):(n+p)),((n+1):(n+p))] <- t(U_hat)%*%U_hat/n+r_3*diag(1,p)
  Hess[((n+p+1):(n+p+K)),(1:n)] <- r_1*t(C2)%*%C1
  Hess[((n+p+1):(n+p+K)),((n+p+1):(n+p+K))] <- r_2*t(D)%*%D+r_1*t(C2)%*%C2

  h_vec <- c(P%*%Y/n + r_1*t(C1)%*%(as.vector(t(u+delta))), t(U_hat)%*%Y/n+r_3*(w+eta), r_1*t(C2)%*%(as.vector(t(u+delta)))+r_2*t(D)%*%(y-v))

  theta <- solve(Hess)%*%h_vec
  return(theta)
}


soft_threshold <- function(x, lambda){
  ifelse(abs(x) <= lambda, 0, sign(x) * (abs(x) - lambda))
}


delta_update <- function(r_1,u,lambda_1,alpha,gamma,grad,C1,C2,n,K){
  u <- as.vector(t(u))
  grad <- as.vector(t(grad))
  vec <- C1%*%alpha+C2%*%gamma-u + lambda_1*grad/r_1
  delta <- soft_threshold(vec,lambda_1/r_1)
  delta <- matrix(delta,n,K,byrow = T)
  return(delta)
}


eta_update <- function(r_3,lambda_2,beta,w){
  eta <- soft_threshold((beta-w),lambda_2/r_3)
  return(eta)
}


y_update <- function(gamma,v,D){
  K <- length(gamma)
  y <- pmin(rep(0,(K-1)),D%*%gamma+v)
  return(y)
}


#' SILFS-Based Subgroup Identification and Variable Selection Optimized by DC-ADMM under the L1 Distance
#'
#' This function employs SILFS method and uses the corresponding Difference of Convex functions-Alternating Direction Method of Multipliers (DC-ADMM) algorithm for optimization to identify subgroup structures and conduct variable selection under the L1 Distance.
#'
#' @import MASS
#' @import glmnet
#' @import Ckmeans.1d.dp
#' @param Y The response vector of length \eqn{n}.
#' @param F_hat The estimated factor matrix of size \eqn{n \times r}.
#' @param U_hat The estimated idiosyncratic factors matrix of size \eqn{n \times p}.
#' @param r_1 The Lagrangian augmentation parameter for constraints of intercepts.
#' @param r_2 The Lagrangian augmentation parameter for constraints of group centers.
#' @param r_3 The Lagrangian augmentation parameter for constraints of coefficients.
#' @param lambda_1 The tuning parameter for Center-Augmented Regularization.
#' @param lambda_2 The tuning parameter for LASSO.
#' @param K The estimated group number.
#' @param alpha_init The initialization of intercept parameter.
#' @param epsilon_1 The user-supplied stopping error for outer loop.
#' @param epsilon_2 The user-supplied stopping error for inner loop.
#' @return A list with the following components:
#' \item{alpha_curr}{The estimated intercept parameter vector of length \eqn{n}.}
#' \item{gamma_curr}{The estimated vector of subgroup centers of length \eqn{K}.}
#' \item{theta_curr}{The estimated regression coefficient vector, matched with common factor terms, with a dimension of \eqn{r}.}
#' \item{beta_curr}{The estimated regression coefficients matched with idiosyncratic factors, with a dimension of \eqn{p}.}
#' @author Yong He, Liu Dong, Fuxin Wang, Mingjuan Zhang, Wenxin Zhou.
#' @references He, Y., Liu, D., Wang, F., Zhang, M., Zhou, W., 2024. High-Dimensional Subgroup Identification under Latent Factor Structures.
#' @export
#' @examples
#' n <- 50
#' p <- 50
#' r <- 3
#' K <- 2
#' alpha <- sample(c(-3,3),n,replace=TRUE,prob=c(1/2,1/2))
#' beta <- c(rep(1,2),rep(0,48))
#' B <- matrix((rnorm(p*r,1,1)),p,r)
#' F_1 <- matrix((rnorm(n*r,0,1)),n,r)
#' U <- matrix(rnorm(p*n,0,0.1),n,p)
#' X <- F_1%*%t(B)+U
#' Y <- alpha + X%*%beta + rnorm(n,0,0.5)
#' alpha_init <- INIT(Y,F_1,0.1)
#' DCADMM_iter_l1(Y,F_1,U,0.5,0.5,0.5,0.01,0.05,K,alpha_init,1,0.3)
DCADMM_iter_l1 <- function(Y,F_hat,U_hat,r_1,r_2,r_3,lambda_1,lambda_2,K,alpha_init,epsilon_1,epsilon_2){
  n <- nrow(U_hat)
  p <- ncol(U_hat)
  C1 <- C1_mat(n,K)
  C2 <- C2_mat(n,K)
  D <- D_mat(K)
  P <- diag(1,n)- 1/n*F_hat%*%t(F_hat)
  alpha_curr <- alpha_init
  gamma_curr <- sort(Ckmedian.1d.dp(alpha_curr,K)$centers)
  delta_curr <- outer(alpha_curr, gamma_curr, "-")
  init_lambda <- cv.glmnet(U_hat,P%*%(Y-alpha_curr),intercept = F)$lambda.min
  beta_curr <- glmnet(U_hat,P%*%(Y-alpha_curr),lambda = init_lambda,intercept = F)$beta
  u_curr <- matrix(0,n,K)
  v_curr <- rep(0,(K-1))
  w_curr <- rep(0,p)
  y_curr <- y_update(gamma_curr,v_curr,D)
  eta_curr <- rep(1,p)
  Res_1 <- list(rep(0,n),rep(1,n))
  #  Res_2 <- list(rep(0,n),rep(1,n))
  count_1 <- 1
  count_2 <- 1
  while(sum(abs(Res_1[[count_1+1]]-Res_1[[count_1]]))>epsilon_1){
    current_grad <- grad_Z2(delta_curr)
    while((sum(abs(beta_curr-eta_curr))+sum(abs(as.vector(t(delta_curr))-C1%*%alpha_curr-C2%*%gamma_curr))+sum(abs(D%*%gamma_curr-y_curr)))>epsilon_2| count_2 < 5){
      eta_curr <-  eta_update(r_3,lambda_2,beta_curr,w_curr)
      delta_curr <- delta_update(r_1,u_curr,lambda_1,alpha_curr,gamma_curr,current_grad,C1,C2,n,K)
      Theta_curr <- Theta_update(Y,U_hat,r_1,r_2,r_3,u_curr,v_curr,w_curr,eta_curr,delta_curr,y_curr,C1,C2,n,K,P,D)
      alpha_curr <- Theta_curr[1:n]
      theta_curr <- 1/n*t(F_hat)%*%(Y-alpha_curr)
      beta_curr <- Theta_curr[(n+1):(n+p)]
      gamma_curr <- Theta_curr[(n+p+1):(n+p+K)]
      y_curr <- y_update(gamma_curr,v_curr,D)

      u_curr_vec <- as.vector(t(u_curr)) + r_1*(as.vector(t(delta_curr))-C1%*%alpha_curr-C2%*%gamma_curr)
      u_curr <- matrix(u_curr_vec,n,K,byrow = T)
      v_curr <- v_curr + r_2*(D%*%gamma_curr-y_curr)
      w_curr <- w_curr + r_3*(eta_curr-beta_curr)
      #      Res_2[[count_2+2]] <- alpha_curr
      count_2 <- count_2 + 1
    }

    Res_1[[count_1+2]] <- alpha_curr
    count_1 <- count_1 + 1
  }
  
  result <- list(c(alpha_curr),c(gamma_curr),c(theta_curr),c(beta_curr))
  return(result)
}





grad_Z2_l2 <- function(delta){
  n <- nrow(delta)
  K <- ncol(delta)
  Grad <- matrix(0,n,K)


  for(i in 1:n){
    Grad[i,K] <- 2*delta[i,K]*ifelse(abs(delta[i,K]) > abs(delta[i,(K-1)]), 1, 0)
    Grad[i,1] <- 2*delta[i,1]*ifelse(abs(delta[i,1]) > abs(delta[i,2]), 1, 0)
  }

  if(K>2){
    for(j in 1:n){
      for(l in 2:(K-1)){
        Grad[j,l] <- 2*delta[j,l]*ifelse(abs(delta[j,l]) > abs(delta[j,(l-1)]), 1, 0) + 2*delta[j,l]*ifelse(abs(delta[j,l]) > abs(delta[j,(l+1)]), 1, 0)
      }
    }
  }


  return(Grad)
}

delta_update_l2 <- function(r_1,u,lambda_1,alpha,gamma,grad,C1,C2,n,K){
  u <- as.vector(t(u))
  grad <- as.vector(t(grad))
  vec <- r_1*C1%*%alpha + r_1*C2%*%gamma + lambda_1*grad - r_1*u
  delta <- (1/(r_1+2*lambda_1))*vec
  delta <- matrix(delta,n,K,byrow = T)
  return(delta)
}


#' SILFS-Based Subgroup Identification and Variable Selection Optimized by DC-ADMM under the L2 Distance
#'
#' This function employs SILFS method and uses the corresponding Difference of Convex functions-Alternating Direction Method of Multipliers (DC-ADMM) algorithm for optimization to identify subgroup structures and conduct variable selection under the L2 Distance.
#'
#' @import MASS
#' @import glmnet
#' @import Ckmeans.1d.dp
#' @param Y The response vector of length \eqn{n}.
#' @param F_hat The estimated factor matrix of size \eqn{n \times r}.
#' @param U_hat The estimated idiosyncratic factors matrix of size \eqn{n \times p}.
#' @param r_1 The Lagrangian augmentation parameter for constraints of intercepts.
#' @param r_2 The Lagrangian augmentation parameter for constraints of group centers.
#' @param r_3 The Lagrangian augmentation parameter for constraints of coefficients.
#' @param lambda_1 The tuning parameter for Center-Augmented Regularization.
#' @param lambda_2 The tuning parameter for LASSO.
#' @param K The estimated group number.
#' @param alpha_init The initialization of intercept parameter.
#' @param epsilon_1 The user-supplied stopping error for outer loop.
#' @param epsilon_2 The user-supplied stopping error for inner loop.
#' @return A list with the following components:
#' \item{alpha_curr}{The estimated intercept parameter vector of length \eqn{n}.}
#' \item{gamma_curr}{The estimated vector of subgroup centers of length \eqn{K}.}
#' \item{theta_curr}{The estimated regression coefficient vector, matched with common factor terms, with a dimension of \eqn{r}.}
#' \item{beta_curr}{The estimated regression coefficients matched with idiosyncratic factors, with a dimension of \eqn{p}.}
#' @author Yong He, Liu Dong, Fuxin Wang, Mingjuan Zhang, Wenxin Zhou.
#' @references He, Y., Liu, D., Wang, F., Zhang, M., Zhou, W., 2024. High-Dimensional Subgroup Identification under Latent Factor Structures.
#' @export
#' @examples
#' n <- 50
#' p <- 50
#' r <- 3
#' K <- 2
#' alpha <- sample(c(-3,3),n,replace=TRUE,prob=c(1/2,1/2))
#' beta <- c(rep(1,2),rep(0,48))
#' B <- matrix((rnorm(p*r,1,1)),p,r)
#' F_1 <- matrix((rnorm(n*r,0,1)),n,r)
#' U <- matrix(rnorm(p*n,0,0.1),n,p)
#' X <- F_1%*%t(B)+U
#' Y <- alpha + X%*%beta + rnorm(n,0,0.5)
#' alpha_init <- INIT(Y,F_1,0.1)
#' DCADMM_iter_l2(Y,F_1,U,0.5,0.5,0.5,0.01,0.05,K,alpha_init,1,0.3)
DCADMM_iter_l2 <- function(Y,F_hat,U_hat,r_1,r_2,r_3,lambda_1,lambda_2,K,alpha_init,epsilon_1,epsilon_2){
  n <- nrow(U_hat)
  p <- ncol(U_hat)
  C1 <- C1_mat(n,K)
  C2 <- C2_mat(n,K)
  D <- D_mat(K)
  P <- diag(1,n)- 1/n*F_hat%*%t(F_hat)
  alpha_curr <- alpha_init
  gamma_curr <- sort(Ckmedian.1d.dp(alpha_curr,K)$centers)
  delta_curr <- outer(alpha_curr, gamma_curr, "-")
  init_lambda <- cv.glmnet(U_hat,P%*%(Y-alpha_curr),intercept = F)$lambda.min
  beta_curr <- glmnet(U_hat,P%*%(Y-alpha_curr),lambda = init_lambda,intercept = F)$beta
  u_curr <- matrix(0,n,K)
  v_curr <- rep(0,(K-1))
  w_curr <- rep(0,p)
  y_curr <- y_update(gamma_curr,v_curr,D)
  eta_curr <- rep(1,p)
  Res_1 <- list(rep(0,n),rep(1,n))
  #  Res_2 <- list(rep(0,n),rep(1,n))
  count_1 <- 1
  count_2 <- 1
  while(sum(abs(Res_1[[count_1+1]]-Res_1[[count_1]]))>epsilon_1){
    current_grad <- grad_Z2_l2(delta_curr)
    while((sum(abs(beta_curr-eta_curr))+sum(abs(as.vector(t(delta_curr))-C1%*%alpha_curr-C2%*%gamma_curr))+sum(abs(D%*%gamma_curr-y_curr)))>epsilon_2| count_2 < 5){
      eta_curr <-  eta_update(r_3,lambda_2,beta_curr,w_curr)
      delta_curr <- delta_update_l2(r_1,u_curr,lambda_1,alpha_curr,gamma_curr,current_grad,C1,C2,n,K)
      Theta_curr <- Theta_update(Y,U_hat,r_1,r_2,r_3,u_curr,v_curr,w_curr,eta_curr,delta_curr,y_curr,C1,C2,n,K,P,D)
      alpha_curr <- Theta_curr[1:n]
      theta_curr <- 1/n*t(F_hat)%*%(Y-alpha_curr)
      beta_curr <- Theta_curr[(n+1):(n+p)]
      gamma_curr <- Theta_curr[(n+p+1):(n+p+K)]
      y_curr <- y_update(gamma_curr,v_curr,D)

      u_curr_vec <- as.vector(t(u_curr)) + r_1*(as.vector(t(delta_curr))-C1%*%alpha_curr-C2%*%gamma_curr)
      u_curr <- matrix(u_curr_vec,n,K,byrow = T)
      v_curr <- v_curr + r_2*(D%*%gamma_curr-y_curr)
      w_curr <- w_curr + r_3*(eta_curr-beta_curr)
      #      Res_2[[count_2+2]] <- alpha_curr
      count_2 <- count_2 + 1
    }

    Res_1[[count_1+2]] <- alpha_curr
    count_1 <- count_1 + 1
  }
  
  result <- list(c(alpha_curr),c(gamma_curr),c(theta_curr),c(beta_curr))
  return(result)
}
