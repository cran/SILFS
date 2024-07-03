GK<-function(Y,Fhat,Uhat,r,alpha_init,epsilon){
  n <- nrow(Uhat)
  p <- nrow(Uhat)
  GKs<-c(rep(0,4))
  Xaug <- cbind(Fhat,Uhat)
  r <- ncol(Fhat)
  for(i in 2:5){
    R <- SILFS(Y,Xaug,r,0.001,0.01,alpha_init,i,epsilon)
    gamma <- unlist(R[2])
    theta1 <- unlist(R[3])
    beta <- unlist(R[4])
    Yhat <- gamma+ Fhat%*%theta1+ Uhat%*%beta
    S <- sum(beta!=0)
    GKs[i-1] <- log(sum((Y-Yhat)^2)/n)+2*log(n*i+p)*(S+i)*log(n)/n
  }
  return(which.min(GKs)+1)
} #BIC of Group number for SILFS


#' Selecting Tuning Parameter for SILFS Method via corresponding BIC
#'
#' This function is to select tuning parameters simultaneously for SILFS method via minimizing the BIC.
#'
#' @param Y The response vector of length \eqn{n}.
#' @param Fhat The estimated common factors matrix of size \eqn{n \times r}.
#' @param Uhat The estimated idiosyncratic factors matrix of size \eqn{n \times p}.
#' @param K The estimated subgroup number.
#' @param alpha_init The initialization of intercept parameter.
#' @param epsilon The user-supplied stopping tolerance.
#' @param lasso_start The user-supplied start search value of the tuning parameters for LASSO.
#' @param lasso_stop The user-supplied stop search value of the tuning parameters for LASSO.
#' @param CAR_start The user-supplied start search value of the tuning parameters for Center-Augmented Regularization.
#' @param CAR_stop The user-supplied stop search value of the tuning parameters for Center-Augmented Regularization.
#' @param grid_1 The user-supplied number of search grid points corresponding to the LASSO tuning parameter.
#' @param grid_2 The user-supplied number of search grid points corresponding to the tuning parameter for Center-Augmented Regularization.
#' @return A list with the following components:
#' \item{lasso}{The tuning parameter of the LASSO penalty selected using BIC.}
#' \item{CAR}{The tuning parameter of the Center Augmented Regularization selected using BIC.}
#' @export
#' @examples
#' n <- 50
#' p <- 50
#' r <- 3
#' K <- 2
#' lasso_start <- sqrt(log(p)/n)*0.01
#' lasso_stop <- sqrt(log(p)/n)*10^(0.5)
#' CAR_start <- 0.001
#' CAR_stop <- 0.1
#' grid_1 <- 5
#' grid_2 <- 5
#' alpha <- sample(c(-3,3),n,replace=TRUE,prob=c(1/2,1/2))
#' beta <- c(rep(1,2),rep(0,48))
#' B <- matrix((rnorm(p*r,1,1)),p,r)
#' F_1 <- matrix((rnorm(n*r,0,1)),n,r)
#' U <- matrix(rnorm(p*n,0,0.1),n,p)
#' X <- F_1%*%t(B)+U
#' Y <- alpha + X%*%beta + rnorm(n,0,0.5)
#' alpha_init <- INIT(Y,F_1,0.1)
#' \donttest{
#' BIC_SILFS(Y,F_1,U,K,alpha_init,lasso_start,lasso_stop,CAR_start,CAR_stop,grid_1,grid_2,0.3)
#' }
BIC_SILFS <- function(Y,Fhat,Uhat,K,alpha_init,lasso_start,lasso_stop,CAR_start,CAR_stop,grid_1,grid_2,epsilon){
  n <- nrow(Fhat)
  p <- ncol(Uhat)
  r <- ncol(Fhat)
  alpha <- alpha_init
  Xaug <- cbind(Fhat,Uhat)
  s <- K
  GCVs<-matrix(0,grid_1,grid_2)
  for(j in 1:grid_1){
    for(k in 1:grid_2){
      Lasso <- seq(lasso_start,lasso_stop,length.out=grid_1)[j]
      CAR <- seq(CAR_start,CAR_stop,length.out=grid_2)[k]
      res <- SILFS(Y,Xaug,r,CAR,Lasso,alpha,s,epsilon)
      gamma <- unlist(res[2])
      theta1 <- res[[3]]
      beta <- unlist(res[4])
      Yhat <- gamma + Fhat%*% theta1+ Uhat%*%beta
      GCVs[j,k]<- sum((Y-Yhat)^2)/(n-length(unique(beta[beta!=0])))^2
    }
  }
  min_positions <- which(GCVs == min(GCVs), arr.ind = TRUE)
  min_row <- min(min_positions[, 1])
  min_col <- min(min_positions[, 2])
  lasso <- seq(lasso_start,lasso_stop,length.out=grid_1)[min_row]
  CAR <- seq(CAR_start,CAR_stop,length.out=grid_2)[min_col]
  return(c(lasso,CAR))
}




G <- function(eta){
  K <- 1
  n <- ncol(eta)
  g1 <- c(rep(0,n))
  group <- list(c(1))
  etasym <- matrix(0,n,n)
  for(s in 1:(n-1)){
    for(l in (s+1):n){
      etasym[l,s] <- eta[s,l]
      etasym[s,l] <- eta[s,l]
    }
  }

  for (s in 2:n){
    for(i in 1:K){
      if (min(abs(etasym[group[[i]],s]))==0){
        group[[i]] <- c(group[[i]],s)
      }
    }
    if (min(abs(etasym[unlist(group),s]))>0){
      group <- append(group,c(s))
    }
    K <- length(group)
  }

  for (m in 1:K){
    g1[group[[m]]] <- m
  }

  return(list(K,g1))
}




#' Selecting Tuning Parameter for Factor Adjusted-Pairwise Fusion Penalty (FA-PFP) Method via corresponding BIC
#'
#' This function is to select tuning parameters simultaneously for FA-PFP method via minimizing the BIC.
#'
#' @param Y The response vector of length \eqn{n}.
#' @param Fhat The estimated common factors matrix of size \eqn{n \times r}.
#' @param Uhat The estimated idiosyncratic factors matrix of size \eqn{n \times p}.
#' @param alpha_init The initialization of intercept parameter.
#' @param lasso_start The user-supplied start search value of the tuning parameters for LASSO.
#' @param lasso_stop The user-supplied stop search value of the tuning parameters for LASSO.
#' @param lam_start The user-supplied start search value of the tuning parameters for Pairwise Fusion Penalty.
#' @param lam_stop The user-supplied stop search value of the tuning parameters for Pairwise Fusion Penalty.
#' @param grid_1 The user-supplied number of search grid points corresponding to the LASSO tuning parameter.
#' @param grid_2 The user-supplied number of search grid points corresponding to the tuning parameter for Pairwise Fusion Penalty.
#' @param epsilon The user-supplied stopping tolerance.
#' @return A list with the following components:
#' \item{lasso}{The tuning parameter of the LASSO penalty selected using BIC.}
#' \item{lambda}{The tuning parameter of the Pairwise Concave Fusion Penalty selected using BIC.}
#' @author Yong He, Liu Dong, Fuxin Wang, Mingjuan Zhang, Wenxin Zhou.
#' @export
#' @examples
#' n <- 50
#' p <- 50
#' r <- 3
#' lasso_start <- sqrt(log(p)/n)*0.1
#' lasso_stop <- sqrt(log(p)/n)
#' lam_start <- 0.3
#' lam_stop <- 1
#' grid_1 <- 5
#' grid_2 <- 5
#' alpha <- sample(c(-3,3),n,replace=TRUE,prob=c(1/2,1/2))
#' beta <- c(rep(1,2),rep(0,48))
#' B <- matrix((rnorm(p*r,1,1)),p,r)
#' F_1 <- matrix((rnorm(n*r,0,1)),n,r)
#' U <- matrix(rnorm(p*n,0,0.1),n,p)
#' X <- F_1%*%t(B)+U
#' Y <- alpha + X%*%beta + rnorm(n,0,0.5)
#' alpha_init <- INIT(Y,F_1,0.1)
#'\donttest{
#' BIC_PFP(Y,F_1,U,alpha_init,lasso_start,lasso_stop,lam_start,lam_stop,grid_1,grid_2,0.3)
#' }
BIC_PFP <- function(Y,Fhat,Uhat,alpha_init,lasso_start,lasso_stop,lam_start,lam_stop,grid_1,grid_2,epsilon){
  n <- nrow(Fhat)
  p <- ncol(Uhat)
  BICs <- matrix(0,grid_1,grid_2)
  for(j in 1:grid_1){
    for(k in 1:grid_2){
      Lasso <- seq(lasso_start,lasso_stop,length.out=grid_1)[j]
      lam <- seq(lam_start,lam_stop,length.out=grid_2)[k]
      res <- FA_PFP(Y,Fhat,Uhat,1,lam,3,alpha_init,Lasso,epsilon)
      alpha_hat <- res[[1]]
      theta_hat <- res[[2]]
      eta <- res[[4]]
      beta_hat <- res[[3]]
      Yhat <- alpha_hat + Fhat%*%theta_hat + Uhat%*%beta_hat
      S <-sum(abs(beta_hat)>0)
      g <- G(eta)[[1]]
      BICs[j,k]<-log(sum((Y-Yhat)^2)/n)+2*log(n*g+p)*(S+g)*log(n)/n
    }
  }
  min_positions <- which(BICs == min(BICs), arr.ind = TRUE)
  min_row <- min(min_positions[, 1])
  min_col <- min(min_positions[, 2])
  results <- c(min_row ,min_col)
  lasso <- seq(lasso_start,lasso_stop,length.out=grid_1)[results[2]]
  lambda <- seq(lam_start,lam_stop,length.out=grid_2)[results[1]]
  return(c(lasso,lambda))
}




Omega <- function(alpha){
  n <- length(alpha)
  O <- matrix(0,n,length(unique(alpha)))
  for (i in (1:length(unique(alpha)))){
    O[,i] <- as.numeric(alpha==(unique(alpha)[i]))
  }
  return(O)
}




BIC_Oracle_fac<-function(Fhat,Uhat,Y,alpha){
  n<-nrow(Fhat)
  p <- ncol(Uhat)
  GCVs<-c()
  for(j in 1:10){
    Lasso <- (sqrt(log(p)/n)*10^(seq(-2,0,length.out=10)[j]))
    R <- Oracle_fac(Y,alpha,Fhat,Uhat,Lasso)
    gamma <- R[[2]]
    theta <- R[[3]]
    beta <- R[[4]]
    Yhat <- gamma + Fhat%*%theta+ Uhat%*%beta
    GCVs[j] <- sum((Y-Yhat)^2)/(n-length(unique(beta[beta!=0])))^2
  }
  min_positions <- min(which(GCVs == min(GCVs)))
  lasso <- sqrt(log(p)/n)*10^(seq(-2,0,length.out=10)[min_positions])
  return(lasso)
}




GK_SCAR<-function(Y,X,n,p,epsilon){
  GKs<-c(rep(0,4))
  for(i in 2:5){
    R <- SCAR(Y,X,5,0.001,0.1,i,epsilon)
    gamma <- unlist(R[1])
    theta1 <- unlist(R[2])
    beta <- unlist(R[3])
    Yhat <- gamma+ X%*%beta
    S <- sum(beta!=0)
    GKs[i-1] <- log(sum((Y-Yhat)^2)/n)+2*log(n*i+p)*(S+i)*log(n)/n
  }
  return(which.min(GKs)+1)
}




BIC_SCAR <- function(X,Y,epsilon){
  n <- nrow(X)
  p <- ncol(X)
  K <- GK_SCAR(Y,X,n,p)
  GCVs<-matrix(0,10,10)
  for(j in 1:10){
    for(k in 1:10){
      Lasso <- (sqrt(log(p)/n)*10^(seq(-1,0.5,length.out=10)[j]))
      CAR <- seq(0.001,0.1,length.out=10)[k]
      R <- SCAR(Y,X,CAR,Lasso,0.1,K,epsilon)
      gamma <- unlist(R[2])
      beta <- unlist(R[3])
      Yhat <- gamma + X%*%beta
      GCVs[j,k]<- sum((Y-Yhat)^2)/(n-length(unique(beta[beta!=0])))^2
    }
  }
  min_positions <- which(GCVs == min(GCVs), arr.ind = TRUE)
  min_row <- min(min_positions[, 1])
  min_col <- min(min_positions[, 2])
  lasso <- (sqrt(log(p)/n)*10^(seq(-1,0.5,length.out=10)[min_row]))
  CAR <- seq(0.001,0.1,length.out=10)[min_col]
  return(c(K,lasso,CAR))
}


