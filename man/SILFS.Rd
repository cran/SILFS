% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Method_function.R
\name{SILFS}
\alias{SILFS}
\title{SILFS-Based Subgroup Identification and Variable Selection Optimized by Coordinate Descent under the L2 Distance}
\usage{
SILFS(Y, X_aug, r, lam_CAR, lam_lasso, alpha_init, K, epsilon)
}
\arguments{
\item{Y}{The response vector of length \eqn{n}.}

\item{X_aug}{The augmented design matrix created by row concatenation of common and idiosyncratic factor matrices, with a size of \eqn{n \times (r+p)}.}

\item{r}{The user supplied number of common factors.}

\item{lam_CAR}{The tuning parameter for Center-Augmented Regularization.}

\item{lam_lasso}{The tuning parameter for LASSO.}

\item{alpha_init}{The initialization of intercept parameter.}

\item{K}{The user-supplied group number.}

\item{epsilon}{The user-supplied stopping tolerance.}
}
\value{
A vector containing the following components:
\item{alpha_m}{The estimated intercept parameter vector of length \eqn{n}.}
\item{gamma}{The estimated vector of subgroup centers of length \eqn{K}.}
\item{theta_m}{The estimated regression coefficient vector, matched with common factor terms, with a dimension of \eqn{r}.}
\item{beta_m}{The estimated regression coefficients matched with idiosyncratic factors, with a dimension of \eqn{p}.}
}
\description{
This function employs SILFS method under L2 distance and uses the Coordinate Descent Algorithm for optimization to effectively identify subgroup structures and perform variable selection.
}
\examples{
n <- 50
p <- 50
r <- 3
K <- 2
alpha <- sample(c(-3,3),n,replace=TRUE,prob=c(1/2,1/2))
beta <- c(rep(1,2),rep(0,48))
B <- matrix((rnorm(p*r,1,1)),p,r)
F_1 <- matrix((rnorm(n*r,0,1)),n,r)
U <- matrix(rnorm(p*n,0,0.1),n,p)
X <- F_1\%*\%t(B)+U
Y <- alpha + X\%*\%beta + rnorm(n,0,0.5)
alpha_init <- INIT(Y,F_1,0.1)
SILFS(Y,cbind(F_1,U),3,0.01,0.05,alpha_init,K,0.3)
}
\references{
He, Y., Liu, D., Wang, F., Zhang, M., Zhou, W., 2024. High-Dimensional Subgroup Identification under Latent Factor Structures.
}
\author{
Yong He, Liu Dong, Fuxin Wang, Mingjuan Zhang, Wenxin Zhou.
}
