#' @title Purity-An external measures
#' @description Purity is a simple and transparent evaluation measure, when the value of purity is close to 1, tell the result of cluster is very good.
#' @param f1 the real classify label 
#' @param f2 the estimated cluster label
#' @param K the number of the communities
#' @param n the number of the nodes
#' @return the proportion of the classified successfully nodes
#' @useDynLib StatComp20025
#' @examples
#' \dontrun{
#' y <- seq(2,50,0.01)
#' plot(y, f(y,2,2),col='blue'.type='l')
#' }
#' @export
Purity=function(f1,f2,K,n){
  L1=list()
  for (k in 1:K){
    Ct=c()
    L1=c(L1,list(Ct))}
  
  for (i in 1:n){L1[[f1[i]]]=c(L1[[f1[i]]],i)}
  
  L2=list()
  for (k in 1:K){
    Ct=c()
    L2=c(L2,list(Ct))}
  for (i in 1:n){L2[[f2[i]]]=c(L2[[f2[i]]],i)}
  
  Ehk=c()
  for (k in 1:K){
    Ekl=c()
    for (l in 1:K){
      e=intersect(L1[[l]],L2[[k]])
      Ekl=c(Ekl,length(e))}
    eo=sort(Ekl, decreasing = TRUE)
    Ehk=c(Ehk,sum(eo[-1]))}
  
  succrate=1-sum(Ehk)/n
  return(succrate)}


#' @title Generate the  theta
#' @description The function generate a matrix for the theta, which will be used in generate data
#' @param clust.size the vector, each element represents the number of a taxon
#' @param K the number of the communities
#' @param n the number of the nodes
#' @return the matrix for theta
#' @useDynLib StatComp20025
#' @examples
#' \dontrun{
#'n <- 200
#'k <- 3 # number of communities
#'m <- 3 # layers 

#'clust.size <- c(rep(floor(n/k), k-1),n-(k-1)*floor(n/k)) # equal cluster size
#'Theta <- Generate.theta(n, k, clust.size)
#' }
#' @export
Generate.theta <- function(n, K, clust.size){
theta <- matrix(0, n, K)
for (k in 1:K){
  if (k==1){id1 <-1} else {id1 <- sum(clust.size[1:(k-1)])+1}
  id2 <- sum(clust.size[1:k])
  theta[id1:id2, k] <- 1
} 
return(theta)
}


#' @title Generate the sample data
#' @description The function generate a \code{m*n*n} multi-layer network sample
#' @param Theta the matrix for theta(membership)
#' @param Btensor the community-wise connectivity parameter 
#' @param self if self==0, the diagonal of the matrix is 1
#' @return the multi-layer network sample 
#' @useDynLib StatComp20025
#' @examples
#' \dontrun{
#'n <- 200
#'k <- 3 # number of communities
#'m <- 3 # layers 

#'clust.size <- c(rep(floor(n/k), k-1),n-(k-1)*floor(n/k)) # equal cluster size
#'Theta <- Generate.theta(n, k, clust.size)
#'Btensor <- array(0, dim = c(3, 3, 3)) # B tensor should be m*k*k

#'Btensor[1, ,] <- matrix(c(0.6,0.4,0.4,0.4,0.2,0.2,0.4,0.2,0.2), 3, 3)
#'Btensor[2, ,] <- matrix(c(0.2,0.4,0.2,0.4,0.6,0.4,0.2,0.4,0.2), 3, 3)
#'Btensor[3, ,] <- matrix(c(0.2,0.2,0.4,0.2,0.2,0.4,0.4,0.4,0.6), 3, 3)
#'A <- Generate.data(Theta, Btensor, self = 0)
#' }
#' @export
Generate.data <- function(Theta, Btensor, self = 0) {
  # generates random adjacency matrix
  p <- dim(Btensor)[1]
  n <- dim(Theta)[1]
  A <- array(0, dim = c(p, n,n))
  for (i in 1:p){ 
    P <- Theta%*%Btensor[i,,]%*%t(Theta)
    lower.tri.ind <- lower.tri(P)
    p.upper <- P[!lower.tri.ind]
    A.upper <- rbinom(n*(n+1)/2, 1, p.upper)
    A0 <- matrix(0, ncol = n, nrow = n)
    A0[!lower.tri.ind] <- A.upper
    tA0 = t(A0)
    A0[lower.tri.ind]  <- tA0[lower.tri.ind]
    if (self == 0){
      diag(A0) <- 0
    }
    A[i,,] = A0
  }
  return(A)
}


#' @title Get the initial estimated Btensor 
#' @description The function obtain an initial estimator of Btensor by simple processing of sample A
#' @param A the sample m*n*n multi-layer network
#' @param idx the initial estiamted membership 
#' @param k the number of the communities
#' @return the initial estimated Btensor
#' @useDynLib StatComp20025
#' @importFrom stats rbinom
#' @export
GetCenter <-function(A, idx, k){
  p <- dim(A)[1]
  m.1 <- dim(A)[2]
  idx.1 <- idx[1:m.1]
  center <- array(0, dim=c(p,k,k))
  for (i in 1:k){
    for (j in i:k){
      if ((sum(idx.1 == i) ==1) & (sum(idx==j)==1) ) {
        center[,i,j] <-0.5
        center[,j,i] <- center[,i,j]
      } else{
        center[,i,j] <- apply(A[,idx.1==i,idx==j],1,mean)
        center[,j,i] <- center[,i,j]
      }
    }
  }
  return(center)
}


#' @title Get distance between the multi-layer network A and the Btensor
#' @description The function calculates the distance between the estimated community-wise connectivity parameter and the real sample, and the return can be used to adjust the estimator
#' @param A the sample m*n*n multi-layer network
#' @param center the m*n*K estimated community-wise connectivity parameter 
#' @return a K*K matrix, each element represents a distance  
#' @useDynLib StatComp20025
#' @export
GetDist <- function(A, center){
  # A: m*n*n
  # center: m*n*K
  k <- dim(center)[3]
  n <- dim(A)[3]
  D <- matrix(0, n, k)
  for (j in 1:n){
    for (i in 1:k){
      D[j,i] <- norm( (A[,,j] - center[,,i]), 'F')^2
    }
  }
  return(D)
}