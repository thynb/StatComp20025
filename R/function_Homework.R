#' @title Pareto(a, b) in the Exercises 3.3, Homework1
#' @description The density function for Pareto distribution with parameter(a,b)
#' @param x the value of the independent variable 
#' @param a the exponent part
#' @param b the Molecular part
#' @return the value when density function with parameter\code{(a,b)} at \code{x}.
#' @useDynLib StatComp20025
#' @examples
#' \dontrun{
#' y <- seq(2,50,0.01)
#' plot(y, f_Pareto(y,2,2),col='blue'.type='l')
#' }
#' @export
f_Pareto=function(x,a,b){
  x1=x^(-a-1)
  c1=a*(b^a)
  return(c1*x1)
}


#' @title Pareto(a, b)^-1 in the Exercises 3.3, Homework1
#' @description The inverse distribution for Pareto with parameter(a,b)
#' @param x the value of the independent variable
#' @param a the exponent part
#' @param b the Molecular part
#' @return the value when inverse distribution with parameter\code{(a,b)} at \code{x}.
#' @importFrom stats runif
#' @examples
#' \dontrun{
#' n <- 1000
#' u <- runif(n)
#' x <- f_inv(u,2,2)
#' hist(x, density = 3,breaks=50, col="red",prob = TRUE)
#' }
#' @export
f_inv=function(x,a,b){
  f1=b
  f2=(1-x)^(1/a)
  return(f1/f2)
}


#' @title rescaled Epanechnikov kernel in Exercise5.13, Homework3
#' @description Devroye and Gyorfi's algorithm for simulation from the rescaled Epanechnikov kernel's distribution, in the Exercise 3.9, homework1.
#' @param n the number of the samples you want to generate 
#' @return the random samples followed the goal distribution
#' @importFrom stats runif 
#' @examples
#' \dontrun{
#' x=G_Rek(1000)
#' hist(x, density = 3,breaks=30, col="red",prob = TRUE) #density histogram of sampl
#' }
#' @export
G_Rek=function(n){
  x=c()
for (i in 1:n) {
  U=runif(3,-1,1)
  if( abs(U[3])>=abs(U[2]) && abs(U[3])>=abs(U[1])){x=c(x,U[2])}
  else{x=c(x,U[3])}
  }
  return(x)}


#' @title The monte carlo in Exercise5.13, Homework3
#' @description a special goal density we need to use the Monte Carlo method to integration
#' @param x the number of the samples you want to generate 
#' @return the value of density function at \code{x}
#' @export
g_goal=function(x){
  return((x^2)*exp(-x^2/2)/sqrt(2*pi))
}


#' @title The important function(1) for monte carlo in Exercise5.13, Homework3
#' @description a important function help us to use the Monte Carlo method to integrate the goal density
#' @param x the value of the independent variable 
#' @return the value of the function at \code{x}
#' @export
f1=function(x){
  return(sqrt(2/pi)*exp(-(x-1)^2/2))
}


#' @title The important function(2) for monte carlo in Exercise5.13, Homework3
#' @description a important function help us to use the Monte Carlo method to integrate the goal density
#' @param x the value of the independent variable 
#' @return the value of the function at \code{x}
#' @export
f2=function(x){
  return(exp(-x+1))
}


#' @title The monte carlo in Exercise5.15, Homework3
#' @description a special goal density we need to use the Monte Carlo piecewise method to integration
#' @param x the value of the independent variable 
#' @return the value of the function at \code{x}
#' @export
g_goal2 <- function(x) {exp(-x - log(1+x^2)) * (x > 0) * (x < 1)}


#' @title The important function in the monte carlo in Exercise5.15, Homework3
#' @description a important function help us to use the Monte Carlo piecewise method to integrate the goal density
#' @param x the value of the independent variable 
#' @param j the order of the piecewise
#' @param k the total number of the piecewise
#' @return the value of the function at \code{x}
#' @export
f_j=function(x,j,k){(exp(-x))/(exp(-(j-1)/k)-exp(-j/k))* (x > 0) * (x < 1)}


#' @title Estimate the power of the skewness test in Exercise 6.7 in Homework4
#' @description The function is used to computes the sample's skewness coeff
#' @param x the sample 
#' @return the value of the sample \code{x}'s skewness
#' @importFrom stats rnorm
#' @examples  
#' \dontrun{
#' X=rnorm(1000,0,2)
#' sk(X)
#' }
#' @export
sk <- function(x) {
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}




#' @title  The power of the CountFive test in Exercise 6.8 in Homework4
#' @description  The two sample “Count Five” test for equality of variance introduced by McGrath and Yeh [193] counts the number of extreme points of each sample relative to the range of the other sample.
#' @param x the sample named X
#' @param y the sample named Y
#' @return the value 1 or 0: reject or accept H0
#' @examples  
#' \dontrun{
#' X=rnorm(1000,0,2)
#' Y=rnorm(1000,0,1)
#' count5test(X,Y)
#' }
#' @export
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx, outy)) > 5))
}


#' @title  Estimate the power of Mardia's multivariate skewness test for in Exercise 6.C in Homework4
#' @description  Mardia proposed tests of multivariate normality based on multivariate generalizations of skewness and kurtosis.
#' @param X the multivariate sample 
#' @return the value of the multivariate skewness statistic 
#' @importFrom stats cov
#' @export
sk_m <- function(X) {
  Xbar=apply(X,2,mean)
  XSigma=cov(X)
  beta=0
 for (i in 1:nrow(X)) {
    for ( j in 1:nrow(X)){
      beta=beta+(matrix((X[i,]-Xbar),nrow = 1)%*%(solve(XSigma))%*%t(matrix(X[j,]-Xbar,nrow = 1)))^3
    }}
  return(beta/((nrow(X))^2))
}


#' @title  Maximum number of extreme points in Exercises 8.3,Homework6
#' @description  The function calculate the maximum number of extreme points that applies when sample sizes are not necessarily equal
#' @param x the sample X
#' @param y the sample y 
#' @return the maximum number of extreme points 
#' @examples 
#' \dontrun{
#' X=rnorm(1000,0,2)
#' Y=rnorm(1000,1,1)
#' maxout(X,Y)
#' }
#' @export
maxout <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return((max(c(outx, outy))))
}


#' @title  Nearest neighbor tests in Discussion, Homework6
#' @description  The nearest neighbor (NN) tests are a type of test based on ordered distances between sample elements, which can be applied when the distributions are continuous.
#' @param z the matrix of sample X and Y
#' @param ix the index of sample
#' @param sizes [1]the number of sample x; [2] the number of sample y
#' @param k the nearest neighbor statistic measures the proportion of first through \code{k}th nearest neighbor coincidences
#' @return the proportion of first through \code{k}th nearest neighbor coincidences 
#' @importFrom RANN nn2
#' @export
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]
  n2 <- sizes[2] 
  n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0)
  z <- z[ix, ]
  NN <- nn2(data=z, k=k+1) 
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5)
  i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}


#' @title  A random walk Metropolis sampler in Exercise 9.4, Homework7
#' @description This is a random walk Metropolis sampler for generating the standard Laplace distribution (see Exercise 3.2). For the increment, simulate from a normal distribution. 
#' @param sigma the standard deviation of a normal distribution
#' @param x0 initial value
#' @param N the number of the goal sample
#' @param f the density function of the sample
#' @return the list of sample subject to the goal distribution
#' @importFrom stats runif rnorm
#' @export
rw.Metropolis <- function(sigma, x0, N,f) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (f(y) / f(x[i-1])))
      x[i] <- y 
    else {
      x[i] <- x[i-1]
      k <- k + 1
    }
  }
  return(list(x=x, k=k))
}


#' @title The Gelman-Rubin method to monitor convergence of the chain, Exercise 9.4+, Homework7
#' @description The Gelman-Rubin method of monitoring convergence of a M-H chain is based on comparing the behavior of several generated chains with respect to the variance of one or more scalar summary statistics. The estimates of the variance of the statistic are analogous to estimates based on between-sample and within-sample mean squared errors in a one-way analysis of variance (ANOVA)
#' @param psi psi[i,j] is the statistic psi(X[i,1:j])
#' @return  the estimated potential scale reduction
#' @importFrom stats var
#' @export
Gelman.Rubin <- function(psi) {
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}


#' @title The Metropolis chain in Exercise 9.4+, Homework7
#' @description Generates a Metropolis chain for Laplace distribution with Normal(X[t], sigma) proposal distribution and starting value X1
#' @param N length of chains
#' @param sigma the standard deviation of a normal distribution
#' @param X1 the initial value of the chains
#' @param f the density function of standard Laplace distribution
#' @return  the chains
#' @importFrom stats var
#' @export
Laplace.chain <- function(sigma, N, X1,f) {
  x <- rep(0, N)
  x[1] <- X1
  u <- runif(N)
  for (i in 2:N) {
    xt <- x[i-1]
    y <- rnorm(1, xt, sigma) #candidate point
    r <- f(y) / f(x[i-1])
    if (u[i] <= r) x[i] <- y else
      x[i] <- xt
  }
  return(x)
}


#' @title The function need to be solved in Exercise 11.4, Homework7
#' @description To find the intersection points for S_k(a) and S_k-1(a), set f_S be the difference between the two equations and then get its root.
#' @param k degrees of freedom
#' @param a the independent variables
#' @return the difference between the two equations
#' @importFrom stats pt
#' @export
f_Sk=function(k,a){
  p1=sqrt(a^2*k/(k+1-a^2))
  p2=sqrt(a^2*(k-1)/(k-a^2))
  return(pt(p1,k)-pt(p2,k-1))}

#' @title lapply variant, in Ex4, Homework8
#' @description The function is a combination of Map() and vapply() to create an lapply() variant that iterates in parallel over all of its inputs and stores its outputs in a vector (or a matrix).
#' @param f The function that operates on data
#' @param n The dimension of data
#' @param type the data type: numeric, character, complex, logical
#' @param ... the input
#' @return the data after the operation
#' @importFrom stats pt
#' @export
F_Mv<-function (f,n,type,...) {  
  NM=unlist(Map(f,...)) 
  if(type=="numeric") return(vapply(NM,cbind,numeric(n))) 
  else if (type=="character") return(vapply(NM,cbind,character(n))) 
  else if (type=="complex") return(vapply(NM,cbind,complex(n))) 
  else if (type=="logical") return(vapply(NM,cbind,logical(n))) 
}

