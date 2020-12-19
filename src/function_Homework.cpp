#include <Rcpp.h>
using namespace Rcpp;

//' @title rw.Metropolis using Rcpp, in Ex1, Homework9
//' @description The function is for random walk Metropolis sampler to generate the standard Laplace distribution, which is written in Rcpp
//' @param sigma the standard deviation of the normal distribution
//' @param x0 the initial value of the chains
//' @param N the length of the chains
//' @return a matrix,include the value of each step, and the reject times
//' @export
// [[Rcpp::export]]
NumericMatrix rwMC(double sigma,double x0,int N) {
  NumericMatrix mat(N, 2);
  mat(0,0)=x0;
  mat(0,1)=0;
  double y=0,u=0;
  for(int i = 2; i < N+1; i++) {
    y=rnorm(1,mat(i-2,0),sigma)[0];
    u=runif(1,0,1)[0];
    if(u<=exp(-abs(y))/exp(-abs(mat(i-2,0)))){
      mat(i-1,0)=y;
      mat(i-1,1)=mat(i-2,1);
    }
    else{
      mat(i-1,0)=mat(i-2,0);
      mat(i-1,1)=mat(i-2,1)+1;
    }
  }
  return(mat);
}

//' @title Runif using Rcpp, in Ex2, Homework9
//' @description The function is for randomly generating sample that are uniformly distributed on (0,1) with Rcpp
//' @param N the number of the sample
//' @return a vector of random sample
//' @export
// [[Rcpp::export]]
NumericVector GetRandomNumber(int N)
{
  NumericVector RandomNumber(N);
  int i;
  int b=10000;
  
  srand((unsigned)time(NULL));
  for(i=0;i<N;i++){
    RandomNumber(i) = rand() % b / (double)(b)  ;
  }
  return RandomNumber;
}


