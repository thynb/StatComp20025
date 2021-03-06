% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_Research.R
\name{GetDist}
\alias{GetDist}
\title{Get distance between the multi-layer network A and the Btensor}
\usage{
GetDist(A, center)
}
\arguments{
\item{A}{the sample m*n*n multi-layer network}

\item{center}{the m*n*K estimated community-wise connectivity parameter}
}
\value{
a K*K matrix, each element represents a distance
}
\description{
The function calculates the distance between the estimated community-wise connectivity parameter and the real sample, and the return can be used to adjust the estimator
}
