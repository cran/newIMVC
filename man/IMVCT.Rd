% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMVC.R
\name{IMVCT}
\alias{IMVCT}
\title{Integrated Mean Variance Correlation Based Hypothesis Test}
\usage{
IMVCT(x, y, K, num_per, NN = 3, type)
}
\arguments{
\item{x}{is the univariate covariate vector}

\item{y}{is the response vector}

\item{K}{is the number of quantile levels}

\item{num_per}{is the number of permutation times}

\item{NN}{is the number of B spline basis, default is 3}

\item{type}{is an indicator for measuring linear or nonlinear correlation, "linear" represents linear correlation and "nonlinear" represents linear or nonlinear correlation using B splines}
}
\value{
The p-value of the corresponding hypothesis test
}
\description{
This function is used to test significance of linear or nonlinear correlation using integrated mean variance correlation
}
\examples{
# linear model
n=100
x=rnorm(n)
y=2*x+rt(n,2)

IMVCT(x,y,K=5,type = "linear")
# nonlinear model
n=100
x=rnorm(n)
y=2*cos(x)+rt(n,2)

IMVCT(x,y,K=5,type = "nonlinear",num_per = 100)
}
