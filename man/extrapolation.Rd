\name{extrapolation}
\alias{extrapolation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extrapolation phase for the SIMEX algorithm on pedigrees
}
\description{
Starting from the simulated values obtained in the previous SIMEX phase, 
it computes averages across all simulations for each error rate 
and fits a function of the estimate in terms of error rate. 
From this function it extrapolates the value corresponding to the null error rate.
The fitted function can be linear, quadratic, cubic or non linear.  
}
\usage{
extrapolation(results, lambda, lambda0, estimate0, fitting.method, B, parameter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results}{
A dataset with $100$ or B values for each error rate 
given by the simulation part of the algorithm.
}
  \item{lambda}{
A numerical vector containing the error rates.
}
  \item{lambda0}{
The initial error rate.
}
  \item{estimate0}{
The initial error prone estimate.
}
  \item{fitting.method}{
A string or a vector of strings containig the fitting methods for the function. 
It can be: 'line', 'quad', 'nonl' or 'cubi'.
}
  \item{B}{
The number of simulation for each error rate.
}
\item{parameter}{
A string containg the parameter of interest P-SIMEX is performed on. 
It can be either 'inbreeding' or 'heritability'.
}
}

\value{

For inbreeding a list: 
 \item{inb_pred}{Inbreeding depression extrapolated error free value}
 \item{inb_pred_se }{Standard error of the error free value: regression component}
 \item{se_pred }{Standard error of the error free value: sampling error component}
 \item{var }{Total variance of the error free value}
 \item{AIC}{The AIC of the fitted function}
 For heritability a list:
  \item{h_pred}{Heritability extrapolated error free value}
  \item{h_pred_se}{Standard error of the error free value: regression component}
  \item{se_pred }{Standard error of the error free value: sampling error component}
  \item{var }{Total variance of the error free value}
  \item{VA_pred}{Extrapolated value for the additive genetic variance} 
  \item{VE_pred}{Extrapolated value for the environmental variance} 
  \item{AIC}{The AIC of the fitted function}
 }
\author{
Erica Ponzi
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ extrapolation }% use one of  RShowDoc("KEYWORDS")
