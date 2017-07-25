\name{fit.nls}
\alias{fit.nls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting a non linear function on the data
}
\description{
It fits a non linear function on a set of data 
}
\usage{
fit.nls(lambda, p.names, estimates)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{
A vector containing the x axis points
}
  \item{p.names}{
Name of the y axis variables to be fitted
}
  \item{estimates}{
Y axis variables to be fitted
}
}

\value{
A list containing the fitted functions of the function in the null value on the x axis.
It contains one function per variable to be estimated.
}

\author{
Erica Ponzi
}




% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ extrapolation }% use one of  RShowDoc("KEYWORDS")
\keyword{nonlinearfit }% __ONLY ONE__ keyword per line
