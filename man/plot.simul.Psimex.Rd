\name{plot.simul.Psimex}
\alias{plot.simul.Psimex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A plot illustrating the simulation phase of the P-SIMEX algorithm
}
\description{
This function produces a plot with the estimate of interest and its 95 \% confidence intervals as a function of the error proportion.
}
\usage{
plot.simul.Psimex(results, lambda, lambda0, estimate0, parameter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results}{
A dataset containg the simulated values. It is the output of simul.replace or simul.na.
}
  \item{lambda}{
A vector of real numbers containing the error proportions to be generated.
}
  \item{lambda0}{
A real number specifying the initial error proportion.
}
  \item{estimate0}{
The initial naive estimate of the parameter
}
\item{parameter}{
A string specifying the parameter of interest. It must be 'inbreeding' or 'heritability'.
}
}


\author{
Erica Ponzi
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~plot }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~simulation }% __ONLY ONE__ keyword per line
