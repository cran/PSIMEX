\name{initial}
\alias{initial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of the initial error prone value for inbreeding depression
}
\description{
This function calculates the starting point for the P-SIMEX 
on inbreeding depression with the initial error in the pedigree.
}
\usage{
initial(model, data, pedigree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
An object specifying the model to calculate inbreeding depression. It can be a lm or glm.
}
  \item{data}{
A dataset containing the phenotypic measurements on the population
and the covariates which are included in the model.
}
  \item{pedigree}{
A dataset containing the initial pedigree structure. 
It must have three columns: id, parent1, parent2.
}
}

\value{
A list:
\item{model}{The correct expression of the model to calculate inbreeding depression}
\item{inb0}{The starting value of inbreeding depression}
\item{se0}{The starting standard error of inbreeding depression}
\item{pval0}{The starting p-value of inbreeding depression}
\item{mean_inb0}{The starting mean value of inbreeding coefficient}
\item{median_inb0}{The starting median value of inbreeding coefficient}
\item{var_inb0}{The starting variance of inbreeding coefficient}


  
  
}

\author{
Erica Ponzi
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~inbreeding }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~startingvalue }% __ONLY ONE__ keyword per line
