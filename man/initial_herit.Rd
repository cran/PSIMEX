\name{initial_herit}
\alias{initial_herit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of the initial error prone value for heritability 
}
\description{
This function calculates the starting point for the P-SIMEX 
on heritability with the initial error in the pedigree.
}
\usage{
initial_herit(model, data, pedigree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
An object specifying the model to calculate inbreeding depression. 
It can be a lm or glm.
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
\item{h0}{The starting value of heritability}
\item{seh0}{The starting standard error of heritability}
\item{VA0}{The starting value of additive genetic variance}
\item{VE0}{The starting value of environmental variance}

}

\author{
Erica Ponzi
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~heritability }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~startingvalue }% __ONLY ONE__ keyword per line
