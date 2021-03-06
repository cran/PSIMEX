\name{simul.replace.uni.herit}
\alias{simul.replace.uni.herit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation phase for P-SIMEX on heritability with misassignment error and uniform replacement.
}
\description{
The function generates increasing amount of errors to the pedigree structure and calculates heritability for each simulation. 
For each level of the error proportion, a number of simulation equal to B are repeated. 
Errors are generated with the same probability across the all pedigree by replacing fathers with other individuals from the same generation. 
}
\usage{
simul.replace.uni.herit(pedigree, lambda, B,  data, 
model, prior, nitt,thin,burnin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pedigree}{
A dataset containing the initial pedigree structure. It must have three columns: id, parent1, parent2.
}
  \item{lambda}{
A vector of real numbers specifying the error proportion to be generated.
}
  \item{B}{
An integer specifying the number of simulations to be run for each error level.
}
  
  \item{data}{
A dataset containing the phenotypic measurements on the population  and the covariates which are included in the model.
The trait shoul be named differently than 'trait' (see MCMCglmm) 
}
  \item{model}{
A object of class MCMCglmm used to calculate heritability.
}
\item{prior}{
Prior distribution for MCMCglmm model
}
\item{nitt}{
Number of iterations for MCMCglmm model
}
\item{thin}{
Thinning interval for MCMCglmm model
}
\item{burnin}{
Burn in period for MCMCglmm model
}
}

\value{
A list:
\item{h}{A dataset containing all the simulated values for heritability}
\item{se_h}{A dataset containing all the simulated standard errors for heritability}
\item{VA}{A dataset containing all the simulated values for additive variance}
\item{VE}{A dataset containing all the simulated values for environmental variance}
}

\author{
Erica Ponzi
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~misassignment }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~heritability }% __ONLY ONE__ keyword per line
