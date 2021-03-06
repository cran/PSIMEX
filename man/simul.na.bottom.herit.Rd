\name{simul.na.bottom.herit}
\alias{simul.na.bottom.herit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation phase for P-SIMEX on heritability with missing paternities and errors in the last generations of the pedigree.
}
\description{
The function generates increasing amount of errors to the pedigree structure and calculates heritability for each simulation. 
For each level of the error proportion, a number of simulation equal to B are repeated. 
Errors are generated within the last nbottom generations of pedigree by eliminating fathers.
}
\usage{
simul.na.bottom.herit(pedigree, pedigree0, lambda, B, 
data, model, nbottom, prior, nitt,thin,burnin)
}
%- maybe also 'usage' for other objects documented here.
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{pedigree}{
A dataset containing the initial pedigree structure. It must have three columns: id, parent1, parent2.
}
 \item{pedigree0}{
A dataset containing the initial pedigree structure. 
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
\item{nbottom}{
An integer specifying the number of the last generations where to generate errors.
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
\keyword{ ~missing }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~heritability }% __ONLY ONE__ keyword per line
\keyword{ ~bottom}