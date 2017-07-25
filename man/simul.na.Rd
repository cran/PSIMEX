\name{simul.na}
\alias{simul.na}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation phase for the P-SIMEX with missing paternities.
}
\description{
This function calls the proper simulation function for the parameter of interest with missing paternities.
}
\usage{
simul.na(pedigree, pedigree0, lambda, lambda0,  B, data, 
model, parameter, way, ntop, nbottom, prior, nitt,thin,burnin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{pedigree}{
A dataset containing the initial pedigree structure. It must have three columns: id, parent1, parent2.
}
\item{pedigree0}{
A dataset containing the initial pedigree structure. It must have five columns: id, parent1, parent2, sex, generation.
}
  \item{lambda}{
A vector of real numbers specifying the error proportion to be generated.
}
  \item{lambda0}{
A real number specifying the initial error rate.
}
  \item{B}{
An integer specifying the number of simulations to be run for each error level.
}

  \item{data}{
A dataset containing the phenotypic measurements on the population  and the covariates which are included in the model.
The trait shoul be named differently than 'trait' (see MCMCglmm) 
}
  \item{model}{
An object specifying the model to calculate the parameter of interest. It can be a lm or glm for inbreeding depression and a MCMCglmm for heritability.
}
\item{parameter}{
A string specifying the parameter of interest. It must be 'inbreeding' or 'heritability'.
}
  \item{way}{
A string specifying how errors are generated. It must be 'uniform' or 'similar'.
}
\item{ntop}{
An integer specifying the number of the first generations where to generate errors, when the parameter 'way' is 'top'.
}

\item{nbottom}{
An integer specifying the number of the last generations where to generate errors,when the parameter 'way' is 'bottom'.
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
For inbreeding a list: 
\item{inb}{A dataset containing all the simulated values for inbreeding depression}
\item{se_inb}{A dataset containing all the simulated standard errors for inbreeding depression}
\item{pval}{A dataset containing all the simulated p-values for inbreeding depression}
\item{mean_inb}{A dataset containing all the simulated values for inbreeding coefficient mean}
\item{median_inb}{A dataset containing all the simulated values for inbreeding coefficient median}
\item{var_inb}{A dataset containing all the simulated values for inbreeding coefficient variance}

For heritability a list:
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
\keyword{ ~simulation}
\keyword{ ~missing}
