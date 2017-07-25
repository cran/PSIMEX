\name{simul.replace.uni}
\alias{simul.replace.uni}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation phase for P-SIMEX on inbreeding depression with misassignment error and uniform replacement.
}
\description{
The function generates increasing amount of errors to the pedigree structure and calculates inbreeding depression for each simulation. 
For each level of the error proportion, a number of simulation equal to B are repeated. 
Errors are generated with the same probability across the all pedigree by replacing fathers with other individuals from the same generation. 
}
\usage{
simul.replace.uni(pedigree, lambda, B, data, model)
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
A object of class model used to calculate inbreeding depression.
}
}

\value{
A list: 
\item{inb}{A dataset containing all the simulated values for inbreeding depression}
\item{se_inb}{A dataset containing all the simulated standard errors for inbreeding depression}
\item{pval}{A dataset containing all the simulated p-values for inbreeding depression}
\item{mean_inb}{A dataset containing all the simulated values for inbreeding coefficient mean}
\item{median_inb}{A dataset containing all the simulated values for inbreeding coefficient median}
\item{var_inb}{A dataset containing all the simulated values for inbreeding coefficient variance}

}

\author{
Erica Ponzi
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~inbreeding }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~misassignment }% __ONLY ONE__ keyword per line
