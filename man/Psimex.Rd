\name{Psimex}
\alias{Psimex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pedigree SIMEX 
}
\description{
This function performs the P-SIMEX on a given dataset and a given pedigree. The parameter of interest can be either heritability or inbreeding depression and the error structure can be chosen between missing or misassigned paternities. In the missing case, the error can be simulated uniformly across the pedigree or just in the first or last generations, chosing the number of them. In the misassignment case, the error can be simulated by replacing fathers with random individuals or with similar individuals. 
After simulation, extrapolation is performed and a correct estimate is given together with its standard error according to the chosen function (linear, quadratic or non linear)
}
\usage{
Psimex(pedigree0, data, lambda, lambda0, B = 100, model, 
parameter = "inbreeding", error = "misassignment", way = "uniform", 
fitting.method = "quadratic", ntop = NA, nbottom = NA, 
prior, nitt, thin, burnin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{pedigree0}{
A dataset containing the initial pedigree structure. It must have five columns: id, parent1, parent2, sex, generation.
}
 \item{data}{
A dataset containing the phenotypic measurements on the population  and the covariates which are included in the model.
The trait shoul be named differently than 'trait' (see MCMCglmm) 
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
  \item{model}{
An object specifying the model which has to be fitted to calculate the parameter of interest. It can be a lm or glm for inbreeding depression and a MCMCglmm for heritability.
}
\item{parameter}{
A string specifying the parameter of interest. It must be 'inbreeding' or 'heritability'.
}
  \item{error}{
A string specifying the type of error. It must be 'missing' or 'misassignment'.
}
  \item{way}{
A string specifying how errors are generated. It must be 'uniform' or 'similar' for the misassignement error or 'uniform', 'top' or 'bottom' for the missing case.
}
  \item{fitting.method}{
A string or a vector of strings specifying the extrapolation functions to be fitted. It must be 'linear', 'quadratic', 'nonlinear' or 'cubic'.
}
\item{ntop}{
An integer specifying the number of the first generations to add errors to. It must be specified when the parameter 'way' is 'top'
}
\item{nbottom}{
An integer specifying the number of the last generations to add errors to. It must be specified when the parameter 'way' is 'bottom'
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
\item{description}{A string describing the analysed case}
\item{error}{A string describing the kind of error}
\item{fitting.method}{A string or a vector of strings describing the extrapolation functions}
\item{way}{A string describing the way of simulating the error}
\item{simul_data}{A list of simulated data, it's the output of one of the simulation functions 'simul.replace' or 'simul.na'}     
\item{extrapolated_data}{A list of extrapolated data, it's the output of the extrapolation function} 
\item{lambda}{The vector of the simulated error proportions}
\item{lambda0}{The initial error proportion}
\item{starting.value}{The initial error prone value of the parameter}
    
}
\examples{
## Seed
set.seed(49494)

# extract data
data(pedigree)
data(data)
pedigree0 <- pedigree

# inbreeding depression case

# fixed error proportions 
lambda <- c(0.2, 0.3, 0.4, 0.5, 0.6)
# initial error proportion
lambda0 <- 0.1
# model used to compute inbreeding depression
model <- lm(y~sex+f_inb, data = data)

# PSIMEX 
results <- Psimex(pedigree0, data, 
                  lambda, lambda0, B = 100, 
                  model, parameter = "inbreeding", 
                  error = "missing", way = "uniform", 
                  fitting.method = c("quadratic", "linear"), 
                  ntop = NA, nbottom = NA, 
                  prior, nitt, thin, burnin)
results$description
results$error
results$fitting.method
results$way

results$extrapolated_data 
results$lambda
results$lambda0 
results$starting.value

\dontrun{
# heritability case
## Seed
set.seed(49494)

# extract data
data(pedigree)
data(data)
pedigree0 <- pedigree

# fixed error proportions 
lambda <- c(0.2, 0.3, 0.4, 0.5, 0.6)
# initial error proportion
lambda0 <- 0.1

# model to compute heritability (MCMCglmm)
# prior specification
prior <- list(G=list(G1=list(V=matrix(1/3),n=1),
                     G2=list(V=matrix(1/3),n=1)),
              R=list(V=matrix(1/3),n=1))

#to fulfill MCMCglmm requirements
pedigree <- pedigree0[ , c(1,2,3)]
names(pedigree) <- c("animal", "dam", "sire")
ord <- orderPed(pedigree)
pedigree <- pedigree[order(ord),]

# model specification
model <- MCMCglmm(y~1+sex, random = ~animal+id, 
                  pedigree = pedigree, data = data, 
                  prior = prior, nitt = 20000, thin = 100, burnin = 1000, 
                  verbose = FALSE)

# PSIMEX
results1 <- Psimex(pedigree0, data, 
                   lambda, lambda0, B = 10, 
                   model, parameter = "heritability", 
                   error = "missing", way = "uniform", 
                   fitting.method = "quadratic", 
                   ntop = NA, nbottom = NA, 
                   prior = prior, nitt = 20000, thin = 100, burnin = 1000)

results1$description
results1$error
results1$fitting.method
results1$way

results1$extrapolated_data 
results1$lambda
results1$lambda0 
results1$starting.value

}
}
\author{
Erica Ponzi
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~pedigree }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~SIMEX }% __ONLY ONE__ keyword per line
