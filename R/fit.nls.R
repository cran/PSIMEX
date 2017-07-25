fit.nls <-
function (lambda, p.names, estimates)
  {
    extrapolation <- list()
    # lambda values for the interpolation
    lambda0 <- c(0, max(lambda) / 2, max(lambda))
    for (d in p.names) {
     
      quad <- lm(estimates[, d] ~ lambda + I(lambda^2))
     
      a.nls <- predict(quad,
                       newdata   =   data.frame(lambda  =  lambda0))
     
      gamma.est.3 <- ((a.nls[2] - a.nls[3]) * lambda0[3] *
                        (lambda0[2] - lambda0[1]) - lambda0[1] *
                        (a.nls[1] - a.nls[2]) * (lambda0[3] - lambda0[2])) /
        ((a.nls[1] - a.nls[2]) * (lambda0[3] - lambda0[2]) -
           (a.nls[2] - a.nls[3]) * (lambda0[2] - lambda0[1]))
      
      gamma.est.2 <- ((a.nls[2] - a.nls[3]) * (gamma.est.3 + lambda0[2]) *
                        (gamma.est.3 + lambda0[3])) / (lambda0[3] - lambda0[2])
      
      gamma.est.1 <- a.nls[1] - (gamma.est.2 / (gamma.est.3 + lambda0[1]))
      
      
      extrapolation[[d]] <-
        nls(estimates[, d] ~ gamma.1 + gamma.2 / (gamma.3 + lambda),
            start  =  list(gamma.1  =  gamma.est.1, gamma.2  =  gamma.est.2,
                           gamma.3  =  gamma.est.3))
    }
    return(extrapolation)
  }
