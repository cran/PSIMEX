Psimex <-
function (
 pedigree0,  #pedigree structure in the form of "id", "dam", "sire", "generation"
 data,       #must contain year of birth in case of misassigned paternities
   #vector of ids
 lambda,     #error proportions
 lambda0,    #initial error proportion
 B = 100, 
 model, 
 parameter = 'inbreeding',
 error = 'misassignment', 
 way = 'uniform', 
 fitting.method = 'quadratic', 
 ntop = NA, 
 nbottom = NA, 
 prior, nitt,thin,burnin) 
  {
  
   for (i in 1:length(fitting.method)) {
    f <- fitting.method[i]
    fitting.method[i] <- substr(f, 1, 4)
    f <- fitting.method[i]
    if (!any(f == c("quad", "line", "nonl", "cubi"))) {
      warning("Fitting Method not implemented. Using: quadratic", call. = FALSE)
      fitting.method[i] <- "quad"
    }
   }
  
    
    if (!any(error == c("misassignment", "missing"))) {
      warning("Error not implemented. Using: misassignment", call. = FALSE)
      error <- "misassignment"
    }
    
    if (error == 'misassignment') {
      
      if (!any(way == c("uniform", "similar"))) {
        warning("Way not implemented. Using: uniform", call. = FALSE)
        way <- "uniform"
      }
      
    }
    
    if (error == 'missing') {
      
      if (!any(way == c("uniform", "top", "bottom"))) {
        warning("Way not implemented. Using: uniform", call. = FALSE)
        way <- "uniform"
      }
      
    }
    
    if (way == 'top' & is.na(ntop)){
      stop("Specify number of generations where to generate NAs",
           call. = FALSE)
    }
    
    
    if (way == 'bottom' & is.na(nbottom)){
      stop("Specify number of generations where to generate NAs",
           call. = FALSE)
    }
    
    
    
        #order pedigree
        pedigree <- pedigree0[ , c(1,2,3)]
        names(pedigree) <- c("animal", "sire", "dam")
        ord <- orderPed(pedigree)
        pedigree <- pedigree[order(ord),]

     
        if(!is.element('year', names(data))){
          data$year <- c()
          for ( i in 1:length(data[ ,1])){
                                k <- which(pedigree$animal == data[i, ]$id)
                                data$year[i] <- pedigree0[k, ]$generation           
        
                               }
        }
    
        if (parameter == 'inbreeding'){
          
        resp <- all.vars(model$call)[1]
    
        #calculate initial values (associated with the initial error rate)    
        estimate0 <- initial(model, data, pedigree)
      
        #calculate right model 
        model <- estimate0$model

}

if (parameter == 'heritability'){
  
  resp <- as.character(model$Fixed$formula[2])
  
  #calculate initial values (associated with the initial error rate)    
  estimate0 <- initial_herit(model, data, pedigree)

  
}

     switch (error, 
              'misassignment' = simulated_data <- simul.replace(pedigree, lambda, lambda0, 
                                                              B, data, model, parameter,
                                                              way, prior, nitt,thin,burnin),
              'missing'= simulated_data <- simul.na(pedigree, pedigree0, lambda, lambda0,  
                                                     B, data, model, parameter, 
                                                     way, ntop, nbottom, prior, 
                                                     nitt,thin,burnin)
             )
      
      extrapolation_results <- extrapolation(simulated_data, lambda, lambda0, estimate0, fitting.method, B, parameter)
      plot.simul.Psimex(simulated_data, lambda, lambda0, estimate0, parameter)      
      plot.Psimex(simulated_data, extrapolation_results, lambda, lambda0, estimate0, parameter, fitting.method)
      ret <- list()
      ret$description <- paste('P SIMEX for', parameter, ':', error, 'paternities in pedigree', sep = " ", collapse = NULL)
      ret$error <- error
      ret$fitting.method <- fitting.method
      ret$way <- way
      ret$simul_data <- simulated_data
      ret$extrapolated_data <- extrapolation_results
      ret$lambda <- lambda
      ret$lambda0 <- lambda0
      ret$starting.value <- estimate0
    
      return(ret)
      

    
    
  }
