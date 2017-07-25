initial_herit <-
function(model, data, pedigree) {
  

    random <- all.vars(model$Random$formula)
    if (length(random) > 1)     VE <- rowSums(model$VCV[ ,random])+model$VCV[ , 'units']
    if (length(random) == 1)     VE <- model$VCV[ ,random]+model$VCV[ ,'units']
  
 
  
 
  h_expr <- model$VCV[,"animal"]/VE
  

  
  init <- list()
  init$h0 <- mean(h_expr)
  init$seh0 <- sd(h_expr)
  init$VA0 <- mean(model$VCV[,"animal"])
  init$VE0 <- mean(VE)
 
  
  return(init)
  
}
